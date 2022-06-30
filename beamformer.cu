#include <assert.h>
#include <cuda.h>
#include <cufft.h>
#include <iostream>

#include "beamformer.h"
#include "cuda_util.h"

using namespace std;

// Helper to calculate a 3d row-major index, ie for:
//   arr[a][b][c]
__host__ __device__ int index3d(int a, int b, int b_end, int c, int c_end) {
  return ((a * b_end) + b) * c_end + c;
}

// Helper to calculate a 4d row-major index, ie for:
//   arr[a][b][c][d]
__host__ __device__ int index4d(int a, int b, int b_end, int c, int c_end, int d, int d_end) {
  return index3d(a, b, b_end, c, c_end) * d_end + d;
}

// Helper to calculate a 5d row-major index, ie for:
//   arr[a][b][c][d][e]
__host__ __device__ int index5d(int a, int b, int b_end, int c, int c_end, int d, int d_end,
                                int e, int e_end) {
  return index4d(a, b, b_end, c, c_end, d, d_end) * e_end + e;
}


/*
  We convert from int8 input with format:
    input[block][antenna][frequency][time-within-block][polarity][real or imag]

  to complex-float output with format:
    fft_buffer[polarity][antenna][frequency][time]

  block and time-within-block combine to form a single time index.
 */
__global__ void convertRaw(const signed char* input, thrust::complex<float>* fft_buffer,
                           int nants, int nblocks, int nchans, int npol, int nsamp,
                           int time_per_block) {
  int time_within_block = threadIdx.x;
  int block = blockIdx.x;
  int antenna = blockIdx.y;
  int chan = blockIdx.z;
  int time = block * time_per_block + time_within_block;

  for (int pol = 0; pol < npol; ++pol) {
    int input_index = 2 * index5d(block, antenna, nants, chan, nchans,
                                  time_within_block, time_per_block, pol, npol);
    int converted_index = index4d(pol, antenna, nants, chan, nchans, time, nsamp);

    fft_buffer[converted_index] = thrust::complex<float>
      (input[input_index] * 1.0, input[input_index + 1] * 1.0);
  }
}

/*
  Beamforming combines the channelized data with format:
    channelized[time][frequency][polarity][antenna]
  with the coefficient data, format:
    coefficients[frequency][beam][polarity][antenna][real or imag]
  to generate output beams with format:
    voltage[time][frequency][beam][polarity]

  We combine channelized with coefficients according to the indices they have in common,
  conjugating the coefficients, and then sum along antenna dimension to reduce.

  TODO: this seems equivalent to a batch matrix multiplication, could we do this with
  a cublas routine?
*/
__global__ void beamform(const thrust::complex<float>* channelized,
                         const float* coefficients,
                         thrust::complex<float>* voltage,
                         int nants, int nbeams, int nchans, int npol) {
  int antenna = threadIdx.x;
  int chan = blockIdx.y;
  int time = blockIdx.x;
  int beam = blockIdx.z;

  const int MAX_ANTS = 64;
  assert(nants <= MAX_ANTS);
  __shared__ thrust::complex<float> reduced[MAX_ANTS];

  for (int pol = 0; pol < npol; ++pol) {
    int channelized_index = index4d(time, chan, nchans, pol, npol, antenna, nants);
    int coeff_index = index4d(chan, beam, nbeams, pol, npol, antenna, nants);
    thrust::complex<float> conjugated = thrust::complex<float>
      (coefficients[2 * coeff_index], -coefficients[2 * coeff_index + 1]);
    reduced[antenna] = channelized[channelized_index] * conjugated;

    __syncthreads();

    for (int k = MAX_ANTS / 2; k > 0; k >>= 1) {
      if (antenna < k && antenna + k < nants) {
        reduced[antenna] += reduced[antenna + k];
      }
      __syncthreads();
    }

    if (antenna == 0) {
      int voltage_index = index4d(time, chan, nchans, beam, nbeams, pol, npol);
      voltage[voltage_index] = reduced[0];
    }
  }
}

/*
  We calculate power and shrink the data at the same time.
  Every polarity and every window of STI adjacent timesteps gets reduced to a single
  power value, by adding the norm of each complex voltage.

  The input voltages have format:
    voltage[time][frequency][beam][polarity]

  and the outpower power has format:
    power[beam][time][frequency]

  where the time dimension has shrunk by a factor of STI, now indexed by [0, nwin).

  TODO: this also seems equivalent to a batch matrix multiplication. could we do this
  with a cublas routine?
 */
__global__ void calculatePower(const thrust::complex<float>* voltage,
                               float* power,
                               int nbeams, int nchans, int npol, int nwindows) {
  int chan = blockIdx.z;
  int beam = blockIdx.y;
  int window = blockIdx.x;

  int subtime = threadIdx.x;
  assert(subtime < STI);
  int time = window * STI + subtime;

  assert(2 == npol);
  int pol0_index = index4d(time, chan, nchans, beam, nbeams, 0, npol);
  int pol1_index = index4d(time, chan, nchans, beam, nbeams, 1, npol);
  int power_index = index3d(beam, window, nwindows, chan, nchans);

  __shared__ float reduced[STI];
  float real0 = voltage[pol0_index].real();
  float imag0 = voltage[pol0_index].imag();
  float real1 = voltage[pol1_index].real();
  float imag1 = voltage[pol1_index].imag();
  reduced[subtime] = real0 * real0 + imag0 * imag0 + real1 * real1 + imag1 * imag1;

  __syncthreads();

  for (int k = STI / 2; k > 0; k >>= 1) {
    if (subtime < k) {
      reduced[subtime] += reduced[subtime + k];
    }
    __syncthreads();
  }

  if (subtime == 0) {
    power[power_index] = reduced[0];
  }
}

/*
  The Beamformer encapsulates the GPU memory allocations we use for beamforming.
  The workflow is to create a beamformer for a particular set of dimensions,
  use it to form many beams, and then destruct it when we want to free the memory.

  TODO: nants and npol are specified twice, one by the recipe file and one by the raw input.
  We should really check to ensure they are the same and handle it cleanly if they aren't.
 */
Beamformer::Beamformer(int fft_size, int nants, int nbeams, int nblocks, int nchans,
                       int npol, int nsamp)
  : fft_size(fft_size), nants(nants), nbeams(nbeams), nblocks(nblocks), nchans(nchans),
    npol(npol), nsamp(nsamp) {
  assert(0 == nsamp % (STI * fft_size));
  assert(0 == nsamp % nblocks);
  
  cudaMallocManaged(&input, 2 * nants * nchans * npol * nsamp * sizeof(signed char));
  checkCuda("Beamformer input malloc");
  
  cudaMallocManaged(&coefficients, 2 * nants * nbeams * nchans * npol * sizeof(float));
  checkCuda("Beamformer coefficients malloc");

  cudaMallocManaged(&fft_buffer, nants * nchans * npol * nsamp * sizeof(thrust::complex<float>));
  checkCuda("Beamformer fft_buffer malloc");
  
  cudaMallocManaged(&channelized, nants * nchans * npol * nsamp * sizeof(thrust::complex<float>));
  checkCuda("Beamformer channelized malloc");

  cudaMallocManaged(&voltage, nbeams * nchans * npol * nsamp * sizeof(thrust::complex<float>));
  checkCuda("Beamformer voltage malloc");

  cudaMallocManaged(&power, nbeams * nchans * nsamp * sizeof(float) / STI);
  checkCuda("Beamformer power malloc");
}

Beamformer::~Beamformer() {
  cudaFree(input);
  cudaFree(coefficients);
  cudaFree(channelized);
  cudaFree(voltage);
  cudaFree(power);
}

char* Beamformer::inputPointer(int block) {
  assert(block < nblocks);

  int total_bytes = nants * nchans * nsamp * npol * 2;
  int bytes_per_block = total_bytes / nblocks;

  return ((char*) input) + (block * bytes_per_block);
}

/*
  The caller should first put the data into *input and *coefficients.
 */
void Beamformer::processInput() {

  int time_per_block = nsamp / nblocks;
  dim3 convert_raw_block(time_per_block, 1, 1);
  dim3 convert_raw_grid(nblocks, nants, nchans);
  convertRaw<<<convert_raw_grid, convert_raw_block>>>
    (input, fft_buffer, nants, nblocks, nchans, npol, nsamp, time_per_block);
  checkCuda("Beamformer convertRaw");
  
  // Run FFTs. TODO: see if there's a faster way
  int num_ffts = nants * npol * nchans * nsamp / fft_size;
  int batch_size = nants * npol;
  int num_batches = num_ffts / batch_size;
  cufftHandle plan;
  cufftPlan1d(&plan, fft_size, CUFFT_C2C, batch_size);
  checkCuda("Beamformer fft planning");
  for (int i = 0; i < num_batches; ++i) {
    cuComplex* pointer = (cuComplex*) fft_buffer + i * batch_size * fft_size;
    cufftExecC2C(plan, pointer, pointer, CUFFT_FORWARD);
  }
  cufftDestroy(plan);
  checkCuda("Beamformer fft operation");

  // TODO: run a shift kernel
  
  dim3 beamform_block(nants, 1, 1);
  dim3 beamform_grid(nsamp, nchans, nbeams);
  beamform<<<beamform_grid, beamform_block>>>(channelized, coefficients, voltage,
                                              nants, nbeams, nchans, npol);
  checkCuda("Beamformer beamform");

  assert(0 == nsamp % STI);
  int nwindows = nsamp / STI;
  dim3 power_block(STI, 1, 1);
  dim3 power_grid(nwindows, nbeams, nchans);
  calculatePower<<<power_grid, power_block>>>(voltage, power, nbeams, nchans, npol, nwindows);
  checkCuda("Beamformer calculatePower");
}

thrust::complex<float> Beamformer::getCoefficient(int antenna, int pol, int beam, int freq) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getCoefficient");
  assert(antenna < nants);
  assert(pol < npol);
  assert(beam < nbeams);
  assert(freq < nchans);
  int i = index4d(freq, beam, nbeams, pol, npol, antenna, nants);
  return thrust::complex<float>(coefficients[2*i], coefficients[2*i+1]);
}

thrust::complex<float> Beamformer::getChannelized(int time, int chan, int pol, int antenna) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getChannelized");
  assert(time < nsamp);
  assert(chan < nchans);
  assert(pol < npol);
  assert(antenna < nants);
  int i = index4d(time, chan, nchans, pol, npol, antenna, nants);
  return channelized[i];
}

thrust::complex<float> Beamformer::getVoltage(int time, int chan, int beam, int pol) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getVoltage");
  assert(time < nsamp);
  assert(chan < nchans);
  assert(beam < nbeams);
  assert(pol < npol);
  int i = index4d(time, chan, nchans, beam, nbeams, pol, npol);
  return voltage[i];
}

float Beamformer::getPower(int beam, int time, int chan) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getPower");
  assert(beam < nbeams);
  assert(time < nsamp / STI);
  assert(chan < nchans);
  int i = index3d(beam, time, nsamp / STI, chan, nchans);
  return power[i];
}
