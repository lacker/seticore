#include <assert.h>
#include <cuda.h>
#include <cufft.h>
#include <iostream>

#include "beamformer.h"
#include "cuda_util.h"
#include "util.h"

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
    input[block][antenna][coarse-channel][time-within-block][polarity][real or imag]

  to complex-float output with format:
    fft_buffer[polarity][antenna][coarse-channel][time]

  block and time-within-block combine to form a single time index.
 */
__global__ void convertRaw(const signed char* input, thrust::complex<float>* fft_buffer,
                           int nants, int nblocks, int num_coarse_channels, int npol, int nsamp,
                           int time_per_block) {
  int time_within_block = threadIdx.x;
  int block = blockIdx.x;
  int antenna = blockIdx.y;
  int chan = blockIdx.z;
  int time = block * time_per_block + time_within_block;

  for (int pol = 0; pol < npol; ++pol) {
    int input_index = 2 * index5d(block, antenna, nants, chan, num_coarse_channels,
                                  time_within_block, time_per_block, pol, npol);
    int converted_index = index4d(pol, antenna, nants, chan, num_coarse_channels, time, nsamp);

    fft_buffer[converted_index] = thrust::complex<float>
      (input[input_index] * 1.0, input[input_index + 1] * 1.0);
  }
}

/*
  shift converts from the post-FFT format with format:
    fft_buffer[polarity][antenna][coarse-channel][time][fine-channel]

  to a format ready for beamforming:
    prebeam[time][channel][polarity][antenna]

  We also toggle the high bit of the frequency fine channel. Hence "shift".
  This is like swapping the low half and the high half of the output of each FFT.
  It would be great for this comment to explain why this shift is necessary, but, I don't
  understand it myself, so I can't explain it.
 */
__global__ void shift(thrust::complex<float>* fft_buffer, thrust::complex<float>* prebeam,
                      int fft_size, int nants, int npol, int num_coarse_channels,
                      int num_timesteps) {
  int antenna = threadIdx.y;
  int pol = threadIdx.z;
  int fine_chan = blockIdx.x;
  int coarse_chan = blockIdx.y;
  int time = blockIdx.z;

  int output_fine_chan = fine_chan ^ (fft_size >> 1);

  int input_index = index5d(pol, antenna, nants, coarse_chan, num_coarse_channels,
                            time, num_timesteps, fine_chan, fft_size);
  int output_index = index5d(time, coarse_chan, num_coarse_channels, output_fine_chan, fft_size,
                             pol, npol, antenna, nants);

  prebeam[output_index] = fft_buffer[input_index];
}

/*
  Beamforming combines the channelized data with format:
    prebeam[time][coarse-channel][fine-channel][polarity][antenna]

  with the coefficient data, format:
    coefficients[coarse-channel][beam][polarity][antenna][real or imag]

  to generate output beams with format:
    voltage[time][coarse-channel][fine-channel][beam][polarity]

  We combine prebeam with coefficients according to the indices they have in common,
  conjugating the coefficients, and then sum along antenna dimension to reduce.

  TODO: this seems equivalent to a batch matrix multiplication, could we do this with
  a cublas routine?
*/
__global__ void beamform(const thrust::complex<float>* prebeam,
                         const float* coefficients,
                         thrust::complex<float>* voltage,
                         int fft_size, int nants, int nbeams, int num_coarse_channels, int npol,
                         int num_timesteps) {
  int antenna = threadIdx.x;
  int fine_chan = blockIdx.x;
  int coarse_chan = blockIdx.y;
  int beam = blockIdx.z;

  const int MAX_ANTS = 64;
  assert(nants <= MAX_ANTS);
  __shared__ thrust::complex<float> reduced[MAX_ANTS];

  for (int pol = 0; pol < npol; ++pol) {
    int coeff_index = index4d(coarse_chan, beam, nbeams, pol, npol, antenna, nants);
    thrust::complex<float> conjugated = thrust::complex<float>
      (coefficients[2 * coeff_index], -coefficients[2 * coeff_index + 1]);
    for (int time = 0; time < num_timesteps; ++time) {
      int prebeam_index = index5d(time, coarse_chan, num_coarse_channels, fine_chan, fft_size,
                                  pol, npol, antenna, nants);
      reduced[antenna] = prebeam[prebeam_index] * conjugated;

      __syncthreads();

      for (int k = MAX_ANTS / 2; k > 0; k >>= 1) {
        if (antenna < k && antenna + k < nants) {
          reduced[antenna] += reduced[antenna + k];
        }
        __syncthreads();
      }

      if (antenna == 0) {
        int voltage_index = index5d(time, coarse_chan, num_coarse_channels, fine_chan, fft_size,
                                    beam, nbeams, pol, npol);
        voltage[voltage_index] = reduced[0];
      }
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
                               int nbeams, int num_channels, int npol, int num_output_timesteps) {
  int chan = blockIdx.x;
  int beam = blockIdx.y;
  int output_timestep = blockIdx.z;

  int fine_timestep = threadIdx.x;
  assert(fine_timestep < STI);
  int time = output_timestep * STI + fine_timestep;

  assert(2 == npol);
  int pol0_index = index4d(time, chan, num_channels, beam, nbeams, 0, npol);
  int pol1_index = index4d(time, chan, num_channels, beam, nbeams, 1, npol);
  int power_index = index3d(beam, output_timestep, num_output_timesteps, chan, num_channels);

  __shared__ float reduced[STI];
  float real0 = voltage[pol0_index].real();
  float imag0 = voltage[pol0_index].imag();
  float real1 = voltage[pol1_index].real();
  float imag1 = voltage[pol1_index].imag();
  reduced[fine_timestep] = real0 * real0 + imag0 * imag0 + real1 * real1 + imag1 * imag1;

  __syncthreads();

  for (int k = STI / 2; k > 0; k >>= 1) {
    if (fine_timestep < k) {
      reduced[fine_timestep] += reduced[fine_timestep + k];
    }
    __syncthreads();
  }

  if (fine_timestep == 0) {
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
Beamformer::Beamformer(int fft_size, int nants, int nbeams, int nblocks, int num_coarse_channels,
                       int npol, int nsamp)
  : fft_size(fft_size), nants(nants), nbeams(nbeams), nblocks(nblocks),
    num_coarse_channels(num_coarse_channels), npol(npol), nsamp(nsamp) {
  assert(0 == nsamp % (STI * fft_size));
  assert(0 == nsamp % nblocks);
  assert(roundUpToPowerOfTwo(fft_size) == fft_size);

  // The metric that is unchanged by the FFT
  int frame_size = num_coarse_channels * nsamp;
  
  cudaMallocManaged(&input, 2 * nants * npol * frame_size * sizeof(signed char));
  checkCuda("Beamformer input malloc");
  
  cudaMallocManaged(&coefficients, 2 * nants * nbeams * num_coarse_channels * npol * sizeof(float));
  checkCuda("Beamformer coefficients malloc");

  cudaMallocManaged(&fft_buffer, nants * npol * frame_size * sizeof(thrust::complex<float>));
  checkCuda("Beamformer fft_buffer malloc");
  
  cudaMallocManaged(&prebeam, nants * npol * frame_size * sizeof(thrust::complex<float>));
  checkCuda("Beamformer prebeam malloc");

  cudaMallocManaged(&voltage, nbeams * npol * frame_size * sizeof(thrust::complex<float>));
  checkCuda("Beamformer voltage malloc");

  cudaMallocManaged(&power, nbeams * frame_size * sizeof(float) / STI);
  checkCuda("Beamformer power malloc");
}

Beamformer::~Beamformer() {
  cudaFree(input);
  cudaFree(coefficients);
  cudaFree(fft_buffer);
  cudaFree(prebeam);
  cudaFree(voltage);
  cudaFree(power);
}

int Beamformer::numOutputChannels() const {
  return num_coarse_channels * fft_size;
}

int Beamformer::numOutputTimesteps() const {
  return nsamp / (fft_size * STI);
}

char* Beamformer::inputPointer(int block) {
  assert(block < nblocks);

  int total_bytes = nants * num_coarse_channels * nsamp * npol * 2;
  int bytes_per_block = total_bytes / nblocks;

  return ((char*) input) + (block * bytes_per_block);
}

/*
  The caller should first put the data into *input and *coefficients.
 */
void Beamformer::processInput() {
  int time_per_block = nsamp / nblocks;
  dim3 convert_raw_block(time_per_block, 1, 1);
  dim3 convert_raw_grid(nblocks, nants, num_coarse_channels);
  convertRaw<<<convert_raw_grid, convert_raw_block>>>
    (input, fft_buffer, nants, nblocks, num_coarse_channels, npol, nsamp, time_per_block);
  checkCuda("Beamformer convertRaw");
  
  // Run FFTs. TODO: see if there's a faster way
  int num_ffts = nants * npol * num_coarse_channels * nsamp / fft_size;
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

  dim3 shift_block(1, nants, npol);
  dim3 shift_grid(fft_size, num_coarse_channels, nsamp / fft_size);
  shift<<<shift_grid, shift_block>>>(fft_buffer, prebeam, fft_size, nants, npol,
                                     num_coarse_channels, nsamp / fft_size);
  checkCuda("Beamformer shift");
  
  dim3 beamform_block(nants, 1, 1);
  dim3 beamform_grid(fft_size, num_coarse_channels, nbeams);
  beamform<<<beamform_grid, beamform_block>>>(prebeam, coefficients, voltage, fft_size,
                                              nants, nbeams, num_coarse_channels, npol,
                                              nsamp / fft_size);
  checkCuda("Beamformer beamform");

  dim3 power_block(STI, 1, 1);
  dim3 power_grid(numOutputChannels(), nbeams, numOutputTimesteps());
  calculatePower<<<power_grid, power_block>>>
    (voltage, power, nbeams, numOutputChannels(), npol, numOutputTimesteps());
  cout << "power block: " << power_block.x << " " << power_block.y << " " << power_block.z << endl;
  cout << "power grid: " << power_grid.x << " " << power_grid.y << " " << power_grid.z << endl;
  checkCuda("Beamformer calculatePower");
}

thrust::complex<float> Beamformer::getCoefficient(int antenna, int pol, int beam,
                                                  int coarse_channel) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getCoefficient");
  assert(antenna < nants);
  assert(pol < npol);
  assert(beam < nbeams);
  assert(coarse_channel < num_coarse_channels);
  int i = index4d(coarse_channel, beam, nbeams, pol, npol, antenna, nants);
  return thrust::complex<float>(coefficients[2*i], coefficients[2*i+1]);
}

thrust::complex<float> Beamformer::getPrebeam(int time, int channel, int pol, int antenna) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getPrebeam");
  assert(time < nsamp);
  assert(channel < num_coarse_channels * fft_size);
  assert(pol < npol);
  assert(antenna < nants);
  int i = index4d(time, channel, num_coarse_channels * fft_size, pol, npol, antenna, nants);
  return prebeam[i];
}

thrust::complex<float> Beamformer::getVoltage(int time, int channel, int beam, int pol) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getVoltage");
  assert(time < nsamp);
  assert(channel < numOutputChannels());
  assert(beam < nbeams);
  assert(pol < npol);
  int i = index4d(time, channel, numOutputChannels(), beam, nbeams, pol, npol);
  return voltage[i];
}

float Beamformer::getPower(int beam, int time, int channel) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getPower");
  assert(beam < nbeams);
  assert(time < numOutputTimesteps());
  assert(channel < numOutputChannels());
  int i = index3d(beam, time, numOutputTimesteps(), channel, numOutputChannels());
  return power[i];
}
