#include <assert.h>
#include "cublas_v2.h"
#include <cuda.h>
#include <cufft.h>
#include <iostream>

#include "beamformer.h"
#include "cuda_util.h"
#include "util.h"

using namespace std;

const cuComplex COMPLEX_ONE = make_cuComplex(1.0, 0.0);
const cuComplex COMPLEX_ZERO = make_cuComplex(0.0, 0.0);

/*
  We convert from int8 input with format:
    input[block][antenna][coarse-channel][time-within-block][polarity][real or imag]

  to complex-float output with format:
    buffer[polarity][antenna][coarse-channel][time]

  block and time-within-block combine to form a single time index.
 */
__global__ void convertRaw(const int8_t* input, int input_size,
                           thrust::complex<float>* buffer, int buffer_size,
                           int nants, int nblocks, int num_coarse_channels, int npol, int nsamp,
                           int time_per_block) {
  int time_within_block = blockIdx.x * CUDA_MAX_THREADS + threadIdx.x;
  if (time_within_block >= time_per_block) {
    return;
  }
  int block = blockIdx.y;
  int antenna = blockIdx.z / num_coarse_channels;
  int chan = blockIdx.z % num_coarse_channels;
  int time = block * time_per_block + time_within_block;
  
  for (int pol = 0; pol < npol; ++pol) {
    int input_index = 2 * index5d(block, antenna, nants, chan, num_coarse_channels,
                                  time_within_block, time_per_block, pol, npol);
    int converted_index = index4d(pol, antenna, nants, chan, num_coarse_channels, time, nsamp);

    assert(input_index + 1 < input_size);
    assert(converted_index < buffer_size);
    
    buffer[converted_index] = thrust::complex<float>
      (input[input_index] * 1.0, input[input_index + 1] * 1.0);
  }
}

/*
  shift converts from the post-FFT format with format:
    buffer[polarity][antenna][coarse-channel][time][fine-channel]

  to a format ready for beamforming:
    prebeam[time][channel][polarity][antenna]

  We also toggle the high bit of the frequency fine channel. Hence "shift".
  This is like swapping the low half and the high half of the output of each FFT.
  It would be great for this comment to explain why this shift is necessary, but, I don't
  understand it myself, so I can't explain it.
 */
__global__ void shift(thrust::complex<float>* buffer, thrust::complex<float>* prebeam,
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

  prebeam[output_index] = buffer[input_index];
}

/*
  Beamforming combines the channelized data with format:
    prebeam[time][coarse-channel][fine-channel][polarity][antenna]

  with the coefficient data, format:
    coefficients[coarse-channel][beam][polarity][antenna]

  to generate output beams with format:
    voltage[time][polarity][coarse-channel][fine-channel][beam]

  We combine prebeam with coefficients according to the indices they have in common,
  not conjugating the coefficients because we expect them to already be in the
  correct conjugation for multiplying, and then sum along antenna dimension to reduce.
*/
__global__ void beamform(const thrust::complex<float>* prebeam,
                         const thrust::complex<float>* coefficients,
                         thrust::complex<float>* voltage,
                         int fft_size, int nants, int nbeams, int num_coarse_channels,
                         int npol, int num_timesteps,
                         int prebeam_size, int voltage_size, int coefficients_size) {
  int antenna = threadIdx.x;
  int fine_chan = blockIdx.x;
  int coarse_chan = blockIdx.y;
  int beam = blockIdx.z;

  const int MAX_ANTS = 64;
  assert(nants <= MAX_ANTS);
  __shared__ thrust::complex<float> reduced[MAX_ANTS];

  for (int pol = 0; pol < npol; ++pol) {
    int coeff_index = index4d(coarse_chan, beam, nbeams, pol, npol, antenna, nants);
    assert(2 * coeff_index + 1 < coefficients_size);
    thrust::complex<float> conjugated = coefficients[coeff_index];
    for (int time = 0; time < num_timesteps; ++time) {
      int prebeam_index = index5d(time, coarse_chan, num_coarse_channels, fine_chan,
                                  fft_size, pol, npol, antenna, nants);
      assert(prebeam_index < prebeam_size);
      assert(antenna < MAX_ANTS);
      reduced[antenna] = prebeam[prebeam_index] * conjugated;

      __syncthreads();

      for (int k = MAX_ANTS / 2; k > 0; k >>= 1) {
        if (antenna < k && antenna + k < nants) {
          assert(antenna + k < MAX_ANTS);
          reduced[antenna] += reduced[antenna + k];
        }
        __syncthreads();
      }

      if (antenna == 0) {
        int voltage_index = index5d(time, pol, npol, coarse_chan, num_coarse_channels,
                                    fine_chan, fft_size, beam, nbeams);
        assert(voltage_index < voltage_size);
        voltage[voltage_index] = reduced[0];
      }
    }
  }
}

/*
  Runs beamforming just for the provided time and polarity, using cublas batch
  matrix multiplication.
  See the comment on the beamform kernel.

  This API is not immediately intuitive. See this blog post:
    https://developer.nvidia.com/blog/cublas-strided-batched-matrix-multiply/
  in particular their explanation of gemmStridedBatched.

  To convert into the notation used by the blog post:

  A = coefficients (which we will transpose)
  B = prebeam
  C = voltage
  m = beam
  n = fine channel
  p = coarse channel
  k = antenna
  alpha = 1.0
  beta = 0.0

  We are converting five-dimensional data plus four-dimensional data into
  five-dimensional output, so we fix time and polarity for each call to cublas.
  We could have fixed fine channel instead of time, or coarse channel instead of
  polarity, but it seems better to fix the smaller dimensions.
  Honestly, there are so many possibilities for how to arrange this data, I have
  not even begun to test all the ways it could work.
 */
void Beamformer::runCublasBeamform(int time, int pol) {
  // Calculate where the matrices start
  int coeff_offset = index4d(0, 0, nbeams, pol, npol, 0, nants);
  auto coeff_start = (const cuComplex*) (coefficients + coeff_offset);
  int prebeam_offset = index5d(time, 0, num_coarse_channels, 0, fft_size, pol, npol,
                               0, nants);
  auto prebeam_start = (const cuComplex*) (prebeam + prebeam_offset);
  int voltage_offset = index5d(time, pol, npol, 0, num_coarse_channels, 0, fft_size,
                               0, nbeams);
  auto voltage_start = (cuComplex*) (buffer + voltage_offset);

  // Calculate strides
  // ldA, the A-m stride (since we are transposing. normally it would be A-k)
  int coeff_beam_stride = index4d(0, 1, nbeams, 0, npol, 0, nants);
  // strideA, the A-p stride
  int coeff_coarse_stride = index4d(1, 0, nbeams, 0, npol, 0, nants);
  // ldB, the B-n stride
  int prebeam_fine_stride = index5d(0, 0, num_coarse_channels, 1, fft_size,
                                    0, npol, 0, nants);
  // strideB, the B-p stride
  int prebeam_coarse_stride = index5d(0, 1, num_coarse_channels, 0, fft_size,
                                      0, npol, 0, nants);
  // ldC, the C-n stride
  int voltage_fine_stride = index5d(0, 0, npol, 0, num_coarse_channels,
                                    1, fft_size, 0, nbeams);
  // strideC, the C-p stride
  int voltage_coarse_stride = index5d(0, 0, npol, 1, num_coarse_channels,
                                      0, fft_size, 0, nbeams);

  cublasCgemm3mStridedBatched
    (cublas_handle, 
     CUBLAS_OP_T, CUBLAS_OP_N,
     nbeams, fft_size, nants,
     &COMPLEX_ONE,
     coeff_start, coeff_beam_stride, coeff_coarse_stride, 
     prebeam_start, prebeam_fine_stride, prebeam_coarse_stride,
     &COMPLEX_ZERO,
     voltage_start, voltage_fine_stride, voltage_coarse_stride,
     num_coarse_channels);
  checkCuda("Beamformer runCublasBeamform");
}

/*
  We calculate power and shrink the data at the same time.
  Every polarity and every window of STI adjacent timesteps gets reduced to a single
  power value, by adding the norm of each complex voltage.

  The input voltages have format: 
    voltage[time][polarity][frequency][beam]

  and the output power has format:
    power[beam][time][frequency]

  where the time dimension has shrunk by a factor of STI, now indexed by [0, nwin).
  We add time_offset to all the times in the output.

  TODO: this also seems equivalent to a batch matrix multiplication. could we do this
  with a cublas routine?
 */
__global__ void calculatePower(const thrust::complex<float>* voltage,
                               float* power,
                               int nbeams, int num_channels, int npol,
                               int num_output_timesteps, int time_offset) {
  int chan = blockIdx.x;
  int beam = blockIdx.y;
  int coarse_timestep = blockIdx.z;
  int output_timestep = coarse_timestep + time_offset;
  
  int fine_timestep = threadIdx.x;
  assert(fine_timestep < STI);
  int time = coarse_timestep * STI + fine_timestep;

  assert(2 == npol);
  int pol0_index = index4d(time, 0, npol, chan, num_channels, beam, nbeams);
  int pol1_index = index4d(time, 1, npol, chan, num_channels, beam, nbeams);
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

  TODO: nants and npol are specified twice, once by the recipe file and once by the input.
  We should check to ensure they are the same and handle it cleanly if they aren't.
 */
Beamformer::Beamformer(int fft_size, int nants, int nbeams, int nblocks,
                       int num_coarse_channels, int npol, int nsamp)
  : fft_size(fft_size), nants(nants), nbeams(nbeams), nblocks(nblocks),
    num_coarse_channels(num_coarse_channels), npol(npol), nsamp(nsamp) {
  assert(0 == nsamp % (STI * fft_size));
  assert(0 == nsamp % nblocks);
  assert(roundUpToPowerOfTwo(fft_size) == fft_size);

  int frame_size = num_coarse_channels * nsamp;
  
  coefficients_size = 2 * nants * nbeams * num_coarse_channels * npol;
  size_t coefficients_bytes = coefficients_size * sizeof(float);
  cudaMallocManaged(&coefficients, coefficients_bytes);
  checkCuda("Beamformer coefficients malloc");
 
  size_t fft_buffer_size = nants * npol * frame_size;
  size_t voltage_size = nbeams * npol * frame_size;
  buffer_size = max(fft_buffer_size, voltage_size);
  size_t buffer_bytes = buffer_size * sizeof(thrust::complex<float>);
  cudaMallocManaged(&buffer, buffer_bytes);
  checkCuda("Beamformer buffer malloc");

  prebeam_size = nants * npol * frame_size;
  size_t prebeam_bytes = prebeam_size * sizeof(thrust::complex<float>);
  cudaMallocManaged(&prebeam, prebeam_bytes);
  checkCuda("Beamformer prebeam malloc");

  power_size = nbeams * frame_size / STI;
  size_t power_bytes = power_size * sizeof(float);
  cudaMallocManaged(&power, power_bytes);
  checkCuda("Beamformer power malloc");

  int batch_size = nants * npol;
  cufftPlan1d(&plan, fft_size, CUFFT_C2C, batch_size);
  checkCuda("Beamformer fft planning");

  cublasCreate(&cublas_handle);
  checkCuda("Beamformer cublas handle");
  
  size_t total_bytes = coefficients_bytes + buffer_bytes + prebeam_bytes + power_bytes;
  cout << "beamformer memory: " << prettyBytes(total_bytes) << endl;
}

Beamformer::~Beamformer() {
  cudaFree(coefficients);
  cudaFree(buffer);
  cudaFree(prebeam);
  cudaFree(power);
  cufftDestroy(plan);
  cublasDestroy(cublas_handle);
}

int Beamformer::numOutputChannels() const {
  return num_coarse_channels * fft_size;
}

int Beamformer::numOutputTimesteps() const {
  return nsamp / (fft_size * STI); 
}

/*
  Power from beamforming the input is written into output, with an offset
  of time_offset.

  The format of the input is row-major:
    input[block][antenna][coarse-channel][time-within-block][polarity][real or imag]

  The format of the output is row-major:
     power[beam][time][channel]
  but its time resolution has been reduced by a factor of (fft_size * STI).
 */
void Beamformer::run(RawBuffer& input, MultibeamBuffer& output, int time_offset) {
  int time_per_block = nsamp / nblocks;
  // Unfortunate overuse of "block"
  int cuda_blocks_per_block = (time_per_block + CUDA_MAX_THREADS - 1) / CUDA_MAX_THREADS;
  dim3 convert_raw_block(CUDA_MAX_THREADS, 1, 1);
  dim3 convert_raw_grid(cuda_blocks_per_block, nblocks, nants * num_coarse_channels);
  convertRaw<<<convert_raw_grid, convert_raw_block>>>
    (input.data, input.data_size,
     buffer, buffer_size,
     nants, nblocks, num_coarse_channels, npol, nsamp, time_per_block);
  checkCuda("Beamformer convertRaw");

  // Run FFTs. TODO: see if there's a faster way
  int num_ffts = nants * npol * num_coarse_channels * nsamp / fft_size;
  int batch_size = nants * npol;
  int num_batches = num_ffts / batch_size;
  for (int i = 0; i < num_batches; ++i) {
    cuComplex* pointer = (cuComplex*) buffer + i * batch_size * fft_size;
    cufftExecC2C(plan, pointer, pointer, CUFFT_FORWARD);
  }
  checkCuda("Beamformer fft operation");

  dim3 shift_block(1, nants, npol);
  dim3 shift_grid(fft_size, num_coarse_channels, nsamp / fft_size);
  shift<<<shift_grid, shift_block>>>(buffer, prebeam, fft_size, nants, npol,
                                     num_coarse_channels, nsamp / fft_size);
  checkCuda("Beamformer shift");
  
  dim3 beamform_block(nants, 1, 1);
  dim3 beamform_grid(fft_size, num_coarse_channels, nbeams);
  beamform<<<beamform_grid, beamform_block>>>(prebeam, coefficients, buffer, fft_size,
                                              nants, nbeams, num_coarse_channels, npol,
                                              nsamp / fft_size, prebeam_size, buffer_size,
                                              coefficients_size);

  dim3 power_block(STI, 1, 1);
  dim3 power_grid(numOutputChannels(), nbeams, numOutputTimesteps());
  calculatePower<<<power_grid, power_block>>>
    (buffer, output.data, nbeams, numOutputChannels(), npol, numOutputTimesteps(),
     time_offset);
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
  return coefficients[i];
}

// The last index can be either time's fine index or the fine channel index, depending
// on whether it's pre-FFT or post-FFT.
thrust::complex<float> Beamformer::getFFTBuffer(int pol, int antenna, int coarse_channel,
                                                int time, int last_index) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getFFTBuffer");
  assert(pol < npol);
  assert(antenna < nants);
  assert(coarse_channel < num_coarse_channels);
  assert(time * fft_size < nsamp);
  assert(last_index < fft_size);
  int i = index5d(pol, antenna, nants, coarse_channel, num_coarse_channels,
                  time, nsamp / fft_size, last_index, fft_size);
  return buffer[i];
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

thrust::complex<float> Beamformer::getVoltage(int time, int pol, int channel, int beam) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getVoltage");
  assert(time < nsamp);
  assert(pol < npol);
  assert(channel < numOutputChannels());
  assert(beam < nbeams);
  int i = index4d(time, pol, npol, channel, numOutputChannels(), beam, nbeams);
  return buffer[i];
}

