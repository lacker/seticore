#include <assert.h>
#include "cuda_util.h"
#include <iostream>
#include <thrust/complex.h>
#include "upchannelizer.h"

using namespace std;


Upchannelizer::Upchannelizer(cudaStream_t stream, int fft_size,
                             int num_input_timesteps, int num_coarse_channels,
                             int num_polarizations, int num_antennas)
  : stream(stream),
    fft_size(fft_size),
    num_input_timesteps(num_input_timesteps),
    num_coarse_channels(num_coarse_channels),
    num_polarizations(num_polarizations),
    num_antennas(num_antennas),
    release_input(true) {
  assert(fft_size > 0);
  assert(num_input_timesteps > 0);
  assert(num_coarse_channels > 0);
  assert(num_polarizations > 0);
  assert(num_antennas > 0);
  assert(num_input_timesteps % fft_size == 0);

  int batch_size = num_antennas * num_polarizations;
  cufftPlan1d(&plan, fft_size, CUFFT_C2C, batch_size);

  checkCuda("Upchannelizer fft planning");
}

Upchannelizer::~Upchannelizer() {
  cufftDestroy(plan);
}

size_t Upchannelizer::requiredInternalBufferSize() const {
  return (size_t) num_antennas * num_coarse_channels * num_polarizations * num_input_timesteps;
}

/*
  We convert from int8 input with format:
    input[block][antenna][coarse-channel][time-within-block][polarization][real or imag]

  to complex-float output with format:
    buffer[polarization][antenna][coarse-channel][time]

  block and time-within-block combine to form a single time index.
 */
__global__ void convertRaw(const int8_t* input, int input_size,
                           thrust::complex<float>* buffer, int buffer_size,
                           int num_antennas, int nblocks, int num_coarse_channels,
                           int num_polarizations, int nsamp, int time_per_block) {
  int time_within_block = blockIdx.x * CUDA_MAX_THREADS + threadIdx.x;
  if (time_within_block >= time_per_block) {
    return;
  }
  int block = blockIdx.y;
  int antenna = blockIdx.z / num_coarse_channels;
  int chan = blockIdx.z % num_coarse_channels;
  int time = block * time_per_block + time_within_block;
  
  for (int pol = 0; pol < num_polarizations; ++pol) {
    long input_index = 2 * index5d(block, antenna, num_antennas, chan, num_coarse_channels,
                                  time_within_block, time_per_block, pol, num_polarizations);
    long converted_index = index4d(pol, antenna, num_antennas, chan, num_coarse_channels, time, nsamp);

    assert(input_index >= 0);
    assert(converted_index >= 0);
    assert(input_index + 1 < input_size);
    assert(converted_index < buffer_size);
    
    buffer[converted_index] = thrust::complex<float>
      (input[input_index] * 1.0, input[input_index + 1] * 1.0);
  }
}

/*
  shift converts from the post-FFT format with format:
    buffer[polarization][antenna][coarse-channel][time][fine-channel]

  to a format ready for beamforming:
    output[time][channel][polarization][antenna]

  We also toggle the high bit of the frequency fine channel. Hence "shift".
  This is like swapping the low half and the high half of the output of each FFT.
  It would be great for this comment to explain why this shift is necessary, but, I don't
  understand it myself, so I can't explain it.
 */
__global__ void shift(thrust::complex<float>* buffer, int buffer_size,
                      thrust::complex<float>* output, int output_size,
                      int fft_size, int num_antennas, int num_polarizations,
                      int num_coarse_channels, int num_timesteps) {
  int antenna = threadIdx.y;
  int pol = threadIdx.z;
  int fine_chan = blockIdx.x;
  int coarse_chan = blockIdx.y;
  int time = blockIdx.z;

  int output_fine_chan = fine_chan ^ (fft_size >> 1);

  long input_index = index5d(pol, antenna, num_antennas, coarse_chan, num_coarse_channels,
			     time, num_timesteps, fine_chan, fft_size);
  long output_index = index5d(time, coarse_chan, num_coarse_channels,
			      output_fine_chan, fft_size,
			      pol, num_polarizations, antenna, num_antennas);

  assert(input_index >= 0);
  assert(output_index >= 0);
  assert(input_index < buffer_size);
  assert(output_index < output_size);
  output[output_index] = buffer[input_index];
}

/*
  Run the upchannelization.
  We release the input buffer when we're done with it.
*/
void Upchannelizer::run(DeviceRawBuffer& input, ComplexBuffer& buffer,
                        MultiantennaBuffer& output) {
  assert(input.num_antennas == num_antennas);
  assert(input.num_coarse_channels == num_coarse_channels);
  assert(input.num_polarizations == num_polarizations);
  assert(input.timesteps_per_block * input.num_blocks == num_input_timesteps);

  assert(buffer.size >= requiredInternalBufferSize());

  assert(output.num_timesteps == num_input_timesteps / fft_size);
  assert(output.num_channels == num_coarse_channels * fft_size);
  assert(output.num_polarizations == num_polarizations);
  assert(output.num_antennas == num_antennas);

  // Unfortunate overuse of "block"
  int cuda_blocks_per_block =
    (input.timesteps_per_block + CUDA_MAX_THREADS - 1) / CUDA_MAX_THREADS;
  dim3 convert_raw_block(CUDA_MAX_THREADS, 1, 1);
  dim3 convert_raw_grid(cuda_blocks_per_block, input.num_blocks,
                        num_antennas * num_coarse_channels);
  convertRaw<<<convert_raw_grid, convert_raw_block, 0, stream>>>
    (input.data, input.size,
     buffer.data, buffer.size,
     num_antennas, input.num_blocks, num_coarse_channels, num_polarizations,
     num_input_timesteps, input.timesteps_per_block);
  checkCuda("Beamformer convertRaw");

  // Release the input buffer when we're done with it
  if (release_input) {
    cudaStreamAddCallback(stream, DeviceRawBuffer::staticRelease, &input, 0);
  }
  
  // Run FFTs. TODO: see if there's a faster way
  int num_ffts = num_antennas * num_polarizations * num_coarse_channels *
    num_input_timesteps / fft_size;
  int batch_size = num_antennas * num_polarizations;
  int num_batches = num_ffts / batch_size;
  for (int i = 0; i < num_batches; ++i) {
    cuComplex* pointer = (cuComplex*) buffer.data + i * batch_size * fft_size;
    cufftExecC2C(plan, pointer, pointer, CUFFT_FORWARD);
  }
  checkCuda("Beamformer fft operation");

  dim3 shift_block(1, num_antennas, num_polarizations);
  dim3 shift_grid(fft_size, num_coarse_channels, num_input_timesteps / fft_size);
  shift<<<shift_grid, shift_block, 0, stream>>>
    (buffer.data, buffer.size, output.data, output.size, fft_size, num_antennas,
     num_polarizations, num_coarse_channels, num_input_timesteps / fft_size);
  checkCuda("Beamformer shift");  
}

int Upchannelizer::numOutputChannels() const {
  return num_coarse_channels * fft_size;
}

int Upchannelizer::numOutputTimesteps() const {
  return num_input_timesteps / fft_size;
}
