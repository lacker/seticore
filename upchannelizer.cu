#include <assert.h>
#include "cuda_util.h"
#include <iostream>
#include <thrust/complex.h>
#include "upchannelizer.h"

using namespace std;


Upchannelizer::Upchannelizer(cudaStream_t stream, int fft_size,
                             int num_input_timesteps, int num_coarse_channels,
                             int num_polarity, int num_antennas)
  : stream(stream),
    fft_size(fft_size),
    num_input_timesteps(num_input_timesteps),
    num_coarse_channels(num_coarse_channels),
    num_polarity(num_polarity),
    num_antennas(num_antennas),
    release_input(true) {
  assert(num_input_timesteps % fft_size == 0);

  int batch_size = num_antennas * num_polarity;
  cufftPlan1d(&plan, fft_size, CUFFT_C2C, batch_size);

  checkCuda("Upchannelizer fft planning");
}

Upchannelizer::~Upchannelizer() {
  cufftDestroy(plan);
}

size_t Upchannelizer::requiredInternalBufferSize() const {
  return (size_t) num_antennas * num_coarse_channels * num_polarity * num_input_timesteps;
}

/*
  We convert from int8 input with format:
    input[block][antenna][coarse-channel][time-within-block][polarity][real or imag]

  to complex-float output with format:
    buffer[polarity][antenna][coarse-channel][time]

  block and time-within-block combine to form a single time index.
 */
__global__ void convertRaw(const int8_t* input, int input_size,
                           thrust::complex<float>* buffer, int buffer_size,
                           int nants, int nblocks, int num_coarse_channels,
                           int npol, int nsamp, int time_per_block) {
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
  Run the upchannelization.
  We release the input buffer when we're done with it.
*/
void Upchannelizer::run(DeviceRawBuffer& input, ComplexBuffer& buffer,
                        MultiantennaBuffer& output) {
  assert(input.num_antennas == num_antennas);
  assert(input.num_coarse_channels == num_coarse_channels);
  assert(input.npol == num_polarity);
  assert(input.timesteps_per_block * input.num_blocks == num_input_timesteps);

  assert(buffer.size >= requiredInternalBufferSize());

  assert(output.num_timesteps == num_input_timesteps / fft_size);
  assert(output.num_channels == num_coarse_channels * fft_size);
  assert(output.num_polarity == num_polarity);
  assert(output.num_antennas == num_antennas);

  // Unfortunate overuse of "block"
  int cuda_blocks_per_block =
    (input.timesteps_per_block + CUDA_MAX_THREADS - 1) / CUDA_MAX_THREADS;
  dim3 convert_raw_block(CUDA_MAX_THREADS, 1, 1);
  dim3 convert_raw_grid(cuda_blocks_per_block, input.num_blocks,
                        num_antennas * num_coarse_channels);
  convertRaw<<<convert_raw_grid, convert_raw_block, 0, stream>>>
    (input.data, input.data_size,
     buffer.data, buffer.size,
     num_antennas, input.num_blocks, num_coarse_channels, num_polarity,
     num_input_timesteps, input.timesteps_per_block);
  checkCuda("Beamformer convertRaw");

  // Release the input buffer when we're done with it
  if (release_input) {
    cudaStreamAddCallback(stream, DeviceRawBuffer::staticRelease, &input, 0);
  }
  
  // Run FFTs. TODO: see if there's a faster way
  int num_ffts = num_antennas * num_polarity * num_coarse_channels *
    num_input_timesteps / fft_size;
  int batch_size = num_antennas * num_polarity;
  int num_batches = num_ffts / batch_size;
  for (int i = 0; i < num_batches; ++i) {
    cuComplex* pointer = (cuComplex*) buffer.data + i * batch_size * fft_size;
    cufftExecC2C(plan, pointer, pointer, CUFFT_FORWARD);
  }
  checkCuda("Beamformer fft operation");

  dim3 shift_block(1, num_antennas, num_polarity);
  dim3 shift_grid(fft_size, num_coarse_channels, num_input_timesteps / fft_size);
  shift<<<shift_grid, shift_block, 0, stream>>>
    (buffer.data, output.data, fft_size, num_antennas, num_polarity, num_coarse_channels,
     num_input_timesteps / fft_size);
  checkCuda("Beamformer shift");  
}
