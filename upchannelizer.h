#pragma once

#include "complex_buffer.h"
#include <cufft.h>
#include "device_raw_buffer.h"
#include "multiantenna_buffer.h"

using namespace std;

/*
  The Upchannelizer converts data from raw buffers to a format with higher frequency
  resolution, lower time resolution, by doing a fast Fourier transform as well as
  some other data shuffling-around.

  It needs an internal buffer to do its processing. This buffer can be shared with
  other operations that are running on the same cuda stream as long as a single thread
  queues up the operations for both the Upchannelizer and the other operations.
 */
class Upchannelizer {
 public:
  // The upchannelizer runs all its operations on one cuda stream.
  const cudaStream_t stream;

  // The factor for upchannelization. The frequency dimension will be expanded by
  // this amount, and the time dimension will be shrunk by this amount.
  const int fft_size;

  // Number of timesteps in the input.
  // This will be reduced by a factor of fft_size.
  const int num_input_timesteps;
  
  // Number of frequency channels in the input.
  // This will be expanded by a multiplicative factor of fft_size.
  const int num_coarse_channels;

  // The same in input and output
  const int num_polarity;
  const int num_antennas;
  
  Upchannelizer(cudaStream_t stream, int fft_size,
                int num_input_timesteps, int num_coarse_channels,
                int num_polarity, int num_antennas);
  ~Upchannelizer();

  // Calculate how large of an internal buffer is needed, measured in
  // amount of complex values
  size_t requiredInternalBufferSize() const;

  void run(DeviceRawBuffer& input, ComplexBuffer& buffer, MultiantennaBuffer& output);

  // Whether to release inputs when we're done with it.
  bool release_input;
  
 private:
  // The plan for the fft.
  cufftHandle plan;
};
