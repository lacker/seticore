#pragma once

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
  // The factor for upchannelization. The frequency dimension will be expanded by
  // this amount, and the time dimension will be shrunk by this amount.
  const int fft_size;

  // Number of antennas
  const int nants;

  // Number of time-blocks that the input is provided in.
  // This only affects the format of the input array.
  const int nblocks;

  // Number of frequency channels in the input.
  // This will be expanded by a multiplicative factor of fft_size.
  const int num_coarse_channels;

  // Number of polarities
  const int npol;

  // Number of timesteps in the input.
  // This will be reduced by a factor of fft_size.
  const int nsamp;
  
  // The upchannelizer runs all its oprations on one cuda stream.
  const cudaStream_t stream;

  Upchannelizer(cudaStream_t stream, int fft_size, int nants, int nblocks,
                int num_coarse_channels, int npol, int nsamp);
  ~Upchannelizer();

  // Calculate how large of an internal buffer is needed, measured in
  // amount of complex values
  size_t requiredInternalBufferSize() const;

  // A buffer of the appropriate size must be provided once before running,
  // and then it can be reused for all subsequent operation.
  setInternalBuffer(shared_ptr<ComplexBuffer> buffer);

  void run(DeviceRawBuffer& input, MultiantennaBuffer& output);
}
