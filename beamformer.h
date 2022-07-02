#pragma once

#include <thrust/complex.h>

using namespace std;

// Factor to reduce the time dimension by in our output
// TODO: make this a runtime parameter
const int STI = 8;

/*
  The beamformer metadata is confusing because it's constantly reshaping the data it's
  operating on.
  In particular, it uses an FFT to exchange time resolution for frequency resolution.
  "Coarse channel" refers to the channel index in the input data.
  Within a coarse channel, there are further "fine channels".
  The combination of a coarse channel index plus a fine channel index can give you
  a global frequency index of
  (coarse_channel * fft_size) + fine_channel
  Whenever data is indexed with adjacent indices of [coarse-channel][fine-channel],
  it's equivalent to a single global index of just [channel].
 */
class Beamformer {
 public:
  // The factor for upchannelization. The frequency dimension will be expanded by this
  // amount, and the time dimension will be shrunk by this amount.
  const int fft_size;
  
  // Number of antennas
  const int nants;

  // Number of beams to form
  const int nbeams;

  // Number of time-blocks that the input is provided in.
  // This only affects the format of the input array.
  const int nblocks;
  
  // Number of frequency channels in the input.
  // This will be expanded by a multiplicative factor of fft_size.
  const int num_coarse_channels;

  // Number of polarities
  const int npol;

  // Number of timesteps in the input.
  // This will be reduced by two multiplicative factors, fft_size and STI.
  const int nsamp;

  Beamformer(int fft_size, int nants, int nbeams, int nblocks, int num_coarse_channels,
             int npol, int nsamp);
  ~Beamformer();

  int numOutputChannels() const;
  int numOutputTimesteps() const;
  
  void processInput();

  // Returns a point for the given block where input can be read to
  char* inputPointer(int block);
  
  // These cause a cuda sync so they are slow, only useful for debugging or testing
  thrust::complex<float> getCoefficient(int antenna, int pol, int beam, int coarse_channel) const;
  thrust::complex<float> getFFTBuffer(int pol, int antenna, int coarse_channel,
                                      int time, int last_index) const;
  thrust::complex<float> getPrebeam(int time, int channel, int pol, int antenna) const;
  thrust::complex<float> getVoltage(int time, int channel, int beam, int pol) const;
  float getPower(int beam, int time, int channel) const;

  // Beamforming coefficients, formatted by row-major:
  //   coefficients[coarse-channel][beam][polarity][antenna][real or imag]
  float *coefficients;

  
 private:

  // The input array is unified to the GPU, used to accept data from the previous step
  // in the data processing pipeline.
  //
  // Its format is row-major:
  //   input[block][antenna][coarse-channel][time-within-block][polarity][real or imag]
  //
  // The signed chars should be interpreted as int8.
  signed char *input;
  int input_size;
  
  // The fft_buffer is where we run FFTs as part of the channelization process.
  // The FFT is in-place so there are two different data formats.
  //
  // The convertRaw kernel populates this buffer with row-major:
  //   fft_buffer[pol][antenna][coarse-channel][time]
  // which is equivalent to
  //   fft_buffer[pol][antenna][coarse-channel][time-coarse-index][time-fine-index]
  //
  // The FFT converts the time fine index to a frequency fine index, leaving this
  // buffer with row-major:
  //   fft_buffer[pol][antenna][coarse-channel][time][fine-channel]
  thrust::complex<float>* fft_buffer;
  int fft_buffer_size;
  
  // The channelized input data, ready for beamforming.
  //
  // Its format is row-major:
  //   prebeam[time][channel][polarity][antenna]
  thrust::complex<float>* prebeam;

  // The beamformed data, as voltages.
  //
  // Its format is row-major:
  //   voltage[time][channel][beam][polarity]
  thrust::complex<float>* voltage;

  // The beamformed data, as power.
  //
  // Its format is row-major:
  //   power[beam][time][channel]
  //
  // but its time resolution has been reduced by a factor of (fft_size * STI).
  float* power;
  
};
