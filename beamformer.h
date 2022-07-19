#pragma once

#include <cufft.h>
#include <thrust/complex.h>

#include "multibeam_buffer.h"
#include "raw_buffer.h"

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
  
  void run(RawBuffer& input, MultibeamBuffer& output, int time_offset);

  // These cause a cuda sync so they are slow, only useful for debugging or testing
  thrust::complex<float> getCoefficient(int antenna, int pol, int beam, int coarse_channel) const;
  thrust::complex<float> getFFTBuffer(int pol, int antenna, int coarse_channel,
                                      int time, int last_index) const;
  thrust::complex<float> getPrebeam(int time, int channel, int pol, int antenna) const;
  thrust::complex<float> getVoltage(int time, int channel, int beam, int pol) const;

  // Beamforming coefficients, formatted by row-major:
  //   coefficients[coarse-channel][beam][polarity][antenna][real or imag]
  thrust::complex<float> *coefficients;
  size_t coefficients_size;

  // The beamformed data, as power.
  float* power;
  size_t power_size;
  
 private:

  // The buffer is reused for a few different stages of the pipeline.
  //
  // The convertRaw kernel populates this buffer with row-major:
  //   buffer[pol][antenna][coarse-channel][time]
  // which is equivalent to
  //   buffer[pol][antenna][coarse-channel][time-coarse-index][time-fine-index]
  //
  // The FFT converts the time fine index in-place to a frequency fine index, leaving this
  // buffer with row-major:
  //   buffer[pol][antenna][coarse-channel][time][fine-channel]
  //
  // After beamforming, the voltage data is also placed here, in row-major:
  //   voltage[time][channel][beam][polarity]
  thrust::complex<float>* buffer;
  size_t buffer_size;
  
  // The channelized input data, ready for beamforming.
  //
  // Its format is row-major:
  //   prebeam[time][channel][polarity][antenna]
  thrust::complex<float>* prebeam;
  size_t prebeam_size;  

  // The plan for the fft.
  cufftHandle plan;

  void runCublasBeamform(int time, int pol);
};
