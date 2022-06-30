#pragma once

#include <thrust/complex.h>

using namespace std;

// Factor to reduce the time dimension by in our output
const int STI = 8;

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
  const int nchans;

  // Number of polarities
  const int npol;

  // Number of timesteps in the input.
  // This will be reduced by two multiplicative factors, fft_size and STI.
  const int nsamp;

  Beamformer(int fft_size, int nants, int nbeams, int nblocks, int nchans, int npol, int nsamp);
  ~Beamformer();

  void processInput();

  // Returns a point for the given block where input can be read to
  char* inputPointer(int block);
  
  // These cause a cuda sync so they are slow, only useful for debugging or testing
  thrust::complex<float> getCoefficient(int antenna, int pol, int beam, int freq) const;
  thrust::complex<float> getPrebeam(int time, int chan, int pol, int antenna) const;
  thrust::complex<float> getVoltage(int time, int chan, int beam, int pol) const;
  float getPower(int beam, int time, int chan) const;

  // Beamforming coefficients, formatted by row-major:
  //   coefficients[frequency][beam][polarity][antenna][real or imag]
  float *coefficients;
  
 private:

  // The input array is unified to the GPU, used to accept data from the previous step
  // in the data processing pipeline.
  //
  // Its format is row-major:
  //   input[block][antenna][frequency][time-within-block][polarity][real or imag]
  //
  // The signed chars should be interpreted as int8.
  signed char *input;

  // The fft_buffer is where we run FFTs as part of the channelization process.
  // The FFT is in-place so there are two different data formats.
  //
  // The convertToFloat kernel populates this buffer with row-major:
  //   fft_buffer[pol][antenna][frequency][time]
  // which is equivalent to
  //   fft_buffer[pol][antenna][frequency][time-coarse-index][time-fine-index]
  //
  // The FFT converts the time fine index to a frequency fine index, leaving this
  // buffer with row-major:
  //   fft_buffer[pol][antenna][frequency-coarse-index][time][frequency-fine-index]
  thrust::complex<float>* fft_buffer;
  
  // The channelized input data, ready for beamforming.
  //
  // Its format is row-major:
  //   prebeam[time][frequency][polarity][antenna]
  thrust::complex<float>* prebeam;

  // The beamformed data, as voltages.
  //
  // Its format is row-major:
  //   voltage[time][frequency][beam][polarity]
  thrust::complex<float>* voltage;

  // The beamformed data, as power.
  //
  // Its format is row-major:
  //   power[beam][time][frequency]
  //
  // but its time resolution has been reduced by a factor of STI.
  float* power;
  
};
