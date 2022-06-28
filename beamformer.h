#pragma once

#include <thrust/complex.h>

using namespace std;

// Factor to reduce the time dimension by in our output
const int STI = 8;

class Beamformer {
 public:

  // Number of antennas
  const int nants;

  // Number of beams to form
  const int nbeams;

  // Number of frequency channels
  const int nchans;

  // Number of polarities
  const int npol;

  // Number of timesteps
  const int nsamp;

  Beamformer(int nants, int nbeams, int nchans, int npol, int nsamp);
  ~Beamformer();

  void processInput();
  
  // The input array is unified to the GPU, used to accept data from the previous step
  // in the data processing pipeline.
  //
  // Its format is row-major:
  //   input[antenna][frequency][time][polarity][real or imag]
  //
  // The signed chars should be interpreted as int8.
  signed char *input;

  // Beamforming coefficients, formatted by row-major:
  //   coefficients[frequency][beam][polarity][antenna][real or imag]
  float *coefficients;

  // The input data transposed and converted to complex.
  //
  // Its format is row-major:
  //   transposed[time][frequency][polarity][antenna]
  thrust::complex<float>* transposed;

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
  
  // These cause a cuda sync so they are slow, only useful for debugging or testing
  thrust::complex<float> getCoefficient(int antenna, int pol, int beam, int freq) const;
  thrust::complex<float> getTransposed(int time, int chan, int pol, int antenna) const;
  thrust::complex<float> getVoltage(int time, int chan, int beam, int pol) const;
  
 private:
  // TODO: put some buffers here
};
