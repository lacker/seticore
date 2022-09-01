#pragma once

#include "complex_buffer.h"

using namespace std;

/*
  A buffer that stores data for multiple antennas.

  Its format is row-major:
    prebeam[time][channel][polarity][antenna]
 */
class MultiantennaBuffer: public ComplexBuffer {
 public:
  MultiantennaBuffer(int num_timesteps, int num_channels, int num_polarity,
                     int num_antennas);

  // No copying
  MultiantennaBuffer(const MultiantennaBuffer&) = delete;
  MultiantennaBuffer& operator=(MultiantennaBuffer&) = delete;

  const int num_timesteps;
  const int num_channels;
  const int num_polarity;
  const int num_antennas;

  // Causes a cuda sync so it's slow. Only useful for debugging or testing
  using ComplexBuffer::get;
  thrust::complex<float> get(int time, int channel, int pol, int antenna) const;
};
