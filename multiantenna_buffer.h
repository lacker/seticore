#pragma once

#include "complex_buffer.h"

using namespace std;

/*
  A buffer that stores data for multiple antennas.

  Its format is row-major:
    prebeam[time][channel][polarization][antenna]
 */
class MultiantennaBuffer: public ComplexBuffer {
 public:
  MultiantennaBuffer(int num_timesteps, int num_channels, int num_polarizations,
                     int num_antennas);

  // No copying
  MultiantennaBuffer(const MultiantennaBuffer&) = delete;
  MultiantennaBuffer& operator=(MultiantennaBuffer&) = delete;

  const int num_timesteps;
  const int num_channels;
  const int num_polarizations;
  const int num_antennas;

  // Causes a cuda sync so it's slow. Only useful for debugging or testing
  using ComplexBuffer::get;
  thrust::complex<float> get(int time, int channel, int pol, int antenna) const;

  // Copies a range of channels to another buffer.
  // Assumes that the range to copy is the entire width of the other buffer.
  // Uses the default cuda stream and copies asynchronously.
  void copyRange(int src_start_channel,
                 MultiantennaBuffer& dest, int dest_start_time) const;
};
