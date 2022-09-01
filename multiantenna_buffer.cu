#include <assert.h>
#include "cuda_util.h"
#include "multiantenna_buffer.h"
#include <iostream>

MultiantennaBuffer::MultiantennaBuffer(int num_timesteps, int num_channels,
                                       int num_polarity, int num_antennas)
  : ComplexBuffer((size_t) num_timesteps * num_channels * num_polarity * num_antennas),
    num_timesteps(num_timesteps), num_channels(num_channels),
    num_polarity(num_polarity), num_antennas(num_antennas) {
}

thrust::complex<float> MultiantennaBuffer::get(int time, int channel,
                                               int polarity, int antenna) const {
  assert(0 <= time && time < num_timesteps);
  assert(0 <= channel && channel < num_channels);
  assert(0 <= polarity && polarity < num_polarity);
  assert(0 <= antenna && antenna < num_antennas);
  int index = index4d(time, channel, num_channels, polarity, num_polarity,
                      antenna, num_antennas);
  return get(index);
}
                                       
