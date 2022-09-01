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
                                       
void MultiantennaBuffer::copyRange(int src_start_channel,
                                   MultiantennaBuffer& dest, int dest_start_time) const {
  assert(src_start_channel >= 0);
  assert(src_start_channel + dest.num_channels <= num_channels);
  assert(dest_start_time >= 0);
  assert(dest_start_time + num_timesteps <= dest.num_timesteps);
  assert(num_polarity == dest.num_polarity);
  assert(num_antennas == dest.num_antennas);

  int src_index = index4d(0, src_start_channel, num_channels,
                             0, num_polarity, 0, num_antennas);
  int dest_index = index4d(dest_start_time, 0, dest.num_channels,
                           0, num_polarity, 0, num_antennas);

  size_t entry_size = sizeof(thrust::complex<float>) * num_polarity * num_antennas;
  size_t src_pitch = entry_size * num_channels;
  size_t dest_pitch = entry_size * dest.num_channels;

  auto src_ptr = data + src_index;
  auto dest_ptr = dest.data + dest_index;
  
  cudaMemcpy2DAsync(dest_ptr, dest_pitch,
                    src_ptr, src_pitch,
                    dest.num_channels, num_timesteps,
                    cudaMemcpyDefault);
  checkCuda("MultiantennaBuffer copyRange");
}
