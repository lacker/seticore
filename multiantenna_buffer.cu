#include <assert.h>
#include "cuda_util.h"
#include "multiantenna_buffer.h"
#include <iostream>

MultiantennaBuffer::MultiantennaBuffer(int num_timesteps, int num_channels,
                                       int num_polarizations, int num_antennas)
  : ComplexBuffer((size_t) num_timesteps * num_channels * num_polarizations * num_antennas),
    num_timesteps(num_timesteps), num_channels(num_channels),
    num_polarizations(num_polarizations), num_antennas(num_antennas) {
}

thrust::complex<float> MultiantennaBuffer::get(int time, int channel,
                                               int polarization, int antenna) const {
  assert(0 <= time && time < num_timesteps);
  assert(0 <= channel && channel < num_channels);
  assert(0 <= polarization && polarization < num_polarizations);
  assert(0 <= antenna && antenna < num_antennas);
  int index = index4d(time, channel, num_channels, polarization, num_polarizations,
                      antenna, num_antennas);
  return get(index);
}
                                       
void MultiantennaBuffer::copyRange(int src_start_channel,
                                   MultiantennaBuffer& dest, int dest_start_time) const {
  assert(src_start_channel >= 0);
  assert(src_start_channel < num_channels);
  assert(src_start_channel + dest.num_channels <= num_channels);
  assert(dest_start_time >= 0);
  assert(dest_start_time + num_timesteps <= dest.num_timesteps);
  assert(num_polarizations == dest.num_polarizations);
  assert(num_antennas == dest.num_antennas);

  int src_index = index4d(0, src_start_channel, num_channels,
                             0, num_polarizations, 0, num_antennas);
  int dest_index = index4d(dest_start_time, 0, dest.num_channels,
                           0, num_polarizations, 0, num_antennas);

  size_t entry_size = sizeof(thrust::complex<float>) * num_polarizations * num_antennas;
  size_t src_pitch = entry_size * num_channels;
  size_t dest_pitch = entry_size * dest.num_channels;

  auto src_ptr = data + src_index;
  auto dest_ptr = dest.data + dest_index;
  
  cudaMemcpy2DAsync(dest_ptr, dest_pitch,
                    src_ptr, src_pitch,
                    dest_pitch, num_timesteps,
                    cudaMemcpyDefault);
  checkCuda("MultiantennaBuffer copyRange");
}
