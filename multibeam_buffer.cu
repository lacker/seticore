#include "multibeam_buffer.h"

#include <assert.h>
#include "cuda_util.h"

using namespace std;

MultibeamBuffer::MultibeamBuffer(int num_beams, int num_timesteps, int num_channels)
  : num_beams(num_beams), num_timesteps(num_timesteps), num_channels(num_channels) {
  cudaMallocManaged(&data, sizeof(float) * num_beams * num_timesteps * num_channels);
  checkCuda("MultibeamBuffer data malloc");
}

MultibeamBuffer::~MultibeamBuffer() {
  cudaFree(data);
}

FilterbankBuffer MultibeamBuffer::getBeam(int beam) {
  assert(0 <= beam && beam < num_beams);
  int beam_size = num_timesteps * num_channels;
  return FilterbankBuffer(num_timesteps, num_channels, data + beam * beam_size);
}

float MultibeamBuffer::getFloat(int beam, int time, int channel) {
  cudaDeviceSynchronize();
  checkCuda("MultibeamBuffer getFloat");
  assert(beam < num_beams);
  assert(time < num_timesteps);
  assert(channel < num_channels);
  int index = (beam * num_timesteps + time) * num_channels + channel;
  return data[index];
}
