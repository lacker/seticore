#include "multibeam_buffer.h"

#include <assert.h>

#include "cuda_util.h"
#include "filterbank_buffer.h"

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

void MultibeamBuffer::set(int beam, int time, int channel, float value) {
  int index = (beam * num_timesteps + time) * num_channels + channel;
  data[index] = value;
}

float MultibeamBuffer::get(int beam, int time, int channel) {
  cudaDeviceSynchronize();
  checkCuda("MultibeamBuffer get");
  assert(beam < num_beams);
  assert(time < num_timesteps);
  assert(channel < num_channels);
  int index = (beam * num_timesteps + time) * num_channels + channel;
  return data[index];
}

void MultibeamBuffer::zeroAsync() {
  size_t data_size = sizeof(float) * num_beams * num_timesteps * num_channels;
  cudaMemsetAsync(data, 0, data_size);
  checkCuda("MultibeamBuffer zeroAsync");
}

void MultibeamBuffer::copyRegionAsync(int beam, int channel_offset,
                                      FilterbankBuffer* output) {
  float* region_start = data + (beam * num_timesteps * num_channels);
  size_t source_pitch = sizeof(float) * num_channels;
  size_t width = sizeof(float) * output->num_channels;
  size_t dest_pitch = width;
  cudaMemcpy2DAsync(output->data, dest_pitch,
                    (void*) region_start, source_pitch,
                    width, num_timesteps,
                    cudaMemcpyDefault);
  checkCuda("MultibeamBuffer copyRegionAsync");
}
