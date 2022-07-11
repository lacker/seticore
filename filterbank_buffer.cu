#include "filterbank_buffer.h"

#include "cuda_util.h"

#include <assert.h>
#include <iostream>

using namespace std;

FilterbankBuffer::FilterbankBuffer(int num_timesteps, int num_channels)
  : num_timesteps(num_timesteps), num_channels(num_channels), managed(true) {
  cudaMallocManaged(&data, sizeof(float) * num_timesteps * num_channels);
  checkCuda("FilterbankBuffer data malloc");
}

FilterbankBuffer::FilterbankBuffer(int num_timesteps, int num_channels, float* data)
  : num_timesteps(num_timesteps), num_channels(num_channels), managed(false), data(data) {
}

FilterbankBuffer::~FilterbankBuffer() {
  if (managed) {
    cudaFree(data);
  }
}

// Set everything to zero
void FilterbankBuffer::zero() {
  memset(data, 0, sizeof(float) * num_timesteps * num_channels);
}

// Inefficient but useful for testing
void FilterbankBuffer::setValue(int time, int channel, float value) {
  assert(0 <= time && time < num_timesteps);
  assert(0 <= channel && channel < num_channels);
  int index = time * num_channels + channel;
  data[index] = value;
}
