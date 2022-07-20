#include "filterbank_buffer.h"

#include "cuda_util.h"

#include <assert.h>
#include <iostream>

using namespace std;

// Creates a buffer that owns its own memory.
FilterbankBuffer::FilterbankBuffer(int num_timesteps, int num_channels)
  : num_timesteps(num_timesteps), num_channels(num_channels), managed(true),
    data_size(num_timesteps * num_channels),
    data_bytes(sizeof(float) * data_size) {
  cudaMallocManaged(&data, data_bytes);
  checkCuda("FilterbankBuffer data malloc");
}

// Creates a buffer that is a view on memory owned by the caller.
FilterbankBuffer::FilterbankBuffer(int num_timesteps, int num_channels, float* data)
  : num_timesteps(num_timesteps), num_channels(num_channels), managed(false),
    data_size(num_timesteps * num_channels),
    data_bytes(sizeof(float) * data_size), data(data) {
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
void FilterbankBuffer::set(int time, int channel, float value) {
  assert(0 <= time && time < num_timesteps);
  assert(0 <= channel && channel < num_channels);
  int index = time * num_channels + channel;
  data[index] = value;
}

float FilterbankBuffer::get(int time, int channel) {
  cudaDeviceSynchronize();
  checkCuda("FilterbankBuffer get");
  int index = time * num_channels + channel;
  return data[index];
}
