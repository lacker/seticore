#include "filterbank_buffer.h"

#include "cuda_util.h"

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

