#include "filterbank_buffer.h"

#include "cuda_util.h"

using namespace std;

FilterbankBuffer::FilterbankBuffer(int num_timesteps, int num_channels)
  : num_timesteps(num_timesteps), num_channels(num_channels) {
  cudaMallocManaged(&data, num_timesteps * num_channels * sizeof(float));
  checkCuda("FilterbankBuffer data malloc");
}

FilterbankBuffer::~FilterbankBuffer() {
  cudaFree(data);
}

