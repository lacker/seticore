#include "raw_buffer.h"

#include <assert.h>
#include "cuda_util.h"

using namespace std;

RawBuffer::RawBuffer(bool gpu, int num_blocks, int num_antennas,
                     int num_coarse_channels,
                     int timesteps_per_block, int npol)
  : gpu(gpu), num_blocks(num_blocks), num_antennas(num_antennas),
    num_coarse_channels(num_coarse_channels),
    timesteps_per_block(timesteps_per_block), npol(npol) {
  data_size = sizeof(int8_t) * num_blocks * num_antennas * num_coarse_channels *
    timesteps_per_block * npol * 2;
  if (gpu) {
    cudaMallocManaged(&data, data_size);
    checkCuda("RawBuffer malloc");
  } else {
    data = (int8_t*) malloc(data_size);
  }
}

RawBuffer::~RawBuffer() {
  if (gpu) {
    cudaFree(data);
  } else {
    free(data);
  }
}

char* RawBuffer::blockPointer(int block) const {
  assert(block < num_blocks);
  size_t bytes_per_block = data_size / num_blocks;
  return ((char*) data) + (block * bytes_per_block);
}

RawBuffer RawBuffer::makeSameSize(bool using_gpu) const {
  return RawBuffer(using_gpu, num_blocks, num_antennas,
                   num_coarse_channels, timesteps_per_block, npol);
}

// Should work for any combination of gpu and non-gpu
void RawBuffer::copyFromAsync(const RawBuffer& other) {
  assert(data_size == other.data_size);
  cudaMemcpyAsync(data, other.data, data_size, cudaMemcpyDefault);
}
