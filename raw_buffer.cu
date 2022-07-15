#include "raw_buffer.h"

#include <assert.h>
#include "cuda_util.h"

using namespace std;

RawBuffer::RawBuffer(int num_blocks, int num_antennas, int num_coarse_channels,
                     int timesteps_per_block, int npol)
  : num_blocks(num_blocks), num_antennas(num_antennas),
    num_coarse_channels(num_coarse_channels),
    timesteps_per_block(timesteps_per_block), npol(npol) {
  data_size = sizeof(int8_t) * num_blocks * num_antennas * num_coarse_channels *
    timesteps_per_block * npol * 2;
  cudaMallocManaged(&data, data_size);
  checkCuda("RawBuffer malloc");
}

RawBuffer::~RawBuffer() {
  cudaFree(data);
}

char* RawBuffer::blockPointer(int block) const {
  assert(block < num_blocks);
  size_t bytes_per_block = data_size / num_blocks;
  return ((char*) data) + (block * bytes_per_block);
}
