#include "raw_buffer.h"

#include <assert.h>
#include "cuda_util.h"

using namespace std;

RawBuffer::RawBuffer(int num_blocks, int num_antennas,
                     int num_coarse_channels,
                     int timesteps_per_block, int npol)
  : num_blocks(num_blocks), num_antennas(num_antennas),
    num_coarse_channels(num_coarse_channels),
    timesteps_per_block(timesteps_per_block), npol(npol) {
  data_size = sizeof(int8_t) * num_blocks * num_antennas * num_coarse_channels *
    timesteps_per_block * npol * 2;
  cudaMallocHost(&data, data_size);
  checkCuda("RawBuffer malloc");
}

RawBuffer::~RawBuffer() {
  cudaFreeHost(data);
}

char* RawBuffer::blockPointer(int block) const {
  assert(block < num_blocks);
  size_t bytes_per_block = data_size / num_blocks;
  return ((char*) data) + (block * bytes_per_block);
}

void RawBuffer::set(int block, int antenna, int coarse_channel,
                    int timestep, int pol, bool imag, int8_t value) {
  int index = 2 * index5d(block, antenna, num_antennas,
                          coarse_channel, num_coarse_channels,
                          timestep, timesteps_per_block, pol, npol);
  if (imag) {
    ++index;
  }
  data[index] = value;
}
