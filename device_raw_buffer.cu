#include "device_raw_buffer.h"

#include <assert.h>
#include "cuda_util.h"

using namespace std;

DeviceRawBuffer::DeviceRawBuffer(int num_blocks, int num_antennas,
                     int num_coarse_channels,
                     int timesteps_per_block, int npol)
  : num_blocks(num_blocks), num_antennas(num_antennas),
    num_coarse_channels(num_coarse_channels),
    timesteps_per_block(timesteps_per_block), npol(npol) {
  data_size = sizeof(int8_t) * num_blocks * num_antennas * num_coarse_channels *
    timesteps_per_block * npol * 2;
  cudaMalloc(&data, data_size);
  checkCuda("DeviceRawBuffer malloc");
}

DeviceRawBuffer::~DeviceRawBuffer() {
  cudaFree(data);
}

// Should work for any combination of gpu and non-gpu
void DeviceRawBuffer::copyFromAsync(const RawBuffer& other) {
  assert(data_size == other.data_size);
  cudaMemcpyAsync(data, other.data, data_size, cudaMemcpyHostToDevice);
}

