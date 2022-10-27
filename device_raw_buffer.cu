#include "device_raw_buffer.h"

#include <assert.h>
#include "cuda_util.h"

using namespace std;


DeviceRawBuffer::DeviceRawBuffer(int num_blocks, int num_antennas,
                                 int num_coarse_channels,
                                 int timesteps_per_block, int num_polarizations)
  : num_blocks(num_blocks), num_antennas(num_antennas),
    num_coarse_channels(num_coarse_channels),
    timesteps_per_block(timesteps_per_block), num_polarizations(num_polarizations),
    state(DeviceRawBufferState::unused) {
  size = sizeof(int8_t) * num_blocks * num_antennas * num_coarse_channels *
    timesteps_per_block * num_polarizations * 2;
  cudaMalloc(&data, size);
  cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking);
  checkCuda("DeviceRawBuffer init");
}

DeviceRawBuffer::~DeviceRawBuffer() {
  cudaFree(data);
  cudaStreamDestroy(stream);
}

// Should only be called from the producer thread
void DeviceRawBuffer::copyFromAsync(const RawBuffer& other) {
  assert(size == other.size);
  waitUntilUnused();

  unique_lock<mutex> lock(m);
  assert(state == DeviceRawBufferState::unused);
  state = DeviceRawBufferState::copying;
  lock.unlock();
  // Nobody waits on copying state, so no need to notify
  
  cudaMemcpyAsync(data, other.data, size, cudaMemcpyHostToDevice, stream);
  cudaStreamAddCallback(stream, DeviceRawBuffer::staticCopyCallback, this, 0);
}

void DeviceRawBuffer::waitUntilReady() {
  unique_lock<mutex> lock(m);
  while (state != DeviceRawBufferState::ready) {
    cv.wait(lock);
  }
}

void DeviceRawBuffer::waitUntilUnused() {
  unique_lock<mutex> lock(m);
  while (state != DeviceRawBufferState::unused) {
    cv.wait(lock);
  }
}

void DeviceRawBuffer::release() {
  unique_lock<mutex> lock(m);
  assert(state == DeviceRawBufferState::ready);
  state = DeviceRawBufferState::unused;
  lock.unlock();
  cv.notify_all();
}

void CUDART_CB DeviceRawBuffer::staticCopyCallback(cudaStream_t stream,
                                                   cudaError_t status,
                                                   void *device_raw_buffer) {
  assert(status == cudaSuccess);
  DeviceRawBuffer* object = (DeviceRawBuffer*) device_raw_buffer;
  object->copyCallback();
}

void CUDART_CB DeviceRawBuffer::staticRelease(cudaStream_t stream,
                                              cudaError_t status,
                                              void *device_raw_buffer) {
  assert(status == cudaSuccess);
  DeviceRawBuffer* object = (DeviceRawBuffer*) device_raw_buffer;
  object->release();
}

void DeviceRawBuffer::copyCallback() {
  // Advance state to "ready"
  unique_lock<mutex> lock(m);
  assert(state == DeviceRawBufferState::copying);
  state = DeviceRawBufferState::ready;
  lock.unlock();
  cv.notify_all();
}
