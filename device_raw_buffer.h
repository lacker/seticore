#pragma once

#include <condition_variable>
#include <cuda.h>
#include <cuda_runtime.h>
#include <mutex>

#include "raw_buffer.h"

using namespace std;

/*
  The DeviceRawBuffer is designed for a single producer thread to produce
  RawBuffers, and a single consumer thread to use the buffer for GPU calculations.

  Its format is row-major:
    input[block][antenna][coarse-channel][time-within-block][polarity][real or imag]

  The states are:
    unused: no producer or consumer is using the data
    copying: data is currently being copied into this buffer from a raw buffer
    ready: the consumer thread is using the data

  Access to state is protected by the mutex.
 */
enum class DeviceRawBufferState { unused, copying, ready };

class DeviceRawBuffer {
 public:
  const int num_blocks;
  const int num_antennas;
  const int num_coarse_channels;
  const int timesteps_per_block;
  const int npol;

  int8_t* data;
  size_t data_size;
  
  DeviceRawBuffer(int num_blocks, int num_antennas, int num_coarse_channels,
                  int timesteps_per_block, int npol);

  ~DeviceRawBuffer();

  void copyFromAsync(const RawBuffer& source);

  // Wait until a copy finishes
  void waitUntilReady();

  // Wait until the consumer thread is done with this buffer
  void waitUntilUnused();
  
  // Return the state from ready to unused
  void release();

  static void CUDART_CB staticCopyCallback(cudaStream_t stream,
                                           cudaError_t status,
                                           void *device_raw_buffer);

  static void CUDART_CB staticRelease(cudaStream_t stream,
                                      cudaError_t status,
                                      void *device_raw_buffer);
  
 private:
  // This stream is just for this raw buffer; the state tracks concurrency
  cudaStream_t stream;

  // Called when a copy completes
  void copyCallback();
  
  DeviceRawBufferState state;
  mutex m;
  condition_variable cv;
};
