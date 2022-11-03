#include "multibeam_buffer.h"

#include <assert.h>

#include "cuda_util.h"
#include "filterbank_buffer.h"
#include "util.h"

using namespace std;

MultibeamBuffer::MultibeamBuffer(int num_beams, int num_timesteps, int num_channels,
                                 int num_write_timesteps)
  : num_beams(num_beams), num_timesteps(num_timesteps), num_channels(num_channels),
    num_write_timesteps(num_write_timesteps) {
  assert(num_write_timesteps <= num_timesteps);
  size_t bytes = sizeof(float) * size();
  cudaMallocManaged(&data, bytes);
  if (bytes > 2000000) {
    cout << "multibeam buffer memory: " << prettyBytes(bytes) << endl;
  }
  checkCudaMalloc("MultibeamBuffer", bytes);

  cudaStreamCreateWithFlags(&prefetch_stream, cudaStreamNonBlocking);
  checkCuda("MultibeamBuffer stream init");
}

MultibeamBuffer::MultibeamBuffer(int num_beams, int num_timesteps, int num_channels)
  : MultibeamBuffer(num_beams, num_timesteps, num_channels, num_timesteps) {}

MultibeamBuffer::~MultibeamBuffer() {
  cudaFree(data);
}

long MultibeamBuffer::size() const {
  return num_beams * num_timesteps * num_channels;
}

FilterbankBuffer MultibeamBuffer::getBeam(int beam) {
  assert(0 <= beam && beam < num_beams);
  int beam_size = num_timesteps * num_channels;
  return FilterbankBuffer(num_timesteps, num_channels, data + beam * beam_size);
}

void MultibeamBuffer::set(int beam, int time, int channel, float value) {
  int index = index3d(beam, time, num_timesteps, channel, num_channels);
  data[index] = value;
}

float MultibeamBuffer::get(int beam, int time, int channel) {
  cudaDeviceSynchronize();
  checkCuda("MultibeamBuffer get");
  assert(beam < num_beams);
  assert(time < num_timesteps);
  assert(channel < num_channels);
  int index = index3d(beam, time, num_timesteps, channel, num_channels);
  return data[index];
}

void MultibeamBuffer::zeroAsync() {
  size_t size = sizeof(float) * num_beams * num_timesteps * num_channels;
  cudaMemsetAsync(data, 0, size);
  checkCuda("MultibeamBuffer zeroAsync");
}

void MultibeamBuffer::copyRegionAsync(int beam, int channel_offset,
                                      FilterbankBuffer* output) {
  float* region_start = data + (beam * num_timesteps * num_channels) + channel_offset;
  size_t source_pitch = sizeof(float) * num_channels;
  size_t width = sizeof(float) * output->num_channels;
  size_t dest_pitch = width;
  cudaMemcpy2DAsync(output->data, dest_pitch,
                    (void*) region_start, source_pitch,
                    width, num_timesteps,
                    cudaMemcpyDefault);
  checkCuda("MultibeamBuffer copyRegionAsync");
}

void MultibeamBuffer::hintWritingTime(int time) {
  prefetchStripes(0, 0, time - num_write_timesteps, time + 2 * num_write_timesteps);
}

void MultibeamBuffer::hintReadingBeam(int beam) {
  prefetchStripes(beam - 1, beam + 1, 0, num_write_timesteps);
}

// Truncates [first_time, last_time]
void MultibeamBuffer::prefetchRange(int beam, int first_time, int last_time,
                                    int destination_device) {
  if (first_time < 0) {
    first_time = 0;
  }
  if (last_time >= num_timesteps) {
    last_time = num_timesteps - 1;
  }
  if (first_time > last_time) {
    return;
  }

  long start_index = index3d(beam, first_time, num_timesteps, 0, num_channels);
  size_t prefetch_size = sizeof(float) * (last_time - first_time + 1) * num_channels;  
  cudaMemPrefetchAsync(data + start_index, prefetch_size, destination_device,
                       prefetch_stream);
  checkCuda("MultibeamBuffer prefetchRange");
}

/*
  Prefetch so that the range of times from [first_time, last_time] and the
  range of beams from [first_beam, last_beam] will be on the GPU.
  The set of values to prefetch is like an intersection of vertical
  and horizontal stripes.

  For (time, beam) pairs that are adjacent to this set in memory, or
  at the first or last times, we don't do any explicit
  prefetching. For (time, beam) pairs that are not adjacent, we
  prefetch them to the CPU. The adjacent regions are like a buffer
  where CUDA can decide that they can go part on the CPU
  and part on the GPU, since the CUDA library rounds prefetch ranges
  to the nearest page size.

  Truncates inputs that are out of range.
 */
void MultibeamBuffer::prefetchStripes(int first_beam, int last_beam,
                                      int first_time, int last_time) {
  if (first_beam < 0) {
    first_beam = 0;
  }
  if (last_beam >= num_beams) {
    last_beam = num_beams - 1;
  }
  if (first_time < 0) {
    first_time = 0;
  }
  if (last_time >= num_timesteps) {
    last_time = num_timesteps - 1;
  }
  
  const int gpu_id = 0;
  for (int beam = 0; beam < num_beams; ++beam) {
    if (first_beam <= beam && beam <= last_beam) {
      // Prefetch this entire beam to the GPU
      prefetchRange(beam, 0, num_timesteps - 1, gpu_id);
      continue;
    }

    // Prefetch the desired time range to the GPU
    prefetchRange(beam, first_time, last_time, gpu_id);

    // I think this could be shrunk to the CPU page size
    int margin = num_write_timesteps;
    
    // Prefetch earlier times to the CPU, with a margin
    prefetchRange(beam, margin, first_time - margin, cudaCpuDeviceId);

    // Prefetch later times to the CPU, with a margin
    prefetchRange(beam, last_time + margin, num_timesteps - margin, cudaCpuDeviceId);
  }
}
