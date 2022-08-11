#include "filterbank_buffer.h"

#include <assert.h>
#include <fmt/core.h>
#include <iostream>

#include "cuda_util.h"
#include "util.h"

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

float FilterbankBuffer::get(int time, int channel) const {
  cudaDeviceSynchronize();
  checkCuda("FilterbankBuffer get");
  int index = time * num_channels + channel;
  return data[index];
}

void FilterbankBuffer::assertEqual(const FilterbankBuffer& other, int drift_block) const {
  assert(num_timesteps == other.num_timesteps);
  assert(num_channels == other.num_channels);
  for (int drift = 0; drift < num_timesteps; ++drift) {
    for (int chan = 0; chan < num_channels; ++chan) {
      int last_chan = chan + (num_timesteps - 1) * drift_block + drift;
      if (last_chan < 0 || last_chan >= num_channels) {
        continue;
      }
      assertFloatEq(get(drift, chan), other.get(drift, chan),
                    fmt::format("data[{}][{}]", drift, chan));
    }
  }
}

// Make a filterbank buffer with a bit of deterministic noise so that
// normalization doesn't make everything infinite SNR.
FilterbankBuffer makeNoisyBuffer(int num_timesteps, int num_channels) {
  FilterbankBuffer buffer(num_timesteps, num_channels);
  buffer.zero();
  for (int chan = 0; chan < buffer.num_channels; ++chan) {
    buffer.set(0, chan, 0.1 * chan / buffer.num_channels);
  }
  return buffer;
}

