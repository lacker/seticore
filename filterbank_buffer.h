#pragma once

using namespace std;

/*
  The FilterbankBuffer stores the contents of a filterbank file in unified memory.
  Just one beam. This can be a single coarse channel, or the entire file.
 */
class FilterbankBuffer {
 public:
  const int num_timesteps;
  const int num_channels;

  // Whether the buffer owns its own memory
  const bool managed;

  const int data_size;
  const size_t data_bytes;
  
  /*
    Row-major indexed by:
      data[time][freq]
   */
  float* data;

  // Create a managed buffer
  FilterbankBuffer(int num_timesteps, int num_channels);

  // Create an unmanaged buffer, essentially a view on a pre-existing buffer
  FilterbankBuffer(int num_timesteps, int num_channels, float* data);
  
  ~FilterbankBuffer();

  // Only implicit moving, no implicit copying
  FilterbankBuffer(const FilterbankBuffer&) = delete;
  FilterbankBuffer& operator=(FilterbankBuffer&) = delete;
  FilterbankBuffer(FilterbankBuffer&& other) = default;
  
  void zero();
  void set(int time, int channel, float value);
  float get(int time, int channel) const;

  // Assert two filterbanks are equal.
  // If drift_block is provided, treat time as a drift and only check if it is valid.
  void assertEqual(const FilterbankBuffer& other, int drift_block) const;
  void assertEqual(const FilterbankBuffer& other) const;
};

// Fill a buffer with meaningless data for testing
FilterbankBuffer makeNoisyBuffer(int num_timesteps, int num_channels);

