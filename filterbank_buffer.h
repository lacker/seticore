#pragma once

using namespace std;

/*
  The FilterbankBuffer stores the contents of a filterbank file in GPU memory.
  Just one beam. This can be a single coarse channel, or the entire file.
 */
class FilterbankBuffer {
 public:
  const int num_timesteps;
  const int num_channels;

  /*
    Row-major indexed by:
      data[time][freq]
   */
  float* data;

  FilterbankBuffer(int num_timesteps, int num_channels);
  ~FilterbankBuffer();
};
