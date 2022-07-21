#pragma once

#include "filterbank_buffer.h"

using namespace std;

/*
  The MultibeamBuffer stores the contents of a filterbank file in unified memory,
  for multiple beams. This can be a single coarse channel, or a number of channels.
*/
class MultibeamBuffer {
 public:
  const int num_beams;
  const int num_timesteps;
  const int num_channels;

  /*
    Row-major indexed by:
      data[beam][time][freq]
   */
  float* data;

  // Create a managed buffer
  MultibeamBuffer(int num_beams, int num_timesteps, int num_channels);

  ~MultibeamBuffer();

  FilterbankBuffer getBeam(int beam);

  void set(int beam, int time, int channel, float value);
  
  float get(int beam, int time, int channel);

  // Zero out all the data as an asynchronous GPU operation
  void zeroAsync();

  // Asynchronously copy out some data to a separate buffer
  void copyRegionAsync(int beam, int channel_offset, FilterbankBuffer* output);
};
