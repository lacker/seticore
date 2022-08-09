#pragma once

#include <cuda.h>
#include <cuda_runtime.h>
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

  // The number of timesteps that will be written in a batch.
  // Used only to optimize managed memory prefetching, so it can be a guess.
  const int num_write_timesteps;
  
  /*
    Row-major indexed by:
      data[beam][time][freq]
   */
  float* data;

  // Create a managed buffer
  MultibeamBuffer(int num_beams, int num_timesteps, int num_channels,
                  int num_write_timesteps);
  MultibeamBuffer(int num_beams, int num_timesteps, int num_channels);  

  ~MultibeamBuffer();

  FilterbankBuffer getBeam(int beam);

  void set(int beam, int time, int channel, float value);
  
  float get(int beam, int time, int channel);

  // Zero out all the data as an asynchronous GPU operation
  void zeroAsync();

  // Asynchronously copy out some data to a separate buffer.
  // Uses default cuda stream.
  void copyRegionAsync(int beam, int channel_offset, FilterbankBuffer* output);

  // Call this when you are writing this time
  void hintWritingTime(int time);

  // Call this when you are reading this beam
  void hintReadingBeam(int beam);

private:
  // This stream is just for prefetching.
  cudaStream_t prefetch_stream;

  void prefetchRange(int beam, int first_time, int last_time, int destinationDevice);
  
  void prefetchStripes(int first_beam, int last_beam, int first_time, int last_time);
};
