#pragma once

#include <cuda.h>

#include "cuda_util.h"

using namespace std;

/*
  A helper function to run the Taylor tree algorithm for one starting
  frequency, calculating the sums of paths of length `path_length` given the
  sums of paths of length `path_length / 2`.

  This function supports source and target of different dimensions. It checks
  boundaries and does not read or write out of bounds.

  Any in-range path will get a value written to it. So you don't have to
  initialize the output buffer, as long as you recognize that
  output[path_offset][start_frequency] is only valid when the last
  frequency of the path is within bounds, i.e.

  0 <= (num_timesteps - 1) * drift_block + path_offset + start_frequency < num_freqs

  target_buffer and num_target_channels can indicate a slice of a larger array,
  because num_target_channels is only used as a stride.

  This code is based on the original kernel by Franklin Antonio, available at
    https://github.com/UCBerkeleySETI/dedopplerperf/blob/main/CudaTaylor5demo.cu
 */
__host__ __device__ inline void
taylorOneStepOneChannel(const float* source_buffer, float* target_buffer,
                        int chan, int num_timesteps, int num_source_channels,
                        int num_target_channels, int path_length, int drift_block) {
  if (chan < 0 || chan >= num_source_channels || chan >= num_target_channels) {
    // We can't calculate any paths for this channel, since it is out of bounds
    return;
  }

  int num_time_blocks = num_timesteps / path_length;
  for (int time_block = 0; time_block < num_time_blocks; ++time_block) {
    for (int path_offset = path_length - 1; path_offset >= 0; path_offset--) {
      // The recursion calculates sums for a target time block based on two
      // different source time blocks.
      // Data for block b comes from blocks 2b and 2b+1.
      // Remember there are twice as many source time blocks, so the
      // indices are different.

      // The recursion adds up two smaller paths, each with offset (path_offset/2).
      // When path_offset is odd, we need to shift the second path
      // over by an extra one, so that the total adds back up to path_offset.
      // freq_shift thus represents the amount we need to shift the
      // second path.
      int half_offset = path_offset / 2;
      int chan_shift = (path_offset + 1) / 2 + drift_block * path_length / 2;

      if (chan + chan_shift < 0 ||
          chan + chan_shift >= num_source_channels) {
        // We can't calculate this path sum, because the last step requires
        // reading out-of-range data.
        continue;
      }

      // With 3d indices this assignment would look like:
      //   target_buffer[time_block][path_offset][freq] =
      //     source_buffer[2*time_block][half_offset][freq] +
      //     source_buffer[2*time_block+1][half_offset][freq + freq_shift]
      //
      // The data is stored as a 1d array in row-major order, so the
      // actual code here multiplies out some indexes and looks more confusing.
      // In row-major order, for an x*y*z array, the element
      // array[i][j][k]
      // is stored at the location
      // array[((i * y) + j) * z + k]

      // Here, the target buffer has dimensions num_time_blocks * path_length * num_freqs
      // The source buffer has dimensions:
      //   (2 * num_time_blocks) * (path_length / 2) * num_freqs
      // so this line of code is just substituting the appropriate
      // dimensions into the above formula.
      target_buffer[(time_block * path_length + path_offset) * num_target_channels + chan] =
        source_buffer[(time_block * path_length + half_offset) * num_source_channels + chan] +
        source_buffer[(time_block * path_length + half_offset + path_length / 2) * num_source_channels + chan + chan_shift];
    }
  }

}

/*
  Maps a parallelogram-shaped region of source_buffer to a rectangular region
  of target_buffer.
  Handles all timesteps for the given channel.

  Both buffers are shaped row-major [time][channel].
  chan_offset is the channel in source_buffer that maps to 0 in target_buffer.
  drift is how much horizontal shift there is for each vertical step.

  Checks range validity for source_buffer, but not for target buffer.
  For paths that are out of range, the contents of target_buffer are left unchanged.
 */
__host__ __device__ inline void
unmapDrift(const float* source_buffer, float* target_buffer, int num_timesteps,
           int chan, int chan_offset, int num_source_channels, int num_target_channels,
           int drift) {
  for (int time = 0; time < num_timesteps; ++time) {
    int source_chan = chan + chan_offset + time * drift;
    if (source_chan < 0 || source_chan >= num_source_channels) {
      // Out of bounds.
      // We might come back in bounds, though, so keep looping.
      continue;
    }

    target_buffer[index2d(time, chan, num_target_channels)] =
      source_buffer[index2d(time, source_chan, num_source_channels)];
  }
}

const float* basicTaylorTree(const float* source_buffer, float* buffer1, float* buffer2,
                             int num_timesteps, int num_freqs, int drift_block);

void tiledTaylorTree(const float* input, float* output, int num_timesteps,
                     int num_channels, int drift_block);
 
