#include <assert.h>
#include <cuda.h>
#include <iostream>

#include "cuda_util.h"
#include "taylor.h"

using namespace std;

/*
  Apologies for the length of this comment, but the Taylor tree algorithm is
  fairly complicated for the number of lines of code it is, so it takes
  a while to explain. My hope is that this code will be comprehensible
  for people that have not seen the Taylor tree algorithm before.

  These paths are diagonal paths through the data, touching one element
  per row, using the nearest cell to a straight line. For example, if
  num_timesteps is 4, the paths for each path_offset look like:

  path_offset:   0        1         2          3
                 ..X..    ..X...    ..X....    ..X.....
                 ..X..    ..X...    ...X...    ...X....
                 ..X..    ...X..    ...X...    ....X...
                 ..X..    ...X..    ....X..    .....X..

  At a high level, you can recursively calculate all the sums of these
  paths in O(n log n) operations by running on the top half, then the
  bottom half, then adding them up appropriately.

  The key to understanding this code is to understand the format of
  the buffers: source_buffer and target_buffer.
  source_buffer and target_buffer store the sum along an
  approximately-linear path through the input data. The best way to
  think of them is as three-dimensional arrays, indexed like

  buffer[time_block][path_offset][start_frequency]

  path_offset is the difference between the frequency of the starting
  point of the path and the ending point of the path. It's in the range:
  [0, path_length)

  start_frequency is the index of the starting frequency, in the range:
  [0, num_freqs)

  time_block is a bit weirder. In our recursion, we don't need to keep
  sums for every possible start time. We can cut the number of start
  times in half, every step through the recursion. After n steps of
  recursion, we only need to keep a sum for every 2^n timesteps. So if
  start_time is the index of the starting time, in [0, num_timesteps),
  time_block obeys the relation

  time_block * path_length = start_time

  time_block thus is in the range:
  [0, num_timesteps / path_length)

  and each time block represents data for sums over path_length
  different start times.

  So each buffer holds (num_timesteps * num_freqs) elements.

  When we read input data, it's normally thought of as a two-dimensional
  array, indexed like

  input[time][frequency]

  Since we just pass the buffers around as one-dimensional arrays, this
  is essentially equivalent to thinking of it as a three-dimensional
  array in the above format, with path_offset always equal to zero.

  When we finish running the Taylor tree algorithm, time_block will
  always be zero. Thus we can think of the final output as a
  two-dimensional array as well, indexed like

  output[path_offset][start_frequency]

  It's really just for understanding the intervening steps that it's
  better to think of these buffers as being three-dimensional arrays.

  There's one more detail: drift blocks. So far we've explained the case
  where drift_block = 0, and we are calculating slopes between vertical
  and one horizontal-step-per-vertical step. You can think of this as
  the drift range [0, 1] when measured in units of
  horizontal-step-per-vertical-step. We can use a similar algorithm to
  calculate sums for all slopes in [drift_block, drift_block+1], if we just shift all
  accesses of the kth timestep by an extra drift_block * k.
*/

/*
  Kernel to run one round of the Taylor tree algorithm on an input array.

  We assume that the caller is using a grid tiling such that
  blockIdx.x * blockDim.x + threadIdx.x
  will cover all frequencies.
*/
__global__ void taylorTreeOneStepKernel(const float* source_buffer, float* target_buffer,
                                        int num_timesteps, int num_freqs, int path_length,
                                        int drift_block) {
  assert(path_length <= num_timesteps);
  int freq = blockIdx.x * blockDim.x + threadIdx.x;
  if (freq < 0 || freq >= num_freqs) {
    return;
  }

  taylorOneStepOneChannel(source_buffer, target_buffer,
                          freq, num_timesteps, num_freqs, num_freqs, path_length,
                          drift_block);
}

/*
  Run all rounds of the Taylor tree algorithm.
  buffer1 and buffer2 are two GPU buffers provided to do work.
  Returns the buffer that the eventual output is in.
 */
const float* fullTaylorTree(const float* input, float* buffer1, float* buffer2,
                            int num_timesteps, int num_channels, int drift_block) {
  // This will create one cuda thread per frequency bin
  int grid_size = (num_channels + CUDA_MAX_THREADS - 1) / CUDA_MAX_THREADS;

  // The dataflow among the buffers looks like:
  // input -> buffer1 -> buffer2 -> buffer1 -> buffer2 -> ...
  // We use the aliases source_buffer and target_buffer to make this simpler.
  // In each pass through the upcoming loop, we are reading from
  // source_buffer and writing to target_buffer.
  const float* source_buffer = input;
  float* target_buffer = buffer1;

  // Each pass through the data calculates the sum of paths that are
  // twice as long as the previous path, until we reach our goal,
  // which is paths of length num_timesteps.
  for (int path_length = 2; path_length <= num_timesteps; path_length *= 2) {

    // Invoke cuda kernel
    taylorTreeOneStepKernel<<<grid_size, CUDA_MAX_THREADS>>>
      (source_buffer, target_buffer, num_timesteps, num_channels,
       path_length, drift_block);
    checkCuda("taylorTreeOneStepKernel");

    // Swap buffer aliases to make the old target the new source
    if (target_buffer == buffer1) {
      source_buffer = buffer1;
      target_buffer = buffer2;
    } else if (target_buffer == buffer2) {
      source_buffer = buffer2;
      target_buffer = buffer1;
    } else {
      cerr << "programmer error; control flow should not reach here\n";
      exit(1);
    }
  }

  // The final path sums are in source_buffer because we did one last
  // alias-swap
  return source_buffer;
}
