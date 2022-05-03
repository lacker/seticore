#include <algorithm>
#include <assert.h>
#include <cuda.h>
#include <fmt/core.h>
#include <functional>
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>

#include "dat_file.h"
#include "h5_file.h"

using namespace std;

static_assert(sizeof(float) == 4, "require 32-bit floats");

const int CUDA_BLOCK_SIZE = 1024;


/*
  Kernel to runs one round of the Taylor tree algorithm, calculating the sums
  of paths of length `path_length`.
  Each thread does the calculations for one frequency.

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
  the buffers, source_buffer and target_buffer.
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

  This kernel is designed to run one thread per frequency in the
  data. If it's trying to calculate the sum of a path that goes out of
  the range, it just won't write to that value. Any in-range path will
  get a value written to it. So you don't have to initialize the
  buffers, as long as you recognize that
  output[path_offset][start_frequency] is only valid when the last
  frequency of the path is within bounds, i.e.

  0 <= (num_timesteps - 1) * drift_block + path_offset + start_frequency < num_freqs

  This code is based on the original kernel by Franklin Antonio, available at
    https://github.com/UCBerkeleySETI/dedopplerperf/blob/main/CudaTaylor5demo.cu
*/
__global__ void taylorTree(const float* source_buffer, float* target_buffer,
                           int num_timesteps, int num_freqs, int path_length, int drift_block) {
  assert(path_length <= num_timesteps);
  int freq = blockIdx.x * blockDim.x + threadIdx.x;
  if (freq < 0 || freq >= num_freqs) {
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
      int freq_shift = (path_offset + 1) / 2 + drift_block * path_length / 2;

      if (freq + freq_shift < 0 ||
          freq + freq_shift >= num_freqs) {
        // We can't calculate this path sum, because it would require
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
      // The source buffer has dimensions (2 * num_time_blocks) * (path_length / 2) * num_freqs
      // so this line of code is just substituting the appropriate
      // dimensions into the above formula.
      target_buffer[(time_block * path_length + path_offset) * num_freqs + freq] =
        source_buffer[(time_block * path_length + half_offset) * num_freqs + freq] +
        source_buffer[(time_block * path_length + half_offset + path_length / 2) * num_freqs + freq + freq_shift];
    }
  }
}


/*
  Gather information about the top hits.

  The eventual goal is for every frequency freq, we want:

  top_path_sums[freq] to contain the largest path sum that starts at freq
  top_drift_blocks[freq] to contain the drift block of that path
  top_path_offsets[freq] to contain the path offset of that path

  path_sums[path_offset][freq] contains one path sum.
  (In row-major order.)
  So we are just taking the max along a column and carrying some
  metadata along as we find it. One thread per freq.

  The function ignores data corresponding to invalid paths. See
  comments in taylorTree for details.
*/
__global__ void findTopPathSums(const float* path_sums, int num_timesteps, int num_freqs,
                                int drift_block, float* top_path_sums,
                                int* top_drift_blocks, int* top_path_offsets) {
  int freq = blockIdx.x * blockDim.x + threadIdx.x;
  if (freq < 0 || freq >= num_freqs) {
    return;
  }

  for (int path_offset = 0; path_offset < num_timesteps; ++path_offset) {
    // Check if the last frequency in this path is out of bounds
    int last_freq = (num_timesteps - 1) * drift_block + path_offset + freq;
    if (last_freq < 0 || last_freq >= num_freqs) {
      // No more of these paths can be valid, either
      return;
    }

    float path_sum = path_sums[num_freqs * path_offset + freq];
    if (path_sum > top_path_sums[freq]) {
      top_path_sums[freq] = path_sum;
      top_drift_blocks[freq] = drift_block;
      top_path_offsets[freq] = path_offset;
    }
  }
}


/*
  Runs a dedoppler algorithm and writes the results to a .dat file.

  input_filename is assumed to be in the .h5 format
  output_filename is where a .dat file will be written with results
  max_drift is the maximum drift we are looking for, in Hz/sec
  snr is the minimum SNR we require to report a signal
 */
void dedoppler(const string& input_filename, const string& output_filename,
               double max_drift, double snr) {
  H5File file(input_filename);
  DatFile output(output_filename, file, max_drift);
  
  // Whether we are calculating drift rates to be compatible with
  // turboseti, as opposed to the way I actually think is correct.
  bool ts_compat = true;
  
  int drift_timesteps = file.num_timesteps - (ts_compat ? 0 : 1);
  double obs_length = drift_timesteps * file.tsamp;
  double drift_rate_resolution = 1e6 * file.foff / obs_length;
  double snr_threshold = 2.0;
  
  // Normalize the max drift in units of "horizontal steps per vertical step"
  double diagonal_drift_rate = drift_rate_resolution * drift_timesteps;
  double normalized_max_drift = max_drift / abs(diagonal_drift_rate);
  int min_drift_block = floor(-normalized_max_drift);
  int max_drift_block = floor(normalized_max_drift);
  
  // For computing Taylor sums, we use three unified memory arrays,
  // each the size of one coarse channel. One array to read the source
  // data, and two buffers to use for Taylor tree calculations.
  float *input, *buffer1, *buffer2;
  cudaMallocManaged(&input, file.coarse_channel_size * file.num_timesteps * sizeof(float));
  cudaMallocManaged(&buffer1, file.coarse_channel_size * file.num_timesteps * sizeof(float));
  cudaMallocManaged(&buffer2, file.coarse_channel_size * file.num_timesteps * sizeof(float));

  // For aggregating top hits, we use three arrays, each the size of
  // one row of the coarse channel. One array to store the largest
  // path sum, one to store which drift block it came from, and one to
  // store its path offset.
  float *top_path_sums;
  int *top_drift_blocks, *top_path_offsets;
  cudaMallocManaged(&top_path_sums, file.coarse_channel_size * sizeof(float));
  cudaMallocManaged(&top_drift_blocks, file.coarse_channel_size * sizeof(int));
  cudaMallocManaged(&top_path_offsets, file.coarse_channel_size * sizeof(int));

  
  // Load and process one coarse channel at a time from the hdf5  
  for (int coarse_channel = 0; coarse_channel < file.num_coarse_channels; ++coarse_channel) {
    file.loadCoarseChannel(coarse_channel, input);

    // For now all cuda operations use the same grid
    int grid_size = (file.coarse_channel_size + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE;

    // Zero out the path sums in between each coarse channel because
    // we pick the top hits separately for each coarse channel
    memset(top_path_sums, 0, file.coarse_channel_size * sizeof(float));

    // We shouldn't need to zero out the index trackers, but the time
    // should be negligible and maybe it helps avoid bugs.
    memset(top_drift_blocks, 0, file.coarse_channel_size * sizeof(int));
    memset(top_path_offsets, 0, file.coarse_channel_size * sizeof(int));
    
    // Calculate distribution statistics
    vector<float> column_sums(file.coarse_channel_size, 0);
    int mid = file.coarse_channel_size / 2;
    for (int row_index = 0; row_index < file.num_timesteps; ++row_index) {
      float* row = input + row_index * file.coarse_channel_size;
      std::transform(column_sums.begin(), column_sums.end(), row,
                     column_sums.begin(), std::plus<float>());

      // Remove the DC spike by making it the average of the adjacent columns
      row[mid] = (row[mid - 1] + row[mid + 1]) / 2.0;
    }
    std::sort(column_sums.begin(), column_sums.end());
    float median;
    if (mid % 2 == 0) {
      median = column_sums[mid];
    } else {
      median = (column_sums[mid - 1] + column_sums[mid]) / 2.0;
    }

    // Use the central 90% to calculate standard deviation
    int begin = ceil(0.05 * column_sums.size());
    int end = floor(0.95 * column_sums.size()) + 1;
    float sum = std::accumulate(column_sums.begin() + begin, column_sums.begin() + end, 0.0);
    float m = sum / (end - begin);
    float accum = 0.0;
    std::for_each(column_sums.begin() + begin, column_sums.begin() + end,
                  [&](const float f) {
                    accum += (f - m) * (f - m);
                  });
    float std_dev = sqrt(accum / (end - begin));

    for (int drift_block = min_drift_block; drift_block <= max_drift_block; ++drift_block) {
    
      // In the Taylor tree algorithm, the dataflow among the buffers looks like:
      // input -> buffer1 -> buffer2 -> buffer1 -> buffer2 -> ...
      // We use the aliases source_buffer and target_buffer to make this simpler.
      // In each pass through the upcoming loop, we are reading from
      // source_buffer and writing to target_buffer.
      float* source_buffer = input;
      float* target_buffer = buffer1;

      // Each pass through the data calculates the sum of paths that are
      // twice as long as the previous path, until we reach our goal,
      // which is paths of length num_timesteps.
      for (int path_length = 2; path_length <= file.num_timesteps; path_length *= 2) {

        // Invoke cuda kernel
        taylorTree<<<grid_size, CUDA_BLOCK_SIZE>>>(source_buffer, target_buffer,
                                                   file.num_timesteps, file.coarse_channel_size,
                                                   path_length, drift_block);

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

      // Invoke cuda kernel
      // The final path sums are in source_buffer because we did one last alias-swap
      findTopPathSums<<<grid_size, CUDA_BLOCK_SIZE>>>(source_buffer, file.num_timesteps,
                                                      file.coarse_channel_size, drift_block,
                                                      top_path_sums, top_drift_blocks,
                                                      top_path_offsets);
    }

    // Now that we have done all the expensive processing for one coarse
    // channel, we can copy the data back out of the GPU.
    cudaDeviceSynchronize();

    // We consider two hits to be duplicates if the distance in their
    // frequency indexes is less than window_size. We only want to
    // output the largest representative of any set of duplicates.
    // window_size is chosen just large enough so that a single bright
    // pixel cannot cause multiple hits.
    // First we break up the data into a set of nonoverlapping
    // windows. Any candidate hit must be the largest within this
    // window.
    float path_sum_threshold = snr_threshold * std_dev + median;
    int window_size = ceil(normalized_max_drift * drift_timesteps);
    for (int i = 0; i * window_size < file.coarse_channel_size; ++i) {
      int candidate_freq = -1;
      float candidate_path_sum = path_sum_threshold;

      for (int j = 0; j < window_size; ++j) {
        int freq = i * window_size + j;
        if (freq >= file.coarse_channel_size) {
          break;
        }
        if (top_path_sums[freq] > candidate_path_sum) {
          // This is the new best candidate of the window
          candidate_freq = freq;
          candidate_path_sum = top_path_sums[freq];
        }
      }
      if (candidate_freq < 0) {
        continue;
      }

      // Check every frequency closer than window_size if we have a candidate
      int window_end = min(file.coarse_channel_size, candidate_freq + window_size);
      bool found_larger_path_sum = false;
      for (int freq = max(0, candidate_freq - window_size + 1); freq < window_end; ++freq) {
        if (top_path_sums[freq] > candidate_path_sum) {
          found_larger_path_sum = true;
          break;
        }
      }
      if (!found_larger_path_sum) {
        // The candidate frequency is the best within its window
        double drift_rate = top_drift_blocks[candidate_freq] * diagonal_drift_rate +
                            top_path_offsets[candidate_freq] * drift_rate_resolution;
        float snr = (candidate_path_sum - median) / std_dev;

        output.reportHit(coarse_channel, candidate_freq, drift_rate, snr);
      }
    }
  }
}
