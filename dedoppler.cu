#include <algorithm>
#include <assert.h>
#include "cublas_v2.h"
#include <cuda.h>
#include <functional>
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>

#include "cuda_util.h"
#include "dedoppler.h"
#include "util.h"

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
                           int num_timesteps, int num_freqs, int path_length,
                           int drift_block) {
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
  Cublas-based equivalent to the taylorTree kernel.

  This is best understood not on its own, but by first understanding taylorTree
  and then understanding how this is a port of that.
 */
void Dedopplerer::runCublasTaylorTree(const float* source_buffer, float* target_buffer,
                                      int path_length, int drift_block) {
  assert(path_length <= rounded_num_timesteps);

  cublasSetStream(cublas_handle, 0);
  
  int num_time_blocks = num_timesteps / path_length;
  for (int path_offset = 0; path_offset < path_length; ++path_offset) {
    int half_offset = path_offset / 2;
    int freq_shift = (path_offset + 1) / 2 + drift_block * path_length / 2;

    // Calculate the range that freq can be in:
    //   [freq_begin, freq_end)
    // freq can be any integer for which freq and freq + freq_shift are both in:
    //   [0, num_channels)
    // Note that due to drift blocks, freq_shift can be positive or negative.
    int freq_begin = max(0, -freq_shift);
    int freq_end = min(num_channels, num_channels - freq_shift);

    // Figure out the starting pointers for the matrices to add.
    // To get the start point, we plug in time_block = 0, freq = freq_begin
    // into the array update line in the taylorTree kernel.
    float* target_start = target_buffer +
      (path_offset * num_channels + freq_begin);
    const float* source_start_1 = source_buffer +
      (half_offset * num_channels + freq_begin);
    const float* source_start_2 = source_buffer +
      ((half_offset + path_length / 2) * num_channels + freq_begin + freq_shift);

    // See:
    //   https://docs.nvidia.com/cuda/cublas/#cublas-lt-t-gt-geam
    // The API is column major so m is frequency, n is time block.
    int stride = path_length * num_channels;
    cublasSgeam(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
                freq_end - freq_begin, num_time_blocks,
                &FLOAT_ONE, source_start_1, stride,
                &FLOAT_ONE, source_start_2, stride,
                target_start, stride);
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
  Sum the columns of a two-dimensional array.
  input is a (num_timesteps x num_freqs) array, stored in row-major order.
  sums is an array of size num_freqs.
 */
__global__ void sumColumns(const float* input, float* sums, int num_timesteps, int num_freqs) {
  int freq = blockIdx.x * blockDim.x + threadIdx.x;
  if (freq < 0 || freq >= num_freqs) {
    return;
  }
  sums[freq] = 0.0;
  for (int i = freq; i < num_timesteps * num_freqs; i += num_freqs) {
    sums[freq] += input[i];
  }
}


/*
  The Dedopplerer encapsulates the logic of dedoppler search. In particular it manages
  the needed GPU memory so that we can reuse the same memory allocation for different searches.
 */
Dedopplerer::Dedopplerer(int num_timesteps, int num_channels, double foff, double tsamp,
                         bool has_dc_spike)
    : num_timesteps(num_timesteps), num_channels(num_channels), foff(foff), tsamp(tsamp),
      has_dc_spike(has_dc_spike), print_hits(false), print_hit_summary(false),
      use_cublas_taylor_tree(false) {
  assert(num_timesteps > 1);
  rounded_num_timesteps = roundUpToPowerOfTwo(num_timesteps);
  drift_timesteps = rounded_num_timesteps - 1;

  drift_rate_resolution = 1e6 * foff / (drift_timesteps * tsamp);
    
  // Allocate everything we need for GPU processing

  cudaMalloc(&buffer1, num_channels * rounded_num_timesteps * sizeof(float));
  checkCuda("Dedopplerer buffer1 malloc");

  cudaMalloc(&buffer2, num_channels * rounded_num_timesteps * sizeof(float));
  checkCuda("Dedopplerer buffer2 malloc");
  
  cudaMalloc(&gpu_column_sums, num_channels * sizeof(float));
  cudaMallocHost(&cpu_column_sums, num_channels * sizeof(float));
  checkCuda("Dedopplerer column_sums malloc");
  
  cudaMalloc(&gpu_top_path_sums, num_channels * sizeof(float));
  cudaMallocHost(&cpu_top_path_sums, num_channels * sizeof(float));
  checkCuda("Dedopplerer top_path_sums malloc");
   
  cudaMalloc(&gpu_top_drift_blocks, num_channels * sizeof(int));
  cudaMallocHost(&cpu_top_drift_blocks, num_channels * sizeof(int));
  checkCuda("Dedopplerer top_drift_blocks malloc");
  
  cudaMalloc(&gpu_top_path_offsets, num_channels * sizeof(int));
  cudaMallocHost(&cpu_top_path_offsets, num_channels * sizeof(int));
  checkCuda("Dedopplerer top_path_offsets malloc");

  cublasCreate(&cublas_handle);
  checkCuda("Dedopplerer cublas handle");
}

Dedopplerer::~Dedopplerer() {
  cudaFree(buffer1);
  cudaFree(buffer2);
  cudaFree(gpu_column_sums);
  cudaFreeHost(cpu_column_sums);
  cudaFree(gpu_top_path_sums);
  cudaFreeHost(cpu_top_path_sums);
  cudaFree(gpu_top_drift_blocks);
  cudaFreeHost(cpu_top_drift_blocks);
  cudaFree(gpu_top_path_offsets);
  cudaFreeHost(cpu_top_path_offsets);
}

// This implementation is an ugly hack
size_t Dedopplerer::memoryUsage() const {
  return num_channels * rounded_num_timesteps * sizeof(float) * 2
    + num_channels * (2 * sizeof(float) + 2 * sizeof(int));
}

/*
  Runs dedoppler search on the input buffer.
  Output is appended to the output vector.
  
  All processing of the input buffer happens on the GPU, so it doesn't need to
  start off with host and device synchronized when search is called, it can still
  have GPU processing pending.
*/
void Dedopplerer::search(const FilterbankBuffer& input,
                         int beam, int coarse_channel,
                         double max_drift, double min_drift, double snr_threshold,
                         vector<DedopplerHit>* output) {
  assert(input.num_timesteps == rounded_num_timesteps);
  assert(input.num_channels == num_channels);

  // Normalize the max drift in units of "horizontal steps per vertical step"
  double diagonal_drift_rate = drift_rate_resolution * drift_timesteps;
  double normalized_max_drift = max_drift / abs(diagonal_drift_rate);
  int min_drift_block = floor(-normalized_max_drift);
  int max_drift_block = floor(normalized_max_drift);

  // For now all cuda operations use the same grid.
  // This will create one cuda thread per frequency bin
  int grid_size = (num_channels + CUDA_MAX_THREADS - 1) / CUDA_MAX_THREADS;

  // Zero out the path sums in between each coarse channel because
  // we pick the top hits separately for each coarse channel
  cudaMemsetAsync(gpu_top_path_sums, 0, num_channels * sizeof(float));

  sumColumns<<<grid_size, CUDA_MAX_THREADS>>>(input.data, gpu_column_sums,
                                              rounded_num_timesteps, num_channels);
  checkCuda("sumColumns");
  
  int mid = num_channels / 2;

  // Do the Taylor tree algorithm
  for (int drift_block = min_drift_block; drift_block <= max_drift_block; ++drift_block) {
    
    // The dataflow among the buffers looks like:
    // input -> buffer1 -> buffer2 -> buffer1 -> buffer2 -> ...
    // We use the aliases source_buffer and target_buffer to make this simpler.
    // In each pass through the upcoming loop, we are reading from
    // source_buffer and writing to target_buffer.
    float* source_buffer = input.data;
    float* target_buffer = buffer1;

    // Each pass through the data calculates the sum of paths that are
    // twice as long as the previous path, until we reach our goal,
    // which is paths of length num_timesteps.
    for (int path_length = 2; path_length <= rounded_num_timesteps; path_length *= 2) {

      // Invoke cuda kernel
      if (use_cublas_taylor_tree) {
        runCublasTaylorTree(source_buffer, target_buffer, path_length, drift_block);
      } else {
        taylorTree<<<grid_size, CUDA_MAX_THREADS>>>(source_buffer, target_buffer,
                                                    rounded_num_timesteps, num_channels,
                                                    path_length, drift_block);
      }
      checkCuda("taylorTree");
      
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
    findTopPathSums<<<grid_size, CUDA_MAX_THREADS>>>(source_buffer, rounded_num_timesteps,
                                                     num_channels, drift_block,
                                                     gpu_top_path_sums,
                                                     gpu_top_drift_blocks,
                                                     gpu_top_path_offsets);
    checkCuda("findTopPathSums");
  }

  // Now that we have done all the GPU processing for one coarse
  // channel, we can copy the data back to host memory.
  // These copies are not async, so they will synchronize to the default stream.
  cudaMemcpy(cpu_column_sums, gpu_column_sums,
             num_channels * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(cpu_top_path_sums, gpu_top_path_sums,
             num_channels * sizeof(float), cudaMemcpyDeviceToHost);
  cudaMemcpy(cpu_top_drift_blocks, gpu_top_drift_blocks,
             num_channels * sizeof(int), cudaMemcpyDeviceToHost);
  cudaMemcpy(cpu_top_path_offsets, gpu_top_path_offsets,
             num_channels * sizeof(int), cudaMemcpyDeviceToHost);
  checkCuda("dedoppler d->h memcpy");
  
  // Use the central 90% of the column sums to calculate standard deviation.
  // We don't need to do a full sort; we can just calculate the 5th,
  // 50th, and 95th percentiles
  auto column_sums_end = cpu_column_sums + num_channels;
  std::nth_element(cpu_column_sums, cpu_column_sums + mid, column_sums_end);
  int first = ceil(0.05 * num_channels);
  int last = floor(0.95 * num_channels);
  std::nth_element(cpu_column_sums, cpu_column_sums + first,
                   cpu_column_sums + mid - 1);
  std::nth_element(cpu_column_sums + mid + 1, cpu_column_sums + last,
                   column_sums_end);
  float median = cpu_column_sums[mid];
    
  float sum = std::accumulate(cpu_column_sums + first, cpu_column_sums + last + 1, 0.0);
  float m = sum / (last + 1 - first);
  float accum = 0.0;
  std::for_each(cpu_column_sums + first, cpu_column_sums + last + 1,
                [&](const float f) {
                  accum += (f - m) * (f - m);
                });
  float std_dev = sqrt(accum / (last + 1 - first));

    
  // We consider two hits to be duplicates if the distance in their
  // frequency indexes is less than window_size. We only want to
  // output the largest representative of any set of duplicates.
  // window_size is chosen just large enough so that a single bright
  // pixel cannot cause multiple hits.
  // First we break up the data into a set of nonoverlapping
  // windows. Any candidate hit must be the largest within this
  // window.
  float path_sum_threshold = snr_threshold * std_dev + median;
  int window_size = 2 * ceil(normalized_max_drift * drift_timesteps);
  for (int i = 0; i * window_size < num_channels; ++i) {
    int candidate_freq = -1;
    float candidate_path_sum = path_sum_threshold;

    for (int j = 0; j < window_size; ++j) {
      int freq = i * window_size + j;
      if (freq >= num_channels) {
        break;
      }
      if (cpu_top_path_sums[freq] > candidate_path_sum) {
        // This is the new best candidate of the window
        candidate_freq = freq;
        candidate_path_sum = cpu_top_path_sums[freq];
      }
    }
    if (candidate_freq < 0) {
      continue;
    }

    // Check every frequency closer than window_size if we have a candidate
    int window_end = min(num_channels, candidate_freq + window_size);
    bool found_larger_path_sum = false;
    for (int freq = max(0, candidate_freq - window_size + 1); freq < window_end; ++freq) {
      if (cpu_top_path_sums[freq] > candidate_path_sum) {
        found_larger_path_sum = true;
        break;
      }
    }
    if (!found_larger_path_sum) {
      // The candidate frequency is the best within its window
      int drift_bins = cpu_top_drift_blocks[candidate_freq] * drift_timesteps +
        cpu_top_path_offsets[candidate_freq];
      double drift_rate = drift_bins * drift_rate_resolution;
      float snr = (candidate_path_sum - median) / std_dev;

      if (abs(drift_rate) >= min_drift) {
        DedopplerHit hit{candidate_freq, drift_bins, drift_rate, snr};
        if (print_hits) {
          cout << "hit: coarse channel = " << coarse_channel << ", "
               << hit.toString() << endl;
        }
        output->push_back(hit);
      }
    }
  }

  if (print_hit_summary && !output->empty()) {
    cout << "found " << pluralize(output->size(), "hit") << " in ";
    if (beam != NO_BEAM) {
      cout << "beam " << beam << ", ";
    }
    cout << "coarse channel " << coarse_channel << endl;
    int display_limit = 5;
    for (int i = 0; i < (int) output->size(); ++i) {
      const DedopplerHit& hit = (*output)[i];
      if (i < display_limit) {
        cout << "  index = " << hit.index << ", drift steps = "
             << hit.drift_steps << ", snr = " << hit.snr << ", drift rate = "
             << hit.drift_rate << endl;
      }
    }
    if ((int) output->size() > display_limit) {
      cout << "  (and " << ((int) output->size() - display_limit) << " more)\n";
    }
  }
}


