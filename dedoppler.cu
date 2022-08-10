#include <algorithm>
#include <assert.h>
#include <cuda.h>
#include <functional>
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>

#include "cuda_util.h"
#include "dedoppler.h"
#include "taylor.h"
#include "util.h"

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
      has_dc_spike(has_dc_spike), print_hits(false), print_hit_summary(false) {
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
      taylorTree<<<grid_size, CUDA_MAX_THREADS>>>(source_buffer, target_buffer,
                                                  rounded_num_timesteps, num_channels,
                                                  path_length, drift_block);
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

