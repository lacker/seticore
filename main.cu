#include <algorithm>
#include <assert.h>
#include <cuda.h>
#include <functional>
#include "hdf5.h"
#include <iostream>
#include <math.h>
#include <numeric>
#include <vector>

#include <fmt/core.h>

using namespace std;

static_assert(sizeof(float) == 4, "require 32-bit floats");


class H5File {
public:
  hid_t file, dataset, dataspace;
  double tsamp, foff;
  int num_timesteps, num_freqs, coarse_channel_size, num_coarse_channels;
  
  H5File(const string& filename) {
    file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file == H5I_INVALID_HID) {
      cerr << "could not open file: " << filename << endl;
      exit(1);
    }
    dataset = H5Dopen2(file, "data", H5P_DEFAULT);
    if (dataset == H5I_INVALID_HID) {
      cerr << "could not open dataset\n";
      exit(1);
    }
    if (!H5Tequal(H5Dget_type(dataset), H5T_NATIVE_FLOAT)) {
      cerr << "dataset is not float\n";
      exit(1);
    }
    
    this->getDoubleAttr("tsamp", &tsamp);
    this->getDoubleAttr("foff", &foff);

    dataspace = H5Dget_space(dataset);
    if (dataspace == H5I_INVALID_HID) {
      cerr << "could not open dataspace\n";
      exit(1);
    }
    if (H5Sget_simple_extent_ndims(dataspace) != 3) {
      cerr << "data is not three-dimensional\n";
      exit(1);
    }
    hsize_t dims[3];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);
    if (dims[1] != 1) {
      cerr << "unexpected second dimension: " << dims[1] << endl;
      exit(1);
    }
    num_timesteps = dims[0];
    num_freqs = dims[2];

    // Guess the coarse channel size
    if (num_timesteps == 16 && num_freqs % 1048576 == 0) {
      // Looks like Green Bank data
      coarse_channel_size = 1048576;
    } else {
      cerr << "unrecognized data dimensions: " << num_timesteps << " x " << num_freqs << endl;
      exit(1);
    }
    num_coarse_channels = num_freqs / coarse_channel_size;
  }

  void getDoubleAttr(const string& name, double* output) {
    auto attr = H5Aopen(dataset, name.c_str(), H5P_DEFAULT);
    if (attr == H5I_INVALID_HID) {
      cerr << "could not access attr " << name << endl;
      exit(1);
    }
    if (H5Aread(attr, H5T_NATIVE_DOUBLE, output) < 0) {
      cerr << "attr " << name << " could not be read as double\n";
      exit(1);
    }
    H5Aclose(attr);
  }

  int sizeOfCoarseChannel() {
    return num_timesteps * coarse_channel_size * sizeof(float);
  }

  // Loads the data in row-major order.
  // output should have size sizeOfCoarseChannel
  void loadCoarseChannel(int i, float* output) {
    // Select a hyperslab containing just the coarse channel we want
    const hsize_t offset[3] = {0, 0, unsigned(i * coarse_channel_size)};
    const hsize_t coarse_channel_dim[3] = {unsigned(num_timesteps), 1,
					   unsigned(coarse_channel_size)};
    if (H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,
			    offset, NULL, coarse_channel_dim, NULL) < 0) {
      cerr << "failed to select coarse channel hyperslab\n";
      exit(1);
    }

    // Define a memory dataspace
    hid_t memspace = H5Screate_simple(3, coarse_channel_dim, NULL);
    if (memspace == H5I_INVALID_HID) {
      cerr << "failed to create memspace\n";
      exit(1);
    }
    
    // Copy from dataspace to memspace
    if (H5Dread(dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, output) < 0) {
      cerr << "h5 read failed\n";
      exit(1);
    }
    
    H5Sclose(memspace);
  }
  
  ~H5File() {
    H5Sclose(dataspace);
    H5Dclose(dataset);
    H5Fclose(file);
  }
};

/*
Kernel to runs one round of the Taylor tree algorithm, calculating the sums
of paths of length `path_length`.
Each thread does the calculations for one frequency.

Apologies if this comment is long, but the Taylor tree algorithm is
fairly complicated for the number of lines of code it is, so it takes
a while to explain.

The key to understanding this algorithm is to understand the format of
the buffers, source_buffer and target_buffer.
source_buffer and target_buffer store the sum along an
approximately-linear path through the input data. The best way to
think of them is as three-dimensional arrays, indexed like

buffer[time_block][path_offset][start_frequency]

First, let's focus on the case where the drift block is zero.

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

TODO: explain drift blocks

This kernel is designed to run one thread per frequency in the
data. If it's trying to calculate the sum of a path that goes out of
the range, it just won't write to that value. Any in-range path will
get a value written to it. So you don't have to initialize the
buffers, as long as you recognize that
output[path_offset][start_frequency] is only valid when

path_offset + start_frequency < num_freqs

TODO: update that statement for drift blocks

This code is based on the original kernel by Franklin Antonio, available at
  https://github.com/UCBerkeleySETI/dedopplerperf/blob/main/CudaTaylor5demo.cu
*/
__global__ void taylorTree(const float* source_buffer, float* target_buffer,
			   int num_timesteps, int num_freqs, int path_length) {
  assert(path_length <= num_timesteps);
  int freq = blockIdx.x * blockDim.x + threadIdx.x;
  bool worker = (freq >= 0) && (freq < num_freqs);
  if (!worker) {
    return;
  }

  int num_time_blocks = num_timesteps / path_length;
  for (int time_block = 0; time_block < num_time_blocks; ++time_block) {
    int j = time_block * path_length;

    for (int path_offset = path_length - 1; path_offset >= 0; path_offset--) {

      // The recursion calculates sums for a target time block based on two
      // different source time blocks.
      // Data for block b comes from blocks 2b and 2b+1.
      // Remember there are twice as many source time blocks, so the
      // indices are different.

      // The recursion adds up two smaller paths, each with offset (path_offset/2).
      // When path_offset is odd, we need to shift the second path
      // over by one, so that the total adds back up to path_offset.
      // freq_shift thus represents the amount we need to shift the
      // second path.
      int half_offset = path_offset / 2;
      int freq_shift = (path_offset + 1) / 2;

      if (freq + freq_shift >= num_freqs) {
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
      int j2 = half_offset + path_length / 2;

      target_buffer[(j + path_offset) * num_freqs + freq] =
	source_buffer[(j + half_offset) * num_freqs + freq] +
	source_buffer[(j + j2) * num_freqs + freq + freq_shift];
    }
  }
}


int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: seticore <h5file>" << endl;
    return 1;
  }

  // Open the file
  string filename = string(argv[1]);
  cout << "argument is: " << filename << endl;
  H5File file(filename);

  double obs_length = file.num_timesteps * file.tsamp;
  double drift_rate_resolution = 1e6 * abs(file.foff) / obs_length;
  double max_drift = 0.4;
  double block_drift = drift_rate_resolution * file.num_timesteps;
  
  // We search through drift blocks in the range [-max_drift_block, max_drift_block].
  // Drift block i represents a search for slopes in the ranges
  // [i, i+1] and [-(i+1), -i]
  // when measured in horizontal pixels per vertical pixel.
  int max_drift_block = floor(max_drift / (drift_rate_resolution * file.num_timesteps));
  
  cout << fmt::format("block_drift: {:.8f}\n", block_drift);
  cout << fmt::format("max_drift_block: {}\n", max_drift_block);

  // We use three unified memory arrays, each the size of one coarse channel. One array to
  // read the source data, and two buffers to use for Taylor tree calculations.
  float *input, *buffer1, *buffer2;
  cudaMallocManaged(&input, file.sizeOfCoarseChannel());
  cudaMallocManaged(&buffer1, file.sizeOfCoarseChannel());
  cudaMallocManaged(&buffer2, file.sizeOfCoarseChannel());

  // Load and process one coarse channel at a time from the hdf5  
  for (int coarse_channel = 0; coarse_channel < file.num_coarse_channels; ++coarse_channel) {
    file.loadCoarseChannel(coarse_channel, input);

    cout << "loaded coarse channel " << coarse_channel << endl;

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
    cout << fmt::format("sample {} median: {:.3f}\n", coarse_channel, median);

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
    float stdev = sqrt(accum / (end - begin));
    cout << fmt::format("sample {} stdev: {:.3f}\n", coarse_channel, stdev);

    // TODO: loop over drift blocks
    
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
      int block_size = 1024;
      int grid_size = (file.coarse_channel_size + block_size - 1) / block_size;
      taylorTree<<<grid_size, block_size>>>(source_buffer, target_buffer,
					    file.num_timesteps, file.coarse_channel_size,
					    path_length);

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

    cudaDeviceSynchronize();

    // The final sums are in source_buffer because we did one last alias-swap
    if (coarse_channel == 0) {
      for (int k = 0; k < file.num_timesteps; ++k) {
	cout << fmt::format("k = {}; sums = {:.3f} {:.3f} {:.3f}\n", k,
			    source_buffer[k * file.coarse_channel_size + 13],
			    source_buffer[k * file.coarse_channel_size + 37],
			    source_buffer[k * file.coarse_channel_size + 123456]);
      }
    }

    if (coarse_channel > 5) {
      return 0;
    }
  }
  return 0;
}
