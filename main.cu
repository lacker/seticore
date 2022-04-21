#include <algorithm>
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
Runs one round of the Taylor tree algorithm, calculating the sums
of paths of length `path_length`. given 

source_buffer contains the sums of paths of length `path_length / 2`
in row-major order.
When path_length is 2, these sums are just single elements, so we
can use the plain data for our base case.

target_buffer will contain the output data in row-major order.
This does not overwrite the elements of the target buffer that do
not correspond to a valid sum, so they may contain garbage data from previous
runs. The caller is responsible for not accessing this garbage
data.

TODO: describe the precise meaning of the elements of source_buffer
and target_buffer, as well as which elements are garbage

num_timesteps is the number of timesteps
num_freqs is the number of frequency bins
So each buffer holds (num_timesteps * num_freqs) elements.
 
Based on the Taylor kernel by Franklin Antonio, available at
  https://github.com/UCBerkeleySETI/dedopplerperf/blob/main/CudaTaylor5demo.cu
which also contains kernel performance testing information.
*/
__global__ void taylorTree(const float* source_buffer, float* target_buffer,
			   int num_timesteps, int num_freqs, int path_length) {

  int k = blockIdx.x * blockDim.x + threadIdx.x;
  bool worker = (k >= 0) && (k < num_freqs) && path_length <= num_timesteps;
  if (!worker) {
    return;
  }
  for (int j = 0; j < num_timesteps; j += path_length) {
    for (int j0 = path_length - 1; j0 >= 0; j0--) {
      int j1 = j0 / 2;
      int j2 = j1 + path_length / 2;
      int j3 = (j0 + 1) / 2;
      if (k + j3 < num_freqs) {
        target_buffer[(j + j0) * num_freqs + k] =
	  source_buffer[(j + j1) * num_freqs + k] +
	  source_buffer[(j + j2) * num_freqs + k + j3];
      }
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
    // cout << fmt::format("sample {} median: {:.3f}\n", coarse_channel, median);

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

    // Cuda parameters
    int block_size = 1024;
    int grid_size = (file.coarse_channel_size + block_size - 1) / block_size;

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
      // XXX launch the kernel

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

    // The final sums are in source_buffer because we did one last alias-swap
  }
  return 0;
}
