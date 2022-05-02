using namespace std;

#include <cstdlib>
#include "hdf5.h"
#include <iostream>

#include "h5_file.h"


/*
  Opens an h5 file for reading, assuming it contains standard radio
  telescope input from one of the known telescopes. If the data is an
  unexpected size or shape we should be conservative and exit.
 */
H5File::H5File(const string& filename) {
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

void H5File::getDoubleAttr(const string& name, double* output) const {
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

// Loads the data in row-major order.
void H5File::loadCoarseChannel(int i, float* output) const {
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
  
H5File::~H5File() {
  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Fclose(file);
}

