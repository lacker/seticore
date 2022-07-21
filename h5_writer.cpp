#include "hdf5.h"
#include <iostream>

#include "h5_writer.h"

using namespace std;

H5Writer::H5Writer(const string& filename, const FilterbankMetadata& metadata)
  : filename(filename), metadata(metadata) {
  // Deletes any already-existing file there
  file = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (file == H5I_INVALID_HID) {
    cerr << "could not open file for writing: " << filename << endl;
    exit(1);
  }

  hsize_t dims[3];
  dims[0] = metadata.num_timesteps;
  dims[1] = 1;
  dims[2] = metadata.num_channels;
  dataspace = H5Screate_simple(3, dims, NULL);
  if (dataspace == H5I_INVALID_HID) {
    cerr << "could not create dataspace\n";
    exit(1);
  }

  // For now, don't set chunking or use bitshuffle.
  dataset = H5Dcreate2(file, "data", H5T_NATIVE_FLOAT, dataspace,
                       H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (dataset == H5I_INVALID_HID) {
    cerr << "could not create dataset\n";
    exit(1);
  }

  setDoubleAttr("fch1", metadata.fch1);
  setDoubleAttr("foff", metadata.foff);
  setDoubleAttr("tstart", metadata.tstart);
  setDoubleAttr("tsamp", metadata.tsamp);
  setDoubleAttr("src_dej", metadata.src_dej);
  setDoubleAttr("src_raj", metadata.src_raj);
  setStringAttr("source_name", metadata.source_name);
  setLongAttr("telescope_id", metadata.telescope_id);
  setLongAttr("nfpc", metadata.coarse_channel_size);
}

H5Writer::~H5Writer() {
  close();
}

void H5Writer::setData(const float* data) {
  auto status = H5Dwrite(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  if (status < 0) {
    cerr << "hdf5 data write failed\n";
    exit(1);
  }
}

void H5Writer::close() {
  if (closed) {
    return;
  }

  H5Sclose(dataspace);
  H5Dclose(dataset);
  H5Fclose(file);
  
  closed = true;
}

void H5Writer::setAttr(const string& name, hid_t type, const void* value) {
  auto attr = H5Aopen(dataset, name.c_str(), H5P_DEFAULT);
  if (attr == H5I_INVALID_HID) {
    cerr << "could not access attr " << name << endl;
    exit(1);
  }
  if (H5Awrite(attr, type, value) < 0) {
    cerr << "could not write to attr " << name << endl;
    exit(1);
  }
  H5Aclose(attr);
}

void H5Writer::setDoubleAttr(const string& name, double value) {
  setAttr(name, H5T_NATIVE_DOUBLE, &value);
}

// Stores the string as fixed-length c-string, to match what rawspec does.
void H5Writer::setStringAttr(const string& name, const string& value) {
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, value.size());
  H5Tset_strpad(type, H5T_STR_NULLTERM);

  setAttr(name, type, value.c_str());
  
  H5Tclose(type);
}

void H5Writer::setLongAttr(const string& name, long value) {
  setAttr(name, H5T_NATIVE_LONG, &value);
}


