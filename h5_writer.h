#pragma once

#include "hdf5.h"

#include "filterbank_metadata.h"

using namespace std;


class H5Writer{
 public:
  const string filename;
  const FilterbankMetadata metadata;
  bool closed;
  
  H5Writer(const string& filename,
           const FilterbankMetadata& metadata);
  ~H5Writer();

  // data must be formatted as row-major:
  //   data[time][freq]
  void setData(const float* data);
  
  void close();
  
 private:
  hid_t file, dataset, dataspace;

  void setAttr(const string& name, hid_t type, const void* value);
  void setDoubleAttr(const string& name, double value);
  void setStringAttr(const string& name, const string& value);
  void setLongAttr(const string& name, long value);
};
