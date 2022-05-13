#pragma once

#include "filterbank_file.h"
#include "hdf5.h"

using namespace std;

/*
  This class contains helper methods for processing .h5 files that are specific to the
  astronomical data format we expect, like extracting particular data from headers.
 */
class H5File: public FilterbankFile {
 private:
  double getDoubleAttr(const string& name) const;
  string getStringAttr(const string& name) const;
  long getLongAttr(const string& name) const;
  bool attrExists(const string& name) const;
  hid_t file, dataset, dataspace;
  int hit_count;
  
 public:
  H5File(const string& filename);
  ~H5File();

  void loadCoarseChannel(int i, float* output) const;
};
