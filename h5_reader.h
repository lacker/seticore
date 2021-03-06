#pragma once

#include "filterbank_file_reader.h"
#include "hdf5.h"

using namespace std;

/*
  This class contains helper methods for processing .h5 files that are specific to the
  astronomical data format we use for filterbanks, like extracting
  particular data from headers.
  Also known as the "FBH5" format.
 */
class H5Reader: public FilterbankFileReader {
 private:
  double getDoubleAttr(const string& name) const;
  string getStringAttr(const string& name) const;
  long getLongAttr(const string& name) const;
  bool attrExists(const string& name) const;
  hid_t file, dataset, dataspace;
  
 public:
  H5Reader(const string& filename);
  ~H5Reader();

  void loadCoarseChannel(int i, FilterbankBuffer* buffer) const;
};
