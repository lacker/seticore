#pragma once

#include "filterbank_buffer.h"
#include "filterbank_metadata.h"

#include <memory>

using namespace std;

/*
  The FilterbankFileReader is a base class for the different file formats that store
  filterbank data.
*/
class FilterbankFileReader: public FilterbankMetadata {
 public:
  FilterbankFileReader(const string& filename) : filename(filename) {}
  
  const string filename;
  
  virtual void loadCoarseChannel(int i, FilterbankBuffer* buffer) const;

  virtual ~FilterbankFileReader() {}  
};

unique_ptr<FilterbankFileReader> loadFilterbankFile(const string& filename);
