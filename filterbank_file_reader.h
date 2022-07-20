#pragma once

#include "filterbank_buffer.h"

#include <memory>

using namespace std;

/*
  The FilterbankFileReader is a base class for the different file formats that store
  filterbank data.
*/
class FilterbankFileReader {
 public:
  FilterbankFileReader(const string& filename) : filename(filename),
                                                 coarse_channel_size(0) {}
  
  const string filename;
  bool has_dc_spike;
  string source_name;
  double fch1, foff, tstart, tsamp, src_dej, src_raj;
  int num_timesteps, num_freqs, coarse_channel_size, num_coarse_channels, telescope_id;
  
  virtual void loadCoarseChannel(int i, FilterbankBuffer* buffer) const;

  virtual ~FilterbankFileReader() {}
  
 protected:
  void inferMetadata();
};

unique_ptr<FilterbankFileReader> loadFilterbankFile(const string& filename);
