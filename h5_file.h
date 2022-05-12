#pragma once

#include "hdf5.h"

using namespace std;

/*
  This class contains helper methods for processing .h5 files that are specific to the
  astronomical data format we expect, like extracting particular data from headers.
 */
class H5File {
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

  const string filename;
  bool has_dc_spike;
  string source_name;
  double fch1, foff, tstart, tsamp, src_dej, src_raj;
  int num_timesteps, num_freqs, coarse_channel_size, num_coarse_channels;
  
  void loadCoarseChannel(int i, float* output) const;
};
