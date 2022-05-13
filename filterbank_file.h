#pragma once

using namespace std;

/*
  The FilterbankFile is a base class for the different file formats that store filterbank data.
*/
class FilterbankFile {
 public:
  FilterbankFile(const string& filename) : filename(filename) {}
  
  const string filename;
  bool has_dc_spike;
  string source_name;
  double fch1, foff, tstart, tsamp, src_dej, src_raj;
  int num_timesteps, num_freqs, coarse_channel_size, num_coarse_channels;
  
  virtual void loadCoarseChannel(int i, float* output) const = 0;
};

FilterbankFile* loadFilterbank(const string& filename);
