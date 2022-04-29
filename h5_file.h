#pragma once
using namespace std;

/*
  This class contains helper methods for processing .h5 files that are specific to the
  astronomical data format we expect, like extracting particular data from headers.
 */
class H5File {
 private:
  void getDoubleAttr(const string &name, double* output);
  hid_t file, dataset, dataspace;

 public:
  H5File(const string& filename);
  ~H5File();
  
  double tsamp, foff;
  int num_timesteps, num_freqs, coarse_channel_size, num_coarse_channels;
  
  void loadCoarseChannel(int i, float* output);
};
