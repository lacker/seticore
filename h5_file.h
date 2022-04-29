using namespace std;

class H5File {
 public:
  H5File(const string& filename);
  ~H5File();
  
  hid_t file, dataset, dataspace;
  double tsamp, foff;
  int num_timesteps, num_freqs, coarse_channel_size, num_coarse_channels;
  
  void getDoubleAttr(const string &name, double* output);
  void loadCoarseChannel(int i, float* output);
};
