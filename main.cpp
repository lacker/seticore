#include <iostream>

#include <highfive/H5File.hpp>

using namespace std;

static_assert(sizeof(float) == 4, "require 32-bit floats");

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: seticore <h5file>" << endl;
    return 1;
  }

  // Open the file
  string filename = string(argv[1]);
  cout << "argument is: " << filename << endl;
  HighFive::File file(filename, HighFive::File::ReadOnly);

  // Find the dimensions of the data
  HighFive::DataSet dataset = file.getDataSet("data");

  string datatype = dataset.getDataType().string();
  if (datatype != "Float32") {
    cout << "this data is " << datatype << " but we can only handle Float32" << endl;
    return 1;
  }
  auto dimensions = dataset.getDimensions();
  for (int i = 0; i < dimensions.size(); ++i) {
    if (i == 0) {
      cout << "hdf5 data dimensions = ";
    } else {
      cout << " x ";
    }
    cout << dimensions[i];
  }
  cout << endl;

  if (dimensions.size() != 3) {
    cout << "expected three data dimensions" << endl;
    return 1;
  }

  if (dimensions[1] != 1) {
    cout << "unexpected second dimension: " << dimensions[1] << endl;
  }

  // Guess the coarse channel size
  int coarse_channel_size;
  if (dimensions[0] == 16 && dimensions[2] % 1048576 == 0) {
    // Looks like Green Bank data
    coarse_channel_size = 1048576;
  } else {
    cout << "unrecognized data dimensions" << endl;
    return 1;
  }
  int num_coarse_channels = dimensions[2] / coarse_channel_size;
  int num_timesteps = dimensions[0];
  
  // Loop through the coarse channels
  for (int coarse_channel = 0; coarse_channel < num_coarse_channels; ++coarse_channel) {
    vector<vector<vector<float> > > data;
    dataset.select({0, 0, coarse_channel * coarse_channel_size},
		   {num_timesteps, 1, coarse_channel_size}).read(data);

    cout << "sample " << coarse_channel << ": " << data[13][0][37] << endl;
  }
  return 0;
}
