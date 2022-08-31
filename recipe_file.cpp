#include <algorithm>
#include <assert.h>
#include <fmt/core.h>
#include "hdf5.h"
#include <iostream>
#include <math.h>
#include <vector>

#include "recipe_file.h"

using namespace std;

// The caller is responsible for closing it
hid_t openFile(const string& filename) {
  auto file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file == H5I_INVALID_HID) {
    cerr << "could not open file: " << filename << endl;
    exit(1);
  }
  return file;
}

long evenlyDivide(long a, long b) {
  assert(0 == a % b);
  return a / b;
}

RecipeFile::RecipeFile(const string& filename) :
  file(openFile(filename)),
  obsid(getStringData("/obsinfo/obsid")),
  src_names(getStringVectorData("/beaminfo/src_names")),
  delays(getVectorData<double>("/delayinfo/delays", H5T_IEEE_F64LE)),
  time_array(getVectorData<double>("/delayinfo/time_array", H5T_IEEE_F64LE)),
  ras(getVectorData<double>("/beaminfo/ras", H5T_IEEE_F64LE)),
  decs(getVectorData<double>("/beaminfo/decs", H5T_IEEE_F64LE)),
  npol(getLongScalarData("/diminfo/npol")),
  nbeams(getLongScalarData("/diminfo/nbeams")),
  cal_all(getComplexVectorData("/calinfo/cal_all")),
  nants(evenlyDivide(delays.size(), nbeams * time_array.size())),
  nchans(evenlyDivide(cal_all.size(), nants * npol)) {
  if ((int) ras.size() < nbeams) {
    cerr << fmt::format("could only load {} ras but nbeams = {}\n", ras.size(), nbeams);
    exit(1);
  }
  if ((int) decs.size() < nbeams) {
    cerr << fmt::format("could only load {} decs but nbeams = {}\n", decs.size(), nbeams);
    exit(1);
  }
}

string makeRecipeFilename(const string& directory, const string& obsid) {
  string fixed_obsid = obsid;
  // Someone changed their mind and all the colons should actually
  // be hyphens in obsid's
  replace(fixed_obsid.begin(), fixed_obsid.end(), ':', '-');
  return fmt::format("{}/{}.bfr5", directory, fixed_obsid);
}

RecipeFile::RecipeFile(const string& directory, const string& obsid)
  : RecipeFile(makeRecipeFilename(directory, obsid)) {}

RecipeFile::~RecipeFile() {
  H5Fclose(file);
}

/*
  This is supposed to work for either UTF8 or ASCII.
 */
string RecipeFile::getStringData(const string& name) const {
  auto dataset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
  auto data_type = H5Dget_type(dataset);
  auto native_type = H5Tget_native_type(data_type, H5T_DIR_DEFAULT);
  int size = H5Tget_size(native_type);
  
  vector<char> buffer(size);
  if (H5Dread(dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &buffer[0]) < 0) {
    cerr << "dataset " << name << " could not be read as string\n";
    exit(1);
  }

  H5Tclose(native_type);
  H5Tclose(data_type);
  H5Dclose(dataset);
  
  string answer(buffer.begin(), buffer.end());
  return answer;
}

/*
  Reads a list of strings into the provided output vector.
  We expect them to have variable length.
  The data is converted to the native character type when we read it.
*/
vector<string> RecipeFile::getStringVectorData(const string& name) const {
  auto dataset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
  auto dataspace = H5Dget_space(dataset);
  int npoints = H5Sget_simple_extent_npoints(dataspace);
  auto native_type = H5Tvlen_create(H5T_NATIVE_CHAR);

  // Intermediately store data as the HDF5-provided "variable length sequence" object
  vector<hvl_t> sequences(npoints);

  if (H5Dread(dataset, native_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &sequences[0]) < 0) {
    cerr << "dataset " << name << " could not be read as string vector\n";
    exit(1);
  }

  vector<string> output;
  for (hvl_t sequence : sequences) {
    char* p = (char*)sequence.p;
    output.push_back(string(p, p + sequence.len));
  }

  H5Tclose(native_type);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  return output;
}

/*
  Reads a variable-length array into a vector.
  You need to provide two types to call this - T is the type that C++ uses, hdf5_type
  is the type that the HDF5 library is using.
  This doesn't work for strings; call getStringVectorData for that.
*/
template <class T>
vector<T> RecipeFile::getVectorData(const string& name, hid_t hdf5_type) const {
  auto dataset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
  auto dataspace = H5Dget_space(dataset);
  int npoints = H5Sget_simple_extent_npoints(dataspace);
  vector<T> output(npoints);
  if (H5Dread(dataset, hdf5_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &output[0]) < 0) {
    cerr << "dataset " << name << " could not be read with getVectorData\n";
    exit(1);
  }
  
  H5Sclose(dataspace);
  H5Dclose(dataset);

  return output;
}

vector<thrust::complex<float> > RecipeFile::getComplexVectorData(const string& name) const {
  hid_t complex_type = H5Tcreate(H5T_COMPOUND, 8);
  H5Tinsert(complex_type, "r", 0, H5T_IEEE_F32LE);
  H5Tinsert(complex_type, "i", 4, H5T_IEEE_F32LE);

  vector<thrust::complex<float> > output =
    getVectorData<thrust::complex<float> >(name, complex_type);
 
  H5Tclose(complex_type);

  return output;
}

// Reads a single long integer.
long RecipeFile::getLongScalarData(const string& name) const {
  auto dataset = H5Dopen(file, name.c_str(), H5P_DEFAULT);
  auto dataspace = H5Dget_space(dataset);
  int npoints = H5Sget_simple_extent_npoints(dataspace);
  if (npoints != 1) {
    cerr << "expected scalar at " << name << " but got " << npoints << " data values\n";
    exit(1);
  }

  long output;
  if (H5Dread(dataset, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &output) < 0) {
    cerr << "dataset " << name << " could not be read with getLongScalarData\n";
    exit(1);
  }
  
  H5Sclose(dataspace);
  H5Dclose(dataset);

  return output;  
}

/*
  Figures out what time array index to use for beamforming, for the given time.
  This is the index for which time_array[index] is closest to the provided time.
 */
int RecipeFile::getTimeArrayIndex(double time) const {  
  // The two candidates are the ones to either side of time
  auto right_iter = lower_bound(time_array.begin(), time_array.end(), time);
  if (time_array.begin() == right_iter) {
    return 0;
  }
  auto left_iter = right_iter - 1;

  // Check which one is closer
  double right_dist = *right_iter - time;
  assert(right_dist >= 0);
  double left_dist = time - *left_iter;
  assert(left_dist >= 0);

  if (right_dist < left_dist) {
    return right_iter - time_array.begin();
  }
  return left_iter - time_array.begin();
}

/*
  Accessor to delays.
  Delays are row-major organized by:
    delays[time_array_index][beam][antenna]
*/
double RecipeFile::getDelay(int time_array_index, int beam, int antenna) const {
  return delays[((time_array_index * nbeams) + beam) * nants + antenna];
}

/*
  Accessor to cal_all.
  cal_all is row-major organized by:
    cal_all[frequency][polarity][antenna]
  TODO: figure out what this should actually be named.
 */
thrust::complex<float> RecipeFile::getCal(int frequency, int polarity, int antenna) const {
  return cal_all[((frequency * npol) + polarity) * nants + antenna];
}

/*
  Generate the beamforming coefficients for the given parameters.

  There are two ways parameters give us different coefficients. They specify the time, and
  the frequency range we are interested in.

  time_array_index corresponds to which entry in time_array to use.

  For the frequency range we need more data to describe the subband we are beamforming for.
  A subband corresponds to a subset of the channels that the recipe file has data for.
  start_channel and num_channels define what subset of the channels in the recipe to use.
  start_channel is the first one of the range, num_channels is the size of the range.
  center_frequency is the center of the subband in MHz. (the OBSFREQ header)
  bandwidth is the width of the subband in MHz, negative for reversed. (the OBSBW header)
 */
void RecipeFile::generateCoefficients(int time_array_index,
                                      int start_channel, int num_channels,
                                      float center_frequency, float bandwidth,
                                      thrust::complex<float>* coefficients) const {
  assert(start_channel + num_channels <= nchans);
  
  int output_index = 0;
  for (int freq = 0; freq < num_channels; ++freq) {

    // Calculate the center of this coarse channel in GHz
    bool use_buggy_logic = true; // for compatibility with hpguppi_proc
    float chan_bandwidth = bandwidth / num_channels * 0.001;
    float center_index = (num_channels - 1.0) / 2.0;
    if (use_buggy_logic) {
      center_index -= 0.5;
    }
    float chan_center = center_frequency * 0.001 + (freq - center_index) * chan_bandwidth;

    for (int beam = 0; beam < nbeams; ++beam) {
      for (int polarity = 0; polarity < npol; ++polarity) {
        for (int antenna = 0; antenna < nants; ++antenna) {
          float tau = getDelay(time_array_index, beam, antenna);
          int global_channel_index = freq + start_channel;
          auto cal = getCal(global_channel_index, polarity, antenna);

          // Figure out how much to rotate
          float angle = 2 * M_PI * chan_center * tau;
          float cos_val = cos(angle);
          float sin_val = sin(angle);

          // Rotate to get the coefficients
          float real = cal.real() * cos_val - cal.imag() * sin_val;
          float imag = cal.real() * sin_val + cal.imag() * cos_val;
          coefficients[output_index++] = thrust::complex<float>(real, imag);
        }
      }
    }
  }
}
