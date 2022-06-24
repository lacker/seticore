#include <assert.h>
#include "hdf5.h"
#include <iostream>
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
  nchans(evenlyDivide(cal_all.size(), nants * npol))
{ }

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
