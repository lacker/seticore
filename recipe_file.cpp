#include "hdf5.h"
#include <iostream>
#include <vector>

#include "recipe_file.h"

using namespace std;

RecipeFile::RecipeFile(const string& filename) {
  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file == H5I_INVALID_HID) {
    cerr << "could not open file: " << filename << endl;
    exit(1);
  }

  cout << "obs id is: " << getStringData("/obsinfo/obsid") << endl;
}

RecipeFile::~RecipeFile() {
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
  This is supposed to work for either UTF8 of ASCII.
  Reads a list of strings into the provided output vector.
*/
void RecipeFile::getStringVectorData(const string& name, vector<string>* output) const {
  auto dataset = H5Dopen(file, name.c_str(), H5P_DEFAULT);

  // TODO: do stuff in here
  
  H5Dclose(dataset);
}
