#include "hdf5.h"
#include <iostream>

#include "recipe_file.h"

using namespace std;

RecipeFile::RecipeFile(const string& filename) {
  file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file == H5I_INVALID_HID) {
    cerr << "could not open file: " << filename << endl;
    exit(1);
  }
}

RecipeFile::~RecipeFile() {
}
