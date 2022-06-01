#pragma once

#include "hdf5.h"
#include <string>

using namespace std;


/*
  This class contains helper methods for processing .h5 files that store information
  you need to run beamforming, also known as "recipe files".
  Often suffixed as ".bfr5" files.
 */
class RecipeFile {
 private:
  hid_t file;

 public:
  RecipeFile(const string& filename);
  ~RecipeFile();
};
