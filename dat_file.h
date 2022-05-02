#pragma once

#include "h5_file.h"

using namespace std;

/*
  This class contains helper methods for the .dat file format, to maintain
  output consistency with turboseti.
*/
class DatFile {
 private:

  
 public:
  DatFile(const string& filename, const H5File& metadata);
  ~DatFile();
};
