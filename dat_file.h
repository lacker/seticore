#pragma once

#include <iostream>
#include <fstream>

#include "h5_file.h"

using namespace std;

/*
  This class contains helper methods for the .dat file format, to maintain
  output consistency with turboseti.
*/
class DatFile {
 private:
  const H5File& metadata;
  ofstream file;
  int hit_count;
  
 public:
  DatFile(const string& filename, const H5File& metadata, double max_drift);
  ~DatFile();

  void reportHit(int coarse_channel, int freq_index, double drift_rate, double snr);
};
