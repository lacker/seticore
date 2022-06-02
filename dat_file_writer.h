#pragma once

#include <iostream>
#include <fstream>

#include "filterbank_file.h"
#include "hit_recorder.h"

using namespace std;

/*
  This class contains helper methods for the .dat file format, to maintain
  output consistency with turboseti.
*/
class DatFileWriter: public HitRecorder {
 private:
  const FilterbankFile& metadata;
  ofstream file;
  int hit_count;
  
 public:
  DatFileWriter(const string& filename, const FilterbankFile& metadata, double max_drift);
  ~DatFileWriter();

  void recordHit(int coarse_channel, int freq_index, int drift_bins, double drift_rate,
                 double snr, const float* input);
};
