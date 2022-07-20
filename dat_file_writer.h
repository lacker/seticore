#pragma once

#include <iostream>
#include <fstream>

#include "filterbank_file_reader.h"
#include "hit_recorder.h"

using namespace std;

/*
  This class contains helper methods for the .dat file format, to maintain
  output consistency with turboseti.
*/
class DatFileWriter: public HitRecorder {
 private:
  const FilterbankFileReader& metadata;
  ofstream file;
  int hit_count;
  
 public:
  DatFileWriter(const string& filename, const FilterbankFileReader& metadata,
                double max_drift);
  ~DatFileWriter();

  void recordHit(DedopplerHit hit, int beam, int coarse_channel, const float* input);
};
