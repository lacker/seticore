#pragma once

#include <iostream>
#include <fstream>

#include "filterbank_file.h"
#include "hit_recorder.h"

using namespace std;

/*
  This class creates a file containing all hit data, similar to the .dat file format,
  but with more useful data.
 */
class HitFileWriter: public HitRecorder {
 private:
  const FilterbankFile& metadata;
  int fd;

 public:
  HitFileWriter(const string& filename, const FilterbankFile& metadata);
  ~HitFileWriter();

  void recordHit(DedopplerHit hit, int beam, int coarse_channel, const float* input);
};
