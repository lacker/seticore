#pragma once

#include <iostream>
#include <fstream>

#include "filterbank_metadata.h"
#include "hit_recorder.h"

using namespace std;

/*
  This class creates a file containing all hit data, similar to the .dat file format,
  but with more useful data.
 */
class HitFileWriter: public HitRecorder {
 private:
  const FilterbankMetadata& metadata;
  int fd;

 public:
  bool verbose;
  
  HitFileWriter(const string& filename, const FilterbankMetadata& metadata);
  ~HitFileWriter();

  void recordHit(DedopplerHit hit, int beam, int coarse_channel, const float* input);
};
