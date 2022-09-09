#pragma once

#include <memory>

#include "dedoppler_hit.h"
#include "filterbank_file_reader.h"

using namespace std;

class HitRecorder {
 public:
  virtual void recordHit(DedopplerHit hit, const float* input) = 0;

  virtual ~HitRecorder() {}
};

unique_ptr<HitRecorder> makeHitRecorder(const string& filename,
                                        const FilterbankFileReader& metadata,
                                        double max_drift);
