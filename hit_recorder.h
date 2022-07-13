#pragma once

#include <memory>

#include "dedoppler_hit.h"
#include "filterbank_file.h"

using namespace std;

class HitRecorder {
 public:
  virtual void recordHit(DedopplerHit hit, int coarse_channel, const float* input) = 0;

  virtual ~HitRecorder() {}
};

unique_ptr<HitRecorder> makeHitRecorder(const string& filename,
                                        const FilterbankFile& metadata,
                                        double max_drift);
