#pragma once

#include <memory>

#include "filterbank_file.h"

using namespace std;

class HitRecorder {
 public:
  virtual void recordHit(int coarse_channel, int freq_index, int drift_bins,
                         double drift_rate, double snr) = 0;
};

unique_ptr<HitRecorder> makeHitRecorder(const string& filename,
                                        const FilterbankFile& metadata, double max_drift);
