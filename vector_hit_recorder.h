#pragma once

#include "dedoppler.h"
#include "hit_recorder.h"

#include <vector>

using namespace std;

/*
  A hit recorder that just accumulates hits in a vector in memory.
*/
class VectorHitRecorder: public HitRecorder {
 public:
  vector<DedopplerHit> hits;
  
  VectorHitRecorder() {}
  ~VectorHitRecorder() {}

  void recordHit(int coarse_channel, int freq_index, int drift_bins,
                 double drift_rate, double snr, const float* input);
};
