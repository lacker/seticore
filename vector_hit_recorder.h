#pragma once

#include "dedoppler.h"
#include "dedoppler_hit.h"
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

  void recordHit(DedopplerHit hit, int beam, int coarse_channel, const float* input);

};
