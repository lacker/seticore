#include "vector_hit_recorder.h"

void VectorHitRecorder::recordHit(DedopplerHit hit, const float* input) {
  hits.push_back(hit);
}
