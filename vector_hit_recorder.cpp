#include "vector_hit_recorder.h"

void VectorHitRecorder::recordHit(DedopplerHit hit, int beam, int coarse_channel,
                                  const float* input) {
  hits.push_back(hit);
}
