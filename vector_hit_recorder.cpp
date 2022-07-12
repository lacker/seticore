#include "vector_hit_recorder.h"

void VectorHitRecorder::recordHit(int coarse_channel, int freq_index, int drift_bins,
                                  double drift_rate, double snr,
                                  const float* input) {
  hits.push_back(DedopplerHit{freq_index, drift_bins, drift_rate, (float) snr});
}
