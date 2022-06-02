#pragma once

using namespace std;

class HitRecorder {
 public:
  virtual void recordHit(int coarse_channel, int freq_index, int drift_bins,
                         double drift_rate, double snr) = 0;
};
