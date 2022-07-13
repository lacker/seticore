#pragma once

#include <string>

using namespace std;

struct DedopplerHit {
  // Which frequency bin the hit starts at
  int index;

  // How many bins the hit drifts over.
  // Like (ending index - starting index), this is positive for rightward drift,
  // negative for leftward drift.
  int drift_steps;

  // The drift rate in Hz/s
  double drift_rate;

  // The signal-to-noise ratio for the hit
  float snr;

  string toString() const;
};



