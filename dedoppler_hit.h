#pragma once

#include <string>

using namespace std;

const int NO_BEAM = -1;

struct DedopplerHit {
  // Which frequency bin the hit starts at, within the coarse channel
  int index;

  // How many bins the hit drifts over.
  // Like (ending index - starting index), this is positive for rightward drift,
  // negative for leftward drift.
  // This is zero for a vertical line.
  int drift_steps;

  // The drift rate in Hz/s
  double drift_rate;

  // The signal-to-noise ratio for the hit
  float snr;

  // Which beam the hit is in. NO_BEAM if there is none, or for the incoherent beam.
  int beam;

  // Which coarse channel the hit is in.
  int coarse_channel;

  string toString() const;

  // Lowest index that contains a bin with this signal
  int lowIndex() const;

  // Highest index that contains a bin with this signal
  int highIndex() const;
};



