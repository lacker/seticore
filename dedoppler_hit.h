#pragma once

#include <string>

#include "filterbank_metadata.h"

using namespace std;

const int NO_BEAM = -1;

class DedopplerHit {
public:
  // The frequency the hit starts at
  double frequency;

  // Which frequency bin the hit starts at, within the coarse channel
  int index;

  // How many bins the hit drifts over.
  // Like (ending index - starting index), this is positive for rightward drift,
  // negative for leftward drift.
  // This is zero for a vertical line.
  // Includes drift over the full rounded-up power-of-two time range.
  int drift_steps;

  // The drift rate in Hz/s
  double drift_rate;

  // The signal-to-noise ratio for the hit
  float snr;

  // Which coarse channel the hit is in.
  int coarse_channel;

  // Which beam the hit is in. NO_BEAM if there is none, or for the incoherent beam.
  int beam;

  // This does *not* use rounded-up-to-a-power-of-two timesteps.
  int num_timesteps;

  // The total power used in the numerator to calculate snr.
  float power;
  
  DedopplerHit(const FilterbankMetadata& metadata, int _index, int _drift_steps,
               double _drift_rate, float _snr, int _beam, int _coarse_channel,
               int _num_timesteps, float _power);

  string toString() const;

  // Lowest index that contains a bin with this signal
  int lowIndex() const;

  // Highest index that contains a bin with this signal
  int highIndex() const;
};


bool operator<(const DedopplerHit& lhs, const DedopplerHit& rhs);
