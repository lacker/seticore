#pragma once

#include <vector>

#include "filterbank_buffer.h"

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


class Dedopplerer {
public:
  // How many timesteps of data are processed by the engine.
  // They may be stored in a larger array.
  const int num_timesteps;

  // How many frequency channels of data are processed by the engine.
  const int num_channels;

  // Frequency difference between adjacent bins, in MHz
  const double foff;

  // Time difference between adjacent bins, in seconds
  const double tsamp;

  // Whether the data we receive has a DC spike
  const bool has_dc_spike;
  
  Dedopplerer(int num_timesteps, int num_channels, double foff, double tsamp, bool has_dc_spike);
  ~Dedopplerer();

  void search(const FilterbankBuffer& input,
              double max_drift, double min_drift, double snr_threshold,
              vector<DedopplerHit>* output);
  
private:
  // For computing Taylor sums, we use three unified memory arrays,
  // each the size of one coarse channel. The input array is to read the source
  // data, and these two buffers are to use for Taylor tree calculations.
  float *buffer1, *buffer2;

  // The size of the larger array that contains the number of timesteps,
  // typically rounded up to the nearest power of two.
  int rounded_num_timesteps;
  
  // For normalization, we sum each column.
  float *column_sums;

  // For aggregating top hits, we use three arrays, each the size of
  // one row of the coarse channel. One array to store the largest
  // path sum, one to store which drift block it came from, and one to
  // store its path offset.
  float *top_path_sums;
  int *top_drift_blocks, *top_path_offsets;

  // How many timesteps the signal drifts in our data
  int drift_timesteps;

  // The difference in adjacent drift rates that we look for, in Hz/s
  double drift_rate_resolution;  
};
