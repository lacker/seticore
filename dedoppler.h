#pragma once

#include <vector>

using namespace std;


struct DedopplerHit {
  // Which frequency bin the hit starts at
  int candidate_freq;

  // How many bins the hit drifts over
  int drift_bins;

  // The drift rate in Hz/s
  double drift_rate;

  // The signal-to-noise ratio for the hit
  float snr;
};


class Dedopplerer {
public:
  // How many timesteps of data are processed by the engine.
  // They may be stored in a larger array.
  int num_timesteps;

  // How many frequency channels of data are processed by the engine.
  int num_channels;

  // The size of the larger array that contains the number of timesteps,
  // typically rounded up to the nearest power of two.
  int rounded_num_timesteps;

  // For computing Taylor sums, we use three unified memory arrays,
  // each the size of one coarse channel. One array to read the source
  // data, and two buffers to use for Taylor tree calculations.
  float *input, *buffer1, *buffer2;

  // For normalization, we sum each column.
  float *column_sums;

  // For aggregating top hits, we use three arrays, each the size of
  // one row of the coarse channel. One array to store the largest
  // path sum, one to store which drift block it came from, and one to
  // store its path offset.
  float *top_path_sums;
  int *top_drift_blocks, *top_path_offsets;

  // Frequency difference between adjacent bins, in MHz
  double foff;

  // Time difference between adjacent bins, in seconds
  double tsamp;
  
  // How many timesteps the signal drifts in our data
  int drift_timesteps;

  // The difference in adjacent drift rates that we look for, in Hz/s
  double drift_rate_resolution;

  // Whether the data we receive has a DC spike
  bool has_dc_spike;

  Dedopplerer(int num_timesteps, int num_channels, double foff, double tsamp, bool has_dc_spike);
  ~Dedopplerer();

  void processInput(double max_drift, double min_drift, double snr_threshold,
                    vector<DedopplerHit>* output);

  
};
