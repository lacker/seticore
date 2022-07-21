#pragma once

#include <string>
#include <vector>

#include "hit_recorder.h"
#include "util.h"

class BeamformingConfig {
public:
  // Required parameters for the run
  const vector<string> raw_files;
  const string output_dir;
  const string recipe_dir;
  const int num_bands;
  const int fft_size;
  const int telescope_id;
  const float snr;
  const float max_drift;
  const float min_drift;

  // Optional additional ways to modify the config

  // If set, we only process a subset of the bands
  int num_bands_to_process;

  // If set, use this to record hits instead of writing a file
  unique_ptr<HitRecorder> hit_recorder;

  // If set, save the beamformed filterbanks as h5 files
  string h5_dir;
  
  BeamformingConfig(const vector<string>& raw_files,
                    const string& output_dir,
                    const string& recipe_dir,
                    int num_bands,
                    int fft_size,
                    int telescope_id,
                    float snr,
                    float max_drift,
                    float min_drift)
    : raw_files(raw_files), output_dir(stripAnyTrailingSlash(output_dir)),
      recipe_dir(recipe_dir), num_bands(num_bands), fft_size(fft_size),
      telescope_id(telescope_id), snr(snr), max_drift(max_drift),
      min_drift(min_drift), num_bands_to_process(num_bands) {
  }

  void run();
};
