#pragma once

#include <string>
#include <vector>

#include "hit_recorder.h"
#include "raw_file_group.h"
#include "util.h"

class BeamformingPipeline {
public:
  // Required parameters for the run
  const vector<string> raw_files;
  const string output_dir;
  const string recipe_filename;
  const int num_bands;
  const int fft_size;
  const int sti;
  const int telescope_id;
  const float snr;
  const float max_drift;

  // Optional additional ways to modify the config

  // If set, we only process a subset of the bands
  int num_bands_to_process;

  // Whether to save hits to a file
  bool record_hits;

  // If set, save the beamformed filterbanks as h5 files
  string h5_dir;

  // recipe_filename can either be a file ending in .bfr5 or a directory
  BeamformingPipeline(const vector<string>& raw_files,
                      const string& output_dir,
                      const string& recipe_filename,
                      int num_bands,
                      int fft_size,
                      int sti,
                      int telescope_id,
                      float snr,
                      float max_drift)
    : raw_files(raw_files), output_dir(stripAnyTrailingSlash(output_dir)),
      recipe_filename(recipe_filename), num_bands(num_bands), fft_size(fft_size),
      sti(sti), telescope_id(telescope_id), snr(snr), max_drift(max_drift),
      num_bands_to_process(num_bands), record_hits(true),
      file_group(raw_files) {
  }

  void findHits();  
  vector<DedopplerHit> hits;

  void makeStamps();
  
private:
  RawFileGroup file_group;
  
};
