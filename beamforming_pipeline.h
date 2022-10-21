#pragma once

#include <fmt/core.h>
#include <string>
#include <vector>

#include "hit_recorder.h"
#include "raw_file_group.h"
#include "util.h"

inline int calculateFFTSize(int num_coarse, int num_fine) {
  if (num_fine <= 0) {
    fatal(fmt::format("invalid number of fine channels: {}", num_fine));
  }

  if (num_fine % num_coarse != 0) {
    fatal(fmt::format("cannot upchannelize from {} channels to {} channels",
                      num_coarse, num_fine));
  }
  
  int fft_size = num_fine / num_coarse;
  if (!isPowerOfTwo(fft_size)) {
    fatal(fmt::format("upchannelizing from {} channels to {} channels gives an "
                      "implicit fft size of {} which we can't handle",
                      num_coarse, num_fine, fft_size));
  }

  return fft_size;
}

class BeamformingPipeline {
public:
  // Required parameters for the run
  const vector<string> raw_files;
  const string output_dir;
  const string recipe_filename;
  const int num_bands;
  const int sti;
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
  // If _fft_size is -1 we calculate from num_fine_channels
  BeamformingPipeline(const vector<string>& raw_files,
                      const string& output_dir,
                      const string& recipe_filename,
                      int num_bands,
                      int sti,
                      int _telescope_id,
                      float snr,
                      float max_drift,
                      int _fft_size,
                      int num_fine_channels)
    : raw_files(raw_files), output_dir(stripAnyTrailingSlash(output_dir)),
      recipe_filename(recipe_filename), num_bands(num_bands), sti(sti), snr(snr),
      max_drift(max_drift), num_bands_to_process(num_bands), record_hits(true),
      file_group(raw_files),
      telescope_id(_telescope_id == NO_TELESCOPE_ID
                   ? file_group.getTelescopeID() : _telescope_id),
      fft_size(_fft_size > 0 ? _fft_size
               : calculateFFTSize(file_group.num_coarse_channels, num_fine_channels)) {
  }

  void findHits();  
  vector<DedopplerHit> hits;

  void makeStamps();
  
private:
  RawFileGroup file_group;

public:
  const int telescope_id;
  const int fft_size;
  
};
