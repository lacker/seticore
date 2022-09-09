#include "beamforming_config.h"
#include "raw_file_group.h"
#include "vector_hit_recorder.h"

#include <iostream>
#include <string>

using namespace std;

/*
  This integration test involves the hardcoded location of about 40G of testing
  files. That makes it implausible to use in any more than one location.
  If you do want to run this, though, you need to create a directory ../data/beamtest
  (it can be a symlink) and download these files into that directory:

  https://bldata.berkeley.edu/pipeline/tmp/MeerKAT-array_1-20220513T043147Z.bfr5
  https://bldata.berkeley.edu/pipeline/tmp/guppi_59712_16307_003760_J1939-6342_0001.0000.raw
  https://bldata.berkeley.edu/pipeline/tmp/guppi_59712_16307_003760_J1939-6342_0001.0001.raw

  TODO: make this based on smaller files which are generally available
 */
int main(int argc, char* argv[]) {
  // Specifying parameters
  string input_dir = "../data/beamtest";
  string output_dir = "mocked";
  string recipe_dir = "../data/beamtest";
  int num_bands = 32;
  int fft_size = 131072;
  int sti = 8;
  int telescope_id = 64;
  float snr = 7.0;
  float max_drift = 0.01;
  float min_drift = 0.01;

  auto groups = scanForRawFileGroups(input_dir);
  assert(groups.size() == 1);
  assert(groups[0].size() == 2);

  BeamformingConfig config(groups[0], output_dir, recipe_dir, num_bands,
                           fft_size, sti, telescope_id, snr, max_drift, min_drift);
  config.num_bands_to_process = 1;
  config.record_hits = false;
  
  config.run();

  assert(114049 == config.hits[0].index);
  assert(-1 == config.hits[0].drift_steps);
  assertFloatEq(7.03593, config.hits[0].snr);
  assertFloatEq(-0.317774, config.hits[0].drift_rate);
  assert(67 == config.hits.size());

  assertStringEq(config.hits[10].toString(),
                 "coarse channel = 0, index = 106914, snr = 10.71197, "
                 "drift rate = -0.31777 (-1 bin)");
  assertStringEq(config.hits[20].toString(),
                 "coarse channel = 1, index = 17237, snr = 7.07149, "
                 "drift rate = -0.31777 (-1 bin)");
  
  return 0;
}
