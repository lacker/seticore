#include "beamforming_config.h"
#include "raw_file_group.h"
#include "vector_hit_recorder.h"

#include <iostream>
#include <string>

using namespace std;

/*
  This integration test involves the hardcoded location of about 40G of testing
  files. That makes it implausible to use in any more than one location.
  TODO: make this based on smaller files which are generally available
 */
int main(int argc, char* argv[]) {
  // Specifying parameters
  string input_dir = "/d/onlytwo";
  string output_dir = "mocked";
  string recipe_dir = "/d/mubf";
  int num_bands = 32;
  int fft_size = 131072;
  int telescope_id = 64;
  float snr = 7.0;
  float max_drift = 0.01;
  float min_drift = 0.01;

  auto groups = scanForRawFileGroups(input_dir);
  assert(groups.size() == 1);
  assert(groups[0].size() == 2);

  BeamformingConfig config(groups[0], output_dir, recipe_dir, num_bands,
                           fft_size, telescope_id, snr, max_drift, min_drift);
  config.num_bands_to_process = 1;
  VectorHitRecorder* recorder = new VectorHitRecorder();
  config.hit_recorder.reset(recorder);

  config.run();

  assert(45 == recorder->hits.size());
  assert(114049 == recorder->hits[0].index);
  assert(-1 == recorder->hits[0].drift_steps);
  assertFloatEq(7.03593, recorder->hits[0].snr);
  assertFloatEq(-0.317774, recorder->hits[0].drift_rate);
  
  return 0;
}
