#include "beamforming_pipeline.h"
#include "raw_file_group.h"

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

  TODO: make this based on smaller files which are generally available, and which
  are actually correct!
 */
int main(int argc, char* argv[]) {
  // Specifying parameters
  string input_dir = "../data/beamtest";
  string output_dir = "mocked";
  string recipe_dir = "../data/beamtest";
  int num_bands = 32;
  // int fft_size = 131072;
  int fine_channels = 8388608;
  int sti = 8;
  int telescope_id = 64;
  float snr = 7.0;
  float max_drift = 0.01;

  auto groups = scanForRawFileGroups(input_dir);
  assert(groups.size() == 1);
  assert(groups[0].size() == 2);

  BeamformingPipeline pipeline(groups[0], output_dir, recipe_dir, num_bands,
                               sti, telescope_id, snr, max_drift, -1, fine_channels);
  pipeline.num_bands_to_process = 1;
  pipeline.record_hits = false;
  
  pipeline.findHits();

  cout << pipeline.hits[0].toString() << endl;
  
  assert(98134 == pipeline.hits[0].index);
  assert(-1 == pipeline.hits[0].drift_steps);
  assertFloatEq(11.93973, pipeline.hits[0].snr);
  assertFloatEq(-0.31777, pipeline.hits[0].drift_rate);

  vector<DedopplerHit> nonzero;
  for (const DedopplerHit& hit : pipeline.hits) {
    if (hit.drift_steps != 0) {
      nonzero.push_back(hit);
    }
  }

  if (73 != nonzero.size()) {
    cerr << "nonzero size was: " << nonzero.size() << endl;
    return 1;
  }
  
  assertStringEq(nonzero[10].toString(),
                 "coarse channel = 0, index = 63010, snr = 8.19147, drift rate = -0.31777 (-1 bin)");
  assertStringEq(nonzero[20].toString(),
                 "coarse channel = 0, index = 50467, snr = 7.30951, drift rate = -0.31777 (-1 bin)");
  
  return 0;
}
