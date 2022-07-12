#include "beamforming_config.h"
#include "raw_file_group.h"

#include <string>

using namespace std;
 
int main(int argc, char* argv[]) {
  // Specifying parameters
  string input_dir = "/d/onlytwo";
  string output_dir = "./data";
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
  config.run();
  
  return 0;
}
