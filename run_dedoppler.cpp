#include <fmt/core.h>
#include <iostream>
#include <vector>

#include "dedoppler.h"
#include "filterbank_file.h"
#include "hit_recorder.h"
#include "run_dedoppler.h"
#include "util.h"

using namespace std;

/*
  Runs a dedoppler algorithm, reading data from an input file, writing the results to a .dat file.

  input_filename is assumed to be in the .h5 format
  output_filename is where a .dat file will be written with results
  max_drift is the maximum drift we are looking for, in Hz/sec
  min_drift is the minimum drift we are looking for. Set to 0 if zero-drift signals are okay
  snr_threshold is the minimum SNR we require to report a signal

  Note that this algorithm does require an input file. In particular, the hit recorder
  copies over some metadata from it. If you wanted to run an algorithm similar to this one
  without any file existing, like with data that was already on the CPU or GPU, the
  best strategy would probably be to keep makeHitRecorder as is, but make some non-file
  object that satisfied the FilterbankFile interface, just to provide metadata.
 */
void runDedoppler(const string& input_filename, const string& output_filename,
                  double max_drift, double min_drift, double snr_threshold) {
  auto file = loadFilterbankFile(input_filename);
  auto recorder = makeHitRecorder(output_filename, *file.get(), max_drift);

  Dedopplerer dedopplerer(file->num_timesteps, file->coarse_channel_size, file->foff,
                          file->tsamp, file->has_dc_spike);
  FilterbankBuffer buffer(roundUpToPowerOfTwo(file->num_timesteps), file->coarse_channel_size);
  
  // Load and process one coarse channel at a time from the filterbank file
  vector<DedopplerHit> hits;
  for (int coarse_channel = 0; coarse_channel < file->num_coarse_channels; ++coarse_channel) {
    file->loadCoarseChannel(coarse_channel, &buffer);
    hits.clear();
    dedopplerer.search(buffer, max_drift, min_drift, snr_threshold, &hits);
    for (DedopplerHit hit : hits) {

      cout << "hit: coarse channel = " << coarse_channel << ", "
           << hit.toString() << endl;
        
      recorder->recordHit(hit, NO_BEAM, coarse_channel, buffer.data);
    }
  }
}
