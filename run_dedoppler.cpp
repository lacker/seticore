#include <vector>

#include "dat_file.h"
#include "dedoppler.h"
#include "filterbank_file.h"
#include "run_dedoppler.h"

using namespace std;

/*
  Runs a dedoppler algorithm, reading data from an input file, writing the results to a .dat file.

  input_filename is assumed to be in the .h5 format
  output_filename is where a .dat file will be written with results
  max_drift is the maximum drift we are looking for, in Hz/sec
  min_drift is the minimum drift we are looking for. Set to 0 if zero-drift signals are okay
  snr_threshold is the minimum SNR we require to report a signal
 */
void runDedoppler(const string& input_filename, const string& output_filename,
                  double max_drift, double min_drift, double snr_threshold) {
  auto file = loadFilterbankFile(input_filename);
  
  DatFile output(output_filename, *file.get(), max_drift);

  Dedopplerer dedopplerer(file->num_timesteps, file->coarse_channel_size, file->foff,
                          file->tsamp, file->has_dc_spike);
  
  // Load and process one coarse channel at a time from the hdf5
  vector<DedopplerHit> hits;
  for (int coarse_channel = 0; coarse_channel < file->num_coarse_channels; ++coarse_channel) {
    file->loadCoarseChannel(coarse_channel, dedopplerer.input);
    hits.clear();
    dedopplerer.processInput(max_drift, min_drift, snr_threshold, &hits);
    for (DedopplerHit hit : hits) {
      output.reportHit(coarse_channel, hit.index, hit.drift_steps, hit.drift_rate, hit.snr);
    }
  }
}
