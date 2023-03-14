#include <assert.h>
#include "dedoppler.h"
#include "filterbank_buffer.h"
#include "filterbank_file_reader.h"
#include "find_events.h"
#include "util.h"

using namespace std;

/*
  Runs a dedoppler algorithm across files in a cadence, assuming they follow an
  "ABACAD" pattern.
  The results that appear in "on" but not "off" data are written to an .events file.

  max_drift is the maximum drift we are looking for, in Hz/sec
  snr_on is the minimum snr we require in the "on" data.
  snr_off is the minimum snr that makes us throw out an event when it
  occurs in the "off" data.
 */
void findEvents(const vector<string>& input_filenames, const string& output_filename,
                double max_drift, double snr_on, double snr_off) {
  vector<unique_ptr<FilterbankFileReader> > files;
  vector<shared_ptr<Dedopplerer> > dedopplerers;
  vector<FilterbankBuffer> buffers;
  
  for (auto& filename : input_filenames) {
    files.push_back(move(loadFilterbankFile(filename)));
    auto& file = files.back();
    shared_ptr<Dedopplerer> dedopplerer(new Dedopplerer(file->num_timesteps,
                                                        file->coarse_channel_size,
                                                        file->foff, file->tsamp,
                                                        file->has_dc_spike));
    dedopplerers.push_back(dedopplerer);
    buffers.emplace_back(roundUpToPowerOfTwo(file->num_timesteps),
                         file->coarse_channel_size);
  }

  // Check the metadata lines up
  int num_timesteps = files[0]->num_timesteps;
  int coarse_channel_size = files[0]->coarse_channel_size;
  double foff = files[0]->foff;
  double tsamp = files[0]->tsamp;  
  bool has_dc_spike = files[0]->has_dc_spike;
  string source_name = files[0]->source_name;

  for (int i = 0; i < (int) files.size(); ++i) {
    assert(files[i]->num_timesteps == num_timesteps);
    assert(files[i]->coarse_channel_size == coarse_channel_size);
    assertFloatEq(files[i]->foff, foff);
    assertFloatEq(files[i]->tsamp, tsamp);
    assert(files[i]->has_dc_spike == has_dc_spike);

    // The targets should work like ABACAD
    if (i % 2 == 0) {
      assert(files[i]->source_name == source_name);
    } else {
      assert(files[i]->source_name != source_name);
    }
  }
  
  // TODO
}
