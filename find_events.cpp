#include "filterbank_file_reader.h"
#include "find_events.h"

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

  // TODO
}
