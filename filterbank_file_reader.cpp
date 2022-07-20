#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <memory>

#include "filterbank_file_reader.h"
#include "fil_reader.h"
#include "h5_reader.h"

using namespace std;

// This is allegedly a SIGPROC standard but the most authoritative source
// I can find is:
//   https://github.com/UCBerkeleySETI/blimpy/blob/master/blimpy/ephemeris/observatory_info.csv
const int NO_TELESCOPE_ID = -1;
const int PARKES = 4;
const int GREEN_BANK = 6;
const int ATA = 9;
const int VLA = 12;
const int MEERKAT = 64;

/*
  Guesses some metadata from other metadata:
    coarse_channel_size
    has_dc_spike
    num_coarse_channels
 */
void FilterbankFileReader::inferMetadata() {
  if (num_timesteps == 16 && num_freqs % 1048576 == 0) {
    // Looks like Green Bank data
    assert(telescope_id == NO_TELESCOPE_ID || telescope_id == GREEN_BANK);
    if (coarse_channel_size == 0) {
      coarse_channel_size = 1048576;
    }
    has_dc_spike = true;
  } else if (num_freqs == 50331648) {
    // Looks like ATA data
    assert(telescope_id == NO_TELESCOPE_ID || telescope_id == ATA);
    if (coarse_channel_size == 0) {
      coarse_channel_size = 262144;
    }
    has_dc_spike = false;
  } else if (coarse_channel_size > 0) {
    // We already have a coarse channel size set, so we only want to infer whether we
    // have a dc spike.
    has_dc_spike = (telescope_id == GREEN_BANK);
  } else {
    cerr << "unable to infer coarse channel size for data with dimensions: " << num_timesteps
         << " x " << num_freqs << ". please set the nfpc header." << endl;
    exit(1);
  }
  num_coarse_channels = num_freqs / coarse_channel_size;
}

void FilterbankFileReader::loadCoarseChannel(int i, FilterbankBuffer* buffer) const {
  cerr << "cannot loadCoarseChannel from a base FilterbankFileReader: "
       << filename << "\n";
  exit(1);
}

unique_ptr<FilterbankFileReader> loadFilterbankFile(const string& filename) {
  if (boost::algorithm::ends_with(filename, ".h5")) {
    return unique_ptr<FilterbankFileReader>(new H5Reader(filename));
  }

  if (boost::algorithm::ends_with(filename, ".fil")) {
    return unique_ptr<FilterbankFileReader>(new FilReader(filename));
  }
  
  cerr << "could not recognize the file type of " << filename << endl;
  exit(1);
}
