#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <memory>

#include "filterbank_file.h"
#include "fil_file.h"
#include "h5_file.h"

using namespace std;

/*
  Guesses some metadata from other metadata:
    coarse_channel_size
    has_dc_spike
    num_coarse_channels
 */
void FilterbankFile::inferMetadata() {
  if (num_timesteps == 16 && num_freqs % 1048576 == 0) {
    // Looks like Green Bank data
    assert(telescope_id == -1 || telescope_id == 6);
    coarse_channel_size = 1048576;
    has_dc_spike = true;
  } else if (num_freqs == 50331648) {
    // Looks like ATA data
    assert(telescope_id == -1 || telescope_id == 9);
    coarse_channel_size = 262144;
    has_dc_spike = false;
  } else {
    cerr << "unrecognized data dimensions: " << num_timesteps << " x " << num_freqs << endl;
    exit(1);
  }
  num_coarse_channels = num_freqs / coarse_channel_size;
}

unique_ptr<FilterbankFile> loadFilterbank(const string& filename) {
  if (boost::algorithm::ends_with(filename, ".h5")) {
    return unique_ptr<FilterbankFile>(new H5File(filename));
  }

  if (boost::algorithm::ends_with(filename, ".fil")) {
    return unique_ptr<FilFile>(new FilFile(filename));
  }
  
  cerr << "could not recognize the file type of " << filename << endl;
  exit(1);
}
