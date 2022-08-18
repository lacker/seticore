#pragma once

#include <string>
#include <vector>

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
  The FilterbankMetadata stores the metadata that is typically stored in conjunction
  with filterbank data.
*/
class FilterbankMetadata {
 public:
  FilterbankMetadata();

  bool has_dc_spike;

  // Represents the boresight source when there are multiple beams
  string source_name;

  // Set when there is one source per beam
  vector<string> source_names;

  double fch1, foff, tstart, tsamp, src_dej, src_raj;
  int num_timesteps, num_channels, coarse_channel_size, num_coarse_channels, telescope_id;
  
  virtual ~FilterbankMetadata() {}

  FilterbankMetadata getSubsetMetadata(int beam, int band, int num_bands);

  // Providing information about beams when metadata is split among beams
  bool isCoherent(int beam) const;
  string getSourceName(int beam) const;
  
 protected:
  void inferMetadata();
};
