#pragma once

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
  string source_name;
  double fch1, foff, tstart, tsamp, src_dej, src_raj;
  int num_timesteps, num_freqs, coarse_channel_size, num_coarse_channels, telescope_id;
  
  virtual ~FilterbankMetadata() {}

 protected:
  void inferMetadata();
};
