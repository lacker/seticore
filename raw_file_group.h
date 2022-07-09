#pragma once

#include <string>
#include <vector>

using namespace std;

/*
  This class provides a unified interface to open a group of
  sequentially-recorded raw files.
  The RawFileGroup also encapsulates missing-block handling, replacing
  missing blocks with zeros.
*/
class RawFileGroup {
 public:
  // Defining a sub-band of the files we are reading in
  const int band;
  const int num_bands;

  // Metadata that should be the same for all blocks
  int nants;
  int num_coarse_channels;
  int npol;
  double obsbw;
  string source_name;
  double tbin;
  int timesteps_per_block;

  // Metadata used for timing
  int synctime;
  int piperblk;  

  // Metadata for the very first block
  int start_pktidx;
  
  // The total number of blocks in this set of raw files.
  // This includes missing blocks. 
  int num_blocks;

  RawFileGroup(const vector<string>& filenames, int band, int num_bands);
  ~RawFileGroup();
};
