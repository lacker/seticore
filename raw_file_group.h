#pragma once

#include <memory>
#include <string>
#include <vector> 

#include "raw/raw.h"

using namespace std;

vector<vector<string> > scanForRawFileGroups(const string& directory);

/*
  This class provides a unified interface to open a group of
  sequentially-recorded raw files.
  The RawFileGroup also encapsulates missing-block handling, replacing
  missing blocks with zeros.
  It is not threadsafe.
*/
class RawFileGroup {
 private:
  // When current_file is -1, these should not be used.
  // Otherwise, header contains the last read from reader.
  raw::Header header;
  unique_ptr<raw::Reader> reader;

  // Index of the current file in filenames that reader is opening.
  // -1 if there is none.
  int current_file;

  // The next pktidx that we will return from read()
  int next_pktidx;

  // Helper to advance the reader and header
  void openNextFile();
  
 public:
  const vector<string> filenames;

  // Has the directory name stripped
  string prefix;
  
  // Defining a sub-band of the files we are reading in
  int band;
  const int num_bands;

  // Metadata that should be the same for all blocks
  int nants;
  int num_coarse_channels;
  int npol;
  int schan;
  double obsbw;
  double obsfreq;
  string obsid;
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

  // How many bytes each read returns 
  int read_size;

  // Starts off reading band zero
  RawFileGroup(const vector<string>& filenames, int num_bands);

  ~RawFileGroup();

  void resetBand(int new_band);

  void read(char* buffer);

  double getStartTime(int block);
};
