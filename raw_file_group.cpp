#include "raw_file_group.h"

#include <assert.h>
#include <iostream>

#include "raw/raw.h"

using namespace std;

RawFileGroup::RawFileGroup(const vector<string>& filenames, int band, int num_bands)
  : band(band), num_bands(num_bands) {
  assert(!filenames.empty());

  // Get metadata from the first file
  raw::Reader first_reader(filenames[0]);
  raw::Header header;
  if (!first_reader.readHeader(&header)) {
    cerr << "error reading raw file " << filenames[0] << endl;
    exit(1);
  }
  
  nants = header.nants;
  num_coarse_channels = header.num_channels;
  npol = header.npol;
  obsbw = header.obsbw;
  source_name = header.src_name;
  start_pktidx = header.pktidx;
  tbin = header.tbin;
  timesteps_per_block = header.num_timesteps;
  
  synctime = header.getUnsignedInt("SYNCTIME", -1);
  assert(synctime > 0);
  piperblk = header.getUnsignedInt("PIPERBLK", -1);
  assert(piperblk > 0);

  // TODO: set num_blocks 
}

RawFileGroup::~RawFileGroup() {}
