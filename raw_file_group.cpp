#include "raw_file_group.h"

#include <assert.h>
#include <iostream>

using namespace std;

RawFileGroup::RawFileGroup(const vector<string>& filenames, int band, int num_bands)
  : current_file(-1), filenames(filenames), band(band), num_bands(num_bands) {
  assert(!filenames.empty());

  // Get metadata from the first file
  raw::Reader first_reader(filenames[0]);
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
  next_pktidx = start_pktidx;
  tbin = header.tbin;
  timesteps_per_block = header.num_timesteps;
  
  synctime = header.getUnsignedInt("SYNCTIME", -1);
  assert(synctime > 0);
  piperblk = header.getUnsignedInt("PIPERBLK", -1);
  assert(piperblk > 0);

  // Calculate the size of each read
  assert(0 == num_coarse_channels % num_bands);
  int channels_per_band = num_coarse_channels / num_bands;
  read_size = nants * channels_per_band * timesteps_per_block * npol * 2;
  
  // Find the last block in the last file
  raw::Reader last_reader(filenames[filenames.size() - 1]);
  while(last_reader.readHeader(&header)) {}
  int pktidx_diff = header.pktidx - start_pktidx;
  assert(0 == pktidx_diff % piperblk);
  num_blocks = (pktidx_diff / piperblk) + 1;
}

RawFileGroup::~RawFileGroup() {}

void RawFileGroup::openNextFile() {
  ++current_file;
  assert(current_file < (int) filenames.size());
  
  reader.reset(new raw::Reader(filenames[current_file]));
  if (!reader->readHeader(&header)) {
    cerr << "error reading first block in " << filenames[current_file] << endl;
    exit(1);
  }

  // Sanity check some values in this header
  assert(nants == header.nants);
  assert(num_coarse_channels == header.num_channels);
  assert(npol == (int) header.npol);  
}

void RawFileGroup::read(char* buffer) {
  if (current_file == -1) {
    openNextFile();
  }

  if (header.pktidx < next_pktidx) {
    // We need to advance the header
    if (reader->readHeader(&header)) {
      if (header.pktidx < next_pktidx) {
        cerr << "error in reading " << reader->filename << " - saw pktidx "
             << header.pktidx << " when expecting at least pktidx "
             << next_pktidx << endl;
        exit(1);
      }
    } else if (reader->error()) {
      cerr << "reader error: " << reader->errorMessage() << endl;
      exit(1);
    } else {
      // This is just the end of a file
      openNextFile();
      if (header.pktidx < next_pktidx) {
        cerr << "first pktidx in " << reader->filename << " is only "
             << header.pktidx << " when we expected at least " << next_pktidx
             << endl;
        exit(1);
      }
    }
  }

  if (header.pktidx == next_pktidx) {
    // We're pointing at the data we want to return
    reader->readBand(header, band, num_bands, buffer);
    next_pktidx += piperblk;
    return;
  }

  // We missed some blocks, so we'll have to return some zeros
  assert(header.pktidx > next_pktidx);
  cout << "missing block with pktidx = " << next_pktidx << endl;
  memset(buffer, 0, read_size);
  next_pktidx += piperblk;
}

// We can't calculate the time from the actual block, because we
// might have missed that block.
double RawFileGroup::getStartTime(int block) {
  assert(block < num_blocks);
  double time_per_block = tbin * timesteps_per_block;
  return synctime + block * time_per_block;
}
