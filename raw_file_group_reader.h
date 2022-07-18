#pragma once

#include <memory>

#include "raw_buffer.h"
#include "raw_file_group.h"

using namespace std;

/*
  This class provides an efficient reader to loop through a RawFileGroup.

  It splits up the raw file into num_bands bands, so that all bands are always handled.

  This reads in "batches" rather than blocks. Each batch is a certain number of blocks,
  defined by blocks_per_batch. We read num_batches batches for each band. If there
  are leftover blocks when we are done reading batches, we just skip them and move
  on to the next band.
 */
class RawFileGroupReader {
 public:
  RawFileGroup& file_group;
  const int num_bands;
  const int num_batches;
  const int blocks_per_batch;

  RawFileGroupReader(RawFileGroup& file_group, int num_batches, int blocks_per_batch);
  ~RawFileGroupReader();

  shared_ptr<RawBuffer> read();
  
 private:
  int current_band;
  int batches_read_in_this_band;
};
