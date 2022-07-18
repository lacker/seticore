#include "raw_file_group_reader.h"

#include <assert.h>

using namespace std;

RawFileGroupReader::RawFileGroupReader(RawFileGroup& file_group,
                                       int num_batches, int blocks_per_batch)
  : file_group(file_group), num_bands(file_group.num_bands),
    num_batches(num_batches), blocks_per_batch(blocks_per_batch) {
  current_band = 0;
  batches_read_in_this_band = 0;
}

RawFileGroupReader::~RawFileGroupReader() { }

shared_ptr<RawBuffer> RawFileGroupReader::makeBuffer(bool gpu) {
  int coarse_channels_per_band = file_group.num_coarse_channels / num_bands;  
  return shared_ptr<RawBuffer>(new RawBuffer(gpu, blocks_per_batch,
                                             file_group.nants,
                                             coarse_channels_per_band,
                                             file_group.timesteps_per_block,
                                             file_group.npol));
}

// Reads into CPU memory
shared_ptr<RawBuffer> RawFileGroupReader::read() {
  if (batches_read_in_this_band == num_batches) {
    ++current_band;
    file_group.resetBand(current_band);
    batches_read_in_this_band = 0;
  }
  assert(current_band < num_bands);

  // Create a buffer that can handle one batch
  auto buffer = makeBuffer(false);

  for (int block = 0; block < buffer->num_blocks; ++block) {
    file_group.read(buffer->blockPointer(block));
  }
  ++batches_read_in_this_band;

  return buffer;
}
