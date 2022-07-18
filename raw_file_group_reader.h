#pragma once

#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>

#include "raw_buffer.h"
#include "raw_file_group.h"

using namespace std;

const int BUFFER_QUEUE_MAX_SIZE = 2;

/*
  This class provides an efficient reader to loop through a RawFileGroup.

  This reads in "batches" rather than blocks. Each batch is a certain number of blocks,
  defined by blocks_per_batch. We read num_batches batches for each band. If there
  are leftover blocks when we are done reading batches, we just skip them and move
  on to the next band.

  If you provide it with a different num_bands than the RawFileGroup has, it will
  only process the provided number of bands.

  The file IO is done in a separate IO thread, reading ahead a few batches,
  synchronizing on buffer_queue so the caller can read() on the user thread.
 */
class RawFileGroupReader {
 public:
  RawFileGroup& file_group;
  const int num_bands;
  const int num_batches;
  const int blocks_per_batch;

  RawFileGroupReader(RawFileGroup& file_group, int num_bands, int num_batches,
                     int blocks_per_batch);
  ~RawFileGroupReader();

  // Just makes a buffer of the correct size
  shared_ptr<RawBuffer> makeBuffer(bool gpu) const;
  
  shared_ptr<RawBuffer> read();
  
 private:
  mutex m;
  condition_variable cv;
  bool destroy;
  thread io_thread;  
  
  void runIOThread();

  // Pushes this buffer onto the output queue, waiting if necessary
  bool push(shared_ptr<RawBuffer> buffer);
  
  // Outputs
  queue<shared_ptr<RawBuffer> > buffer_queue;
};
