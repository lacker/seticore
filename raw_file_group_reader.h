#pragma once

#include <condition_variable>
#include <memory>
#include <mutex>
#include <queue>
#include <thread>

#include "device_raw_buffer.h"
#include "raw_buffer.h"
#include "raw_file_group.h"

using namespace std;

const int BUFFER_QUEUE_MAX_SIZE = 4;

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
  read() should only be called from a single thread.
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

  // Makes a buffer of the correct size, reusing our pool of extra buffers if possible
  unique_ptr<RawBuffer> makeBuffer();

  // Makes a buffer of the correct size in GPU memory
  shared_ptr<DeviceRawBuffer> makeDeviceBuffer();
  
  unique_ptr<RawBuffer> read();

  // The caller can return ownership of extra buffers to avoid future mallocs
  void returnBuffer(unique_ptr<RawBuffer> buffer);

 private:
  mutex m;
  condition_variable cv;
  bool destroy;
  thread io_thread;  
  
  void runIOThread();

  // Pushes this buffer onto the output queue, waiting if necessary
  bool push(unique_ptr<RawBuffer> buffer);
  
  // Outputs
  queue<unique_ptr<RawBuffer> > buffer_queue;

  queue<unique_ptr<RawBuffer> > extra_buffers;
};
