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

/*
  The RawFileGroupReader is the most efficient way to read data from a RawFileGroup.
  It only supports one particular access pattern, where we read one frequency band
  of the data at a time, looping all the way through all the raw files for each band.

  This reads in "batches" rather than blocks. Each batch is a certain number of blocks,
  defined by blocks_per_batch. We read num_batches batches for each band. If there
  are leftover blocks when we are done reading batches, we just skip them and move
  on to the next band.

  If you provide it with a different num_bands than the RawFileGroup has, it will
  only process the provided number of bands.

  Externally, the RawFileGroupReader should be used from a single thread. Call
  read() to get the next RawBuffer, and when you're done with it, call
  returnBuffer on it to reuse the pinned memory. You probably only want to
  create one at a time, since it will attempt to use most of the available
  memory for this buffer.

  Internally, there are several threads doing things.
  runIOThread is reading buffers and passing them to buffer_queue so that they
  are ready when the client thread calls read(). This is possible since the access
  pattern is defined when the RawFileGroupReader is created.
  The io thread itself creates multiple helper threads to do the file reading. In
  testing this does seem to help a significant amount on SSDs, even though it makes the
  reads happen out of order.

  Client code should be able to ignore the multithreading and just treat the
  RawFileGroupReader like a single-threaded reader.
 */
class RawFileGroupReader {
 public:
  RawFileGroup& file_group;
  const int num_bands;
  const int num_batches;
  const int blocks_per_batch;
  const int coarse_channels_per_band;
  
  RawFileGroupReader(RawFileGroup& file_group, int num_bands, int num_batches,
                     int blocks_per_batch);
  ~RawFileGroupReader();

  // Makes a buffer of the same size the read() returns, in GPU memory
  shared_ptr<DeviceRawBuffer> makeDeviceBuffer();
  
  unique_ptr<RawBuffer> read();

  // The caller can return ownership of extra buffers to avoid future mallocs
  void returnBuffer(unique_ptr<RawBuffer> buffer);

 private:
  int buffer_queue_max_size;
  
  mutex m;
  condition_variable cv;
  bool destroy;
  thread io_thread;  

  // Makes a buffer of the correct size, reusing our pool of extra buffers if possible
  unique_ptr<RawBuffer> makeBuffer();
  
  void runIOThread();

  // Pushes this buffer onto the output queue, waiting if necessary
  bool push(unique_ptr<RawBuffer> buffer);
  
  // Outputs
  queue<unique_ptr<RawBuffer> > buffer_queue;

  queue<unique_ptr<RawBuffer> > extra_buffers;
};
