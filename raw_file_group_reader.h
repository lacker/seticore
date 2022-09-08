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
  The RawFileGroupReader is the most efficient way to read data from a RawFileGroup
  onto the GPU. It uses two sets of buffers internally, one in pinned memory on the
  CPU, and another on the GPU.
  It only supports one particular access pattern, where we read one frequency band
  of the data at a time, looping all the way through all the raw files for each band.  

  This reads in "batches" rather than blocks. Each batch is a certain number of blocks,
  defined by blocks_per_batch. We read num_batches batches for each band. If there
  are leftover blocks when we are done reading batches, we just skip them and move
  on to the next band.

  num_bands doesn't have to match the underlying RawFileGroup. You can process only
  a subrange by specifying start_band and num_bands.

  Externally, the RawFileGroupReader should be used from a single thread. Call
  readToDevice() to get the next DeviceRawBuffer, and when you're done with it, call
  release() on it to reuse the pinned memory. You probably only want to
  create one at a time, since it will attempt to use most of the available
  memory for this buffer.

  Alternatively, if you don't want data to go onto the GPU, call readToHost() to get
  the next RawBuffer, and when you're done with it, call reader.returnBuffer.

  Internally, there are several threads doing things.
  runInputThread is reading buffers and passing them to buffer_queue so that they
  are ready when the client thread calls read(). This is possible since the access
  pattern is defined when the RawFileGroupReader is created.
  The input thread itself creates multiple helper threads to do the file reading. In
  testing this does seem to help a significant amount on SSDs, even though it makes the
  reads happen out of order.
  Finally, the DeviceRawBuffer itself maintains a cuda stream to be used just for
  the CPU -> GPU copy. Synchronization there happens with locks rather than cuda
  stream synchronization.

  Client code should be able to ignore the multithreading and just treat the
  RawFileGroupReader like a single-threaded reader, as long as it calls release() when
  done using the DeviceRawBuffer.
 */
class RawFileGroupReader {
 public:
  RawFileGroup& file_group;
  const int num_bands;
  const int first_band;
  const int last_band;
  const int num_batches;
  const int blocks_per_batch;
  const int coarse_channels_per_band;
  
  RawFileGroupReader(RawFileGroup& file_group, int num_bands,
                     int first_band, int last_band,
                     int num_batches, int blocks_per_batch);
  ~RawFileGroupReader();

  unique_ptr<RawBuffer> readToHost();

  // The caller can return ownership of extra buffers to avoid future mallocs
  // For convenience, returning nullptr is a no-op
  void returnBuffer(unique_ptr<RawBuffer> buffer);

  shared_ptr<DeviceRawBuffer> readToDevice();
  
 private:
  int buffer_queue_max_size;
  
  mutex m;
  condition_variable cv;
  bool destroy;
  thread io_thread;  

  // Makes a buffer of the same size the read() returns, in GPU memory
  shared_ptr<DeviceRawBuffer> makeDeviceBuffer();  

  // Makes a buffer of the correct size, reusing our pool of extra buffers if possible
  unique_ptr<RawBuffer> makeBuffer();
  
  void runInputThread();

  // Pushes this buffer onto the output queue, waiting if necessary
  bool push(unique_ptr<RawBuffer> buffer);

  // The read buffer containing data ready for transfer to the GPU
  unique_ptr<RawBuffer> read_buffer;

  // The buffer containing data on the GPU, for client code to use
  shared_ptr<DeviceRawBuffer> device_raw_buffer;
  
  // Buffers that contain data we will need in the future
  queue<unique_ptr<RawBuffer> > buffer_queue;

  // Buffers that contain nothing useful
  queue<unique_ptr<RawBuffer> > extra_buffers;
};
