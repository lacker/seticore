#include "raw_file_group_reader.h"

#include <assert.h>
#include <iostream>

using namespace std;


RawFileGroupReader::RawFileGroupReader(RawFileGroup& file_group, int num_bands,
                                       int num_batches, int blocks_per_batch)
  : file_group(file_group), num_bands(num_bands), num_batches(num_batches),
    blocks_per_batch(blocks_per_batch), destroy(false)  {
  io_thread = thread(&RawFileGroupReader::runIOThread, this);
}

RawFileGroupReader::~RawFileGroupReader() {
  unique_lock<mutex> lock(m);
  destroy = true;
  lock.unlock();
  cv.notify_all();

  if (io_thread.joinable()) {
    io_thread.join();
  }
}

shared_ptr<RawBuffer> RawFileGroupReader::makeBuffer(bool gpu) const {
  int coarse_channels_per_band = file_group.num_coarse_channels / file_group.num_bands;  
  return shared_ptr<RawBuffer>(new RawBuffer(gpu, blocks_per_batch,
                                             file_group.nants,
                                             coarse_channels_per_band,
                                             file_group.timesteps_per_block,
                                             file_group.npol));
}

shared_ptr<RawBuffer> RawFileGroupReader::read() {
  unique_lock<mutex> lock(m);
  while (buffer_queue.empty()) {
    cv.wait(lock);
  }
  auto buffer = buffer_queue.front();
  buffer_queue.pop();
  lock.unlock();
  cv.notify_one();
  return buffer;
}

// Returns false if the reader gets destroyed before a new item is pushed
bool RawFileGroupReader::push(shared_ptr<RawBuffer> buffer) {
  unique_lock<mutex> lock(m);
  while (!destroy && buffer_queue.size() >= BUFFER_QUEUE_MAX_SIZE) {
    cv.wait(lock);
  }
  if (destroy) {
    return false;
  }
  buffer_queue.push(buffer);
  lock.unlock();
  cv.notify_one();
  return true;
}

// Reads all the input and passes it to the buffer_queue
void RawFileGroupReader::runIOThread() {
  for (int band = 0; band < num_bands; ++band) {
    file_group.resetBand(band);
    for (int batch = 0; batch < num_batches; ++batch) {
      auto buffer = makeBuffer(false);

      for (int block = 0; block < buffer->num_blocks; ++block) {
        file_group.read(buffer->blockPointer(block));
      }

      if (!push(buffer)) {
        return;
      }
    }
  }
}
