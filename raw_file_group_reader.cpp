#include "raw_file_group_reader.h"

#include <assert.h>
#include <iostream>

#include "thread_util.h"

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

unique_ptr<RawBuffer> RawFileGroupReader::makeBuffer() {
  unique_lock<mutex> lock(m);
  if (!extra_buffers.empty()) {
    auto buffer = move(extra_buffers.front());
    extra_buffers.pop();
    return buffer;
  }
  lock.unlock();

  int coarse_channels_per_band = file_group.num_coarse_channels / file_group.num_bands;  
  return make_unique<RawBuffer>(blocks_per_batch,
                                file_group.nants,
                                coarse_channels_per_band,
                                file_group.timesteps_per_block,
                                file_group.npol);
}

shared_ptr<DeviceRawBuffer> RawFileGroupReader::makeDeviceBuffer() {
  int coarse_channels_per_band = file_group.num_coarse_channels / file_group.num_bands;  
  return make_shared<DeviceRawBuffer>(blocks_per_batch,
                                      file_group.nants,
                                      coarse_channels_per_band,
                                      file_group.timesteps_per_block,
                                      file_group.npol);
}

unique_ptr<RawBuffer> RawFileGroupReader::read() {
  unique_lock<mutex> lock(m);
  while (buffer_queue.empty()) {
    cv.wait(lock);
  }
  auto buffer = move(buffer_queue.front());
  buffer_queue.pop();
  lock.unlock();
  cv.notify_one();
  return buffer;
}

void RawFileGroupReader::returnBuffer(unique_ptr<RawBuffer> buffer) {
  // We might want to check it's the right size
  if (buffer.get() == nullptr) {
    return;
  }
  unique_lock<mutex> lock(m);
  extra_buffers.push(move(buffer));
}

// Returns false if the reader gets destroyed before a new item is pushed
bool RawFileGroupReader::push(unique_ptr<RawBuffer> buffer) {
  unique_lock<mutex> lock(m);
  while (!destroy && buffer_queue.size() >= BUFFER_QUEUE_MAX_SIZE) {
    cv.wait(lock);
  }
  if (destroy) {
    return false;
  }
  buffer_queue.push(move(buffer));
  lock.unlock();
  cv.notify_one();
  return true;
}

// Reads all the input and passes it to the buffer_queue
void RawFileGroupReader::runIOThread() {
  for (int band = 0; band < num_bands; ++band) {
    file_group.resetBand(band);
    for (int batch = 0; batch < num_batches; ++batch) {
      auto buffer = makeBuffer();

      vector<function<bool()> > tasks;
      for (int block = 0; block < buffer->num_blocks; ++block) {
        file_group.readTasks(buffer->blockPointer(block), &tasks);
      }

      // Testing on meerkat, any more than 4 threads doesn't help
      int num_threads = 4;
      runInParallel(move(tasks), num_threads);

      if (!push(move(buffer))) {
        return;
      }
    }
  }
}
