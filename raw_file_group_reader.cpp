#include "raw_file_group_reader.h"

#include <assert.h>
#include <fmt/core.h>
#include <iostream>
#include <sys/sysinfo.h>
#include "thread_util.h"
#include "util.h"

using namespace std;


/*
  Reads [first_band, last_band], inclusive, out of num_bands total.
 */
RawFileGroupReader::RawFileGroupReader(RawFileGroup& file_group, int num_bands,
                                       int first_band, int last_band,
                                       int num_batches, int blocks_per_batch)
  : file_group(file_group), num_bands(num_bands), first_band(first_band),
    last_band(last_band), num_batches(num_batches), blocks_per_batch(blocks_per_batch),
    coarse_channels_per_band(file_group.num_coarse_channels / num_bands),
    stopped(false)  {

  // Limit queue size depending on total memory.
  struct sysinfo info;
  sysinfo(&info);
  size_t mb = 1024 * 1024;
  size_t gb = mb * 1024;
  if ((size_t) info.totalram > 50 * gb) {
    // Looks like a production machine.
    size_t buffer_size = rawBufferSize(blocks_per_batch, file_group.nants,
                                       coarse_channels_per_band,
                                       file_group.timesteps_per_block,
                                       file_group.npol);
    buffer_queue_max_size = (int) (0.25 * info.totalram / buffer_size);
    cout << fmt::format("limiting raw file input buffer memory to {:.1f} GB\n",
                        1.0 * buffer_queue_max_size * buffer_size / gb);
  } else {
    // Looks like a dev machine.
    buffer_queue_max_size = 4;
  }

  io_thread = thread(&RawFileGroupReader::runInputThread, this);

  device_raw_buffer = makeDeviceBuffer();
  cout << "raw buffer memory: " << prettyBytes(device_raw_buffer->size) << endl;
}

void RawFileGroupReader::stop() {
  unique_lock<mutex> lock(m);
  stopped = true;
  lock.unlock();
  cv.notify_all();
}

RawFileGroupReader::~RawFileGroupReader() {
  stop();
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

  return make_unique<RawBuffer>(blocks_per_batch,
                                file_group.nants,
                                coarse_channels_per_band,
                                file_group.timesteps_per_block,
                                file_group.npol);
}

shared_ptr<DeviceRawBuffer> RawFileGroupReader::makeDeviceBuffer() {
  return make_shared<DeviceRawBuffer>(blocks_per_batch,
                                      file_group.nants,
                                      coarse_channels_per_band,
                                      file_group.timesteps_per_block,
                                      file_group.npol);
}

unique_ptr<RawBuffer> RawFileGroupReader::readToHost() {
  unique_lock<mutex> lock(m);
  while (!stopped && buffer_queue.empty()) {
    cv.wait(lock);
  }

  if (stopped) {
    fatal("RawFileGroupReader stopped");
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

// Returns false if the reader gets stopped before a new item is pushed
bool RawFileGroupReader::push(unique_ptr<RawBuffer> buffer) {
  unique_lock<mutex> lock(m);
  while (!stopped && (int) buffer_queue.size() >= buffer_queue_max_size) {
    cv.wait(lock);
  }

  if (stopped) {
    return false;
  }
  buffer_queue.push(move(buffer));
  lock.unlock();
  cv.notify_one();
  return true;
}

// Reads all the input and passes it to the buffer_queue
void RawFileGroupReader::runInputThread() {
  setThreadName("input");
  for (int band = first_band; band <= last_band; ++band) {
    file_group.resetBand(band, num_bands);
    for (int batch = 0; batch < num_batches; ++batch) {
      if (stopped) {
        return;
      }

      auto buffer = makeBuffer();

      vector<function<bool()> > tasks;
      for (int block = 0; block < buffer->num_blocks; ++block) {
        file_group.readTasks(buffer->blockPointer(block), &tasks);
      }

      // Testing on meerkat, any more than 4 threads doesn't help
      int num_threads = 4;
      if (!runInParallel(move(tasks), num_threads)) {
        stop();
        return;
      }

      if (!push(move(buffer))) {
        return;
      }
    }
  }
}

shared_ptr<DeviceRawBuffer> RawFileGroupReader::readToDevice() {
  // Client code could still be using the raw buffers.
  // Wait for it to finish.
  device_raw_buffer->waitUntilUnused();

  returnBuffer(move(read_buffer));
  read_buffer = readToHost();
  
  // Useful for debugging the output of raw file reading
  // cerr << "read_buffer[0] = " << (int) read_buffer->data[0] << endl;

  device_raw_buffer->copyFromAsync(*read_buffer);
  device_raw_buffer->waitUntilReady();

  return device_raw_buffer;
}
