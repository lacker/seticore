#include <fmt/core.h>
#include <iostream>
#include <time.h>

#include "raw_file_group.h"
#include "raw_file_group_reader.h"

using namespace std;

/*
  This benchmark reads the first bands for the first raw
  file group in the ../benchmark directory, or the provided argument.

  You may have to drop disk caches first for this test to be meaningful:

  echo 3 | sudo tee /proc/sys/vm/drop_caches
 */
int main(int argc, char* argv[]) {
  // Specifying parameters
  string dir;
  if (argc < 2) {
    dir = "../benchmark";
  } else {
    dir = argv[1];
  }
  auto file_lists = scanForRawFileGroups(dir);
  assert(file_lists.size() == 1);

  int tstart = time(NULL);
  
  int num_bands = 16;
  RawFileGroup file_group(file_lists[0], num_bands);

  int blocks_per_batch = 128;
  int num_batches = file_group.num_blocks / blocks_per_batch;
  long bytes_read = 0;
  
  // Only process some bands
  int num_bands_to_process = 1;
  RawFileGroupReader reader(file_group, num_bands_to_process, num_batches,
                            blocks_per_batch);

  for (int band = 0; band < num_bands_to_process; ++band) {
    for (int batch = 0; batch < num_batches; ++batch) {
      auto buffer = reader.read();
      bytes_read += buffer->data_size;
      reader.returnBuffer(move(buffer));
      cerr << "done band " << band << " batch " << batch << endl;
    }
  }

  int tstop = time(NULL);
  int elapsed = tstop - tstart;
  cerr << fmt::format("file io benchmark elapsed time {:d}s\n", elapsed);
  float giga = 1024.0 * 1024.0 * 1024.0;
  float gb = bytes_read / giga;
  float gbps = gb / elapsed;
  cerr << fmt::format("{:.1f} GB read at a rate of {:.2f} GB/s\n", gb, gbps);
}
