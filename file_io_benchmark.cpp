#include <fmt/core.h>
#include <iostream>
#include "util.h"

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

  int num_bands = 1;
  RawFileGroup file_group(file_lists[0]);

  int blocks_per_batch = 128;
  int num_batches = file_group.num_blocks / blocks_per_batch;
  
  // Only process some bands
  int num_bands_to_process = 1;
  RawFileGroupReader reader(file_group, num_bands, 0, num_bands_to_process - 1,
                            num_batches, blocks_per_batch);

  for (int band = 0; band < num_bands_to_process; ++band) {
    long tstart = timeInMS();
    long bytes_read = 0;
  
    for (int batch = 0; batch < num_batches; ++batch) {
      auto buffer = reader.readToHost();
      bytes_read += buffer->size;
      reader.returnBuffer(move(buffer));
      cerr << "done band " << band << " batch " << batch << endl;
    }

    long tstop = timeInMS();
    long elapsed_ms = tstop - tstart;
    float elapsed_s = elapsed_ms / 1000.0;
    cerr << "band " << band << " stats:\n";
    cerr << fmt::format("file io benchmark elapsed time {:.3f}s\n", elapsed_s);
    float giga = 1024.0 * 1024.0 * 1024.0;
    float gb = bytes_read / giga;
    float gbps = gb / elapsed_s;
    cerr << fmt::format("{:.1f} GB read at a rate of {:.2f} GB/s\n", gb, gbps);
    
  }

}
