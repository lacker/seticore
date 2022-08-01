#include <fmt/core.h>
#include <iostream>
#include <time.h>

#include "raw_file_group.h"
#include "raw_file_group_reader.h"

using namespace std;

/*
  This benchmark reads the first out of 16 bands for the first raw
  file group in the ../benchmark directory.

  You have to drop disk caches first for this test to be meaningful:

  echo 3 | sudo tee /proc/sys/vm/drop_caches
 */
int main(int argc, char* argv[]) {
  // Specifying parameters
  string dir = "../benchmark";
  auto file_lists = scanForRawFileGroups(dir);
  assert(file_lists.size() == 1);

  int tstart = time(NULL);
  
  int num_bands = 16;
  RawFileGroup file_group(file_lists[0], num_bands);

  int blocks_per_batch = 128;
  int num_batches = file_group.num_blocks / blocks_per_batch;

  // Only process one band
  int num_bands_to_process = 1;
  RawFileGroupReader reader(file_group, num_bands_to_process, num_batches,
                            blocks_per_batch);

  for (int band = 0; band < num_bands_to_process; ++band) {
    for (int batch = 0; batch < num_batches; ++batch) {
      auto buffer = reader.read();
      reader.returnBuffer(move(buffer));
      cerr << "done band " << band << " batch " << batch << endl;
    }
  }

  int tstop = time(NULL);
  cerr << fmt::format("file io benchmark elapsed time {:d}s\n", tstop - tstart);  
}
