#include <fmt/core.h>
#include <iostream>

#include "raw_file_group.h"
#include "util.h"

using namespace std;

/*
  Displays information about a directory full of raw files.

  Usage:
    rawls <directory>
 */
int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "usage: rawls <directory>\n";
    exit(1);
  }
  string dir(argv[1]);

  auto file_lists = scanForRawFileGroups(dir);
  cout << dir << " has " << pluralize(file_lists.size(), "raw file group") << ":\n";

  for (auto& file_list : file_lists) {
    RawFileGroup group(file_list, 1);
    cout << endl;
    cout << pluralize(group.filenames.size(), "file") << " match " << group.prefix
         << ".*.raw\n";
    cout << pluralize(group.nants, "antenna") << endl;
    cout << fmt::format("{:.1f} MHz bandwidth\n", group.obsbw);
    cout << fmt::format("{:.1f} MHz center frequency\n", group.obsfreq);
    cout << pluralize(group.num_coarse_channels, "coarse channel") << endl;
    cout << pluralize(group.timesteps_per_block, "timestep") << " per block\n";
    cout << (group.read_size / 1024 / 1024) << " MB block size\n";
    cout << pluralize(group.num_blocks, "total block") << endl;
    cout << fmt::format("{:.1f} GB total data\n", group.totalDataGB());
    cout << fmt::format("{:.1f}s total time\n", group.totalTime());
  }
}
