#include <fmt/core.h>
#include <iostream>

#include "raw_file.h"

using namespace std;

/*
  Displays block-by-block information about a raw file.

  Usage:
    blockls <filename>
 */
int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "usage: blockls <filename>\n";
    exit(1);
  }
  string fname(argv[1]);

  RawFile f(fname);
  for (auto& header : f.headers()) {
    vector<char> buffer(header.blocsize);
    assert(f.reader().readBand(header, 0, 1, &buffer[0]));
    cout << "pktidx " << header.pktidx << ", data[:3] = "
         << (int)(buffer[0]) << ", "
         << (int)(buffer[1]) << ", "
         << (int)(buffer[2]) << endl;
  }
}
