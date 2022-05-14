using namespace std;

#include <fstream>
#include <iostream>

#include "fil_file.h"

/*
  Opens a sigproc filterbank file for reading.
*/
FilFile::FilFile(const string& filename) : FilterbankFile(filename),
                                           file(filename, ifstream::binary) {
  for (int i = 0; i < 4; ++i) {
    char c = file.get();
    printf("%d\n", c);
  }
  cerr << "TODO: finish implementing FilFile constructor\n";
  exit(1);
}

// Loads the data in row-major order.
void FilFile::loadCoarseChannel(int i, float* output) const {
  cerr << "TODO: implement loadCoarseChannel\n";
  exit(1);
}

FilFile::~FilFile() {
  // TODO: clean up file handle
}
