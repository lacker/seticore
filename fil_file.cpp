using namespace std;

#include <fstream>
#include <iostream>

#include "fil_file.h"

/*
  Opens a sigproc filterbank file for reading.
*/
FilFile::FilFile(const string& filename) : FilterbankFile(filename),
                                           file(filename, ifstream::binary) {
  // Read the headers
  bool seen_header_start = false;
  while (true) {
    uint32_t num_bytes;
    file.read((char*)&num_bytes, sizeof(num_bytes));
    cerr << num_bytes << endl;
    break;
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
