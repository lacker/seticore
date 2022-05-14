using namespace std;

#include <assert.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "fil_file.h"

/*
  Opens a sigproc filterbank file for reading.
*/
FilFile::FilFile(const string& filename) : FilterbankFile(filename),
                                           file(filename, ifstream::binary) {
  // Read the headers
  // Note: this code will fail on big-endian systems.
  // If this is not working, you may want to compare it to the code at:
  //   https://github.com/UCBerkeleySETI/blimpy/blob/master/blimpy/io/sigproc.py
  string header_start = readString();
  if (header_start != "HEADER_START") {
    cerr << "The file " << filename << " did not start with HEADER_START. "
         << "Is it really a .fil file?\n";
    exit(1);
  }

  while (true) {
    string attr_name = readString();
    cerr << "attr name: " << attr_name << endl;
    if (attr_name == "telescope_id") {
      telescope_id = readBasic<int>();
    } else if (attr_name == "machine_id") {
      readBasic<int>();
    } else if (attr_name == "data_type") {
      readBasic<int>();
    } else if (attr_name == "barycentric") {
      readBasic<int>();
    } else if (attr_name == "pulsarcentric") {
      readBasic<int>();
    } else if (attr_name == "nbits") {
      readBasic<int>();
    } else if (attr_name == "nsamples") {
      num_timesteps = readBasic<int>();
    } else if (attr_name == "nchans") {
      readBasic<int>();
    } else if (attr_name == "nifs") {
      readBasic<int>();
    } else if (attr_name == "nbeams") {
      readBasic<int>();
    } else if (attr_name == "ibeam") {
      readBasic<int>();
    } else if (attr_name == "rawdatafile") {
      readString();
    } else if (attr_name == "source_name") {
      source_name = readString();
    } else if (attr_name == "az_start") {
      readBasic<double>();
    } else if (attr_name == "za_start") {
      readBasic<double>();
    } else if (attr_name == "tstart") {
      tstart = readBasic<double>();
    } else if (attr_name == "tsamp") {
      tsamp = readBasic<double>();
    } else if (attr_name == "fch1") {
      fch1 = readBasic<double>();
    } else if (attr_name == "foff") {
      foff = readBasic<double>();
    } else if (attr_name == "refdm") {
      readBasic<double>();
    } else if (attr_name == "period") {
      readBasic<double>();
    } else if (attr_name == "src_raj") {
      src_raj = readBasic<double>();
    } else if (attr_name == "src_dej") {
      src_dej = readBasic<double>();
    } else if (attr_name == "HEADER_END") {
      break;
    } else {
      cerr << "unhandled attr name: " << attr_name << endl;
      exit(1);
    }
  }
  
  cerr << "TODO: finish implementing FilFile constructor\n";
  exit(1);
}

template <class T> T FilFile::readBasic() {
  T answer;
  file.read((char*)&answer, sizeof(answer));
  return answer;
}

// Strings are encoded in filterbank headers as first a uint32 containing the string length,
// then the string data
string FilFile::readString() {
  uint32_t num_bytes = readBasic<uint32_t>();
  assert(num_bytes < 256);
  vector<char> buffer(num_bytes);
  file.read(&buffer[0], buffer.size());
  string answer(buffer.begin(), buffer.end());
  return answer;
}

// Loads the data in row-major order.
void FilFile::loadCoarseChannel(int i, float* output) const {
  cerr << "TODO: implement loadCoarseChannel\n";
  exit(1);
}

FilFile::~FilFile() {
  // TODO: clean up file handle
}
