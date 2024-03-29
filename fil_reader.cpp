#include <assert.h>
#include <fstream>
#include <iostream>
#include <math.h>
#include <vector>

#include "fil_reader.h"

using namespace std;

/*
  ra and dec are stored in fil files as:
  HHMMSS.S
  This converts to the more-sane "number of hours with a fraction" format.
 */
double convertFromSigprocRaOrDec(double sigproc) {
  bool negative = sigproc < 0;
  double abs_val = fabs(sigproc);
  int hours_or_degrees = floor(abs_val / 10000.0);
  double remnant = abs_val - (hours_or_degrees * 10000.0);
  int minutes = floor(remnant / 100.0);
  double seconds = remnant - (minutes * 100.0);
  double abs_answer = hours_or_degrees + (minutes / 60.0) + (seconds / 3600.0);
  return negative ? -abs_answer : abs_answer;
}

/*
  Opens a sigproc filterbank file for reading.
  This class is not threadsafe.
*/
FilReader::FilReader(const string& filename) : FilterbankFileReader(filename),
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
      int nbits = readBasic<int>();
      // We only handle 32-bit float data
      assert(nbits == 32);
    } else if (attr_name == "nsamples") {
      num_timesteps = readBasic<int>();
    } else if (attr_name == "nchans") {
      num_channels = readBasic<int>();
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
      src_raj = convertFromSigprocRaOrDec(readBasic<double>());
    } else if (attr_name == "src_dej") {
      src_dej = convertFromSigprocRaOrDec(readBasic<double>());
    } else if (attr_name == "HEADER_END") {
      break;
    } else {
      cerr << "unhandled attr name: " << attr_name << endl;
      exit(1);
    }
  }

  // We have finally figured out where the actual data starts.
  data_start = file.tellg();
  
  // Sometimes the number of timesteps isn't in the header, because the process writing
  // the file didn't know how long it was going to be when it wrote the headers.
  // So figure out the amount of data based on file size.
  file.seekg(0, file.end);
  long num_data_bytes = file.tellg() - data_start;
  if (num_data_bytes % sizeof(float) != 0) {
    cerr << "indivisible amount of data is " << num_data_bytes << " bytes\n";
    exit(1);
  }
  long num_floats = num_data_bytes / sizeof(float);
  if (num_floats % num_channels != 0) {
    cerr << "we have " << num_floats << " which does not divide into " << num_channels
         << " frequencies\n";
    exit(1);
  }
  long inferred_num_timesteps = num_floats / num_channels;
  if (num_timesteps == 0) {
    num_timesteps = inferred_num_timesteps;
  } else if (num_timesteps != inferred_num_timesteps) {
    cerr << "header num timesteps is " << num_timesteps
         << " but inferred num timesteps is " << inferred_num_timesteps << endl;
    exit(1);
  }
  
  inferMetadata();
}

template <class T> T FilReader::readBasic() {
  T answer;
  file.read((char*)&answer, sizeof(answer));
  return answer;
}

// Strings are encoded in filterbank headers as first a uint32 containing the string length,
// then the string data
string FilReader::readString() {
  uint32_t num_bytes = readBasic<uint32_t>();
  assert(num_bytes < 256);
  vector<char> buffer(num_bytes);
  file.read(&buffer[0], buffer.size());
  string answer(buffer.begin(), buffer.end());
  return answer;
}

// Loads the data in row-major order.
void FilReader::loadCoarseChannel(int i, FilterbankBuffer* buffer) const {
  cerr << "TODO: implement FilReader::loadCoarseChannel\n";
  exit(1);
}

FilReader::~FilReader() {}
