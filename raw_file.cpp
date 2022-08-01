#include "raw_file.h"

using namespace std;

RawFile::RawFile(string filename, int num_bands)
  : _reader(filename), filename(filename), num_bands(num_bands) {
}

// Idempotent
void RawFile::readAllHeaders() {
  if (!_headers.empty()) {
    return;
  }
  
  while (true) {
    raw::Header header;
    if (!_reader.readHeader(&header)) {
      break;
    }
    _headers.push_back(move(header));
  }
}

const vector<raw::Header>& RawFile::headers() {
  readAllHeaders();
  return _headers;
}

const raw::Reader& RawFile::reader() {
  return _reader;
}

