#include "raw_file.h"

using namespace std;

RawFile::RawFile(string filename, int num_bands)
  : _reader(filename), filename(filename), num_bands(num_bands) {

  while (true) {
    raw::Header header;
    if (!_reader.readHeader(&header)) {
      break;
    }
    _headers.push_back(move(header));
  }

}

RawFile::~RawFile() {}

const vector<raw::Header>& RawFile::headers() const {
  return _headers;
}

const raw::Reader& RawFile::reader() const {
  return _reader;
}

