#include "raw_file.h"

#include <iostream>

using namespace std;

RawFile::RawFile(string filename)
  : _reader(filename), filename(filename) {

  while (true) {
    raw::Header header;
    if (!_reader.readHeader(&header)) {
      break;
    }
    _headers.push_back(move(header));
  }

  if (_reader.error()) {
    cerr << "error reading " << filename << endl;
    cerr << _reader.errorMessage() << endl;
    exit(1);
  }

  if (_headers.empty()) {
    cerr << "no headers found in " << filename << endl;
    exit(1);
  }
}

RawFile::~RawFile() {}

const vector<raw::Header>& RawFile::headers() const {
  return _headers;
}

const raw::Reader& RawFile::reader() const {
  return _reader;
}

