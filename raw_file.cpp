#include "raw_file.h"

#include <fmt/core.h>
#include <iostream>
#include "util.h"

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
    fatal(fmt::format("error reading {} : {}", filename, _reader.errorMessage()));
  }

  if (_headers.empty()) {
    fatal("no headers found in", filename);
  }
}

RawFile::~RawFile() {}

const vector<raw::Header>& RawFile::headers() const {
  return _headers;
}

const raw::Reader& RawFile::reader() const {
  return _reader;
}

