#pragma once

#include <vector>

#include "raw/raw.h"

using namespace std;

/*
  The RawFile reads raw files while caching all the header information.
  It is designed to be faster when reading raw files one band at a time.
 */
class RawFile {
 private:
  vector<raw::Header> _headers;
  raw::Reader _reader;
  
 public:
  const string filename;

  RawFile(string filename);
  ~RawFile();

  RawFile(const RawFile&) = delete;
  RawFile& operator=(RawFile&) = delete;

  const vector<raw::Header>& headers() const;
  const raw::Reader& reader() const;
};
