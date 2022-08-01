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
  
  void readAllHeaders();
  
 public:
  const string filename;
  const int num_bands;

  RawFile(string filename, int num_bands);
  ~RawFile();

  RawFile(const RawFile&) = delete;
  RawFile& operator=(RawFile&) = delete;
  
  const vector<raw::Header>& headers();
  const raw::Reader& reader();
};
