#pragma once

#include <string>
#include <vector>

#include "dedoppler_hit.h"

using namespace std;

class EventFileWriter {
 private:
  const vector<FilterbankMetadata> metadatas;
  int fd;

  const string tmp_filename;
  const string final_filename;
  
 public:
  EventFileWriter(const string& filename, const vector<FilterbankMetadata>& metadatas);
  ~EventFileWriter();

  void write(const vector<DedopplerHit*>& hits, const vector<float*> inputs);
};
