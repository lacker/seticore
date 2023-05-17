#pragma once

#include <string>
#include <vector>

#include "dedoppler_hit.h"
#include "filterbank_buffer.h"
#include "filterbank_file_reader.h"

using namespace std;

class EventFileWriter {
 private:
  const vector<shared_ptr<FilterbankFileReader> > metadatas;
  int fd;

  const string tmp_filename;
  const string final_filename;
  
 public:
  // The number of extra columns on each side of the filterbank-data to store
  int channel_padding;

  EventFileWriter(const string& filename,
                  const vector<shared_ptr<FilterbankFileReader> >& metadatas);
  ~EventFileWriter();

  void write(const vector<DedopplerHit*>& hits,
             const vector<shared_ptr<FilterbankBuffer> >& buffers);
};
