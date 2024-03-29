#pragma once

#include <capnp/message.h>
#include "dedoppler_hit.h"
#include "raw_file_group.h"


using namespace std;

class StampExtractor {
 private:
  RawFileGroup& file_group;
  const int fft_size;
  const int telescope_id;
  const string tmp_filename;
  const string final_filename;
  int fd;
  bool opened;

  // Open the output file.
  // No-op if it's already open.
  // We do this lazily so that extracting doesn't create a file if extract()
  // is never called.
  void openOutputFile();
  
 public:
  StampExtractor(RawFileGroup& file_group, int fft_size, int telescope_id,
                 const string& output_filename);
  ~StampExtractor();

  // No copying
  StampExtractor(const StampExtractor&) = delete;
  StampExtractor& operator=(StampExtractor&) = delete;

  // hit may be nullptr, in which case no hit information is added.
  void extract(const DedopplerHit* hit, int coarse_channel,
               int start_channel, int num_channels);
};
