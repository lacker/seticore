#pragma once

#include <capnp/message.h>
#include "raw_file_group.h"


using namespace std;

class StampExtractor {
 private:
  RawFileGroup& file_group;
  const int fft_size;
  const int telescope_id;
  int fd;
  
 public:
  StampExtractor(RawFileGroup& file_group, int fft_size, int telescope_id,
                 const string& output_filename);
  ~StampExtractor();

  // No copying
  StampExtractor(const StampExtractor&) = delete;
  StampExtractor& operator=(StampExtractor&) = delete;

  void extract(int coarse_channel, int start_channel, int num_channels);
};
