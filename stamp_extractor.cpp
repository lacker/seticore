#include <capnp/serialize-packed.h>
#include <fmt/core.h>
#include <iostream>
#include "raw_file_group_reader.h"
#include "stamp.capnp.h"
#include "upchannelizer.h"
#include "util.h"

using namespace std;

StampExtractor(RawFileGroup& file_group, int fft_size, int telescope_id)
  : file_group(file_group), fft_size(fft_size), telescope_id(telescope_id) {}

void StampExtractor::extract(int coarse_channel, int start_channel, int num_channels,
                             int fd) {

}
