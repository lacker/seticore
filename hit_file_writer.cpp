#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <errno.h>
#include <fmt/core.h>
#include "hit.capnp.h"
#include "hit_file_writer.h"
#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include "util.h"

using namespace std;

// The number of extra columns on each side of the hit to store
const int EXTRA_COLUMNS = 40;

HitFileWriter::HitFileWriter(const string& filename,
                             const FilterbankMetadata& metadata)
  : metadata(metadata), verbose(true) {
  fd = open(filename.c_str(), O_WRONLY | O_CREAT, 0664);
  if (fd < 0) {
    int err = errno;
    fatal(fmt::format("could not open {} for writing. errno = {}", filename, err));
  }
}

HitFileWriter::~HitFileWriter() {
  close(fd);
}

void HitFileWriter::recordHit(DedopplerHit dedoppler_hit, const float* input) {
  ::capnp::MallocMessageBuilder message;

  Hit::Builder hit = message.initRoot<Hit>();

  dedoppler_hit.buildSignal(hit.getSignal());

  // Most metadata is copied from some input
  Filterbank::Builder filterbank = hit.getFilterbank();
  filterbank.setSourceName(metadata.getBeamSourceName(dedoppler_hit.beam));
  filterbank.setRa(metadata.getBeamRA(dedoppler_hit.beam));
  filterbank.setDec(metadata.getBeamDec(dedoppler_hit.beam));
  filterbank.setTelescopeId(metadata.telescope_id);
  filterbank.setFoff(metadata.foff);
  filterbank.setTsamp(metadata.tsamp);
  filterbank.setTstart(metadata.tstart);
  filterbank.setNumTimesteps(metadata.num_timesteps);
  filterbank.setCoarseChannel(dedoppler_hit.coarse_channel);
  filterbank.setBeam(dedoppler_hit.beam);
  
  // Extract the subset of columns near the hit
  // final_index is the index of the signal at the last point in time we dedopplered for
  // This could be extra if we padded the signal; if that looks weird then fix it
  int final_index = dedoppler_hit.index + dedoppler_hit.drift_steps;
  int leftmost_index, rightmost_index;
  if (final_index < dedoppler_hit.index) {
    // The hit is drifting left
    leftmost_index = final_index;
    rightmost_index = dedoppler_hit.index;
  } else {
    // The hit is drifting right
    leftmost_index = dedoppler_hit.index;
    rightmost_index = final_index;
  }

  // Find the interval [begin_index, end_index) that we actually want to copy over
  // Pad with extra columns but don't run off the edge
  int begin_index = max(leftmost_index - EXTRA_COLUMNS, 0);
  int end_index = min(rightmost_index + EXTRA_COLUMNS, metadata.coarse_channel_size) - 1;
  int num_channels = end_index - begin_index;
  int coarse_offset = dedoppler_hit.coarse_channel * metadata.coarse_channel_size;
  filterbank.setNumChannels(num_channels);
  filterbank.setFch1(metadata.fch1 + (coarse_offset + begin_index) * metadata.foff);
  filterbank.setStartChannel(begin_index);
  filterbank.initData(metadata.num_timesteps * num_channels);
  auto data = filterbank.getData();
  
  // Copy out the subset which would be data[:][begin_index:end_index] in numpy
  for (int tstep = 0; tstep < metadata.num_timesteps; ++tstep) {
    for (int chan = 0; chan < num_channels; ++chan) {
      data.set(tstep * num_channels + chan,
               input[tstep * metadata.coarse_channel_size + begin_index + chan]);
    }
  }

  writeMessageToFd(fd, message);
}
