#include <boost/filesystem.hpp>
#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include <errno.h>
#include "event_file_writer.h"
#include <fmt/core.h>
#include "hit.capnp.h"
#include "hit_file_writer.h"
#include <fcntl.h>
#include <iostream>
#include <unistd.h>
#include "util.h"

using namespace std;


EventFileWriter::EventFileWriter(const string& filename,
                 const vector<shared_ptr<FilterbankFileReader> >& metadatas)
  : metadatas(metadatas), tmp_filename(filename + ".tmp"), final_filename(filename), channel_padding(40)
{
  fd = open(tmp_filename.c_str(), O_WRONLY | O_CREAT, 0664);
  if (fd < 0) {
    int err = errno;
    fatal(fmt::format("could not open {} for writing. errno = {}", tmp_filename, err));
  }
}

EventFileWriter::~EventFileWriter() {
  close(fd);
  boost::filesystem::rename(tmp_filename, final_filename);
}

// hits are null when there is no hit for a region.
// Both of the input vectors should be parallel to the input files.
void EventFileWriter::write(const vector<DedopplerHit*>& hits,
                            const vector<shared_ptr<FilterbankBuffer> >& buffers) {
  assert(hits.size() == metadatas.size());
  assert(buffers.size() == metadatas.size());

  ::capnp::MallocMessageBuilder message;
  Event::Builder event = message.initRoot<Event>();
  auto hits_builder = event.initHits(hits.size());

  assert(hits[0] != NULL);

  int low_index = hits[0]->lowIndex();
  int high_index = hits[0]->highIndex();
  for (int i = 1; i < (int) hits.size(); ++i) {
    if (hits[i] == NULL) {
      continue;
    }
    low_index = min(low_index, hits[i]->lowIndex());
    high_index = max(high_index, hits[i]->highIndex());
  }
  
  low_index -= channel_padding;
  high_index += channel_padding;
  
  int beam = hits[0]->beam;
  int coarse_channel = hits[0]->coarse_channel;

  for (int i = 0; i < (int) hits.size(); ++i) {
    Filterbank::Builder filterbank = hits_builder[i].getFilterbank();
    buildFilterbank(*metadatas[i], beam, coarse_channel, low_index, high_index,
                    buffers[i]->data, filterbank);
    
    if (hits[i] != NULL) {
      buildSignal(*hits[i], hits_builder[i].getSignal());
    }
  }

  writeMessageToFd(fd, message);
}
