#include <capnp/message.h>
#include <capnp/serialize-packed.h>
#include "hit.capnp.h"
#include "hit_file_writer.h"
#include <fcntl.h>
#include <unistd.h>

using namespace std;

HitFileWriter::HitFileWriter(const string& filename, const FilterbankFile& metadata)
  : metadata(metadata) {
  fd = open(filename.c_str(), O_WRONLY | O_CREAT);
}

HitFileWriter::~HitFileWriter() {
  close(fd);
}

void HitFileWriter::recordHit(int coarse_channel, int freq_index, int drift_bins,
                              double drift_rate, double snr) {

  ::capnp::MallocMessageBuilder message;

  Hit::Builder hit = message.initRoot<Hit>();

}
