#include "hit_file_writer.h"

using namespace std;

HitFileWriter::HitFileWriter(const string& filename, const FilterbankFile& metadata)
  : metadata(metadata) {
  // TODO: open file
}

HitFileWriter::~HitFileWriter() {
  // TODO: close file
}

void HitFileWriter::recordHit(int coarse_channel, int freq_index, int drift_bins,
                              double drift_rate, double snr) {
  // TODO: make a Hit object, write it to the file
}
