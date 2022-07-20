#include <boost/algorithm/string/predicate.hpp>
#include <memory>

#include "dat_file_writer.h"
#include "filterbank_file_reader.h"
#include "hit_file_writer.h"

using namespace std;

unique_ptr<HitRecorder> makeHitRecorder(const string& filename,
                                        const FilterbankFileReader& metadata,
                                        double max_drift) {
  if (boost::algorithm::ends_with(filename, ".dat")) {
    return unique_ptr<DatFileWriter>(new DatFileWriter(filename, metadata, max_drift));
  }

  if (boost::algorithm::ends_with(filename, ".hits")) {
    return unique_ptr<HitFileWriter>(new HitFileWriter(filename, metadata));
  }

  cerr << "could not recognize the file type of " << filename << endl;
  exit(1);
}
