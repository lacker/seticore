#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <memory>

#include "filterbank_file_reader.h"
#include "fil_reader.h"
#include "h5_reader.h"

using namespace std;

void FilterbankFileReader::loadCoarseChannel(int i, FilterbankBuffer* buffer) const {
  cerr << "cannot loadCoarseChannel from a base FilterbankFileReader: "
       << filename << "\n";
  exit(1);
}

unique_ptr<FilterbankFileReader> loadFilterbankFile(const string& filename) {
  if (boost::algorithm::ends_with(filename, ".h5")) {
    return unique_ptr<FilterbankFileReader>(new H5Reader(filename));
  }

  if (boost::algorithm::ends_with(filename, ".fil")) {
    return unique_ptr<FilterbankFileReader>(new FilReader(filename));
  }
  
  cerr << "could not recognize the file type of " << filename << endl;
  exit(1);
}
