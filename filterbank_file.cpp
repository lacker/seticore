#include <boost/algorithm/string/predicate.hpp>
#include <iostream>

#include "filterbank_file.h"
#include "h5_file.h"

using namespace std;

FilterbankFile* loadFilterbank(const string& filename) {
  if (boost::algorithm::ends_with(filename, ".h5")) {
    return new H5File(filename);
  }

  cerr << "could not recognize the file type of " << filename << endl;
  exit(1);
}
