#include <boost/algorithm/string/predicate.hpp>
#include <iostream>
#include <memory>

#include "filterbank_file.h"
#include "fil_file.h"
#include "h5_file.h"

using namespace std;

unique_ptr<FilterbankFile> loadFilterbank(const string& filename) {
  if (boost::algorithm::ends_with(filename, ".h5")) {
    return unique_ptr<FilterbankFile>(new H5File(filename));
  }

  if (boost::algorithm::ends_with(filename, ".fil")) {
    return unique_ptr<FilFile>(new FilFile(filename));
  }
  
  cerr << "could not recognize the file type of " << filename << endl;
  exit(1);
}
