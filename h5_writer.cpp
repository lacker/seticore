#include <iostream>

#include "h5_writer.h"

/*
  data must be formatted as row-major:
    data[time][freq]
 */
void writeH5(const string& filename,
             const FilterbankMetadata& metadata,
             const float* data) {
  cerr << "TODO: implement writeH5\n";
  exit(1);
}
