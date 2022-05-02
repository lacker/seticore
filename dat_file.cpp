#include <string>

#include "dat_file.h"
#include "h5_file.h"

using namespace std;

/*
  Opens a dat file for writing, to log hits that we find.
  Some metadata for headers is copied out from the h5 file.
 */
DatFile::DatFile(const string& filename, const H5File& metadata) {
  file.open(filename);
}


DatFile::~DatFile() {
  file.close();
}
