#include <boost/algorithm/string.hpp>
#include <string>
#include <vector>

#include "dat_file.h"
#include "h5_file.h"

using namespace std;

/*
  Opens a dat file for writing, to log hits that we find.
  Some metadata for headers is copied out from the h5 file.
 */
DatFile::DatFile(const string& filename, const H5File& metadata) {
  file.open(filename);

  vector<string> parts;
  boost::split(parts, filename, boost::is_any_of("/"));
  file << "# -------------------------- o --------------------------\n";
  file << "# File ID: " << parts[parts.size() - 1] << " \n";
  file << "# -------------------------- o --------------------------\n";
  file << "# Source:" << metadata.source_name << "\n";
  file << flush;
}


DatFile::~DatFile() {
  file.close();
}


void DatFile::reportHit() {
}
