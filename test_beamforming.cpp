#include <assert.h>
#include <iostream>

#include "raw/raw.h"
#include "recipe_file.h"

using namespace std;

int main(int argc, char* argv[]) {
  RecipeFile recipe("/d/minput/guppi_59599_76288_000062_J03323-28075_0001.bfr5");
  cout << "from recipe file:\n";
  cout << "nants: " << recipe.nants << endl;
  cout << "nbeams: " << recipe.nbeams << endl;
  cout << "nchans: " << recipe.nchans << endl;
  cout << "ndelays: " << recipe.delays.size() << endl;

  // Read one block from the raw file
  raw::Reader reader("/d/minput/guppi_59599_76288_000062_J03323-28075_0001.0000.raw");
  raw::Header header;
  if (!reader.readHeader(&header)) {
    cout << "raw error: " << reader.errorMessage() << endl;
    return 1;
  }
  vector<char> data(header.blocsize);
  reader.readData(data.data());
  
  cout << "\nfrom raw file:\n";
  cout << "nants: " << header.nants << endl;
  cout << "nbeam: " << header.nbeam << endl;
  cout << "nchans: " << header.num_channels << endl;

  cout << "OK\n";
  return 0;
}
