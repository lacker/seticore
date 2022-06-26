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
  cout << "nchans: " << header.num_channels << endl;

  double mid_time = header.getMidTime();
  int time_array_index = recipe.getTimeArrayIndex(mid_time);
  int schan = header.getInt("SCHAN", 0);
  vector<float> coefficients(header.num_channels * recipe.nbeams * recipe.npol * recipe.nants * 2);
  
  recipe.generateCoefficients(time_array_index, schan, header.num_channels,
                              header.obsfreq, header.obsbw, coefficients.data());
  
  cout << "mid_time: " << mid_time << endl;
  cout << "tai: " << time_array_index << endl;
  cout << "schan: " << schan << endl;
  cout << "generated " << coefficients.size() << " coefficients\n";
  
  cout << "OK\n";
  return 0;
}
