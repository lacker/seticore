#include <assert.h>
#include <iostream>

#include "beamformer.h"
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

  Beamformer beamformer(header.nants, recipe.nbeams, header.num_channels,
                        recipe.npol, header.num_timesteps);

  reader.readData((char*) beamformer.input);
  
  cout << "\nfrom raw file:\n";
  cout << "nants: " << header.nants << endl;
  cout << "nchans: " << header.num_channels << endl;

  double mid_time = header.getMidTime();
  int time_array_index = recipe.getTimeArrayIndex(mid_time);
  int schan = header.getInt("SCHAN", 0);
  recipe.generateCoefficients(time_array_index, schan, header.num_channels,
                              header.obsfreq, header.obsbw, beamformer.coefficients);
  
  cout << "mid_time: " << mid_time << endl;
  cout << "tai: " << time_array_index << endl;
  cout << "schan: " << schan << endl;

  beamformer.debugCoefficients(0, 0, 0, 0);
  beamformer.debugCoefficients(1, 1, 1, 1);
  beamformer.debugCoefficients(9, 0, 9, 9);
  beamformer.debugCoefficients(32, 0, 32, 256);
  beamformer.debugCoefficients(48, 1, 48, 384);
    
  cout << "OK\n";
  return 0;
}
