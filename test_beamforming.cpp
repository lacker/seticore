#include <assert.h>
#include <iostream>

#include "beamformer.h"
#include "raw/raw.h"
#include "recipe_file.h"

using namespace std;

void assertComplexEq(thrust::complex<float> c, float real, float imag) {
  if (abs(c.real() - real) > 0.001) {
    cerr << "c.real() = " << c.real() << " but real = " << real << endl;
    exit(1);
  }
  if (abs(c.imag() - imag) > 0.001) {
    cerr << "c.imag() = " << c.imag() << " but imag = " << imag << endl;
    exit(1);
  }
}

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

  cout << "reading " << header.blocsize << " bytes\n";
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

  assertComplexEq(beamformer.getCoefficient(0, 0, 0, 0), 20.0037, -94.5755);
  assertComplexEq(beamformer.getCoefficient(1, 1, 1, 1), -41.4678, -83.6412);
  assertComplexEq(beamformer.getCoefficient(9, 0, 9, 9), -91.5771, 1.49745);
  assertComplexEq(beamformer.getCoefficient(32, 0, 32, 256), -35.1951, -102.186);
  assertComplexEq(beamformer.getCoefficient(48, 1, 48, 384), 97.9193, -3.7384);

  beamformer.beamform();

  assertComplexEq(beamformer.getTransposed(9, 7, 1, 3), 8.0, 7.0);
  
  cout << "OK\n";
  return 0;
}
