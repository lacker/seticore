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

void assertFloatEq(float a, float b) {
  while (a > 100) {
    a /= 2.0;
    b /= 2.0;
  }
  if (abs(a - b) > 0.001) {
    cerr << a << " != " << b << endl;
    exit(1);
  }
}

const string& RAW_FILE = "/d/minput/guppi_59599_76288_000062_J03323-28075_0001.0000.raw";
const string& RECIPE_FILE = "/d/minput/guppi_59599_76288_000062_J03323-28075_0001.bfr5";

int testWithoutChannelizing() {
  RecipeFile recipe(RECIPE_FILE);
  cout << "from recipe file:\n";
  cout << "nants: " << recipe.nants << endl;
  cout << "nbeams: " << recipe.nbeams << endl;
  cout << "nchans: " << recipe.nchans << endl;
  cout << "ndelays: " << recipe.delays.size() << endl;

  // Read one block from the raw file
  raw::Reader reader(RAW_FILE);
  raw::Header header;
  if (!reader.readHeader(&header)) {
    cout << "raw error: " << reader.errorMessage() << endl;
    return 1;
  }

  Beamformer beamformer(1, header.nants, recipe.nbeams, 1, header.num_channels,
                        recipe.npol, header.num_timesteps);

  cout << "reading " << header.blocsize << " bytes\n";
  reader.readData(beamformer.inputPointer(0));
  
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

  beamformer.processInput();

  assertComplexEq(beamformer.getPrebeam(9, 7, 1, 3), 8.0, 7.0);

  assertComplexEq(beamformer.getVoltage(0, 0, 0, 0), -4715.979492, 2549.703125);
  assertComplexEq(beamformer.getVoltage(4, 3, 2, 1), 562.406616, -10480.619141);

  assertFloatEq(beamformer.getPower(9, 8, 7), 820444800.0);
  
  cout << "OK\n";
  return 0;
}

int main(int argc, char* argv[]) {
  RecipeFile recipe(RECIPE_FILE);
  raw::Reader reader(RAW_FILE);
  raw::Header header;

  if (!reader.readHeader(&header)) {
    cout << "raw error: " << reader.errorMessage() << endl;
    return 1;
  }
  
  int nbands = 4;
  int fft_size = 1024;
  int nblocks = 16;
  int nsamp = header.num_timesteps * nblocks;
  
  Beamformer beamformer(fft_size, header.nants, recipe.nbeams, nblocks,
                        header.num_channels / nbands, recipe.npol, nsamp);

  cout << "\none block from the raw file:\n";
  cout << "nants: " << header.nants << endl;
  cout << "nchans: " << header.num_channels << endl;
  cout << "npol: " << header.npol << endl;
  cout << "num_timesteps: " << header.num_timesteps << endl;

  int block = 0;
  while (true) {
    // TODO: shift this by block
    reader.readBand(header, 2, nbands, beamformer.inputPointer(block));
    ++block;
    cout << "read block " << block << endl;
    if (block == 16) {
      break;
    }
    
    assert(reader.readHeader(&header));
  }

  cout << "OK\n";
  return 0;
}
