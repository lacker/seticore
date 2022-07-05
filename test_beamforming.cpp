#include <assert.h>
#include <iostream>

#include "beamformer.h"
#include "raw/raw.h"
#include "recipe_file.h"
#include "util.h"

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

const string& RAW_FILE_0 = "/d/mubf/guppi_59712_16307_003760_J1939-6342_0001.0000.raw";
const string& RAW_FILE_1 = "/d/mubf/guppi_59712_16307_003760_J1939-6342_0001.0001.raw";
const string& RECIPE_FILE = "/d/mubf/MeerKAT-array_1-20220513T043147Z.bfr5";

int main(int argc, char* argv[]) {
  RecipeFile recipe(RECIPE_FILE);
  raw::Reader reader(RAW_FILE_0);
  raw::Header header;

  if (!reader.readHeader(&header)) {
    cout << "raw error: " << reader.errorMessage() << endl;
    return 1;
  }

  cout << "\none block from " << RAW_FILE_0 << endl;
  cout << "nants: " << header.nants << endl;
  cout << "nchans: " << header.num_channels << endl;
  cout << "npol: " << header.npol << endl;
  cout << "num_timesteps: " << header.num_timesteps << endl;

  int timesteps_per_block = 8192;
  
  int fft_size = 131072;
  int nants = 61;
  int nbands = 16;
  int nbeams = 64;
  int nblocks = 128;
  int num_coarse_channels = 4;
  int npol = 2;
  int nsamp = timesteps_per_block * nblocks;
  
  Beamformer beamformer(fft_size, nants, nbeams, nblocks,
                        num_coarse_channels, npol, nsamp);

  int block = 0;
  int band = 0;
  while (block < nblocks) {
    reader.readBand(header, band, nbands, beamformer.inputPointer(block));
    ++block;

    if (block % 8 == 0) {
      cout << "read " << block << " blocks, band " << band << endl;
    }
    
    if (!reader.readHeader(&header)) {
      break;
    }

    if (block == nblocks / 2) {
      // Start time for this block is mid time for the beamformer
      double mid_time = header.getStartTime();
      int time_array_index = recipe.getTimeArrayIndex(mid_time);
      int schan = header.getInt("SCHAN", 0);
      recipe.generateCoefficients(time_array_index, schan, beamformer.num_coarse_channels,
                                  header.obsfreq, header.obsbw, beamformer.coefficients);
    }
  }

  beamformer.processInput();

  auto val = beamformer.getPrebeam(0, 0, 0, 0);
  assertComplexEq(val, 483.0, -3011.0);
  cout << "prebeam[0]: " << cToS(val) << endl;

  val = beamformer.getVoltage(0, 0, 0, 0);
  assertComplexEq(val, -1938482.875, 8989387.0);
  cout << "voltage[0]: " << cToS(val) << endl;

  float power = beamformer.getPower(0, 0, 0);
  cout << "power[0]: " << power << endl;
  
  return 0;
}
