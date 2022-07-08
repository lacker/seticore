#include <assert.h>
#include <iostream>

#include "beamformer.h"
#include "dedoppler.h"
#include "hit_file_writer.h"
#include "filterbank_buffer.h"
#include "filterbank_file.h"
#include "multibeam_buffer.h"
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
const string& OUTPUT_HITS = "./data/beamforming.hits";

int main(int argc, char* argv[]) {
  RecipeFile recipe(RECIPE_FILE);
  raw::Reader reader(RAW_FILE_0);
  raw::Header header;

  if (!reader.readHeader(&header)) {
    cout << "raw error: " << reader.errorMessage() << endl;
    return 1;
  }

  cout << "\nreading raw file: " << RAW_FILE_0 << endl;
  cout << "gathering metadata from first block:\n";
  cout << "nants: " << header.nants << endl;
  cout << "nchans: " << header.num_channels << endl;
  cout << "npol: " << header.npol << endl;
  cout << "num_timesteps: " << header.num_timesteps << endl;
  cout << "tbin: " << header.tbin << endl;

  // Can be passed as arguments
  int fft_size = 131072;
  int telescope_id = 64;

  // This we need to detect from existing files
  int nblocks = 128;

  // We need "planner" type functionality to figure out how thin to slice a band
  int nbands = 16;
  
  int timesteps_per_block = header.num_timesteps;
  int nants = header.nants;
  int nbeams = recipe.nbeams;
  int num_coarse_channels = header.num_channels / nbands;
  int npol = header.npol;
  int nsamp = timesteps_per_block * nblocks;

  assert(0 == header.num_channels % nbands);
  
  // double obsfreq = header.obsfreq;
  double obsbw = header.obsbw;
  double tbin = header.tbin;
  
  Beamformer beamformer(fft_size, nants, nbeams, nblocks,
                        num_coarse_channels, npol, nsamp);

  // Calculate metadata for output
  FilterbankFile metadata("");
  metadata.has_dc_spike = false;
  metadata.source_name = string(header.src_name);
  double output_bandwidth = obsbw / nbands;
  metadata.fch1 = header.obsfreq - 0.5 * obsbw;
  metadata.foff = output_bandwidth / beamformer.numOutputChannels();
  metadata.tsamp = tbin * nblocks * timesteps_per_block / beamformer.numOutputTimesteps();
  metadata.num_timesteps = beamformer.numOutputTimesteps();
  metadata.num_freqs = beamformer.numOutputChannels();
  metadata.coarse_channel_size = metadata.num_freqs;
  metadata.telescope_id = telescope_id;

  // Right now the dedoppler treats everything in the subband as a single
  // coarse channel. That is probably not the right way to do it, but it
  // only affects the SNR normalization logic.
  metadata.num_coarse_channels = 1;

  cout << "dedopplering: " << metadata.num_freqs << " channels, "
       << metadata.num_timesteps << " timesteps\n";

  cout << "foff: " << metadata.foff << endl;
  cout << "tsamp: " << metadata.tsamp << endl;
  
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

  MultibeamBuffer multibeam(beamformer.nbeams, beamformer.numOutputTimesteps(),
                            beamformer.numOutputChannels());
  
  cout << "\nbeamforming...\n";
  beamformer.processInput(multibeam, 0);
  cout << "spot checking data.\n";
  
  auto val = beamformer.getPrebeam(0, 0, 0, 0);
  assertComplexEq(val, 483.0, -3011.0);
  cout << "prebeam[0]: " << cToS(val) << endl;

  val = beamformer.getVoltage(0, 0, 0, 0);
  assertComplexEq(val, -1938482.875, 8989387.0);
  cout << "voltage[0]: " << cToS(val) << endl;

  float power = multibeam.getFloat(0, 0, 0);
  cout << "power[0]: " << power << endl;

  // Beam zero
  cout << "\nstarting dedoppler search.\n";
  cout << "searching beam 0...\n";
  FilterbankBuffer buffer = multibeam.getBeam(0);
  Dedopplerer dedopplerer(beamformer.numOutputTimesteps(), beamformer.numOutputChannels(),
                          metadata.foff, metadata.tsamp, false);
  vector<DedopplerHit> hits;
  float snr = 8.0;
  cout << "SNR threshold: " << snr << endl;
  dedopplerer.search(buffer, 0.01, 0.01, snr, &hits);

  // Write hits to output
  auto recorder = HitFileWriter(OUTPUT_HITS, metadata);
  for (DedopplerHit hit : hits) {
    cout << "  index = " << hit.index << ", drift steps = " << hit.drift_steps
         << ", snr = " << hit.snr << ", drift rate = " << hit.drift_rate << endl;
    recorder.recordHit(0, hit.index, hit.drift_steps, hit.drift_rate, hit.snr, beamformer.power);
  }
  cout << "wrote " << hits.size() << " hits to " << OUTPUT_HITS << endl;  
  return 0;
}
