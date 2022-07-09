#include <assert.h>
#include <iostream>

#include "beamformer.h"
#include "dedoppler.h"
#include "hit_file_writer.h"
#include "filterbank_buffer.h"
#include "filterbank_file.h"
#include "multibeam_buffer.h"
#include "raw/raw.h"
#include "raw_file_group.h"
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


// Just fill up the beamformer once.
void readRawBand(RawFileGroup& file_group,
                 Beamformer& beamformer,
                 const RecipeFile& recipe,
                 int band) {
  file_group.resetBand(band);
  assert(beamformer.nblocks <= file_group.num_blocks);
  
  for (int block = 0; block < beamformer.nblocks; ++block) {
    file_group.read(beamformer.inputPointer(block));
    if ((block + 1) % 32 == 0) {
      cout << "read " << block + 1 << " blocks, band " << band << endl;
    }
  }

  double mid_time = file_group.getStartTime(beamformer.nblocks / 2);
  int time_array_index = recipe.getTimeArrayIndex(mid_time);
  recipe.generateCoefficients(time_array_index, file_group.schan,
                              beamformer.num_coarse_channels,
                              file_group.obsfreq, file_group.obsbw,
                              beamformer.coefficients);
}


// Construct metadata for the data that would be created by a given
// file group and beamformer, running the beamform as many times as
// it can be filled by the file group.
FilterbankFile combineMetadata(const RawFileGroup& file_group,
                               const Beamformer& beamformer,
                               int band, int telescope_id) {
  FilterbankFile metadata("");
  metadata.has_dc_spike = false;
  metadata.source_name = file_group.source_name;
  double output_bandwidth = file_group.obsbw / file_group.num_bands;
  double file_group_fch1 = file_group.obsfreq - 0.5 * file_group.obsbw;
  metadata.fch1 = band * output_bandwidth + file_group_fch1;
  metadata.foff = output_bandwidth / beamformer.numOutputChannels();
  int beamformer_runs = file_group.num_blocks / beamformer.nblocks;
  double time_per_block = file_group.tbin * file_group.timesteps_per_block;
  double time_per_beamform = time_per_block * beamformer.nblocks;
  metadata.tsamp = time_per_beamform / beamformer.numOutputTimesteps();
  metadata.num_timesteps = beamformer.numOutputTimesteps() * beamformer_runs;
  metadata.num_freqs = beamformer.numOutputChannels();
  metadata.telescope_id = telescope_id;

  // Currently, we treat every band as a single coarse channel for the
  // purposes of seti search. This might not be the right way to do it.
  // TODO: figure out if we should split up seti searching
  metadata.coarse_channel_size = metadata.num_freqs;
  metadata.num_coarse_channels = 1;

  return metadata;
}


int main(int argc, char* argv[]) {
  // TODO: How should we figure out how thin to slice a band?
  int nbands = 32;

  vector<string> filenames;
  filenames.push_back(RAW_FILE_0);
  RawFileGroup file_group(filenames, nbands);

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

  // TODO: infer this
  int nblocks = 128;

  int timesteps_per_block = header.num_timesteps;
  int nants = header.nants;
  int nbeams = recipe.nbeams;
  int num_coarse_channels = header.num_channels / nbands;
  int npol = header.npol;
  int nsamp = timesteps_per_block * nblocks;

  Beamformer beamformer(fft_size, nants, nbeams, nblocks,
                        num_coarse_channels, npol, nsamp);

  FilterbankFile metadata = combineMetadata(file_group, beamformer, 0, telescope_id);
  
  readRawBand(file_group, beamformer, recipe, 0);
  
  // Create a buffer that, for now, is large enough to hold one beamformer output
  MultibeamBuffer multibeam(beamformer.nbeams, beamformer.numOutputTimesteps(),
                            beamformer.numOutputChannels());
  
  cout << "\nbeamforming...\n";
  beamformer.processInput(multibeam, 0);

  // Spot check the beamformed power
  float power = multibeam.getFloat(0, 0, 0);
  cout << "power[0]: " << power << endl;
  assertFloatEq(power / 1.0e14, 1.87708);
  
  // Search beam zero
  cout << "\ndedoppler searching " << metadata.num_freqs << " channels, "
       << metadata.num_timesteps << " timesteps\n";
  
  FilterbankBuffer buffer = multibeam.getBeam(0);
  Dedopplerer dedopplerer(beamformer.numOutputTimesteps(),
                          beamformer.numOutputChannels(),
                          metadata.foff, metadata.tsamp, false);
  vector<DedopplerHit> hits;
  float snr = 8.0;
  cout << "SNR threshold: " << snr << endl;
  dedopplerer.search(buffer, 0.01, 0.01, snr, &hits);

  // Spot check the hits
  assert(69884 == hits[0].index);
  assert(-2 == hits[1].drift_steps);
  assertFloatEq(8.20315, hits[0].snr);
  assertFloatEq(-0.847396, hits[1].drift_rate);
  
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
