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
  // Specifying the non-file parameters
  int nbands = 32;
  int fft_size = 131072;
  int telescope_id = 64;
  int nblocks = 128;
  float snr = 7.0;
  float max_drift = 0.01;
  float min_drift = 0.01;

  // Specifying the input files
  vector<string> filenames;
  filenames.push_back(RAW_FILE_0);
  filenames.push_back(RAW_FILE_1);
  RawFileGroup file_group(filenames, nbands);
  RecipeFile recipe(RECIPE_FILE);

  int num_coarse_channels = file_group.num_coarse_channels / nbands;
  int nsamp = file_group.timesteps_per_block * nblocks;

  Beamformer beamformer(fft_size, file_group.nants, recipe.nbeams, nblocks,
                        num_coarse_channels, file_group.npol, nsamp);

  // Create a buffer large enough to hold all beamformer runs for one band  
  int beamformer_runs = file_group.num_blocks / beamformer.nblocks;
  int valid_num_timesteps = beamformer.numOutputTimesteps() * beamformer_runs;
  int rounded_num_timesteps = roundUpToPowerOfTwo(valid_num_timesteps);
  MultibeamBuffer multibeam(beamformer.nbeams,
                            rounded_num_timesteps,
                            beamformer.numOutputChannels());

  if (valid_num_timesteps != rounded_num_timesteps) {
    // Zero out the buffer, because it has some extra area
    multibeam.zeroAsync();
  }
  
  // Hardcoded for now
  int band = 0;
  
  FilterbankFile metadata = combineMetadata(file_group, beamformer, band,
                                            telescope_id);
  file_group.resetBand(band);
  cout << endl;
  
  for (int run = 0; run < beamformer_runs; ++run) {
    if (run > 0) {
      // We need to sync so that we don't overwrite on reads
      cudaDeviceSynchronize();
    }

    for (int block = 0; block < beamformer.nblocks; ++block) {
      file_group.read(beamformer.inputPointer(block));
    } 
  
    cout << "beamforming band " << band << " run " << run << "...\n";
    
    int block_right_after_mid = run * beamformer.nblocks + beamformer.nblocks / 2;
    double mid_time = file_group.getStartTime(block_right_after_mid);
    int time_array_index = recipe.getTimeArrayIndex(mid_time);
    recipe.generateCoefficients(time_array_index, file_group.schan,
                                beamformer.num_coarse_channels,
                                file_group.obsfreq, file_group.obsbw,
                                beamformer.coefficients);
    int time_offset = beamformer.numOutputTimesteps() * run;
    beamformer.processInput(multibeam, time_offset);
  }
  
  // Spot check the beamformed power
  float power = multibeam.getFloat(0, 0, 0);
  cout << "power[0]: " << power << endl;
  assertFloatEq(power / 1.0e14, 1.87708);
  
  // Search beam zero
  cout << "\ndedoppler searching " << metadata.num_freqs << " channels, "
       << metadata.num_timesteps << " timesteps\n";
  
  FilterbankBuffer buffer = multibeam.getBeam(0);
  Dedopplerer dedopplerer(valid_num_timesteps,
                          beamformer.numOutputChannels(),
                          metadata.foff, metadata.tsamp, false);
  vector<DedopplerHit> hits;
  dedopplerer.search(buffer, max_drift, min_drift, snr, &hits);

  // Spot check the hits
  assert(3 == hits.size());
  assert(61685 == hits[0].index);
  assert(-1 == hits[1].drift_steps);
  assertFloatEq(7.22453, hits[2].snr);
  assertFloatEq(1.08951, hits[0].drift_rate);
  
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
