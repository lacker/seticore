#include <assert.h>
#include <iostream>

#include "beamformer.h"
#include "dedoppler.h"
#include "hit_file_writer.h"
#include "filterbank_buffer.h"
#include "filterbank_file.h"
#include <fmt/core.h>
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

// Construct metadata for the data that would be created by a given
// file group and beamformer, running the beamform as many times as
// it can be filled by the file group.
FilterbankFile combineMetadata(const RawFileGroup& file_group,
                               const Beamformer& beamformer,
                               int telescope_id) {
  FilterbankFile metadata("");
  metadata.has_dc_spike = false;
  metadata.source_name = file_group.source_name;
  metadata.fch1 = file_group.obsfreq - 0.5 * file_group.obsbw;
  double output_bandwidth = file_group.obsbw / file_group.num_bands;
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
  metadata.num_coarse_channels = file_group.num_bands;

  return metadata;
}

class BeamformingConfig {
public:
  // Required parameters for the run
  const vector<string> raw_files;
  const string output_dir;
  const string recipe_dir;
  const int num_bands;
  const int fft_size;
  const int telescope_id;
  const float snr;
  const float max_drift;
  const float min_drift;

  // Optional additional ways to modify the config

  // If set, we only process a subset of the bands
  int num_bands_to_process;

  // If set, use this to record hits instead of writing a file
  unique_ptr<HitRecorder> hit_recorder;
  
  BeamformingConfig(const vector<string>& raw_files,
                    const string& output_dir,
                    const string& recipe_dir,
                    int num_bands,
                    int fft_size,
                    int telescope_id,
                    float snr,
                    float max_drift,
                    float min_drift)
    : raw_files(raw_files), output_dir(stripAnyTrailingSlash(output_dir)),
      recipe_dir(recipe_dir), num_bands(num_bands), fft_size(fft_size),
      telescope_id(telescope_id), snr(snr), max_drift(max_drift),
      min_drift(min_drift), num_bands_to_process(num_bands) {
  }

  // Runs the beamforming pipeline from raw files to hits files, based on
  // the current config
  void run() {
    RawFileGroup file_group(raw_files, num_bands);
    RecipeFile recipe(recipe_dir, file_group.obsid);

    // Do enough blocks per beamformer run to handle one STI block
    assert((STI * fft_size) % file_group.timesteps_per_block == 0);
    int blocks_per_run = (STI * fft_size) / file_group.timesteps_per_block;
  
    int num_coarse_channels = file_group.num_coarse_channels / num_bands;
    int nsamp = file_group.timesteps_per_block * blocks_per_run;

    Beamformer beamformer(fft_size, file_group.nants, recipe.nbeams, blocks_per_run,
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

    FilterbankFile metadata = combineMetadata(file_group, beamformer, telescope_id);

    if (hit_recorder == NULL) {
      string output_filename = fmt::format("{}/{}.hits", output_dir,
                                           file_group.prefix);
      cout << "recording hits to " << output_filename << endl;
      hit_recorder.reset(new HitFileWriter(output_filename, metadata));
    }

    Dedopplerer dedopplerer(valid_num_timesteps,
                            beamformer.numOutputChannels(),
                            metadata.foff, metadata.tsamp, false);
  
    for (int band = 0; band < num_bands_to_process; ++band) {
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
  
      // Search beam zero
      cout << "\nband " << band << ": dedoppler searching " << metadata.num_freqs
           << " channels, " << metadata.num_timesteps << " timesteps\n";
  
      FilterbankBuffer buffer = multibeam.getBeam(0);
      vector<DedopplerHit> hits;
      dedopplerer.search(buffer, max_drift, min_drift, snr, &hits);

      if (band == 0) {
        // Spot check the beamformed power
        float power = multibeam.getFloat(0, 0, 0);
        cout << "power[0]: " << power << endl;
        assertFloatEq(power / 1.0e14, 4.82461);
  
        // Spot check the hits
        assert(1 == hits.size());
        assert(114049 == hits[0].index);
        assert(-1 == hits[0].drift_steps);
        assertFloatEq(7.05269, hits[0].snr);
        assertFloatEq(-0.317774, hits[0].drift_rate);
      }
  
      // Write hits to output
      for (DedopplerHit hit : hits) {
        cout << "  index = " << hit.index << ", drift steps = " << hit.drift_steps
             << ", snr = " << hit.snr << ", drift rate = " << hit.drift_rate << endl;
        hit_recorder->recordHit(band, hit.index, hit.drift_steps, hit.drift_rate,
                                hit.snr, beamformer.power);
      }
      cout << "recorded " << hits.size() << " hits" << endl;
 
    }
  }
};


int main(int argc, char* argv[]) {
  // Specifying parameters
  string input_dir = "/d/onlytwo";
  string output_dir = "./data";
  string recipe_dir = "/d/mubf";
  int num_bands = 32;
  int fft_size = 131072;
  int telescope_id = 64;
  float snr = 7.0;
  float max_drift = 0.01;
  float min_drift = 0.01;

  auto groups = scanForRawFileGroups(input_dir);
  assert(groups.size() == 1);
  assert(groups[0].size() == 2);

  BeamformingConfig config(groups[0], output_dir, recipe_dir, num_bands,
                           fft_size, telescope_id, snr, max_drift, min_drift);
  config.num_bands_to_process = 1;
  config.run();
  
  return 0;
}
