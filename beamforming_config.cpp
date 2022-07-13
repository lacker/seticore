#include "beamforming_config.h"

#include <assert.h>
#include <iostream>

#include "beamformer.h"
#include "dedoppler.h"
#include "hit_file_writer.h"
#include "hit_recorder.h"
#include "filterbank_buffer.h"
#include "filterbank_file.h"
#include <fmt/core.h>
#include "multibeam_buffer.h"
#include "raw_file_group.h"
#include "raw/raw.h"
#include "recipe_file.h"
#include "util.h"

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


// Runs the beamforming pipeline from raw files to hits files, based on
// the current config
void BeamformingConfig::run() {
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

    // Write hits to output
    for (DedopplerHit hit : hits) {
      cout << "  index = " << hit.index << ", drift steps = " << hit.drift_steps
           << ", snr = " << hit.snr << ", drift rate = " << hit.drift_rate << endl;
      hit_recorder->recordHit(hit, band, beamformer.power);
    }
    cout << "recorded " << hits.size() << " hits" << endl;
 
  }
}
