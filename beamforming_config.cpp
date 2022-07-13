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
  int beamformer_batches = file_group.num_blocks / beamformer.nblocks;
  double time_per_block = file_group.tbin * file_group.timesteps_per_block;
  double time_per_beamform = time_per_block * beamformer.nblocks;
  metadata.tsamp = time_per_beamform / beamformer.numOutputTimesteps();
  metadata.num_timesteps = beamformer.numOutputTimesteps() * beamformer_batches;
  metadata.num_freqs = beamformer.numOutputChannels();
  metadata.telescope_id = telescope_id;

  // Currently, we treat every band as a single coarse channel for the
  // purposes of seti search. This might not be the right way to do it.
  // TODO: figure out if we should split up seti searching
  metadata.coarse_channel_size = metadata.num_freqs;
  metadata.num_coarse_channels = file_group.num_bands;

  return metadata;
}

/*
  Runs the beamforming pipeline from raw files to hits files, based on
  the current config.

  Roughly, the input is a group of raw files. The data in them is logically
  the same sort of data, it's just broken up among a number of raw files to
  keep any one from getting too large.

  The beamformer has a limited size capacity, so we don't beamform all the
  data at once. We break it up in two ways. First, we divide the input
  data by frequency, into bands of equal width. We process bands one at a time.

  Second, within a band, we do multiple batches of beamforming. In each batch, we
  load as many timesteps as we can into the beamformer, then beamform, and
  append the output to a multibeam buffer. When we finish beamforming all timesteps
  within a band, we run a dedoppler search on the power values in the accumulated
  buffer.
*/
void BeamformingConfig::run() {
  RawFileGroup file_group(raw_files, num_bands);
  RecipeFile recipe(recipe_dir, file_group.obsid);

  // Do enough blocks per beamformer batch to handle one STI block
  assert((STI * fft_size) % file_group.timesteps_per_block == 0);
  int blocks_per_batch = (STI * fft_size) / file_group.timesteps_per_block;
  
  int num_coarse_channels = file_group.num_coarse_channels / num_bands;
  int nsamp = file_group.timesteps_per_block * blocks_per_batch;

  Beamformer beamformer(fft_size, file_group.nants, recipe.nbeams, blocks_per_batch,
                        num_coarse_channels, file_group.npol, nsamp);

  // Create a buffer large enough to hold all beamformer batches for one band  
  int beamformer_batches = file_group.num_blocks / beamformer.nblocks;
  int num_multibeam_timesteps = beamformer.numOutputTimesteps() * beamformer_batches;
  MultibeamBuffer multibeam(beamformer.nbeams,
                            num_multibeam_timesteps,
                            beamformer.numOutputChannels());

  // Create a buffer for dedoppler input, padding with zeros
  FilterbankBuffer buffer(roundUpToPowerOfTwo(multibeam.num_timesteps),
                          multibeam.num_channels);
  buffer.zero();
  
  FilterbankFile metadata = combineMetadata(file_group, beamformer, telescope_id);

  if (hit_recorder == NULL) {
    string output_filename = fmt::format("{}/{}.hits", output_dir,
                                         file_group.prefix);
    cout << "recording hits to " << output_filename << endl;
    hit_recorder.reset(new HitFileWriter(output_filename, metadata));
  }

  Dedopplerer dedopplerer(multibeam.num_timesteps,
                          buffer.num_channels,
                          metadata.foff, metadata.tsamp, false);
  
  for (int band = 0; band < num_bands_to_process; ++band) {
    file_group.resetBand(band);
    cout << endl;
  
    for (int batch = 0; batch < beamformer_batches; ++batch) {
      if (batch > 0) {
        // We need to sync so that we don't overwrite on reads
        cudaDeviceSynchronize();
      }

      for (int block = 0; block < beamformer.nblocks; ++block) {
        file_group.read(beamformer.inputPointer(block));
      } 
  
      cout << "beamforming band " << band << ", batch " << batch << "...\n";
    
      int block_right_after_mid = batch * beamformer.nblocks + beamformer.nblocks / 2;
      double mid_time = file_group.getStartTime(block_right_after_mid);
      int time_array_index = recipe.getTimeArrayIndex(mid_time);
      recipe.generateCoefficients(time_array_index, file_group.schan,
                                  beamformer.num_coarse_channels,
                                  file_group.obsfreq, file_group.obsbw,
                                  beamformer.coefficients);
      int time_offset = beamformer.numOutputTimesteps() * batch;
      beamformer.processInput(multibeam, time_offset);
    }

    cout << "band " << band << ": dedoppler searching "
         << beamformer.nbeams << " beams, " << metadata.num_freqs << " channels, "
         << metadata.num_timesteps << " timesteps\n";
    
    int total_hits = 0;
    for (int beam = 0; beam < beamformer.nbeams; ++beam) {
      // TODO: this is setting coarse channel to zero
      multibeam.copyRegionAsync(beam, 0, &buffer);
      vector<DedopplerHit> hits;
      dedopplerer.search(buffer, max_drift, min_drift, snr, &hits);
      
      // Write hits to output
      if (hits.empty()) {
        continue;
      }
      total_hits += hits.size();
      cout << "found " << pluralize(hits.size(), "hit") << " in beam " << beam << endl;
      int display_limit = 5;
      for (int i = 0; i < (int) hits.size(); ++i) {
        const DedopplerHit& hit = hits[i];
        if (i < display_limit) {
          cout << "  index = " << hit.index << ", drift steps = " << hit.drift_steps
               << ", snr = " << hit.snr << ", drift rate = " << hit.drift_rate << endl;
        }
        hit_recorder->recordHit(hit, beam, band, beamformer.power);
      }
      if ((int) hits.size() > display_limit) {
        cout << "  (and " << ((int) hits.size() - display_limit) << " more)\n";
      }
    }
    cout << "recorded " << total_hits << " hits in band " << band << endl;
  }
}
