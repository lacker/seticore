#include "beamforming_pipeline.h"

#include <assert.h>
#include <iostream>

#include "beamformer.h"
#include "cuda_util.h"
#include "dedoppler.h"
#include "dedoppler_hit.h"
#include "dedoppler_hit_group.h"
#include "h5_writer.h"
#include "hit_file_writer.h"
#include "hit_recorder.h"
#include "filterbank_buffer.h"
#include "filterbank_file_reader.h"
#include <fmt/core.h>
#include "multibeam_buffer.h"
#include "raw_file_group.h"
#include "raw_file_group_reader.h"
#include "raw_buffer.h"
#include "recipe_file.h"
#include "stamp_extractor.h"
#include "util.h"

// Construct metadata for the data created by a RawFileGroup and Beamformer.
// This metadata should apply to the entire span of channels, not just one band.
FilterbankMetadata combineMetadata(const RawFileGroup& file_group,
                                   const Beamformer& beamformer,
                                   const RecipeFile& recipe,
                                   int telescope_id) {
  FilterbankMetadata metadata;
  metadata.source_name = file_group.source_name;
  metadata.fch1 = file_group.getFch1(beamformer.fft_size);
  double output_bandwidth = file_group.obsbw / file_group.num_bands;
  metadata.foff = output_bandwidth / beamformer.numOutputChannels();
  int beamformer_batches = file_group.num_blocks / beamformer.num_blocks;
  double time_per_block = file_group.tbin * file_group.timesteps_per_block;
  double time_per_beamform = time_per_block * beamformer.num_blocks;
  metadata.tsamp = time_per_beamform / beamformer.numOutputTimesteps();
  metadata.tstart = unixTimeToMJD(file_group.getStartTime(0));
  metadata.src_raj = file_group.ra;
  metadata.src_dej = file_group.dec;
  metadata.num_timesteps = beamformer.numOutputTimesteps() * beamformer_batches;
  metadata.num_channels = beamformer.numOutputChannels() * file_group.num_bands;
  metadata.telescope_id = telescope_id;
  metadata.coarse_channel_size = beamformer.fft_size;
  if (0 != metadata.num_channels % metadata.coarse_channel_size) {
    fatal(fmt::format("error in combineMetadata: num fine channels = {} but "
                      "coarse channel size = {}", metadata.num_channels,
                      metadata.coarse_channel_size));
  }
  metadata.num_coarse_channels = metadata.num_channels / metadata.coarse_channel_size;
  metadata.source_names = recipe.src_names;
  metadata.ras = recipe.getRAsInHours();
  metadata.decs = recipe.getDecsInDegrees();
  
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
void BeamformingPipeline::findHits() {
  cout << fmt::format("processing {:.1f}s of data from {}.*.raw\n",
                      file_group.totalTime(), file_group.prefix);  
  cout << "using beamforming recipe from " << recipe_filename << endl;

  RecipeFile recipe(recipe_filename, file_group.obsid);
  recipe.validateRawRange(file_group.schan, file_group.num_coarse_channels);
  
  // Do enough blocks per beamformer batch to handle one STI block
  if (0 != (sti * fft_size) % file_group.timesteps_per_block) {
    fatal(fmt::format("invalid parameters: sti = {}, fft_size = {}, "
                      "timesteps per block = {}", sti, fft_size,
                      file_group.timesteps_per_block));
  }
  int blocks_per_batch = (sti * fft_size) / file_group.timesteps_per_block;

  if (file_group.num_coarse_channels % num_bands != 0) {
    fatal(fmt::format("{}.*.raw has {} coarse channels, so we cannot "
                      "divide it into {} bands", file_group.prefix,
                      file_group.num_coarse_channels, num_bands));
  }

  int coarse_channels_per_band = file_group.num_coarse_channels / num_bands;
  int nsamp = file_group.timesteps_per_block * blocks_per_batch;

  Beamformer beamformer(0, fft_size, file_group.nants, recipe.nbeams,
                        blocks_per_batch, coarse_channels_per_band, file_group.npol,
                        nsamp, sti);

  // Create a buffer large enough to hold all beamformer batches for one band  
  int num_batches = file_group.num_blocks / beamformer.num_blocks;
  int num_multibeam_timesteps = beamformer.numOutputTimesteps() * num_batches;
  if (num_multibeam_timesteps < 2) {
    cout << "this recording is too short to process. the output would only have "
         << num_multibeam_timesteps << " timestep"
         << (num_multibeam_timesteps == 1 ? "s" : "") << "." << endl;
    return;
  }
  MultibeamBuffer multibeam(beamformer.num_beams + 1,
                            num_multibeam_timesteps,
                            beamformer.numOutputChannels(),
                            beamformer.numOutputTimesteps());
  
  RawFileGroupReader reader(file_group, num_bands, 0, num_bands_to_process - 1,
                            num_batches, blocks_per_batch);
  
  // Create a buffer for dedopplering a single coarse channel, padding
  // timesteps with zeros.
  FilterbankBuffer fb_buffer(roundUpToPowerOfTwo(multibeam.num_timesteps), fft_size);
  fb_buffer.zero();
  cout << "filterbank buffer memory: " << prettyBytes(fb_buffer.bytes) << endl;
  
  FilterbankMetadata metadata = combineMetadata(file_group, beamformer, recipe,
                                                telescope_id);

  Dedopplerer dedopplerer(multibeam.num_timesteps,
                          fb_buffer.num_channels,
                          metadata.foff, metadata.tsamp, false);
  dedopplerer.print_hit_summary = true;
  cout << "dedoppler memory: " << prettyBytes(dedopplerer.memoryUsage()) << endl;

  unique_ptr<HitFileWriter> hit_recorder;
  if (record_hits) {
    string output_filename = fmt::format("{}/{}.hits", output_dir,
                                         file_group.prefix);
    cout << "recording hits to " << output_filename << endl;
    auto hfw = new HitFileWriter(output_filename, metadata);
    hfw->verbose = false;
    hit_recorder.reset(hfw);
  }
  
  cout << "processing " << pluralize(beamformer.num_beams, "beam") << " and "
       << pluralize(num_bands_to_process, "band") << endl;
  cout << "each band has "
       << pluralize(beamformer.num_coarse_channels, "coarse channel")
       << ", for a total of " << file_group.num_coarse_channels << endl;
  cout << "dedoppler input is " << fb_buffer.num_timesteps << " timesteps x "
       << fb_buffer.num_channels << " fine channels\n";
  cout << fmt::format("dedoppler resolution is {:.1f} s, {:.1f} hz\n",
                      file_group.tbin * fft_size * sti,
                      file_group.coarseChannelBandwidth() / fft_size * 1'000'000);

  for (int band = 0; band < num_bands_to_process; ++band) {
    // We read data into the read buffer, and copy data to the GPU from the
    // work buffer.
    // At the start of this loop, neither buffer is being used, because
    // dedoppler analysis for any previous loop synchronized cuda devices.
    cout << "beamforming band " << band << "...\n";
    for (int batch = 0; batch < num_batches; ++batch) {
    
      int block_after_mid = batch * beamformer.num_blocks + beamformer.num_blocks / 2;
      double mid_time = file_group.getStartTime(block_after_mid);
      int time_array_index = recipe.getTimeArrayIndex(mid_time);
      recipe.generateCoefficients(time_array_index,
				  file_group.schan,
				  file_group.num_coarse_channels,
				  file_group.obsfreq,
				  file_group.obsbw,
				  coarse_channels_per_band * band,
                                  coarse_channels_per_band,
                                  beamformer.coefficients);

      int time_offset = beamformer.numOutputTimesteps() * batch;

      shared_ptr<DeviceRawBuffer> device_raw_buffer = reader.readToDevice();

      // At this point, the beamformer could still be processing the
      // previous batch, but that's okay.
      multibeam.hintWritingTime(time_offset);
      beamformer.run(*device_raw_buffer, multibeam, time_offset);
    }

    for (int beam = 0; beam < multibeam.num_beams; ++beam) {
      multibeam.hintReadingBeam(beam);
      
      if (!h5_dir.empty()) {
        // Write out data for this band and beam to a file
        cudaDeviceSynchronize();
        string beam_name = metadata.isCoherentBeam(beam) ?
          fmt::format("beam{}", zeroPad(beam, numDigits(beamformer.num_beams))) :
          "incoherent";
        string h5_filename =
          fmt::format("{}/{}.band{}.{}.h5",
                      h5_dir,
                      file_group.prefix,
                      zeroPad(band, numDigits(num_bands)),
                      beam_name);
        FilterbankMetadata band_metadata = metadata.getSubsetMetadata(beam, band,
                                                                      num_bands);
        FilterbankBuffer output(multibeam.getBeam(beam));
        H5Writer writer(h5_filename, band_metadata);
        writer.setData(output.data);
        writer.close();
      }
      
      // local_coarse_channel is the index of the coarse channel within the band
      for (int local_coarse_channel = 0;
           local_coarse_channel < coarse_channels_per_band;
           ++local_coarse_channel) {
        int coarse_channel = band * coarse_channels_per_band + local_coarse_channel;

        multibeam.copyRegionAsync(beam, local_coarse_channel * fft_size, &fb_buffer);

        vector<DedopplerHit> local_hits;
        dedopplerer.search(fb_buffer, metadata, beam, coarse_channel, max_drift,
                           0.0, snr, &local_hits);
        for (DedopplerHit hit : local_hits) {
          if (record_hits) {
            hit_recorder->recordHit(hit, fb_buffer.data);
          }
          hits.push_back(hit);
        }
      }
    }
  }
}

/*
  Runs the beamforming pipeline from hits to stamps.
 */
void BeamformingPipeline::makeStamps() {
  int margin = 30;
  int max_stamps = 3;

  string output_filename = fmt::format("{}/{}.stamps", output_dir,
                                       file_group.prefix);
  StampExtractor extractor(file_group, fft_size, telescope_id, output_filename);

  int stamps_created = 0;
  vector<DedopplerHitGroup> groups = makeHitGroups(hits, margin);
  for (const DedopplerHitGroup& group : groups) {
    if (stamps_created >= max_stamps) {
      break;
    }
    const DedopplerHit& top_hit = group.topHit();
    
    if (top_hit.drift_steps == 0) {
      // This is a vertical line. No drift = terrestrial. Skip it
      continue;
    }

    // Extract the stamp
    int first_channel = max(0, top_hit.lowIndex() - margin);
    int last_channel = min(fft_size - 1, top_hit.highIndex() + margin);
    cout << "top hit: " << top_hit.toString() << endl;
    cout << fmt::format("extracting fine channels {} to {} from coarse channel {}\n",
                        first_channel, last_channel, top_hit.coarse_channel);
    extractor.extract(&top_hit,
                      top_hit.coarse_channel,
                      first_channel,
                      last_channel - first_channel + 1);
    
    stamps_created++;
  }

  cout << "saved " << pluralize(stamps_created, "stamp") << " to " << output_filename
       << endl;
}
