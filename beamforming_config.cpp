#include "beamforming_config.h"

#include <assert.h>
#include <iostream>

#include "beamformer.h"
#include "cuda_util.h"
#include "dedoppler.h"
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
#include "util.h"

// Construct metadata for the data created by a RawFileGroup and Beamformer.
// This metadata should apply to the entire span of channels, not just one band.
FilterbankMetadata combineMetadata(const RawFileGroup& file_group,
                                   const Beamformer& beamformer,
                                   const RecipeFile& recipe,
                                   int telescope_id) {
  FilterbankMetadata metadata;
  metadata.source_name = file_group.source_name;
  metadata.fch1 = file_group.obsfreq - 0.5 * file_group.obsbw;
  double output_bandwidth = file_group.obsbw / file_group.num_bands;
  metadata.foff = output_bandwidth / beamformer.numOutputChannels();
  int beamformer_batches = file_group.num_blocks / beamformer.nblocks;
  double time_per_block = file_group.tbin * file_group.timesteps_per_block;
  double time_per_beamform = time_per_block * beamformer.nblocks;
  metadata.tsamp = time_per_beamform / beamformer.numOutputTimesteps();
  metadata.tstart = unixTimeToMJD(file_group.getStartTime(0));
  metadata.src_raj = file_group.ra;
  metadata.src_dej = file_group.dec;
  metadata.num_timesteps = beamformer.numOutputTimesteps() * beamformer_batches;
  metadata.num_channels = beamformer.numOutputChannels() * file_group.num_bands;
  metadata.telescope_id = telescope_id;
  metadata.coarse_channel_size = beamformer.fft_size;
  metadata.num_coarse_channels = metadata.num_channels / metadata.coarse_channel_size;
  metadata.source_names = recipe.src_names;
  metadata.ras = recipe.ras;
  metadata.decs = recipe.decs;
  
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

  cout << fmt::format("processing {:.1f}s of data from {}.*.raw\n",
                      file_group.totalTime(), file_group.prefix);
  
  RecipeFile recipe(recipe_dir, file_group.obsid);

  // Do enough blocks per beamformer batch to handle one STI block
  assert((sti * fft_size) % file_group.timesteps_per_block == 0);
  int blocks_per_batch = (sti * fft_size) / file_group.timesteps_per_block;
  
  int coarse_channels_per_band = file_group.num_coarse_channels / num_bands;
  int nsamp = file_group.timesteps_per_block * blocks_per_batch;

  Beamformer beamformer(0, fft_size, file_group.nants, recipe.nbeams,
                        blocks_per_batch, coarse_channels_per_band, file_group.npol,
                        nsamp, sti);

  // Create a buffer large enough to hold all beamformer batches for one band  
  int num_batches = file_group.num_blocks / beamformer.nblocks;
  int num_multibeam_timesteps = beamformer.numOutputTimesteps() * num_batches;
  MultibeamBuffer multibeam(beamformer.nbeams + 1,
                            num_multibeam_timesteps,
                            beamformer.numOutputChannels(),
                            beamformer.numOutputTimesteps());
  
  RawFileGroupReader reader(file_group, num_bands_to_process, num_batches,
                            blocks_per_batch);
  
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
  
  // These buffers hold data we have read from raw input while we are working on them
  unique_ptr<RawBuffer> read_buffer;
  shared_ptr<DeviceRawBuffer> device_raw_buffer = reader.makeDeviceBuffer();
  cout << "raw buffer memory: " << prettyBytes(device_raw_buffer->size) << endl;

  
  if (hit_recorder == nullptr) {
    string output_filename = fmt::format("{}/{}.hits", output_dir,
                                         file_group.prefix);
    cout << "recording hits to " << output_filename << endl;
    auto hfw = new HitFileWriter(output_filename, metadata);
    hfw->verbose = false;
    hit_recorder.reset(hfw);
  }
  
  cout << "processing " << pluralize(beamformer.nbeams, "beam") << " and "
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
    cout << endl;

    // We read data into the read buffer, and copy data to the GPU from the
    // work buffer.
    // At the start of this loop, neither buffer is being used, because
    // dedoppler analysis for any previous loop synchronized cuda devices.
    for (int batch = 0; batch < num_batches; ++batch) {
      cout << "beamforming band " << band << ", batch " << batch << "...\n";
    
      int block_after_mid = batch * beamformer.nblocks + beamformer.nblocks / 2;
      double mid_time = file_group.getStartTime(block_after_mid);
      int time_array_index = recipe.getTimeArrayIndex(mid_time);
      recipe.generateCoefficients(time_array_index, file_group.schan,
                                  beamformer.num_coarse_channels,
                                  file_group.obsfreq, file_group.obsbw,
                                  beamformer.coefficients);
      int time_offset = beamformer.numOutputTimesteps() * batch;

      // The beamformer could still be using the raw buffers.
      // So we have to wait for it to finish.
      device_raw_buffer->waitUntilUnused();

      reader.returnBuffer(move(read_buffer));
      read_buffer = reader.read();
      device_raw_buffer->copyFromAsync(*read_buffer);
      device_raw_buffer->waitUntilReady();

      // At this point, the beamformer could still be processing the
      // previous batch, but that's okay.
      multibeam.hintWritingTime(time_offset);
      beamformer.run(*device_raw_buffer, multibeam, time_offset);
    }

    cout << endl;
    for (int beam = 0; beam < multibeam.num_beams; ++beam) {
      multibeam.hintReadingBeam(beam);
      
      if (!h5_dir.empty()) {
        // Write out data for this band and beam to a file
        cudaDeviceSynchronize();
        string beam_name = metadata.isCoherentBeam(beam) ?
          fmt::format("beam{}", zeroPad(beam, numDigits(beamformer.nbeams))) :
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

        vector<DedopplerHit> hits;
        dedopplerer.search(fb_buffer, beam, coarse_channel, max_drift, min_drift, snr,
                           &hits);
        for (DedopplerHit hit : hits) {
          hit_recorder->recordHit(hit, beam, coarse_channel, fb_buffer.data);
        }
      }
    }
  }
}
