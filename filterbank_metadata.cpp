#include <assert.h>
#include <fmt/core.h>
#include <iostream>
#include "util.h"

#include "filterbank_metadata.h"

using namespace std;

FilterbankMetadata::FilterbankMetadata(): has_dc_spike(false), coarse_channel_size(0) {}

/*
  Gives the slightly altered metadata that you get if you split this file into one beam
  and one subband.
  The coarse channels have to divide evenly into bands.
 */
FilterbankMetadata FilterbankMetadata::getSubsetMetadata(int beam, int band,
                                                         int num_bands) const {
  assert(num_bands > 0);
  assert(0 <= band && band < num_bands);
  if (0 != num_coarse_channels % num_bands) {
    fatal(fmt::format("in getSubsetMetadata, num_coarse_channels = {} but "
                      "num_bands = {}", num_coarse_channels, num_bands));
  }

  FilterbankMetadata subset;
  subset.has_dc_spike = has_dc_spike;
  subset.source_name = getBeamSourceName(beam);
  subset.source_names = source_names;
  subset.ras = ras;
  subset.decs = decs;

  subset.num_timesteps = num_timesteps;
  subset.num_channels = num_channels / num_bands;
  subset.coarse_channel_size = coarse_channel_size;
  subset.num_coarse_channels = num_coarse_channels / num_bands;

  assert(subset.num_channels > 0);

  subset.fch1 = fch1 + (foff * subset.num_channels * band);
  subset.foff = foff;
  subset.tstart = tstart;
  subset.tsamp = tsamp;

  subset.src_raj = getBeamRA(beam);
  subset.src_dej = getBeamDec(beam);
  subset.telescope_id = telescope_id;
  
  return subset;
}

/*
  Guesses some metadata from other metadata:
    coarse_channel_size
    has_dc_spike
    num_coarse_channels
 */
void FilterbankMetadata::inferMetadata() {
  if ((num_timesteps == 16 || num_timesteps == 32) && num_channels % 1048576 == 0) {
    // Looks like Green Bank data
    assert(telescope_id == NO_TELESCOPE_ID || telescope_id == GREEN_BANK);
    if (coarse_channel_size == 0) {
      coarse_channel_size = 1048576;
    }
    has_dc_spike = true;
  } else if (num_channels == 50331648) {
    // Looks like ATA data
    assert(telescope_id == NO_TELESCOPE_ID || telescope_id == ATA);
    if (coarse_channel_size == 0) {
      coarse_channel_size = 262144;
    }
    has_dc_spike = false;
  } else if (coarse_channel_size > 0) {
    // We already have a coarse channel size set, so we only want to infer whether we
    // have a dc spike.
    has_dc_spike = (telescope_id == GREEN_BANK);
  } else {
    cerr << "unable to infer coarse channel size for data with dimensions: "
         << num_timesteps << " x " << num_channels << ". please set the nfpc header."
         << endl;
    exit(1);
  }
  num_coarse_channels = num_channels / coarse_channel_size;
}

string FilterbankMetadata::getBeamSourceName(int beam) const {
  if (!isCoherentBeam(beam)) {
    return source_name;
  }
  return source_names[beam];
}

double FilterbankMetadata::getBeamRA(int beam) const {
  if (!isCoherentBeam(beam)) {
    return src_raj;
  }
  return ras[beam];
}

double FilterbankMetadata::getBeamDec(int beam) const {
  if (!isCoherentBeam(beam)) {
    return src_dej;
  }
  return decs[beam];
}

bool FilterbankMetadata::isCoherentBeam(int beam) const {
  return 0 <= beam && beam < (int) source_names.size();
}
