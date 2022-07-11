#include "catch/catch.hpp"
#include <iostream>
#include <vector>

#include "dedoppler.h"
#include "filterbank_buffer.h"

// Make a filterbank buffer with a bit of deterministic noise so that
// normalization doesn't make everything infinite SNR.
FilterbankBuffer noise(int num_timesteps, int num_channels) {
  FilterbankBuffer buffer(num_timesteps, num_channels);
  buffer.zero();
  for (int chan = 0; chan < buffer.num_channels; ++chan) {
    buffer.setValue(0, chan, 0.1 * chan / buffer.num_channels);
  }
  return buffer;
}

TEST_CASE("basic functionality", "[dedoppler]") {
  int num_timesteps = 8;
  int num_channels = 1000;
  FilterbankBuffer buffer = noise(num_timesteps, num_channels);

  // Draw a line drifting a few timesteps to the right
  buffer.setValue(0, 70, 1.0);
  buffer.setValue(1, 70, 1.0);
  buffer.setValue(2, 71, 1.0);
  buffer.setValue(3, 71, 1.0);
  buffer.setValue(4, 72, 1.0);
  buffer.setValue(5, 72, 1.0);
  buffer.setValue(6, 73, 1.0);
  buffer.setValue(7, 73, 1.0);

  Dedopplerer dedopplerer(num_timesteps, num_channels, 1.0, 1.0, false);
  vector<DedopplerHit> hits;
  dedopplerer.search(buffer, 0.01, 0.01, 200.0, &hits);
  REQUIRE(hits.size() == 1);
  REQUIRE(hits[0].index == 70);
  REQUIRE(hits[0].drift_steps == 3);
}
