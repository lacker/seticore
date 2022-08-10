#include "catch/catch.hpp"
#include <iostream>
#include <vector>

#include "dedoppler.h"
#include "filterbank_buffer.h"

TEST_CASE("basic functionality", "[dedoppler]") {
  int num_timesteps = 8;
  int num_channels = 1000;
  FilterbankBuffer buffer(makeNoisyBuffer(num_timesteps, num_channels));

  // Draw a line drifting a few timesteps to the right
  buffer.set(0, 70, 1.0);
  buffer.set(1, 70, 1.0);
  buffer.set(2, 71, 1.0);
  buffer.set(3, 71, 1.0);
  buffer.set(4, 72, 1.0);
  buffer.set(5, 72, 1.0);
  buffer.set(6, 73, 1.0);
  buffer.set(7, 73, 1.0);

  Dedopplerer dedopplerer(num_timesteps, num_channels, 1.0, 1.0, false);
  vector<DedopplerHit> hits;
  dedopplerer.search(buffer, NO_BEAM, 555, 0.01, 0.01, 200.0, &hits);
  REQUIRE(hits.size() == 1);
  REQUIRE(hits[0].index == 70);
  REQUIRE(hits[0].drift_steps == 3);
}

