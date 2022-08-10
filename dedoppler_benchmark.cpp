#include <assert.h>
#include <fmt/core.h>
#include <iostream>
#include <time.h>
#include "util.h"

#include "dedoppler.h"

/*
  Performance testing on dedoppler.
 */
int main(int argc, char* argv[]) {
  FilterbankBuffer buffer(makeNoisyBuffer(16, 1 << 24));

  // Draw a line drifting a few timesteps to the right
  buffer.set(0, 70, 1.0);
  buffer.set(1, 70, 1.0);
  buffer.set(2, 71, 1.0);
  buffer.set(3, 71, 1.0);
  buffer.set(4, 72, 1.0);
  buffer.set(5, 72, 1.0);
  buffer.set(6, 73, 1.0);
  buffer.set(7, 73, 1.0);

  Dedopplerer dedopplerer(buffer.num_timesteps, buffer.num_channels, 1.0, 1.0, false);
  vector<DedopplerHit> hits;

  long start = timeInMS();
  dedopplerer.search(buffer, NO_BEAM, 555, 0.01, 0.01, 200.0, &hits);
  long end = timeInMS();
  cout << fmt::format("elapsed time: {:.3f}s\n", (end - start) / 1000.0);

  assert(hits.size() == 1);
  assert(hits[0].index == 70);
}
