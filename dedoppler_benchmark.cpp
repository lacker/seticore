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

  // Draw a line
  buffer.set(0, 70, 1.0);
  buffer.set(1, 70, 1.0);
  buffer.set(2, 70, 1.0);
  buffer.set(3, 70, 1.0);
  buffer.set(4, 70, 1.0);
  buffer.set(5, 70, 1.0);
  buffer.set(6, 70, 1.0);
  buffer.set(7, 70, 1.0);
  buffer.set(8, 71, 1.0);
  buffer.set(9, 71, 1.0);
  buffer.set(10, 71, 1.0);
  buffer.set(11, 71, 1.0);
  buffer.set(12, 71, 1.0);
  buffer.set(13, 71, 1.0);
  buffer.set(14, 71, 1.0);
  buffer.set(15, 71, 1.0);

  Dedopplerer dedopplerer(buffer.num_timesteps, buffer.num_channels, 1.0, 1.0, false);
  vector<DedopplerHit> hits;

  long start = timeInMS();
  dedopplerer.search(buffer, NO_BEAM, 555, 0.01, 0.01, 200.0, &hits);
  long end = timeInMS();
  cout << fmt::format("elapsed time: {:.3f}s\n", (end - start) / 1000.0);

  cout << "hits size: " << hits.size() << endl;
  for (auto hit : hits) {
    cout << "index: " << hit.index << ", steps: " << hit.drift_steps << endl;
  }
}
