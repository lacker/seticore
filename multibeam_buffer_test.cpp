#include "catch/catch.hpp"

#include "filterbank_buffer.h"
#include "multibeam_buffer.h"
#include "util.h"

TEST_CASE("copyRegionAsync", "[multibeam_buffer]") {
  auto mb = MultibeamBuffer(10, 10, 10);
  auto fb = FilterbankBuffer(10, 5);
  for (int beam = 0; beam < 10; ++beam) {
    for (int time = 0; time < 10; ++time) {
      for (int chan = 0; chan < 10; ++chan) {
        mb.set(beam, time, chan, 100.0 * beam + 10.0 * time + 1.0 * chan);
      }
    }
  }

  // Copy with an offset of 3 should take [7][7] to [7][4]
  mb.copyRegionAsync(5, 3, &fb);
  float value = fb.get(7, 4);
  REQUIRE(value == Approx(574.0));
}
