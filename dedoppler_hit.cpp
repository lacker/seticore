#include "dedoppler_hit.h"

#include <fmt/core.h>
#include <string>

#include "util.h"

using namespace std;

string DedopplerHit::toString() const {
  return fmt::format("coarse channel = {}, index = {}, snr = {:.5f}, "
                     "drift rate = {:.5f} ({})",
                     coarse_channel, index, snr, drift_rate,
                     pluralize(drift_steps, "bin"));
}

int DedopplerHit::lowIndex() const {
  return min(index, index + drift_steps);
}

int DedopplerHit::highIndex() const {
  return max(index, index + drift_steps);
}
