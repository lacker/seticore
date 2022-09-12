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

// Sort first by coarse channel, then by low index, then by high index
bool operator<(const DedopplerHit& lhs, const DedopplerHit& rhs) {
  if (lhs.coarse_channel != rhs.coarse_channel) {
    return lhs.coarse_channel < rhs.coarse_channel;
  }
  if (lhs.lowIndex() != rhs.lowIndex()) {
    return lhs.lowIndex() < rhs.lowIndex();
  }
  return lhs.highIndex() < rhs.highIndex();
}
