#include "dedoppler_hit.h"

#include <fmt/core.h>
#include <string>

#include "util.h"

using namespace std;

string DedopplerHit::toString() const {
  return fmt::format("index = {}, snr = {:.5f}, drift rate = {:.5f} ({})",
                     index, snr, drift_rate, pluralize(drift_steps, "bin"));
}
