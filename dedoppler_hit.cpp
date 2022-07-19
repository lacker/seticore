#include "dedoppler_hit.h"

#include <fmt/core.h>
#include <string>

using namespace std;

string DedopplerHit::toString() const {
  return fmt::format("index = {}, snr = {:.5f}, drift rate = {:.5f} ({} bins)",
                     index, snr, drift_rate, drift_steps);
}
