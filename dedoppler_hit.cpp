#include "dedoppler_hit.h"

#include <fmt/core.h>
#include <string>

using namespace std;

string DedopplerHit::toString() const {
  return fmt::format("index = {}, snr = {:.6f}, drift rate = {:.6f} ({} bins)",
                     index, snr, drift_rate, drift_steps);
}
