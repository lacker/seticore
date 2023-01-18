#include "dedoppler_hit.h"

#include <assert.h>
#include <fmt/core.h>
#include <string>

#include "util.h"

using namespace std;

DedopplerHit::DedopplerHit(const FilterbankMetadata& metadata, int _index,
                           int _drift_steps, double _drift_rate,
                           float _snr, int _beam, int _coarse_channel,
                           int _num_timesteps, float _power)
  : index(_index), drift_steps(_drift_steps), drift_rate(_drift_rate),
    snr(_snr), coarse_channel(_coarse_channel),
    beam(metadata.isCoherentBeam(_beam) ? _beam : NO_BEAM),
    num_timesteps(_num_timesteps), power(_power), incoherent_power(0.0) {

  int coarse_offset = coarse_channel * metadata.coarse_channel_size;
  int global_index = coarse_offset + index;
  frequency = metadata.fch1 + global_index * metadata.foff;
}

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

// If we have incoherent power, we want to use the ratio of power to incoherent power.
// Otherwise, we just want to sort by SNR.
float DedopplerHit::score() const {
  if (incoherent_power > 0) {
    return power / incoherent_power;
  } else {
    assert(snr >= 0.0);
    return -1.0 / (1.0 + snr);
  }
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


// Alternate sort comparator that just compares drift steps
bool driftStepsLessThan(const DedopplerHit& lhs, const DedopplerHit& rhs) {
  return lhs.drift_steps < rhs.drift_steps;
}

