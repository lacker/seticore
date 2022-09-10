#include <assert.h>
#include "dedoppler_hit_group.h"

using namespace std;

DedopplerHitGroup::DedopplerHitGroup(DedopplerHit hit)
  : coarse_channel(hit.coarse_channel) {
  hits.push_back(hit);
  low_index = hit.lowIndex();
  high_index = hit.highIndex();
}

void DedopplerHitGroup::add(DedopplerHit hit) {
  assert(hit.coarse_channel == coarse_channel);
  low_index = min(low_index, hit.lowIndex());
  high_index = max(high_index, hit.highIndex());
  hits.push_back(hit);
}
