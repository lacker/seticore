#include <algorithm>
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

/*
  Group together any hits that do not have 'margin' columns between them.
  This does modify the input hits by sorting them.
*/
vector<DedopplerHitGroup> makeHitGroups(vector<DedopplerHit>& hits, int margin) {
  vector<DedopplerHitGroup> output;

  sort(hits.begin(), hits.end());

  for (const DedopplerHit& hit : hits) {
    if (output.empty()) {
      // Start the first group with this hit
      output.emplace_back(hit);
      continue;
    }

    if (hit.coarse_channel == output.back().coarse_channel &&
        hit.lowIndex() <= output.back().high_index + margin) {
      // Merge this hit into the last group
      output.back().add(hit);
      continue;
    }

    output.emplace_back(hit);
  }

  return output;
}
