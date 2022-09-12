#pragma once

#include "dedoppler_hit.h"
#include <vector>

using namespace std;

class DedopplerHitGroup {
 public:
  int coarse_channel;
  
  vector<DedopplerHit> hits;

  // Lowest lowIndex across all hits
  int low_index;

  // Highest highIndex across all hits
  int high_index;

  // Index of the hit with highest snr across all hits
  int top_hit_index;
  
  DedopplerHitGroup(DedopplerHit hit);
  void add(DedopplerHit hit);
  ~DedopplerHitGroup() {}

  const DedopplerHit& topHit() const;

  // No copying
  DedopplerHitGroup(const DedopplerHitGroup& that) = delete;

  // Moving is okay
  DedopplerHitGroup(DedopplerHitGroup&& that) = default;
  DedopplerHitGroup& operator=(DedopplerHitGroup&&) = default;
};

vector<DedopplerHitGroup> makeHitGroups(vector<DedopplerHit>& hits, int margin);

bool operator<(const DedopplerHitGroup& lhs, const DedopplerHitGroup& rhs);
