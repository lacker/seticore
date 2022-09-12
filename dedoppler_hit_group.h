#pragma once

#include "dedoppler_hit.h"
#include <vector>

using namespace std;

class DedopplerHitGroup {
 public:
  const int coarse_channel;
  
  vector<DedopplerHit> hits;

  // Lowest lowIndex across all hits
  int low_index;

  // Highest highIndex across all hits
  int high_index;

  DedopplerHitGroup(DedopplerHit hit);
  void add(DedopplerHit hit);
  ~DedopplerHitGroup() {}
};

vector<DedopplerHitGroup> makeHitGroups(vector<DedopplerHit>& hits, int margin);
