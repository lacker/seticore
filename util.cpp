#include "util.h"

using namespace std;

int roundUpToPowerOfTwo(int n) {
  int rounded = 1;
  while (rounded < n) {
    rounded *= 2;
  }
  return rounded;
}
