#include "util.h"

#include <fmt/core.h>

using namespace std;

int roundUpToPowerOfTwo(int n) {
  int rounded = 1;
  while (rounded < n) {
    rounded *= 2;
  }
  return rounded;
}

string cToS(thrust::complex<float> c) {
  return fmt::format("{:.6f} + {:.6f} i", c.real(), c.imag());
}

string stripAnyTrailingSlash(const string& s) {
  if (s.empty()) {
    return s;
  }
  if (s[s.size() - 1] == '/') {
    return s.substr(0, s.size() - 1);
  }
  return s;
}
