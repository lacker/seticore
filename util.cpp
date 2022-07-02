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
  if (c.imag() >= 0) {
    return fmt::format("{:.6f} + {:.6f} i", c.real(), c.imag());
  }
  return fmt::format("{:.6f} - {:.6f} i", c.real(), -c.imag());
}
