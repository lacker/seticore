#include "util.h"

#include <assert.h>
#include <fmt/core.h>
#include <iostream>

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

void assertComplexEq(thrust::complex<float> c, float real, float imag) {
  if (abs(c.real() - real) > 0.001) {
    cerr << "c.real() = " << c.real() << " but real = " << real << endl;
    exit(1);
  }
  if (abs(c.imag() - imag) > 0.001) {
    cerr << "c.imag() = " << c.imag() << " but imag = " << imag << endl;
    exit(1);
  }
}

void assertFloatEq(float a, float b) {
  while (a > 100) {
    a /= 2.0;
    b /= 2.0;
  }
  if (abs(a - b) > 0.001) {
    cerr << a << " != " << b << endl;
    exit(1);
  }
}

string pluralize(int n, const string& noun) {
  if (n == 1) {
    return fmt::format("1 {}", noun);
  }
  return fmt::format("{} {}s", n, noun);
}
