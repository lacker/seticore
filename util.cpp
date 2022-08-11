#include "util.h"

#include <assert.h>
#include <chrono>
#include <fmt/core.h>
#include <iostream>
#include <pthread.h>

using namespace std;

int roundUpToPowerOfTwo(int n) {
  int rounded = 1;
  while (rounded < n) {
    rounded *= 2;
  }
  return rounded;
}

int numDigits(int n) {
  if (n < 10) {
    return 1;
  }
  return 1 + numDigits(n / 10);
}

string zeroPad(int n, int size) {
  string formatter = fmt::format("{{:0{}d}}", size);
  return fmt::format(formatter, n);
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
  assertFloatEq(a, b, "");
}

void assertFloatEq(float a, float b, const string& tag) {
  float initial_a = a;
  float initial_b = b;
  while (a > 100) {
    a /= 2.0;
    b /= 2.0;
  }
  if (abs(a - b) > 0.001) {
    if (!tag.empty()) {
      cerr << tag << ": ";
    }
    cerr << fmt::format("{:.3f} != {:.3f}\n", initial_a, initial_b);
    exit(1);
  }
}

void assertStringEq(const string& lhs, const string& rhs, const string& tag) {
  if (lhs == rhs) {
    return;
  }
  if (tag.empty()) {
    cerr << "assertStringEq failed:\n";
  } else {
    cerr << "assertStringEq failed for " << tag << ":\n";
  }
  cerr << "lhs: " << lhs << endl;
  cerr << "rhs: " << rhs << endl;
  exit(1);
}

void assertStringEq(const string& lhs, const string& rhs) {
  assertStringEq(lhs, rhs, "");
}

string pluralize(int n, const string& noun) {
  if (n == 1 || n == -1) {
    return fmt::format("{} {}", n, noun);
  }
  return fmt::format("{} {}s", n, noun);
}

string prettyBytes(size_t n) {
  size_t mb = 1024 * 1024;
  size_t gb = 1024 * mb;
  if (n >= gb) {
    return fmt::format("{:.1f} GB", ((float) n) / gb);
  }
  if (n >= 10 * mb) {
    return fmt::format("{} MB", n / mb);
  }
  if (n >= mb) {
    return fmt::format("{:.1f} MB", ((float) n) / mb);
  }  
  return fmt::format("{} bytes", n);
}

long timeInMS() {
  return chrono::duration_cast<chrono::milliseconds>
    (chrono::system_clock::now().time_since_epoch()).count();
}
