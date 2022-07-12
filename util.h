#pragma once

#include <thrust/complex.h>
#include <string>

using namespace std;

int roundUpToPowerOfTwo(int n);
string cToS(thrust::complex<float> c);
string stripAnyTrailingSlash(const string& s);
