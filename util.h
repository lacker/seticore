#pragma once

#include <thrust/complex.h>
#include <string>

using namespace std;

int roundUpToPowerOfTwo(int n);
string cToS(thrust::complex<float> c);
string stripAnyTrailingSlash(const string& s);
void assertComplexEq(thrust::complex<float> c, float real, float imag);
void assertFloatEq(float a, float b);
string pluralize(int n, const string& noun);
