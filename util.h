#pragma once

#include <thrust/complex.h>
#include <string>

using namespace std;

const float FLOAT_ONE = 1.0f;

int roundUpToPowerOfTwo(int n);
int numDigits(int n);
string zeroPad(int n, int size);
string cToS(thrust::complex<float> c);
string stripAnyTrailingSlash(const string& s);
void assertComplexEq(thrust::complex<float> c, float real, float imag);
void assertFloatEq(float a, float b);
void assertStringEq(const string& lhs, const string& rhs);
string pluralize(int n, const string& noun);
string prettyBytes(size_t n);
