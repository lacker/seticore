#pragma once

#include <thrust/complex.h>
#include <string>

using namespace std;

int roundUpToPowerOfTwo(int n);
bool isPowerOfTwo(int n);
int numDigits(int n);
string zeroPad(int n, int size);
string cToS(thrust::complex<float> c);
string stripAnyTrailingSlash(const string& s);
void assertComplexEq(thrust::complex<float> c, float real, float imag);
void assertFloatEq(float a, float b);
void assertFloatEq(float a, float b, const string& tag);
void assertStringEq(const string& lhs, const string& rhs);
void assertStringEq(const string& lhs, const string& rhs, const string& tag);
string pluralize(int n, const string& noun);
string prettyBytes(size_t n);
long timeInMS();
double unixTimeToMJD(double unix_time);
double hoursToRadians(double hours);
double degreesToRadians(double degrees);
