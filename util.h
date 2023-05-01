#pragma once

#include <cuda_runtime.h>
#include <thrust/complex.h>
#include <string>

using namespace std;

const string VERSION = "1.0.2";

// This is allegedly a SIGPROC standard but the most authoritative source
// I can find is:
//   https://github.com/UCBerkeleySETI/blimpy/blob/master/blimpy/ephemeris/observatory_info.csv
const int NO_TELESCOPE_ID = -1;
const int PARKES = 4;
const int GREEN_BANK = 6;
const int ATA = 9;
const int VLA = 12;
const int MEERKAT = 64;

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
double radiansToHours(double radians);
double degreesToRadians(double degrees);
double radiansToDegrees(double radians);
void logErrorTimestamp();
void logError(const string& message);
void fatal(const string& message);
void fatal(const string& message1, const string& message2);
int telescopeID(const string& telescope);
