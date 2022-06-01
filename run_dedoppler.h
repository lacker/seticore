#pragma once

#include <string>

using namespace std;

void runDedoppler(const string& input_filename, const string& output_filename,
                  double max_drift, double snr, double min_drift);

