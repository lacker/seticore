#pragma once

#include <string>
#include <vector>

using namespace std;

void findEvents(const vector<string>& input_filenames, const string& output_filename,
                double max_drift, double snr_on, double snr_off);
