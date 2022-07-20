#pragma once

#include "hdf5.h"

#include "filterbank_metadata.h"

using namespace std;

void writeH5(const string& filename,
             const FilterbankMetadata& metadata,
             const float* data);
