#pragma once

#include <cuda.h>
#include <iostream>

using namespace std;

static_assert(sizeof(float) == 4, "require 32-bit floats");

const int CUDA_MAX_THREADS = 1024;

// Helpers to nicely display cuda errors
void checkCuda(const string& tag);

