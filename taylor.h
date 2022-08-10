#pragma once

#include <cuda.h>

using namespace std;

__global__ void taylorTree(const float* source_buffer, float* target_buffer,
                           int num_timesteps, int num_freqs, int path_length,
                           int drift_block);
