#pragma once

#include <cuda.h>

using namespace std;

__global__ void taylorTreeOneStepKernel(const float* source_buffer, float* target_buffer,
                                        int num_timesteps, int num_freqs, int path_length,
                                        int drift_block);

const float* fullTaylorTree(const float* source_buffer, float* buffer1, float* buffer2,
                            int num_timesteps, int num_freqs, int drift_block);
