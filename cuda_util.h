#pragma once

#include <cuda.h>
#include <iostream>

using namespace std;

static_assert(sizeof(float) == 4, "require 32-bit floats");

const int CUDA_MAX_THREADS = 1024;

// Helpers to nicely display cuda errors
void checkCuda(const string& tag);

// Helper to check errors and clean up
class Stream {
public:
  cudaStream_t stream;
  Stream();
  ~Stream();
};

// Helper to calculate a 3d row-major index, ie for:
//   arr[a][b][c]
__host__ __device__ inline int index3d(int a, int b, int b_end, int c, int c_end) {
  return ((a * b_end) + b) * c_end + c;
}

// Helper to calculate a 4d row-major index, ie for:
//   arr[a][b][c][d]
__host__ __device__ inline int index4d(int a, int b, int b_end, int c, int c_end, int d, int d_end) {
  return index3d(a, b, b_end, c, c_end) * d_end + d;
}

// Helper to calculate a 5d row-major index, ie for:
//   arr[a][b][c][d][e]
__host__ __device__ inline int index5d(int a, int b, int b_end, int c, int c_end, int d, int d_end,
                                int e, int e_end) {
  return index4d(a, b, b_end, c, c_end, d, d_end) * e_end + e;
}


