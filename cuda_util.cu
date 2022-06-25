#include <iostream>

#include "cuda_util.h"

using namespace std;

void checkCuda(cudaError_t err) {
  if (err != 0) {
    cerr << "cuda error " << err << ": " << cudaGetErrorString(err) << endl;
    exit(1);
  }
}

