#include <iostream>

#include "cuda_util.h"

using namespace std;

void checkCuda(const string& tag) {
  auto err = cudaGetLastError();
  if (err != cudaSuccess) {
    cerr << tag << ": cuda error " << err << ": " << cudaGetErrorString(err) << endl;
    exit(1);
  }
}


 
