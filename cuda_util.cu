#include <iostream>

#include "cuda_util.h"

using namespace std;

void checkCuda(const string& tag) {
  auto err = cudaGetLastError();
  if (err != cudaSuccess) {
    cerr << tag << ": cuda error " << err << ": " << cudaGetErrorString(err) << endl;
    exit(2);
  }
}

void checkCudaMalloc(const string& tag, size_t bytes) {
  auto err = cudaGetLastError();
  if (err != cudaSuccess) {
    cerr << tag << ": cuda error " << err << " allocating " <<
      bytes << " bytes: " << cudaGetErrorString(err) << endl;
    exit(2);
  }
}

Stream::Stream() {
  cudaStreamCreate(&stream);
  checkCuda("cudaStreamCreate");
}

Stream::~Stream() {
  cudaStreamDestroy(stream);
  checkCuda("cudaStreamDestroy");
}
