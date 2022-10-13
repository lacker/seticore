#include <assert.h>
#include "complex_buffer.h"
#include "cuda_util.h"

using namespace std;

ComplexBuffer::ComplexBuffer(size_t size) :
  size(size), bytes(sizeof(thrust::complex<float>) * size) {

  cudaMallocManaged(&data, bytes);
  checkCudaMalloc("ComplexBuffer", bytes);
}

ComplexBuffer::~ComplexBuffer() {
  cudaFree(data);
}

thrust::complex<float> ComplexBuffer::get(int index) const {
  assert(index >= 0);
  assert(index < size);
  cudaDeviceSynchronize();
  return data[index];
}
