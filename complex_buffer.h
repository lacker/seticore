#pragma once

#include <thrust/complex.h>

using namespace std;

/*
  A buffer that holds complex numbers on the GPU. This is generic so that it
  can be reused by multiple stages of the pipeline.
 */
class ComplexBuffer {
 public:
  ComplexBuffer(size_t size);
  virtual ~ComplexBuffer();

  // No copying
  ComplexBuffer(const ComplexBuffer&) = delete;
  ComplexBuffer& operator=(ComplexBuffer&) = delete;
  
  // The number of complex entries
  const size_t size;

  // The number of bytes allocated
  const size_t bytes;

  // The data itself
  thrust::complex<float>* data;

  // Causes a cuda sync so it's slow. Only useful for debugging or testing
  thrust::complex<float> get(int index) const;
};
