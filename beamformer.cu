#include <cuda.h>
#include <iostream>
#include <thrust/complex.h>

#include "beamformer.h"
#include "cuda_util.h"

using namespace std;

/*
  The input encodes each value as two bytes, so it has size (2*num_values).
  You can think of it as row-major:
    input[index][0 for real, 1 for imag]

  We convert from signed char to float, treating the signed char as a signed 8-bit int.

  TODO: see if we can use thrust::complex<int8> instead of signed char, for input
 */
__global__ void convertToFloat(const signed char* input, thrust::complex<float>* output,
                               int num_values) {
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < 0 || i >= num_values) {
    return;
  }
  output[i].real(input[2*i] * 1.0f);
  output[i].imag(input[2*i+1] * 1.0f);
}

/*
  The Beamformer encapsulates the GPU memory allocations we use for beamforming.
  The workflow is to create a beamformer for a particular set of dimensions,
  use it to form many beams, and then destruct it when we want to free the memory.

  TODO: nants and npol are specified twice, one by the recipe file and one by the raw input.
  We should really check to ensure they are the same and handle it cleanly if they aren't.
 */
Beamformer::Beamformer(int nants, int nbeams, int nchans, int npol, int nsamp)
  : nants(nants), nbeams(nbeams), nchans(nchans), npol(npol), nsamp(nsamp) {
  checkCuda(cudaMallocManaged(&input, nants * nchans * npol * nsamp * sizeof(signed char)));
  checkCuda(cudaMallocManaged(&coefficients, 2 * nants * nbeams * nchans * npol * sizeof(float)));
}

Beamformer::~Beamformer() {
  cudaFree(input);
  cudaFree(coefficients);
}

/*
  The caller should first set *input and *coefficients.
 */
void Beamformer::beamform() {
  // TODO: implement
}

void Beamformer::debugCoefficients(int antenna, int pol, int beam, int freq) const {
  int i = (((freq * nbeams) + beam) * npol + pol) * nants + antenna;
  cout << "coefficient[" << antenna << "][" << pol << "][" << beam << "][" << freq << "] = "
       << coefficients[2*i] << " + " << coefficients[2*i+1] << "i\n";
}
