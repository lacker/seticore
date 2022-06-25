#include <cuda.h>

#include "beamformer.h"
#include "cuda_util.h"

using namespace std;

/*
  The Beamformer encapsulates the GPU memory allocations we use for beamforming.
  The workflow is to create a beamformer for a particular set of dimensions,
  use it to form many beams, and then destruct it when we want to free the memory.
 */
Beamformer::Beamformer(int nants, int nbeams, int nchans, int npol, int nsamp)
  : nants(nants), nbeams(nbeams), nchans(nchans), npol(npol), nsamp(nsamp) {
  checkCuda(cudaMallocManaged(&input, nants * nchans * npol * nsamp * sizeof(signed char)));
  checkCuda(cudaMallocManaged(&coefficients, nants * nbeams * nchans * npol * sizeof(float)));
}

Beamformer::~Beamformer() {
  cudaFree(input);
  cudaFree(coefficients);
}
