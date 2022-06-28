#include <assert.h>
#include <cuda.h>
#include <iostream>

#include "beamformer.h"
#include "cuda_util.h"

using namespace std;

// Helper to calculate a 4d row-major index, ie for:
//   arr[a][b][c][d]
__host__ __device__ int index4d(int a, int b, int b_end, int c, int c_end, int d, int d_end) {
  return ((((a * b_end) + b) * c_end) + c) * d_end + d;
}

/*
  We convert from int8 input with format:
    input[antenna][frequency][time][polarity][real or imag]

  to complex-float output with format:
    transposed[time][frequency][polarity][antenna]
 */
__global__ void transpose(const signed char* input, thrust::complex<float>* transposed,
                          int nants, int nchans, int npol, int nsamp) {
  int antenna = threadIdx.x;
  int pol = threadIdx.y;
  int chan = blockIdx.y;
  int time = blockIdx.x;

  int input_index = 2 * index4d(antenna, chan, nchans, time, nsamp, pol, npol);
  int transposed_index = index4d(time, chan, nchans, pol, npol, antenna, nants);

  transposed[transposed_index] = thrust::complex<float>
    (input[input_index] * 1.0, input[input_index + 1] * 1.0);
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
  cudaMallocManaged(&input, 2 * nants * nchans * npol * nsamp * sizeof(signed char));
  checkCuda("beamformer input malloc");
  
  cudaMallocManaged(&coefficients, 2 * nants * nbeams * nchans * npol * sizeof(float));
  checkCuda("beamformer coefficients malloc");

  cudaMallocManaged(&transposed, 2 * nants * nchans * npol * nsamp * sizeof(float));
  checkCuda("beamformer transposed malloc");
}

Beamformer::~Beamformer() {
  cudaFree(input);
  cudaFree(coefficients);
  cudaFree(transposed);
}

/*
  The caller should first put the data into *input and *coefficients.
 */
void Beamformer::beamform() {
  dim3 transpose_block(nants, npol, 1);
  dim3 transpose_grid(nsamp, nchans, 1);
  transpose<<<transpose_grid, transpose_block>>>(input, transposed, nants, nchans, npol, nsamp);
  checkCuda("transpose cuda kernel");
}

thrust::complex<float> Beamformer::getCoefficient(int antenna, int pol, int beam, int freq) const {
  cudaDeviceSynchronize();
  checkCuda("getCoefficient");
  assert(antenna < nants);
  assert(pol < npol);
  assert(beam < nbeams);
  assert(freq < nchans);
  int i = index4d(freq, beam, nbeams, pol, npol, antenna, nants);
  return thrust::complex<float>(coefficients[2*i], coefficients[2*i+1]);
}

thrust::complex<float> Beamformer::getTransposed(int time, int chan, int pol, int antenna) const {
  cudaDeviceSynchronize();
  checkCuda("getTransposed");
  assert(time < nsamp);
  assert(chan < nchans);
  assert(pol < npol);
  assert(antenna < nants);
  int i = index4d(time, chan, nchans, pol, npol, antenna, nants);
  return transposed[i];
}
