#include <assert.h>
#include <fmt/core.h>
#include <iostream>

#include "filterbank_buffer.h"
#include "taylor.h"
#include "util.h"

const int nk_timesteps = 16;
const int nk_width = 256;

// NOTE: this screws up the edges, assumes everything divides evenly, etc
// It's about four times faster than the old kernel though.
// TODO: fix up the edge cases
__global__ void buggyKernel(const float* global_input, float* global_output,
                            int num_timesteps, int num_channels, int drift_block) {
  assert(num_timesteps == nk_timesteps);

  int tx = threadIdx.x;
  int block_start = blockIdx.x * nk_width;
  __shared__ float buffer1[nk_timesteps * nk_width];
  __shared__ float buffer2[nk_timesteps * nk_width];

  // This thread block operates on shifted subarrays
  const float* input = global_input + block_start;
  float* output = global_output + block_start;
  
  // Path length goes 2, 4, 8, 16
  // input -> buffer1 -> buffer2 -> output
  taylorOneStepOneChannel(input, buffer1, tx, 16, num_channels, nk_width, 2, drift_block);
  __syncthreads();
  taylorOneStepOneChannel(buffer1, buffer2, tx, 16, nk_width, nk_width, 4, drift_block);
  __syncthreads();
  taylorOneStepOneChannel(buffer2, buffer1, tx, 16, nk_width, nk_width, 8, drift_block);
  __syncthreads();
  taylorOneStepOneChannel(buffer1, output, tx, 16, nk_width, num_channels, 16, drift_block);
}



/*
  Performance testing the taylor tree inner loop.
 */
int main(int argc, char* argv[]) {
  int num_timesteps = nk_timesteps;
  int num_channels = 1 << 24;
  FilterbankBuffer input(makeNoisyBuffer(num_timesteps, num_channels));

  // Draw a line
  input.set(0, 70, 1.0);
  input.set(1, 70, 1.0);
  input.set(2, 70, 1.0);
  input.set(3, 70, 1.0);
  /*
  input.set(4, 70, 1.0);
  input.set(5, 70, 1.0);
  input.set(6, 70, 1.0);
  input.set(7, 70, 1.0);
  input.set(8, 71, 1.0);
  input.set(9, 71, 1.0);
  input.set(10, 71, 1.0);
  input.set(11, 71, 1.0);
  input.set(12, 71, 1.0);
  input.set(13, 71, 1.0);
  input.set(14, 71, 1.0);
  input.set(15, 71, 1.0);
  */

  FilterbankBuffer buffer1(num_timesteps, num_channels);
  FilterbankBuffer buffer2(num_timesteps, num_channels);

  for (int drift_block = -2; drift_block <= 2; ++drift_block) {
    cout << "\ndrift block " << drift_block << endl;
    long start = timeInMS();
    auto out = fullTaylorTree(input.data, buffer1.data, buffer2.data,
                              num_timesteps, num_channels, drift_block);
    assert(out == buffer2.data);
    cudaDeviceSynchronize();
    long end = timeInMS();
    cout << fmt::format("old algorithm: elapsed time {:.3f}s\n",
                        (end - start) / 1000.0);

    start = timeInMS();
    buggyKernel<<<num_channels / nk_width, nk_width>>>(input.data, buffer1.data,
                                                       num_timesteps, num_channels,
                                                       drift_block);
    cudaDeviceSynchronize();
    end = timeInMS();
    cout << fmt::format("new algorithm: elapsed time {:.3f}s\n",
                        (end - start) / 1000.0);
    
  }

}
