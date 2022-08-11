#include <assert.h>
#include <fmt/core.h>
#include <iostream>

#include "filterbank_buffer.h"
#include "taylor.h"
#include "util.h"

/*
  Performance testing on dedoppler.
 */
int main(int argc, char* argv[]) {
  int num_timesteps = 16;
  int num_channels = 1 << 24;
  FilterbankBuffer input(makeNoisyBuffer(num_timesteps, num_channels));

  // Draw a line
  input.set(0, 70, 1.0);
  input.set(1, 70, 1.0);
  input.set(2, 70, 1.0);
  input.set(3, 70, 1.0);
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

  FilterbankBuffer buffer1(num_timesteps, num_channels);
  FilterbankBuffer buffer2(num_timesteps, num_channels);

  for (int block = -2; block <= 2; ++block) {
    long start = timeInMS();
    fullTaylorTree(input.data, buffer1.data, buffer2.data,
                   num_timesteps, num_channels, block);
    cudaDeviceSynchronize();
    long end = timeInMS();
    cout << fmt::format("elapsed time {:.3f}s for block {}\n",
                        (end - start) / 1000.0, block);
  }

}
