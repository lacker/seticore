#include <assert.h>
#include <fmt/core.h>
#include <iostream>

#include "filterbank_buffer.h"
#include "taylor.h"
#include "util.h"

/*
  Performance testing the taylor tree inner loop.
 */
int main(int argc, char* argv[]) {
  int num_timesteps = 16;
  int num_channels = 1 << 24;
  FilterbankBuffer input(makeNoisyBuffer(num_timesteps, num_channels));

  FilterbankBuffer buffer1(num_timesteps, num_channels);
  FilterbankBuffer buffer2(num_timesteps, num_channels);

  for (int drift_block = -2; drift_block <= 2; ++drift_block) {
    cout << "\ndrift block " << drift_block << endl;
    long start = timeInMS();
    basicTaylorTree(input.data, buffer1.data, buffer2.data,
                    num_timesteps, num_channels, drift_block);
    cudaDeviceSynchronize();
    long end = timeInMS();
    cout << fmt::format("old algorithm: elapsed time {:.3f}s\n",
                        (end - start) / 1000.0);

    start = timeInMS();
    tiledTaylorTree(input.data, buffer1.data, num_timesteps, num_channels, drift_block);
    cudaDeviceSynchronize();
    end = timeInMS();
    cout << fmt::format("new algorithm: elapsed time {:.3f}s\n",
                        (end - start) / 1000.0);
    
  }

}
