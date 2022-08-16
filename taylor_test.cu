#include "catch/catch.hpp"
#include <fmt/core.h>
#include <iostream>

#include "filterbank_buffer.h"
#include "taylor.h"

TEST_CASE("tiled taylor outputs match basic algorithm", "[taylor]") {
  for (int num_timesteps = 4; num_timesteps <= 32; num_timesteps *= 2) {
    int num_channels = 2000;
    FilterbankBuffer input(num_timesteps, num_channels);
    for (int time = 0; time < num_timesteps; ++time) {
      for (int chan = 0; chan < num_channels; ++chan) {
        input.set(time, chan, 1.0 * (1 + time + chan)); 
      }
    }
  
    FilterbankBuffer buffer1(num_timesteps, num_channels);
    FilterbankBuffer buffer2(num_timesteps, num_channels);
    FilterbankBuffer tiled(num_timesteps, num_channels);

    for (int drift_block = -2; drift_block <= 2; ++drift_block) {
      // Avoid confusion while debugging
      buffer1.zero();
      buffer2.zero();
      tiled.zero();
      cudaDeviceSynchronize();
    
      const float* out_ptr = basicTaylorTree(input.data, buffer1.data, buffer2.data,
                                             num_timesteps, num_channels, drift_block);
      const FilterbankBuffer& basic = (out_ptr == buffer1.data) ? buffer1 : buffer2;

      tiledTaylorTree(input.data, tiled.data, num_timesteps, num_channels, drift_block);

      cudaDeviceSynchronize();

      basic.assertEqual(tiled, drift_block);
    }
  }
}

TEST_CASE("two-pass taylor outputs match basic algorithm", "[taylor]") {
  int num_timesteps = 64;
  int num_channels = 1000;
  FilterbankBuffer input(num_timesteps, num_channels);
  for (int time = 0; time < num_timesteps; ++time) {
    for (int chan = 0; chan < num_channels; ++chan) {
      input.set(time, chan, 1.0 * (1 + time + chan)); 
    }
  }
  
  FilterbankBuffer buffer1(num_timesteps, num_channels);
  FilterbankBuffer buffer2(num_timesteps, num_channels);

  FilterbankBuffer two_pass_buffer(num_timesteps, num_channels);
  FilterbankBuffer two_pass_output(num_timesteps, num_channels);

  for (int drift_block = 0; drift_block <= 2; ++drift_block) {
    // Avoid confusion while debugging
    buffer1.zero();
    buffer2.zero();
    two_pass_buffer.zero();
    two_pass_output.zero();
    cudaDeviceSynchronize();
    
    const float* out_ptr = basicTaylorTree(input.data, buffer1.data, buffer2.data,
                                           num_timesteps, num_channels, drift_block);
    const FilterbankBuffer& basic = (out_ptr == buffer1.data) ? buffer1 : buffer2;
    // const FilterbankBuffer& other = (out_ptr == buffer1.data) ? buffer2 : buffer1;

    twoPassTaylorTree(input.data, two_pass_buffer.data, two_pass_output.data,
                      num_timesteps, num_channels, drift_block);

    cudaDeviceSynchronize();

    basic.assertEqual(two_pass_output, drift_block);
  }
  
}
