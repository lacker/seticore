#include "catch/catch.hpp"
#include <fmt/core.h>
#include <iostream>

#include "filterbank_buffer.h"
#include "taylor.h"

TEST_CASE("optimized taylor outputs match basic algorithm", "[taylor]") {
  for (int num_timesteps = 4; num_timesteps <= 128; num_timesteps *= 2) {
    int num_channels = 1000;
    FilterbankBuffer input(num_timesteps, num_channels);
    for (int time = 0; time < num_timesteps; ++time) {
      for (int chan = 0; chan < num_channels; ++chan) {
        input.set(time, chan, 1.0 * (1 + time + chan)); 
      }
    }
  
    FilterbankBuffer basic_buffer1(num_timesteps, num_channels);
    FilterbankBuffer basic_buffer2(num_timesteps, num_channels);

    FilterbankBuffer opt_buffer1(num_timesteps, num_channels);
    FilterbankBuffer opt_buffer2(num_timesteps, num_channels);

    for (int drift_block = -2; drift_block <= 2; ++drift_block) {
      // Avoid confusion while debugging
      basic_buffer1.zero();
      basic_buffer2.zero();
      opt_buffer1.zero();
      opt_buffer2.zero();
      cudaDeviceSynchronize();
    
      const float* basic_ptr = basicTaylorTree(input.data,
                                               basic_buffer1.data, basic_buffer2.data,
                                               num_timesteps, num_channels, drift_block);
      const FilterbankBuffer& basic = (basic_ptr == basic_buffer1.data)
        ? basic_buffer1 : basic_buffer2;

      const float* opt_ptr = optimizedTaylorTree(input.data,
                                                 opt_buffer1.data, opt_buffer2.data,
                                                 num_timesteps, num_channels,
                                                 drift_block);
      const FilterbankBuffer& opt = (opt_ptr == opt_buffer1.data)
        ? opt_buffer1 : opt_buffer2;
      
      cudaDeviceSynchronize();

      basic.assertEqual(opt, drift_block);
    }
  
  }
}

