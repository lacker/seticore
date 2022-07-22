#include "catch/catch.hpp"

#include "beamformer.h"

TEST_CASE("cublasBeamform", "[beamformer]") {
  int nants = 8;
  int nbeams = 8;
  int nblocks = 8;
  int fft_size = 8;
  int num_coarse_channels = 8;
  int nsamp = 512;
  Beamformer beamformer(0, fft_size, nants, nbeams, nblocks, num_coarse_channels,
                        2, nsamp);
  
}
