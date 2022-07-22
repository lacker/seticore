#include "catch/catch.hpp"

#include "beamformer.h"

TEST_CASE("cublasBeamform", "[beamformer]") {
  int nants = 8;
  int nbeams = 8;
  int nblocks = 8;
  int fft_size = 8;
  int num_coarse_channels = 8;
  int npol = 2;
  int nsamp = 512;
  Beamformer beamformer(0, fft_size, nants, nbeams, nblocks, num_coarse_channels,
                        npol, nsamp);

  RawBuffer raw(nblocks, nants, num_coarse_channels, nsamp / nblocks, npol);
  raw.set(1, 2, 3, 4, 1, false, 100);
  DeviceRawBuffer input(nblocks, nants, num_coarse_channels, nsamp / nblocks, npol);
  input.copyFromAsync(raw);
  input.waitUntilReady();

  // Try both ways
  MultibeamBuffer output1(nbeams, beamformer.numOutputTimesteps(),
                          beamformer.numOutputChannels());
  MultibeamBuffer output2(nbeams, beamformer.numOutputTimesteps(),
                          beamformer.numOutputChannels());
  beamformer.use_cublas_beamform = true;
  beamformer.release_input = false;
  beamformer.run(input, output1, 0);
  beamformer.use_cublas_beamform = false;
  beamformer.run(input, output2, 0);

  float value1 = output1.get(1, 2, 3);
  float value2 = output2.get(1, 2, 3);
  REQUIRE(value1 == Approx(value2));
}
