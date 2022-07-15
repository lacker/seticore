#pragma once

using namespace std;

/*
  The RawBuffer stores data from a raw file, or perhaps just a sub-band of it,
  in GPU memory.

  Its format is row-major:
    input[block][antenna][coarse-channel][time-within-block][polarity][real or imag]
*/

class RawBuffer {
 public:
  const int num_blocks;
  const int num_antennas;
  const int num_coarse_channels;
  const int timesteps_per_block;
  const int npol;

  int8_t* data;
  size_t data_size;
  
  RawBuffer(int num_blocks, int num_antennas, int num_coarse_channels,
            int timesteps_per_block, int npol);

  ~RawBuffer();

  // Returns a pointer to where data for the given block should start
  char* blockPointer(int block) const;
};
