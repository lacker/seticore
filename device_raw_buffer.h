#pragma once

#include "raw_buffer.h"

using namespace std;

/*
  The DeviceRawBuffer stores data from a raw file in GPU memory.

  Its format is row-major:
    input[block][antenna][coarse-channel][time-within-block][polarity][real or imag]
 */

class DeviceRawBuffer {
 public:
  const int num_blocks;
  const int num_antennas;
  const int num_coarse_channels;
  const int timesteps_per_block;
  const int npol;

  int8_t* data;
  size_t data_size;
  
  DeviceRawBuffer(int num_blocks, int num_antennas, int num_coarse_channels,
                  int timesteps_per_block, int npol);

  ~DeviceRawBuffer();

  void copyFromAsync(const RawBuffer& source);

};
