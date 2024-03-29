#pragma once

#include "cublas_v2.h"
#include <thrust/complex.h>

#include "complex_buffer.h"
#include "device_raw_buffer.h"
#include "multiantenna_buffer.h"
#include "multibeam_buffer.h"
#include "upchannelizer.h"

using namespace std;

/*
  The beamformer metadata is confusing because it's constantly reshaping the data it's
  operating on.
  In particular, it upchannelizes to convert time resolution for frequency resolution.
  "Coarse channel" refers to the channel index in the input data.
  Within a coarse channel, there are further "fine channels".
  The combination of a coarse channel index plus a fine channel index can give you
  a global frequency index of
  (coarse_channel * fft_size) + fine_channel
  Whenever data is indexed with adjacent indices of [coarse-channel][fine-channel],
  it's equivalent to a single global index of just [channel].

  The overall data flow is thus:

  upchannelizer:
    raw input -> prebeam

  incoherent:
    prebeam -> output

  beamform:
    prebeam -> buffer

  calculatePower:
    buffer -> output

  Note that "buffer" is also used internally by the upchannelizer to save memory.
 */
class Beamformer {
 public:
  // The factor for upchannelization. The frequency dimension will be expanded by this
  // amount, and the time dimension will be shrunk by this amount.
  const int fft_size;
  
  // Number of antennas in the input data
  const int num_antennas;

  // Number of beams to form
  const int num_beams;

  // Number of time-blocks that the input is provided in.
  // This only affects the format of the input array.
  const int num_blocks;
  
  // Number of frequency channels in the input.
  // This will be expanded by a multiplicative factor of fft_size.
  const int num_coarse_channels;

  // Number of polarizations
  const int num_polarizations;

  // Number of timesteps in the input.
  // This will be reduced by two multiplicative factors, fft_size and STI.
  const int num_input_timesteps;

  // Number of timesteps for short time integration.
  const int sti;
  
  // The beamformer runs all its operations on one stream.
  const cudaStream_t stream;
  
  Beamformer(cudaStream_t stream, int fft_size, int num_antennas, int num_beams,
             int num_blocks, int num_coarse_channels, int num_polarizations,
             int num_input_timesteps, int sti);
  ~Beamformer();

  int numOutputChannels() const;
  int numOutputTimesteps() const;

  void setCoefficient(int chan, int beam, int pol, int antenna, float real, float imag);
  
  void run(DeviceRawBuffer& input, MultibeamBuffer& output, int time_offset);

  // These cause a cuda sync so they are slow, only useful for debugging or testing
  thrust::complex<float> getCoefficient(int antenna, int pol, int beam, int coarse_channel) const;
  thrust::complex<float> getFFTBuffer(int pol, int antenna, int coarse_channel,
                                      int time, int last_index) const;
  thrust::complex<float> getPrebeam(int time, int channel, int pol, int antenna) const;
  thrust::complex<float> getVoltage(int time, int pol, int channel, int beam) const;
  float getPower(int beam, int time, int channel) const;
  
  // Beamforming coefficients, formatted by row-major:
  //   coefficients[coarse-channel][beam][polarization][antenna][real or imag]
  thrust::complex<float> *coefficients;
  size_t coefficients_size;

  // The square magnitude of the beamforming coefficients, formatted by row-major:
  //   square_magnitudes[coarse-channel][polarization][antenna]
  float *square_magnitudes;
  size_t square_magnitudes_size;
  
  // The beamformed data, as power. Not owned by the beamformer.
  // Its format is row-major:
  //   power[beam][time][frequency]
  float* power;
  size_t power_size;

  // Selects which of two alternatives for the beamform kernel to use
  bool use_cublas_beamform;

  // Selects whether we weight the incoherent beam
  bool weight_incoherent_beam;
  
  // Hack for unit testing
  void setReleaseInput(bool flag);
  
 private:
  unique_ptr<Upchannelizer> upchannelizer;
  
  // The buffer is reused for a few different stages of the pipeline.
  //
  // The convertRaw kernel populates this buffer with row-major:
  //   buffer[pol][antenna][coarse-channel][time]
  // which is equivalent to
  //   buffer[pol][antenna][coarse-channel][time-coarse-index][time-fine-index]
  //
  // The FFT converts the time fine index in-place to a frequency fine index, leaving this
  // buffer with row-major:
  //   buffer[pol][antenna][coarse-channel][time][fine-channel]
  //
  // After beamforming, the voltage data is also placed here, in row-major:
  //   voltage[time][polarization][channel][beam]
  unique_ptr<ComplexBuffer> buffer;
  
  // The channelized input data, ready for beamforming.
  //
  // Its format is row-major:
  //   prebeam[time][channel][polarization][antenna]
  unique_ptr<MultiantennaBuffer> prebeam;

  cublasHandle_t cublas_handle;

  void unweightedIncoherentBeam(float* output);
  void runCublasBeamform(int time, int pol);
};
