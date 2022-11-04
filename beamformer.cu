#include <assert.h>
#include "cublas_v2.h"
#include <cuda.h>
#include <cufft.h>
#include <fmt/core.h>
#include <iostream>

#include "beamformer.h"
#include "cuda_util.h"
#include "util.h"

using namespace std;

const int MAX_STI = 8;

const cuComplex COMPLEX_ONE = make_cuComplex(1.0, 0.0);
const cuComplex COMPLEX_ZERO = make_cuComplex(0.0, 0.0);
const float ONE = 1.0f;
const float ZERO = 0.0f;

/*
  Beamforming combines the channelized data with format:
    prebeam[time][coarse-channel][fine-channel][polarization][antenna]

  with the coefficient data, format:
    coefficients[coarse-channel][beam][polarization][antenna]

  to generate output beams with format:
    voltage[time][polarization][coarse-channel][fine-channel][beam]

  We combine prebeam with coefficients according to the indices they have in common,
  not conjugating the coefficients because we expect them to already be in the
  correct conjugation for multiplying, and then sum along antenna dimension to reduce.
*/
__global__ void beamform(const thrust::complex<float>* prebeam,
                         const thrust::complex<float>* coefficients,
                         thrust::complex<float>* voltage,
                         int fft_size, int num_antennas, int num_beams,
                         int num_coarse_channels, int num_polarizations, int num_timesteps,
                         int prebeam_size, int voltage_size, int coefficients_size) {
  int antenna = threadIdx.x;
  int fine_chan = blockIdx.x;
  int coarse_chan = blockIdx.y;
  int beam = blockIdx.z;

  const int MAX_ANTS = 64;
  assert(num_antennas <= MAX_ANTS);
  __shared__ thrust::complex<float> reduced[MAX_ANTS];

  for (int pol = 0; pol < num_polarizations; ++pol) {
    int coeff_index = index4d(coarse_chan, beam, num_beams, pol, num_polarizations,
                              antenna, num_antennas);
    assert(0 <= coeff_index);
    assert(coeff_index < coefficients_size);
    thrust::complex<float> conjugated = thrust::conj(coefficients[coeff_index]);
    for (int time = 0; time < num_timesteps; ++time) {
      int prebeam_index = index5d(time, coarse_chan, num_coarse_channels, fine_chan,
                                  fft_size, pol, num_polarizations, antenna, num_antennas);
      assert(0 <= prebeam_index);
      assert(prebeam_index < prebeam_size);
      assert(0 <= antenna);
      assert(antenna < MAX_ANTS);
      reduced[antenna] = prebeam[prebeam_index] * conjugated;

      __syncthreads();

      for (int k = MAX_ANTS / 2; k > 0; k >>= 1) {
        if (antenna < k && antenna + k < num_antennas) {
          assert(antenna + k < MAX_ANTS);
          reduced[antenna] += reduced[antenna + k];
        }
        __syncthreads();
      }

      if (antenna == 0) {
        long voltage_index = index5d(time, pol, num_polarizations,
                                     coarse_chan, num_coarse_channels,
                                     fine_chan, fft_size, beam, num_beams);
        assert(0 <= voltage_index);
        assert(voltage_index < voltage_size);
        voltage[voltage_index] = reduced[0];
      }
    }
  }
}

/*
  Runs beamforming just for the provided time and polarization, using cublas batch
  matrix multiplication.
  See the comment on the beamform kernel.

  This API is not immediately intuitive. See this blog post:
    https://developer.nvidia.com/blog/cublas-strided-batched-matrix-multiply/
  in particular their explanation of gemmStridedBatched.

  To convert into the notation used by the blog post:

  A = coefficients (which we will conjugate and transpose)
  B = prebeam
  C = voltage
  m = beam
  n = fine channel
  p = coarse channel
  k = antenna
  alpha = 1.0
  beta = 0.0

  We are converting five-dimensional data plus four-dimensional data into
  five-dimensional output, so we fix time and polarization for each call to cublas.
  We could have fixed fine channel instead of time, or coarse channel instead of
  polarization, but it seems better to fix the smaller dimensions.
  Honestly, there are so many possibilities for how to arrange this data, I have
  not even begun to test all the ways it could work.
 */
void Beamformer::runCublasBeamform(int time, int pol) {
  // Calculate where the matrices start
  int coeff_offset = index4d(0, 0, num_beams, pol, num_polarizations, 0, num_antennas);
  auto coeff_start = (const cuComplex*) (coefficients + coeff_offset);
  int prebeam_offset = index5d(time, 0, num_coarse_channels, 0, fft_size,
                               pol, num_polarizations, 0, num_antennas);
  auto prebeam_start = (const cuComplex*) (prebeam->data + prebeam_offset);
  int voltage_offset = index5d(time, pol, num_polarizations, 0, num_coarse_channels,
                               0, fft_size, 0, num_beams);
  auto voltage_start = (cuComplex*) (buffer->data + voltage_offset);

  // Calculate strides
  // ldA, the A-m stride (since we are transposing. normally it would be A-k)
  int coeff_beam_stride = index4d(0, 1, num_beams, 0, num_polarizations, 0, num_antennas);
  // strideA, the A-p stride
  int coeff_coarse_stride = index4d(1, 0, num_beams, 0, num_polarizations, 0, num_antennas);
  // ldB, the B-n stride
  int prebeam_fine_stride = index5d(0, 0, num_coarse_channels, 1, fft_size,
                                    0, num_polarizations, 0, num_antennas);
  // strideB, the B-p stride
  int prebeam_coarse_stride = index5d(0, 1, num_coarse_channels, 0, fft_size,
                                      0, num_polarizations, 0, num_antennas);
  // ldC, the C-n stride
  int voltage_fine_stride = index5d(0, 0, num_polarizations, 0, num_coarse_channels,
                                    1, fft_size, 0, num_beams);
  // strideC, the C-p stride
  int voltage_coarse_stride = index5d(0, 0, num_polarizations, 1, num_coarse_channels,
                                      0, fft_size, 0, num_beams);

  cublasSetStream(cublas_handle, stream);
  cublasCgemm3mStridedBatched
    (cublas_handle, 
     CUBLAS_OP_C, CUBLAS_OP_N,
     num_beams, fft_size, num_antennas,
     &COMPLEX_ONE,
     coeff_start, coeff_beam_stride, coeff_coarse_stride, 
     prebeam_start, prebeam_fine_stride, prebeam_coarse_stride,
     &COMPLEX_ZERO,
     voltage_start, voltage_fine_stride, voltage_coarse_stride,
     num_coarse_channels);
  checkCuda("Beamformer runCublasBeamform");
}

/*
  The prebeamforming data, interpreted as complex, has format:
    prebeam[time][channel][pol][antenna]

  We want to combine by STI, so break down this time into "big_timestep"
  and "little_timestep". Also, polarization, antenna, and real-vs-imaginary all get
  handled the same way; they all get squared and added into the output value.
  So we can interpret the data as:
    prebeam[big_timestep][little_timestep][channel][combined_index]

  The output should just have format:
    output[big_timestep][channel]

  This method should write numOutputChannels() * numOutputTimesteps() floats,
  starting at output.

  We want to compute x^T x, ie the norm of x. To do this we treat x as a
  1 x num_combined matrix and use matrix-matrix multiplication.

  Note: once we upgrade to cuda 11.6.2 or later everywhere, we can implement this
  with a gemv routine, ie a matrix-vector multiplication, which might be faster.
  That might not be necessary, because it seems like some cuda versions, like
  11.5 at least, are smart enough to do this already.
 */
void Beamformer::formIncoherentBeam(float* output) {
  float* input = (float*) prebeam->data;
  int num_combined = num_polarizations * num_antennas * 2;
  int output_timesteps = numOutputTimesteps();

  cublasSetStream(cublas_handle, stream);
  
  for (int big_timestep = 0; big_timestep < output_timesteps; ++big_timestep) {
    // Add to the big_timestep row in the output
    float* output_start = output + big_timestep * numOutputChannels();
    
    for (int little_timestep = 0; little_timestep < sti; ++little_timestep) {
      int input_timestep = big_timestep * sti + little_timestep;
      float* input_start = input + index3d(input_timestep,
                                           0, numOutputChannels(),
                                           0, num_combined);

      cublasSgemmStridedBatched
        (cublas_handle, CUBLAS_OP_N, CUBLAS_OP_T,
         1, 1, num_combined,
         &ONE,
         input_start, 1, num_combined,
         input_start, 1, num_combined,
         (little_timestep == 0) ? &ZERO : &ONE,
         output_start, 1, 1,
         numOutputChannels());
      
       checkCuda("Beamformer formIncoherentBeam");
    }
  }
}

/*
  We calculate power and shrink the data at the same time.
  Every polarization and every window of STI adjacent timesteps gets reduced to a single
  power value, by adding the norm of each complex voltage.

  The input voltages have format: 
    voltage[time][polarization][frequency][beam]

  and the output power has format:
    power[beam][time][frequency]

  integrated_timestep refers to an index after integration. So a voltage at time t has:
    integrated_timestep = t / STI

  The power array is typically a much larger array in which we are only populating a
  subset of times.
  num_power_timesteps refers to the size of the time dimension of this array, and
  power_time_offset is the time at which to start writing to it.

  TODO: this also seems equivalent to a batch matrix multiplication. could we do this
  with a cublas routine?
 */
__global__ void calculatePower(const thrust::complex<float>* voltage, long voltage_size,
                               float* power, long power_size,
                               int num_beams, int num_channels, int num_polarizations, int sti,
			       int num_power_timesteps, int power_time_offset) {
  int chan = blockIdx.x;
  int beam = blockIdx.y;
  int integrated_timestep = blockIdx.z;
  int output_timestep = integrated_timestep + power_time_offset;
  
  int subintegration_timestep = threadIdx.x;
  assert(0 <= subintegration_timestep && subintegration_timestep < sti);
  int time = integrated_timestep * sti + subintegration_timestep;

  assert(2 == num_polarizations);
  long pol0_index = index4d(time, 0, num_polarizations, chan, num_channels,
                            beam, num_beams);
  assert(0 <= pol0_index && pol0_index < voltage_size);
  long pol1_index = index4d(time, 1, num_polarizations, chan, num_channels,
                            beam, num_beams);
  assert(0 <= pol1_index && pol1_index < voltage_size);
  long power_index = index3d(beam, output_timestep, num_power_timesteps,
			     chan, num_channels);
  assert(0 <= power_index && power_index < power_size);

  __shared__ float reduced[MAX_STI];
  float real0 = voltage[pol0_index].real();
  float imag0 = voltage[pol0_index].imag();
  float real1 = voltage[pol1_index].real();
  float imag1 = voltage[pol1_index].imag();
  reduced[subintegration_timestep] =
    real0 * real0 + imag0 * imag0 + real1 * real1 + imag1 * imag1;

  __syncthreads();

  for (int k = sti / 2; k > 0; k >>= 1) {
    if (subintegration_timestep < k) {
      reduced[subintegration_timestep] += reduced[subintegration_timestep + k];
    }
    __syncthreads();
  }

  assert(power_index >= 0);
  if (subintegration_timestep == 0) {
    power[power_index] = reduced[0];
  }
}

/*
  The Beamformer encapsulates the GPU memory allocations we use for beamforming.
  The workflow is to create a beamformer for a particular set of dimensions,
  use it to form many beams, and then destroy it when we want to free the memory.
 */
Beamformer::Beamformer(cudaStream_t stream, int fft_size, int num_antennas, int num_beams,
                       int num_blocks, int num_coarse_channels, int num_polarizations,
                       int num_input_timesteps, int sti)
  : fft_size(fft_size), num_antennas(num_antennas), num_beams(num_beams),
    num_blocks(num_blocks), num_coarse_channels(num_coarse_channels),
    num_polarizations(num_polarizations), num_input_timesteps(num_input_timesteps), sti(sti),
    stream(stream), use_cublas_beamform(true) {
  
  assert(0 == num_input_timesteps % (sti * fft_size));
  assert(0 == num_input_timesteps % num_blocks);
  assert(roundUpToPowerOfTwo(fft_size) == fft_size);
  assert(num_coarse_channels > 0);
  assert(fft_size > 0);
  
  upchannelizer = make_unique<Upchannelizer>(stream, fft_size,
                                             num_input_timesteps, num_coarse_channels,
                                             num_polarizations, num_antennas);
  
  size_t frame_size = num_coarse_channels * num_input_timesteps;
  
  coefficients_size = num_antennas * num_beams * num_coarse_channels * num_polarizations;
  size_t coefficients_bytes = coefficients_size * sizeof(thrust::complex<float>);
  cudaMallocManaged(&coefficients, coefficients_bytes);
  checkCudaMalloc("Beamformer coefficients", coefficients_bytes);

  // Sanity checking parameters for better debug messages
  size_t upper_bound = (size_t) 16 * 1024 * 1024 * 1024;
  
  size_t fft_buffer_size = num_antennas * num_polarizations * frame_size;

  if (fft_buffer_size > upper_bound) {
    fatal(fmt::format("fft buffer size is too large. it's a product of: "
                      "num_antennas = {}, "
                      "num_polarizations = {}, "
                      "num_coarse_channels = {}, "
                      "num_input_timesteps = {}",
                      num_antennas, num_polarizations,
                      num_coarse_channels, num_input_timesteps));
  }

  size_t voltage_size = num_beams * num_polarizations * frame_size;
  size_t buffer_size = max(fft_buffer_size, voltage_size);

  buffer = make_unique<ComplexBuffer>(buffer_size);
  prebeam = make_unique<MultiantennaBuffer>(num_input_timesteps / fft_size,
                                            num_coarse_channels * fft_size,
                                            num_polarizations, num_antennas);

  cublasCreate(&cublas_handle);
  checkCuda("Beamformer cublas handle");
  
  size_t total_bytes = coefficients_bytes + buffer->bytes + prebeam->bytes;
  if (total_bytes > 2000000) {
    cout << "beamformer memory: " << prettyBytes(total_bytes) << endl;
  }
}

Beamformer::~Beamformer() {
  cudaFree(coefficients);
  cublasDestroy(cublas_handle);
}

void Beamformer::setReleaseInput(bool flag) {
  upchannelizer->release_input = flag;
}

int Beamformer::numOutputChannels() const {
  return num_coarse_channels * fft_size;
}

int Beamformer::numOutputTimesteps() const {
  return num_input_timesteps / (fft_size * sti); 
}

/*
  Power from beamforming the input is written into output, with an offset
  of time_offset.

  The format of the input is row-major:
    input[block][antenna][coarse-channel][time-within-block][polarization][real or imag]

  The format of the output is row-major:
     power[beam][time][channel]
  but its time resolution has been reduced by a factor of (fft_size * STI), and we are
  only writing into a subrange of the time, starting at power_time_offset.

  The input must be ready to go when run is called.

  It is okay to call run before the previous run completes, as long as you
  only call it from one thread.
 */
void Beamformer::run(DeviceRawBuffer& input, MultibeamBuffer& output,
                     int power_time_offset) {
  assert(input.num_blocks == num_blocks);
  assert(input.num_antennas == num_antennas);
  assert(input.num_coarse_channels == num_coarse_channels);
  assert(input.timesteps_per_block * input.num_blocks == num_input_timesteps);
  assert(input.num_polarizations == num_polarizations);
  assert(output.num_beams == num_beams || output.num_beams == num_beams + 1);
  assert(output.num_channels == numOutputChannels());

  assert(sti <= MAX_STI);
  
  // If the output has an extra beam, fill it with incoherent beamforming
  bool incoherent = (output.num_beams > num_beams);
  
  upchannelizer->run(input, *buffer, *prebeam);

  if (incoherent) {
    // The incoherent beam goes into the beam numbered num_beams in the output
    long data_offset = index3d(num_beams,
			       power_time_offset, output.num_timesteps,
			       0, output.num_channels);
    int num_output_floats = numOutputChannels() * numOutputTimesteps();
    assert(data_offset >= 0);
    assert(data_offset + num_output_floats <= output.size());
    formIncoherentBeam(output.data + data_offset);
  }
  
  if (use_cublas_beamform) {
    for (int time = 0; time < num_input_timesteps / fft_size; ++time) {
      for (int pol = 0; pol < num_polarizations; ++pol) {
        runCublasBeamform(time, pol);
      }
    }
  } else {
    dim3 beamform_block(num_antennas, 1, 1);
    dim3 beamform_grid(fft_size, num_coarse_channels, num_beams);
    beamform<<<beamform_grid, beamform_block, 0, stream>>>
      (prebeam->data, coefficients, buffer->data, fft_size, num_antennas, num_beams,
       num_coarse_channels, num_polarizations, num_input_timesteps / fft_size, prebeam->size,
       buffer->size, coefficients_size);
  }
  
  dim3 power_block(sti, 1, 1);
  dim3 power_grid(numOutputChannels(), num_beams, numOutputTimesteps());
  calculatePower<<<power_grid, power_block, 0, stream>>>
    (buffer->data, buffer->size, output.data, output.size(), num_beams,
     numOutputChannels(), num_polarizations, sti, output.num_timesteps,
     power_time_offset);
  checkCuda("Beamformer calculatePower");
}

thrust::complex<float> Beamformer::getCoefficient(int antenna, int pol, int beam,
                                                  int coarse_channel) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getCoefficient");
  assert(antenna < num_antennas);
  assert(pol < num_polarizations);
  assert(beam < num_beams);
  assert(coarse_channel < num_coarse_channels);
  int i = index4d(coarse_channel, beam, num_beams, pol, num_polarizations,
                  antenna, num_antennas);
  assert(i < coefficients_size);
  return coefficients[i];
}

// The last index can be either time's fine index or the fine channel index, depending
// on whether it's pre-FFT or post-FFT.
thrust::complex<float> Beamformer::getFFTBuffer(int pol, int antenna, int coarse_channel,
                                                int time, int last_index) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getFFTBuffer");
  assert(pol < num_polarizations);
  assert(antenna < num_antennas);
  assert(coarse_channel < num_coarse_channels);
  assert(time * fft_size < num_input_timesteps);
  assert(last_index < fft_size);
  int i = index5d(pol, antenna, num_antennas, coarse_channel, num_coarse_channels,
                  time, num_input_timesteps / fft_size, last_index, fft_size);
  return buffer->data[i];
}

thrust::complex<float> Beamformer::getPrebeam(int time, int channel, int pol,
                                              int antenna) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getPrebeam");
  assert(time < num_input_timesteps);
  assert(channel < num_coarse_channels * fft_size);
  assert(pol < num_polarizations);
  assert(antenna < num_antennas);
  int i = index4d(time, channel, num_coarse_channels * fft_size, pol, num_polarizations,
                  antenna, num_antennas);
  return prebeam->data[i];
}

thrust::complex<float> Beamformer::getVoltage(int time, int pol, int channel,
                                              int beam) const {
  cudaDeviceSynchronize();
  checkCuda("Beamformer getVoltage");
  assert(time < num_input_timesteps);
  assert(pol < num_polarizations);
  assert(channel < numOutputChannels());
  assert(beam < num_beams);
  int i = index4d(time, pol, num_polarizations, channel, numOutputChannels(),
                  beam, num_beams);
  return buffer->data[i];
}

