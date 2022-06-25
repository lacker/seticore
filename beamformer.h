#pragma once

using namespace std;

class Beamformer {
 public:

  // Number of antennas
  const int nants;

  // Number of beams to form
  const int nbeams;

  // Number of frequency channels
  const int nchans;

  // Number of polarities
  const int npol;

  // Number of timesteps
  const int nsamp;

  Beamformer(int nants, int nbeams, int nchans, int npol, int nsamp);
  ~Beamformer();

  // The input array is unified to the GPU, used to accept data from the previous step
  // in the data processing pipeline.
  //
  // Its format is row-major:
  //   input[antenna][frequency][time][polarity]
  //
  // TODO: explain what each signed char of input actually is. Some weird raw-file-format thing
  signed char *input;

  // Beamforming coefficients, formatted by row-major:
  //   coefficients[frequency][beam][polarity][antenna]
  float *coefficients;
  
 private:
  // TODO: put more buffers here
};
