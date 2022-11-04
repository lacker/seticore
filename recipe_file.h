#pragma once

#include "beamformer.h"
#include "hdf5.h"
#include <string>
#include <thrust/complex.h>
#include <vector>

using namespace std;


/*
  This class contains helper methods for processing .h5 files that store information
  you need to run beamforming, also known as "recipe files".
  Often suffixed as ".bfr5" files.
 */
class RecipeFile {
 private:
  const hid_t file;
  
  string getStringData(const string& name) const;
  vector<string> getStringVectorData(const string& name) const;
  template <class T> vector<T> getVectorData(const string& name, hid_t hdf5_type) const;
  vector<thrust::complex<float> > getComplexVectorData(const string& name) const;

  // Note that these are in radians, and in most places they are stored in hours/degrees.
  // They are private to keep people from accidentally screwing up and using the
  // radians directly.
  const vector<double> ras;
  const vector<double> decs;
  
 public:
  long getLongScalarData(const string& name) const;

  const string filename;
  const string obsid;
  const vector<string> src_names;
  const vector<double> delays;
  const vector<double> time_array;

  const long nants;
  const long nbeams;
  const long nchans;
  const long npol;
  const vector<thrust::complex<float> > cal_all;

  
  RecipeFile(const string& filename);

  // filename can be a directory, in which case we open the right file within the directory
  RecipeFile(const string& filename, const string& obsid);

  ~RecipeFile();

  int getTimeArrayIndex(double time) const;
  
  double getDelay(int time_array_index, int beam, int antenna) const;

  thrust::complex<float> getCal(int frequency, int polarization, int antenna) const;

  // Logs information and exits if this recipe file does not match the given raw
  // file parameters.
  void validateRawRange(int schan, int num_coarse_channels) const;
  
  void generateCoefficients(int time_array_index,
			    int raw_start_channel,
			    int raw_num_channels,
			    float raw_center_mhz,
			    float raw_bandwidth_mhz,
			    int subband_start,
			    int subband_size,
                            Beamformer* beamformer) const;

  vector<double> getRAsInHours() const;
  vector<double> getDecsInDegrees() const;
};
