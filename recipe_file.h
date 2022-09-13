#pragma once

#include <thrust/complex.h>
#include "hdf5.h"
#include <string>
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
  long getLongScalarData(const string& name) const;
  vector<thrust::complex<float> > getComplexVectorData(const string& name) const;

  // Note that these are in radians, and in most places they are stored in hours/degrees.
  // They are private to keep people from accidentally screwing up and using the
  // radians directly.
  const vector<double> ras;
  const vector<double> decs;
  
 public:
  const string obsid;
  const vector<string> src_names;
  const vector<double> delays;
  const vector<double> time_array;

  const long npol;
  const long nbeams;
  const vector<thrust::complex<float> > cal_all;

  // This data isn't provided in the recipe file, but we can infer it
  const int nants;
  const int nchans;
  
  RecipeFile(const string& filename);
  RecipeFile(const string& directory, const string& obsid);
  ~RecipeFile();

  int getTimeArrayIndex(double time) const;
  
  double getDelay(int time_array_index, int beam, int antenna) const;

  thrust::complex<float> getCal(int frequency, int polarity, int antenna) const;

  void generateCoefficients(int time_array_index,
                            int start_channel, int num_channels,
                            float center_frequency, float bandwidth,
                            thrust::complex<float>* coefficients) const;

  vector<double> getRAsInHours() const;
  vector<double> getDecsInDegrees() const;
};
