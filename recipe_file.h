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
  
 public:
  const string obsid;
  const vector<string> src_names;
  const vector<double> delays;
  const vector<double> time_array;
  const vector<double> ras;
  const vector<double> decs;
  const long npol;
  const long nbeams;
  const vector<thrust::complex<float> > cal_all;
  
  RecipeFile(const string& filename);
  ~RecipeFile();
};
