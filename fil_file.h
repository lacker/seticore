#pragma once

#include <fstream>
#include <iostream>

#include "filterbank_file.h"
#include "fil_file.h"

using namespace std;

/*
  This class reads in sigproc filterbank files, typically ending in the .fil suffix.
 */
class FilFile: public FilterbankFile {
 private:
  ifstream file;

  streampos data_start;
  
  template <class T> T readBasic();
  string readString();
  
 public:
  FilFile(const string& filename);
  ~FilFile();

  void loadCoarseChannel(int i, float* output) const;
};
