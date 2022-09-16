#pragma once

#include <fstream>
#include <iostream>

#include "filterbank_file_reader.h"

using namespace std;

double convertFromSigprocRaOrDec(double sigproc);

/*
  This class reads in sigproc filterbank files, typically ending in the .fil suffix.
 */
class FilReader: public FilterbankFileReader {
 private:
  ifstream file;

  streampos data_start;
  
  template <class T> T readBasic();
  string readString();
  
 public:
  FilReader(const string& filename);
  ~FilReader();

  void loadCoarseChannel(int i, FilterbankBuffer* buffer) const;
};
