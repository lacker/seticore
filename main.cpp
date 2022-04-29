#include <iostream>
#include <string>

#include "dedoppler.h"

using namespace std;

int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "usage: seticore <h5file>" << endl;
    return 1;
  }

  // Open the file
  string filename = string(argv[1]);
  cout << "argument is: " << filename << endl;

  dedoppler(filename);
}

