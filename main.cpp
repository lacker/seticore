#include <iostream>

using namespace std;

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: seticore <h5file>" << endl;
    return 1;
  }

  auto filename = string(argv[1]);
  cout << "argument is: " << filename << endl;
  return 0;
}
