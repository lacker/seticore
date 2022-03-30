#include <iostream>

extern "C" {
  int hello_world(int n) {
    std::cout << "hello python-wrapping-c++ world! " << n << std::endl;
    return 0;
  }
}
