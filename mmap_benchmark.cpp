#include <assert.h>
#include <chrono>
#include <fcntl.h>
#include <iostream>
#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

// A speed test to see how fast mmapping can be.
// run: mmap <filename>
int main(int argc, char* argv[]) {
  assert(argc >= 2);
  string filename = argv[1];
  auto begin = chrono::system_clock::now();
  int fd = open(filename.c_str(), O_RDONLY);
  auto file_size = lseek(fd, 0, SEEK_END);
  lseek(fd, 0, SEEK_SET);
  void *mmapped = mmap(NULL, file_size, PROT_READ, MAP_PRIVATE | MAP_POPULATE, fd, 0);
  int8_t *everything = reinterpret_cast<int8_t*>(mmapped);

  // Check we actually got the data
  int8_t foo = 0;
  for (long j = 0; j < file_size; j += 1000) {
    foo ^= everything[j];
  }

  auto duration = chrono::system_clock::now() - begin;
  auto ms = chrono::duration_cast<chrono::milliseconds>(duration).count();
  
  munmap(mmapped, file_size);  
  
  cout << "file size: " << file_size << endl;
  cout << "checksum: " << (int) foo << endl;
  cout << "time elapsed: " << ms << " ms" << endl;
  double gb = (double) file_size / 1024.0 / 1024.0 / 1024.0;
  double s = ms / 1000.0;
  cout << "GB/s: " << (gb/s) << endl;
}
