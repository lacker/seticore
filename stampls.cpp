#include <capnp/message.h>
#include <capnp/serialize.h>
#include <errno.h>
#include <fcntl.h>
#include <fmt/core.h>
#include <iostream>
#include "stamp.capnp.h"
#include <time.h>
#include <unistd.h>
#include "util.h"

using namespace std;

/*
  Displays information about a stamp file.

  Usage:
    stampls <filename>
 */
int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "usage: stampls <filename>\n";
    exit(1);
  }

  string filename(argv[1]);

  int fd = open(filename.c_str(), O_RDONLY);
  if (fd < 0) {
    int err = errno;
    cerr << "could not open " << filename << " for reading. errno = "
         << err << endl;
    exit(1);
  }

  int count = 0;
  kj::FdInputStream fd_stream(fd);
  kj::BufferedInputStreamWrapper buffered_stream(fd_stream);
  while (buffered_stream.tryGetReadBuffer() != nullptr) {
    capnp::InputStreamMessageReader message(buffered_stream);
    Stamp::Reader stamp = message.getRoot<Stamp>();
    cout << "stamp " << count << ":\n";
    cout << "version: " << stamp.getSeticoreVersion().cStr() << endl;
    cout << "source: " << stamp.getSourceName().cStr() << endl;
    cout << "ra: " << stamp.getRa() << endl;
    cout << "dec: " << stamp.getDec() << endl;
    cout << "fch1: " << stamp.getFch1() << endl;
    cout << "foff: " << stamp.getFoff() << endl;

    double tstart = stamp.getTstart();
    long num_seconds = (long) tstart;
    double remainder = tstart - num_seconds;
    long microseconds = floor(remainder * 1000000);
    time_t time = num_seconds;
    int n = 100;
    char date_string[n];
    strftime(date_string, n, "%F %T", gmtime(&time));
    cout << "tstart: " << date_string << fmt::format(".{:06d}", microseconds) << endl;

    cout << "tsamp: " << stamp.getTsamp() << endl;
    cout << "telescope: " << stamp.getTelescopeId() << endl;
    cout << "timesteps: " << stamp.getNumTimesteps() << endl;
    cout << "channels: " << stamp.getNumChannels() << endl;
    cout << "polarizations: " << stamp.getNumPolarizations() << endl;
    cout << "antennas: " << stamp.getNumAntennas() << endl;
    cout << "coarse channel: " << stamp.getCoarseChannel() << endl;
    cout << "fft size: " << stamp.getFftSize() << endl;
    cout << "start channel: " << stamp.getStartChannel() << endl;
    cout << endl;
    ++count;
  }
  
  close(fd);
  cout << pluralize(count, "stamp") << " in the file\n";
}
