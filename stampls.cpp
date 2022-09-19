#include <capnp/message.h>
#include <capnp/serialize.h>
#include <errno.h>
#include <fcntl.h>
#include <fmt/core.h>
#include <iostream>
#include "stamp.capnp.h"
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
    cerr << "usage: stamps <filename>\n";
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
    cout << "tstart: " << stamp.getTstart() << endl;
    cout << "tsamp: " << stamp.getTsamp() << endl;
    cout << "telescope: " << stamp.getTelescopeId() << endl;
    cout << "timesteps: " << stamp.getNumTimesteps() << endl;
    cout << "channels: " << stamp.getNumChannels() << endl;
    cout << "polarities: " << stamp.getNumPolarities() << endl;
    cout << "antennas: " << stamp.getNumAntennas() << endl;
    cout << endl;
    ++count;
  }
  
  close(fd);
  cout << pluralize(count, "stamp") << " in the file\n";
}
