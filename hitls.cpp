#include <capnp/message.h>
#include <capnp/serialize.h>
#include <errno.h>
#include <fcntl.h>
#include <fmt/core.h>
#include <iostream>
#include "hit.capnp.h"
#include <unistd.h>
#include "util.h"

using namespace std;

/*
  Displays information about a hits file.

  Usage:
    hitls <filename>
*/
int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "usage: hitls <filename>\n";
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
    Hit::Reader hit = message.getRoot<Hit>();
    cout << "hit " << count << ":\n";
    cout << "snr " << hit.getSignal().getSnr() << endl;
    cout << "coarse channel " << hit.getFilterbank().getCoarseChannel() << endl;
    int beam = hit.getFilterbank().getBeam();
    if (beam == -1) {
      cout << "incoherent beam\n";
    } else {
      cout << "beam " << beam << endl;
    }
    cout << "starting at fine channel " << hit.getSignal().getIndex() << endl;
    cout << "drift steps: " << hit.getSignal().getDriftSteps() << endl;
    cout << endl;
    ++count;
  }

  close(fd);
  cout << pluralize(count, "hit") << " in the file\n";
}
