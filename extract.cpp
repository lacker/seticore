#include <boost/program_options.hpp>
#include <fmt/core.h>
#include <iostream>
#include "raw_file_group.h"
#include "upchannelizer.h"

using namespace std;

namespace po = boost::program_options;

/*
  Extract a single "postage stamp" from a group of raw files.
 */
int main(int argc, char* argv[]) {
  po::options_description desc("extract options");
  desc.add_options()

    ("raw_prefix", po::value<string>(),
     "prefix for the raw files to extract a stamp from")

    ("output", po::value<string>(),
     "location to write the stamp file")

    ("band", po::value<int>(),
     "parameter for raw-file-group reading")

    ("num_bands", po::value<int>(),
     "parameter for raw-file-group reading")

    ("coarse_channel", po::value<int>(),
     "coarse channel to extract data from (within the band)")

    ("fft_size", po::value<int>(),
     "size of the fft to use for upchannelization")

    ("start_channel", po::value<int>(),
     "channel to start extraction after fft")

    ("num_channels", po::value<int>(),
     "how many post-fft channels to extract")

    ;

  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  string raw_prefix = vm["raw_prefix"].as<string>();
  string output = vm["output"].as<string>();
  // int band = vm["band"].as<int>();
  int num_bands = vm["num_bands"].as<int>();
  // int coarse_channel = vm["coarse_channel"].as<int>();
  int fft_size = vm["fft_size"].as<int>();
  // int start_channel = vm["start_channel"].as<int>();
  // int num_channels = vm["num_channels"].as<int>();
  
  cout << fmt::format("extracting stamp from {}.*.raw\n", raw_prefix);
  vector<string> raw_files = getFilesMatchingPrefix(raw_prefix);
  cout << fmt::format("reading {} raw files:\n", raw_files.size());
  for (auto f : raw_files) {
    cout << "  " << f << endl;
  }

  RawFileGroup file_group(raw_files, num_bands);

  // One of the fft size and block duration should divide the other
  int blocks_per_batch;
  if (file_group.timesteps_per_block >= fft_size) {
    assert(file_group.timesteps_per_block % fft_size == 0);
    blocks_per_batch = 1;
  } else {
    assert(fft_size % file_group.timesteps_per_block == 0);
    blocks_per_batch = fft_size / file_group.timesteps_per_block;
  }

  Upchannelizer upchannelizer(0, fft_size,
                              file_group.timesteps_per_block * blocks_per_batch,
                              file_group.num_coarse_channels,
                              file_group.npol,
                              file_group.nants);

  ComplexBuffer internal(upchannelizer.requiredInternalBufferSize());

  MultiantennaBuffer multiantenna(upchannelizer.numOutputTimesteps(),
                                  upchannelizer.numOutputChannels(),
                                  file_group.npol,
                                  file_group.nants);

  
}
