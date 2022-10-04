#include <boost/program_options.hpp>
#include <fmt/core.h>
#include <iostream>
#include "raw_file_group.h"
#include "stamp_extractor.h"
#include "util.h"

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

    ("coarse_channel", po::value<int>(),
     "coarse channel to extract data from")

    ("fft_size", po::value<int>(),
     "size of the fft to use for upchannelization")

    ("start_channel", po::value<int>(),
     "fine channel at which to start extraction")

    ("num_channels", po::value<int>(),
     "how many fine channels to extract")

    ("telescope_id", po::value<int>(),
     "telescope id")
    
    ;

  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).run(), vm);
  po::notify(vm);

  if (!vm.count("raw_prefix") || !vm.count("output")) {
    cerr << "usage: all flags are mandatory. sorry for the inconvenience\n";
    cerr << desc << endl;
    return 1;
  }
  
  string raw_prefix = vm["raw_prefix"].as<string>();
  string output_filename = vm["output"].as<string>();
  int coarse_channel = vm["coarse_channel"].as<int>();
  int fft_size = vm["fft_size"].as<int>();
  int start_channel = vm["start_channel"].as<int>();
  int num_channels = vm["num_channels"].as<int>();
  int telescope_id = vm["telescope_id"].as<int>();
  
  cout << fmt::format("extracting stamp from {}.*.raw\n", raw_prefix);
  vector<string> raw_files = getFilesMatchingPrefix(raw_prefix);
  cout << fmt::format("reading {} raw files:\n", raw_files.size());
  for (auto f : raw_files) {
    cout << "  " << f << endl;
  }

  RawFileGroup file_group(raw_files);

  StampExtractor extractor(file_group, fft_size, telescope_id, output_filename);
  extractor.extract(coarse_channel, start_channel, num_channels);
}
