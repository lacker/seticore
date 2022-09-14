#include <assert.h>
#include "beamforming_pipeline.h"
#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options.hpp>
#include <fmt/core.h>
#include <iostream>
#include "raw_file_group.h"
#include "run_dedoppler.h"
#include <string>
#include "thread_util.h"
#include <time.h>
#include "util.h"

using namespace std;

namespace po = boost::program_options;

const string VERSION = "0.2.0";

int beamformingMode(const po::variables_map& vm) {
  cout << "running in beamforming mode.\n";

  string input_dir = vm["input"].as<string>();
  string output_dir = vm["output"].as<string>();
  string recipe_dir = vm["recipe_dir"].as<string>();
  int num_bands = vm["num_bands"].as<int>();
  int fft_size = vm["fft_size"].as<int>();
  int sti = vm["sti"].as<int>();
  int telescope_id = vm["telescope_id"].as<int>();
  float snr = vm["snr"].as<double>();
  float max_drift = vm["max_drift"].as<double>();
  float min_drift = vm["min_drift"].as<double>();

  auto groups = scanForRawFileGroups(input_dir);
  cout << "found " << pluralize(groups.size(), "group") << " of raw files.\n";
  for (auto group : groups) {
    BeamformingPipeline pipeline(group, output_dir, recipe_dir, num_bands,
                                 fft_size, sti, telescope_id, snr, max_drift, min_drift);

    if (vm.count("h5_dir")) {
      pipeline.h5_dir = vm["h5_dir"].as<string>();
    }
    int tstart = time(NULL);
    pipeline.findHits();
    int tmid = time(NULL);
    cerr << fmt::format("time to find hits: {:d}s\n", tmid - tstart);
    pipeline.makeStamps();
    int tstop = time(NULL);
    cerr << fmt::format("time to make stamps: {:d}s\n", tstop - tmid);
  }
  return 0;
}

int dedopplerMode(const po::variables_map& vm) {
  cout << "running in dedoppler mode.\n";
  string input = vm["input"].as<string>();
  string output;
  if (!vm.count("output")) {
    // By default, output to a .dat file in the same location as the input file
    auto index = input.find_last_of(".");
    if (index >= input.size()) {
      cerr << "unrecognized input filename: " << input << endl;
      return 1;
    }
    output = input.substr(0, index) + ".dat";
  } else {
    output = vm["output"].as<string>();
  }

  double max_drift = vm["max_drift"].as<double>();
  double snr = vm["snr"].as<double>();
  double min_drift = vm["min_drift"].as<double>();

  cout << "loading input from " << input << endl;
  cout << fmt::format("dedoppler parameters: max_drift={:.2f} min_drift={:.4f} "
                      "snr={:.2f}\n",
                      max_drift, min_drift, snr);
  cout << "writing output to " << output << endl;
  int tstart = time(NULL);
  runDedoppler(input, output, max_drift, min_drift, snr);
  int tstop = time(NULL);
  cerr << fmt::format("dedoppler elapsed time: {:d}s\n", tstop - tstart);
  return 0;
}
  
// This method just handles command line parsing, and the real work is done
// via the dedoppler function.
int main(int argc, char* argv[]) {
  setThreadName("main");
  
  po::options_description desc("seticore options");
  desc.add_options()
    ("help,h", "produce help message")

    ("input", po::value<string>(),
     "alternate way of setting the input file or directory")

    ("output", po::value<string>(),
     "the output as .dat file, .hits file, or directory. defaults to <input>.dat")

    ("max_drift,M", po::value<double>()->default_value(10.0),
     "maximum drift in Hz/sec")

    ("min_drift,m", po::value<double>()->default_value(0.0001),
     "minimum drift in Hz/sec")

    ("snr,s", po::value<double>()->default_value(25.0),
     "minimum SNR to report a hit")

    ("recipe_dir", po::value<string>(),
     "the directory to find beamforming recipes in. set this to beamform.")

    ("h5_dir", po::value<string>(),
     "optional directory to save .h5 files containing post-beamform data")
    
    ("num_bands", po::value<int>()->default_value(1),
     "number of bands to break input into")

    ("fft_size", po::value<int>(),
     "size of the fft for upchannelization")

    ("telescope_id", po::value<int>(),
     "SIGPROC standard id for the telescope that provided this data")

    ("sti", po::value<int>()->default_value(8),
     "duration of the Short Time Integration to compress post-beamforming data")
    ;

  po::positional_options_description p;
  p.add("input", -1);
  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);
  
  if (!vm.count("input") || vm.count("help")) {
    cerr << "usage: seticore [input]\n";
    cerr << "seticore version: " << VERSION << endl;
    cerr << desc << "\n";
    return 1;
  }

  cout << "welcome to seticore, version " << VERSION << endl;

  if (vm.count("recipe_dir")) {
    return beamformingMode(vm);
  }
  
  return dedopplerMode(vm);
}

