#include <assert.h>
#include "beamforming_pipeline.h"
#include <boost/algorithm/string.hpp>
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

int beamformingMode(const po::variables_map& vm) {
  string input_dir = vm["input"].as<string>();
  string output_dir = vm["output"].as<string>();

  string recipe_filename;
  if (vm.count("recipe")) {
    if (vm.count("recipe_dir")) {
      fatal("you cannot specify both --recipe and --recipe_dir");
    }
    recipe_filename = vm["recipe"].as<string>();
    if (!boost::algorithm::ends_with(recipe_filename, ".bfr5")) {
      fatal("expected --recipe to end with .bfr5 but it is", recipe_filename);
    }
  } else {
    if (!vm.count("recipe_dir")) {
      fatal("you must specify --recipe or --recipe_dir to beamform");
    }
    recipe_filename = vm["recipe_dir"].as<string>();
    if (boost::algorithm::ends_with(recipe_filename, ".bfr5")) {
      fatal("expected a directory for --recipe_dir but got", recipe_filename);
    }
  }

  int num_bands = vm["num_bands"].as<int>();
  int sti = vm["sti"].as<int>();
  int telescope_id = vm["telescope_id"].as<int>();
  float snr = vm["snr"].as<double>();
  float max_drift = vm["max_drift"].as<double>();
  int fft_size = vm["fft_size"].as<int>();
  int num_fine_channels = vm["fine_channels"].as<int>();
  
  if (vm.count("min_drift")) {
    cout << "the min_drift flag is ignored in beamforming mode.\n";
  }

  auto groups = scanForRawFileGroups(input_dir);
  cout << "found " << pluralize(groups.size(), "group") << " of raw files.\n";
  for (auto group : groups) {
    BeamformingPipeline pipeline(group, output_dir, recipe_filename, num_bands,
                                 sti, telescope_id, snr, max_drift, fft_size,
                                 num_fine_channels);
    if (vm.count("h5_dir")) {
      pipeline.h5_dir = vm["h5_dir"].as<string>();
    }
    int tstart = time(NULL);
    pipeline.findHits();
    int tmid = time(NULL);
    cout << fmt::format("time to find hits: {:d}s\n", tmid - tstart);
    pipeline.makeStamps();
    int tstop = time(NULL);
    cout << fmt::format("time to make stamps: {:d}s\n", tstop - tmid);
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
      fatal(fmt::format("unrecognized input filename: {}", input));
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

    ("recipe", po::value<string>(),
     "the beamforming recipe file. set this to beamform.")
    
    ("h5_dir", po::value<string>(),
     "optional directory to save .h5 files containing post-beamform data")
    
    ("num_bands", po::value<int>()->default_value(1),
     "number of bands to break input into")

    ("fft_size", po::value<int>()->default_value(-1),
     "size of the fft for upchannelization. -1 to calculate from fine_channels")

    ("fine_channels", po::value<int>()->default_value(-1),
     "how many channels to upchannelize to. -1 to calculate from fft_size")
    
    ("telescope_id", po::value<int>()->default_value(NO_TELESCOPE_ID),
     "SIGPROC standard id for the telescope. If not provided, try to infer from input.")

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

  if (vm.count("recipe_dir") || vm.count("recipe")) {
    vector<string> parts;
    parts.push_back("running seticore version");
    parts.push_back(VERSION);
    parts.push_back("with:");
    for (int i = 0; i < argc; ++i) {
      parts.push_back(string(argv[i]));
    }
    logError(boost::algorithm::join(parts, " "));
    return beamformingMode(vm);
  }
  
  cout << "welcome to seticore, version " << VERSION << endl;  
  return dedopplerMode(vm);
}

