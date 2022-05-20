#include <assert.h>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/program_options.hpp>
#include <fmt/core.h>
#include <iostream>
#include <string>
#include <time.h>

#include "dedoppler.h"

using namespace std;

namespace po = boost::program_options;

const string VERSION = "0.0.4";

// This method just handles command line parsing, and the real work is done
// via the dedoppler function.
int main(int argc, char* argv[]) {
  po::options_description desc("seticore options");
  desc.add_options()
    ("help,h", "produce help message")
    ("input", po::value<string>(), "alternate way of setting the input .h5 file")
    ("output", po::value<string>(), "the output .dat file. if not provided, uses the input filename but replaces its suffix with .dat")
    ("max_drift,M", po::value<double>()->default_value(10.0), "maximum drift in Hz/sec")
    ("min_drift,m", po::value<double>()->default_value(0.0001), "minimum drift in Hz/sec")
    ("snr,s", po::value<double>()->default_value(25.0), "minimum SNR to report a hit")
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
    if (!boost::algorithm::ends_with(output, ".dat")) {
      cerr << "the output " << output << " is expected to be a .dat file\n";
      return 1;
    }
  }

  double max_drift = vm["max_drift"].as<double>();
  double snr = vm["snr"].as<double>();
  double min_drift = vm["min_drift"].as<double>();

  cout << "welcome to seticore, version " << VERSION << endl;
  cout << "loading input from " << input << endl;
  cout << fmt::format("dedoppler parameters: max_drift={:.2f} min_drift={:.4f} snr={:.2f}\n",
                      max_drift, min_drift, snr);
  cout << "writing output to " << output << endl;
  int tstart = time(NULL);
  dedoppler(input, output, max_drift, min_drift, snr);
  int tstop = time(NULL);
  cout << fmt::format("dedoppler elapsed time {:d}s\n", tstop - tstart);
}

