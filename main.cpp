#include <boost/program_options.hpp>
#include <iostream>
#include <string>

#include "dedoppler.h"

using namespace std;

namespace po = boost::program_options;

int main(int argc, char* argv[]) {
  po::options_description desc("seticore options");
  desc.add_options()
    ("input", po::value<string>(), "the input .h5 file")
    ;

  po::positional_options_description p;
  p.add("input", -1);
  
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
  po::notify(vm);
  
  if (!vm.count("input")) {
    cerr << desc << "\n";
    return 1;
  }

  string filename = vm["input"].as<string>();
  cout << "argument is: " << filename << endl;
  dedoppler(filename);
}

