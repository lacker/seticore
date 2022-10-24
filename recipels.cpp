#include <fmt/core.h>
#include <iostream>

#include "recipe_file.h"
#include "util.h"

using namespace std;

/*
  Displays information about a recipe file.

  Usage:
    recipels <filename>
*/
int main(int argc, char* argv[]) {
  if (argc != 2) {
    cerr << "usage: recipels <filename>\n";
    exit(1);
  }
  string filename(argv[1]);

  cout << "loading recipe from " << filename << endl;
  RecipeFile recipe(filename);
  cout << "obsid: " << recipe.obsid << endl;
  cout << "npol: " << recipe.npol << endl;
  cout << "nbeams: " << recipe.nbeams << endl;
  cout << "delays size: " << recipe.delays.size() << endl;
  cout << "cal_all size: " << recipe.cal_all.size() << endl;
  cout << "time_array size: " << recipe.time_array.size() << endl;
  cout << "inferred nants: " << recipe.nants << endl;
  cout << "inferred nchans: " << recipe.nchans << endl;
  cout << "/diminfo/nchans = " << recipe.getLongScalarData("/diminfo/nchans") << endl;
}
