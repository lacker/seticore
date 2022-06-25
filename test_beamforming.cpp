#include <iostream>

#include "raw/raw.h"
#include "recipe_file.h"

using namespace std;

int main(int argc, char* argv[]) {
  RecipeFile recipe("data/golden_synthesized_input.bfr5");
  cout << "nants: " << recipe.nants << endl;
  cout << "nbeams: " << recipe.nbeams << endl;
  cout << "nchans: " << recipe.nchans << endl;
  cout << "ndelays: " << recipe.delays.size() << endl;
  cout << "OK\n";
  return 0;
}
