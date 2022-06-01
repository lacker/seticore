#include <iostream>

#include "recipe_file.h"

using namespace std;

int main(int argc, char* argv[]) {
  RecipeFile recipe("data/golden_synthesized_input.bfr5");
  cout << "OK\n";
  return 0;
}
