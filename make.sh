#!/bin/bash -e
# This script builds seticore in a "normal" environment, where cmake and nvcc are on the
# system path, so that they can be autodetected.

cd build

# Configure
cmake ..

# Actually compile & link
cmake --build .
