#!/bin/bash -e

cd build

# Configure
cmake ..

# Actually compile & link
cmake --build .
