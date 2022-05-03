#!/bin/bash -e
CMAKE_BIN=/home/obs/cmake/bin/cmake

cd build

# Configure
CUDACXX=/usr/local/cuda/bin/nvcc $CMAKE_BIN ..

# Compile & link
$CMAKE_BIN --build .
