#!/bin/bash -e
# This script builds seticore in the "Breakthrough Listen" environment, where
# cmake and nvcc are not on the system path, so we explicitly locate them.

# To make this work you have to install your own cmake in ~/cmake - the default BL one is
# too old. I used 3.10.2. Run the script you get from
#   wget https://github.com/Kitware/CMake/releases/download/v3.10.2/cmake-3.10.2-Linux-x86_64.sh

CMAKE_BIN=~/cmake/bin/cmake

if [ ! -f "$CMAKE_BIN" ]; then
    echo "this script expects you to put your own version of cmake in $CMAKE_BIN"
    echo "read bl_make.sh for more details"
    exit 1
fi

cd build

# Configure
CUDACXX=/usr/local/cuda/bin/nvcc $CMAKE_BIN ..

# Compile & link
$CMAKE_BIN --build .
