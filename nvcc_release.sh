#!/bin/bash -e
# Prints the nvcc release as major.minor
nvcc --version | grep release | sed 's:.*release ::' | sed 's:,.*::'
