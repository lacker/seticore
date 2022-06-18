#!/bin/bash -e

H5=blc17_guppi_59544_62191_HIP99317_0059.rawspec.0000.h5
DATA_URL="https://bldata.berkeley.edu/pipeline/AGBT21B_999_31/blc17_blp17/$H5"

if [ ! -f "data/$H5" ]; then
    echo downloading h5 data for regression testing...
    wget $DATA_URL -P data/
fi
