#!/bin/bash -e

# Extracting a hit from blpn129

OUTPUT=data/voyager.events

INPUT=\
data/single_coarse_guppi_59046_80036_DIAG_VOYAGER-1_0011.rawspec.0000.h5,\
data/single_coarse_guppi_59046_80354_DIAG_VOYAGER-1_0012.rawspec.0000.h5,\
data/single_coarse_guppi_59046_80672_DIAG_VOYAGER-1_0013.rawspec.0000.h5,\
data/single_coarse_guppi_59046_80989_DIAG_VOYAGER-1_0014.rawspec.0000.h5,\
data/single_coarse_guppi_59046_81310_DIAG_VOYAGER-1_0015.rawspec.0000.h5,\
data/single_coarse_guppi_59046_81628_DIAG_VOYAGER-1_0016.rawspec.0000.h5

rm -f $OUTPUT

./build/seticore \
    --input=$INPUT \
    --output=$OUTPUT \
    --snr=10.0

