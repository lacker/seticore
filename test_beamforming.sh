#!/bin/bash

RAW=golden_synthesized_input.0000.raw
BFR=golden_synthesized_input.bfr5

if [ ! -d "data/$RAW" ]; then
    echo downloading raw data for regression testing...
    wget https://bldata.berkeley.edu/pipeline/testdata/$RAW -P data/
fi

if [ ! -d "data/$BFR" ]; then
    echo downloading beamformer recipe data for regression testing...
    wget https://bldata.berkeley.edu/pipeline/testdata/$BFR -P data/
fi

echo TODO: actually test beamforming
