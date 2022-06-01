#!/bin/bash -e

H5=blc17_guppi_59544_62191_HIP99317_0059.rawspec.0000.h5
DATA_URL="https://bldata.berkeley.edu/pipeline/AGBT21B_999_31/blc17_blp17/$H5"

if [ ! -f "data/$H5" ]; then
    echo downloading h5 data for regression testing...
    wget $DATA_URL -P data/
fi

./build/seticore data/$H5 --max_drift=0.4 --snr=10 --min_drift=0 | tee data/output.txt
echo diffing against expected output.
diff data/output.txt data/golden.txt | tee data/diff.txt

if [ -s data/diff.txt ]; then
    echo output did not match. either fix the bug, or if the new behavior is correct,
    echo update the golden file with:
    echo   cp data/output.txt data/golden.txt
else
    echo output matches. regression test looks good
fi
