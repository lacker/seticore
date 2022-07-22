#!/bin/bash -e
#
# Runs most of the unit tests and regression tests.
# This doesn't run tests that require too much data to download everywhere.

echo unit testing...

./build/tests

H5=blc17_guppi_59544_62191_HIP99317_0059.rawspec.0000.h5

./download_test_data.sh

echo regression testing...

./build/seticore data/$H5 --max_drift=0.4 --snr=10 --min_drift=0 --output=data/testout.hits \
    | tee data/output.txt
echo diffing against expected output.
diff data/output.txt data/golden.txt | tee data/diff.txt

if [ -s data/diff.txt ]; then
    echo output did not match. either fix the bug, or if the new behavior is correct,
    echo update the golden file with:
    echo   cp data/output.txt data/golden.txt
    exit 1
else
    echo output matches. regression test looks good
fi
