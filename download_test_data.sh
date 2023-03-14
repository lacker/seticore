#!/bin/bash -e

for URL in `cat ./test_data_files.txt`; do
    H5=`basename $URL`
    if [ ! -f "data/$H5" ]; then
        wget $URL -P data/
    fi
done

