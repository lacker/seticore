#!/bin/bash -e

INPUT_DIR="/buf0/20220512/0009/Unknown/GUPPI/"
OUTPUT_DIR="/mydatag/test/lacker/"
RECIPE_DIR="/home/davidm/bfr5.test/"

/home/lacker/seticore/build/seticore \
    --input $INPUT_DIR \
    --output $OUTPUT_DIR \
    --max_drift 10.0 \
    --snr 10.0 \
    --recipe_dir $RECIPE_DIR \
    --num_bands 16 \
    --fft_size 131072 \
    --telescope_id 64
    
