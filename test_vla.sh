#!/bin/bash -e

FFT_SIZE=16384
INPUT_DIR="/home/lacker/test"
RECIPE_DIR="/home/lacker/test"
NUM_BANDS=1
OUTPUT_DIR="/home/lacker/output"
EXTRA_FLAGS=

/home/lacker/seticore/build/seticore \
    --input $INPUT_DIR \
    --output $OUTPUT_DIR \
    --max_drift 10.0 \
    --min_drift 0.0 \
    --snr 10.0 \
    --recipe_dir $RECIPE_DIR \
    --num_bands $NUM_BANDS \
    --fft_size $FFT_SIZE \
    --telescope_id 12 \
    $EXTRA_FLAGS
    
