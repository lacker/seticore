#!/bin/bash -e

FFT_SIZE=524288
INPUT_DIR="/home/lacker/mars"
RECIPE_DIR="/home/lacker/mars"
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
    --telescope_id 9 \
    $EXTRA_FLAGS
    
