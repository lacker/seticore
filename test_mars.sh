#!/bin/bash -e

FFT_SIZE=524288
INPUT_DIR="/home/lacker/mars"
RECIPE="/home/lacker/mars/guppi_59856_57567_103605834_mars_0001.bfr5"
NUM_BANDS=1
OUTPUT_DIR="/home/lacker/output"
EXTRA_FLAGS=

/home/lacker/seticore/build/seticore \
    --input $INPUT_DIR \
    --output $OUTPUT_DIR \
    --max_drift 10.0 \
    --snr 10.0 \
    --recipe $RECIPE \
    --num_bands $NUM_BANDS \
    --fft_size $FFT_SIZE \
    --telescope_id 9 \
    --sti 1 \
    $EXTRA_FLAGS
    
