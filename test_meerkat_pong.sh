#!/bin/bash -e
# Gathers a performance profile.
#
# Make sure the machine is set up for profiling. Run:
#   nsys status --environment
# and it should say the environment is "OK".
#
# In particular, this should be 1:
#   cat /proc/sys/kernel/perf_event_paranoid
# If not, run:
#   sudo sh -c 'echo 1 >/proc/sys/kernel/perf_event_paranoid'

# The fft sizes for different modes
FFT1K=524288
FFT4K=131072
FFT32K=16384

OUTPUT_DIR="/mydatag/test/lacker/"

# Testing the Sep 28 dataset on blpn130
if true; then
    INPUT_DIR="/buf0ro/lacker/small"
    RECIPE_DIR="/home/obs/bfr5/"
    FFT_SIZE=$FFT4K
    NUM_BANDS=16
    OUTPUT_DIR="/scratch/data/lacker/seticore/"
    EXTRA_FLAGS="--h5_dir $OUTPUT_DIR"
fi

# Testing the August 1 dataset on blpn131
if false; then
    INPUT_DIR="/buf0ro/20220801/0055/Unknown/GUPPI/"
    RECIPE_DIR="/home/obs/bfr5/"
    FFT_SIZE=$FFT4K
    NUM_BANDS=16
fi

# Testing the August 16 dataset on blpn129
if false; then
    INPUT_DIR="/buf0/20220816/"
    RECIPE_DIR="/home/obs/bfr5/"
    FFT_SIZE=$FFT4K
    NUM_BANDS=16
fi
    


# Remove old reports
rm -f *qdrep

# Use nvidia profiling tools
# nsys profile -t nvtx,cuda \

/home/lacker/seticore/build/seticore \
    --input $INPUT_DIR \
    --output $OUTPUT_DIR \
    --max_drift 10.0 \
    --min_drift 0.0 \
    --snr 10.0 \
    --recipe_dir $RECIPE_DIR \
    --num_bands $NUM_BANDS \
    --fft_size $FFT_SIZE \
    --telescope_id 64 \
    $EXTRA_FLAGS
    
