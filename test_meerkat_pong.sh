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

# Testing the July 19 dataset on blpn130
# INPUT_DIR="/buf0ro/20220719/0026/Unknown/GUPPI/"
# RECIPE_DIR="/home/obs/bfr5/"
# FFT_SIZE=$FFT1K

# Testing the August 1 dataset on blpn131
INPUT_DIR="/buf0ro/20220801/0055/Unknown/GUPPI/"
RECIPE_DIR="/home/obs/bfr5/"
FFT_SIZE=$FFT4K

OUTPUT_DIR="/mydatag/test/lacker/"


# Remove old reports
rm -f *qdrep

# Use nvidia profiling tools
nsys profile -t nvtx,cuda \
/home/lacker/seticore/build/seticore \
    --input $INPUT_DIR \
    --output $OUTPUT_DIR \
    --max_drift 10.0 \
    --snr 10.0 \
    --recipe_dir $RECIPE_DIR \
    --num_bands 16 \
    --fft_size $FFT_SIZE \
    --telescope_id 64
    
