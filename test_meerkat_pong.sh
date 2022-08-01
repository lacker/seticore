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

# Testing the smaller May 12 dataset on blpn129
# INPUT_DIR="/home/lacker/benchmark"
# RECIPE_DIR="/home/davidm/bfr5.test/"

# Testing the July 19 dataset on blpn130
INPUT_DIR="/buf0ro/20220719/0026/Unknown/GUPPI/"
RECIPE_DIR="/home/obs/bfr5/"

OUTPUT_DIR="/mydatag/test/lacker/"


# Remove old reports
rm -f *qdrep

# Use nvidia profiling tools
nsys profile \
/home/lacker/seticore/build/seticore \
    --input $INPUT_DIR \
    --output $OUTPUT_DIR \
    --max_drift 10.0 \
    --snr 10.0 \
    --recipe_dir $RECIPE_DIR \
    --num_bands 16 \
    --fft_size 131072 \
    --telescope_id 64
    
