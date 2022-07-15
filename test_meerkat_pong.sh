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

# The full five minutes
# INPUT_DIR="/buf0/20220512/0009/Unknown/GUPPI/"

# Just one target
INPUT_DIR="/home/lacker/benchmark"

OUTPUT_DIR="/mydatag/test/lacker/"
RECIPE_DIR="/home/davidm/bfr5.test/"

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
    
