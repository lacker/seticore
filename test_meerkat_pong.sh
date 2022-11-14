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

OUTPUT_DIR="/mydatag/test/lacker/"
RECIPE_DIR="/home/obs/bfr5/"
NUM_BANDS=16

# Testing the setCoefficient failure on blpn130, Nov 14
if true; then
    INPUT_DIR="/scratch/data/danielc/20221020/0031/Unknown/GUPPI"
fi

# Testing on blpn131
# 20221103/0012 dataset, cuda illegal memory access, Nov 3
if false; then
    INPUT_DIR="/buf0ro/lacker/GUPPI/"
fi

# Testing the 20221014/0014 dataset (32k mode crash) on blpn129
if false; then
    INPUT_DIR="/buf0/lacker/20221014/"
fi
    


# Remove old reports
rm -f *qdrep

# Use nvidia profiling tools
# nsys profile -t nvtx,cuda

/home/lacker/seticore/build/seticore \
    --input $INPUT_DIR \
    --output $OUTPUT_DIR \
    --max_drift 10.0 \
    --snr 6 \
    --recipe_dir $RECIPE_DIR \
    --num_bands $NUM_BANDS \
    --fine_channels 8388608 \
    --telescope_id 64 \
    $EXTRA_FLAGS
    
