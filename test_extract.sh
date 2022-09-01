#!/bin/bash -e

# Try to extract this hit:
#
# found 2 hits in beam 42, coarse channel 0
#   index = 86216, drift steps = 1, snr = 12.3893, drift rate = 0.317774
#   <snipped>

./build/extract \
    --raw_prefix=../data/beamtest/guppi_59712_16307_003760_J1939-6342_0001 \
    --output=data/test.stamp \
    --band=0 \
    --num_bands=32 \
    --coarse_channel=0 \
    --fft_size=131072 \
    --start_channel=86100 \
    --num_channels=300


echo done
