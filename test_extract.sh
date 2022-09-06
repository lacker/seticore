#!/bin/bash -e

# Extracting a hit from blpn129


./build/extract \
    --raw_prefix=/buf0ro/20220816/guppi_59806_84779_003380_J1939-6342_0001 \
    --output=/mydatag/test/lacker/20220816.stamp \
    --band=12 \
    --num_bands=16 \
    --coarse_channel=0 \
    --fft_size=131072 \
    --start_channel=38417 \
    --num_channels=80 \
    --telescope_id=64

echo done
