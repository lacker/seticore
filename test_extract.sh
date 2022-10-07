#!/bin/bash -e

# Extracting a hit from blpn129

STAMP=data/J1939.stamp

rm -f $STAMP

./build/extract \
    --raw_prefix=../data/beamtest/guppi_59712_16307_003760_J1939-6342_0001 \
    --output=$STAMP \
    --coarse_channel=1 \
    --fft_size=131072 \
    --start_channel=86100 \
    --num_channels=200 \
    --telescope_id=64

if [[ `./build/stampls data/J1939.stamp | grep 1163` == "" ]]; then
    echo bad stamp file
    exit 1
fi

echo OK
