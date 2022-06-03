#!/bin/bash -e

cd capnproto/c++
autoreconf -i
./configure
sudo make install
