# This dockerfile relies on the Nvidia container toolkit to use your GPU from Docker.
# See:
#   https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html

FROM nvidia/cuda:11.7.0-devel-ubuntu22.04

RUN apt-get update -y && apt-get install -y \
        bitshuffle \
        cmake \
        git \
        libboost-all-dev \
        libhdf5-dev \
        ninja-build \
        pkg-config \
        python3-pip

RUN pip install meson

COPY . /seticore

WORKDIR /seticore
RUN meson setup build
RUN meson compile -C build
