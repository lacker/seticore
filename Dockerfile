# This dockerfile relies on the Nvidia container toolkit to use your GPU from Docker.
# See:
#   https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html
#
# It would be nice to target Ubuntu 22, but first we need the standard host environment
# to have at least cuda 11.7. See here for options:
#   https://hub.docker.com/r/nvidia/cuda

FROM nvidia/cuda:11.0.3-devel-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y \
        cmake \
        git \
        libboost-all-dev \
        libhdf5-dev \
        ninja-build \
        pkg-config \
        python3-pip

RUN pip install \
        hdf5plugin \
        meson

# This is a hack, but I don't know a better way to get the hdf5 plugins installed globally.
RUN mkdir /usr/lib/x86_64-linux-gnu/hdf5/plugins
RUN cp /usr/local/lib/python3.8/dist-packages/hdf5plugin/plugins/* /usr/lib/x86_64-linux-gnu/hdf5/plugins/

COPY . /seticore

WORKDIR /seticore
RUN meson setup build
RUN meson compile -C build
