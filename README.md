# seticore
A high-performance implementation of some core SETI algorithms, in particular turboSETI's
dedoppler algorithm that creates .dat files from .h5 files, and a beamforming algorithm
to produce filterbanks from .raw and .bfr5 files.

## Quickstart

Install dependencies:

```
sudo apt-get install cmake g++ libboost-all-dev libhdf5-dev pkg-config
```

You will also need a Python environment for your build tools:

```
python3 -m pip install meson
python3 -m pip install ninja
```

Git clone this repo, then update submodules, then run the make scripts:

```
git submodule init
git submodule update
meson setup build
cd build
meson compile
```

then run on an h5 file:

```
./seticore /path/to/your.h5
```

It will behave roughly like turboseti. If something doesn't immediately work, try the more
detailed instructions below.

## Fixing hdf5 plugin errors

Depending on how you installed hdf5, you may not have the plugins that you need, in particular
the bitshuffle plugin. If seticore
tells you that you need to make sure that plugin files are in the plugin directory, check if
that directory exists. If not, we'll have to install them.

This is not straightforward because the hdf5 group does not support plugins on Ubuntu, and
the bitshuffle plugin developer only supports Python. So we use the Python install and copy
its files over.

* Create the plugin directory that it expects
* `pip install hdf5plugin`
* Find the plugin file location: `pip show hdf5plugin | grep Location`
* Look around for a directory full of `.so` files. Maybe it's in `$LOCATION/hdf5plugin/plugins`
* Copy all of those `.so` files to the directory you created

We don't need the `hdf5plugin` install any more once its files have been copied out.

## Testing on your machine

You can run the unit tests with:

```
./unit_tests.sh
```

This will automatically download a few gigabytes of test data from Berkeley and ensure
you get the expected output.

The main supported target environment is defined by [this
Dockerfile](https://github.com/lacker/seticore/blob/master/Dockerfile). If
you have the [Nvidia container
toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
set up to let Docker access the GPU, you can test the code
using this Docker environment with:

```
./test_docker.sh
```

## Debugging

To get a binary with debug symbols, so you can run gdb:

```
meson configure --buildtype=debug build
```

When you're done, configure it back:

```
meson configure --buildtype=release build
```

## Selecting a Specific GPU Device ID

Reference: https://shawnliu.me/post/nvidia-gpu-id-enumeration-in-linux/

It's easy to screw this up because `nvidia-smi` doesn't give GPUs the same numbering as the CUDA
library does by default. To use the nth GPU, where you've determined `n` by looking at `nvidia-smi`:

```
CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=$n seticore <args>
```

This makes only the nth GPU visible to seticore. This same procedure works for any binary that
uses CUDA, not just seticore.

## Profiling

Install some dependencies:

```
sudo apt-get install perf
pushd ~
git clone git@github.com:brendangregg/FlameGraph.git
popd
```

Then do a profiling run with your arguments of choice:

```
perf record -g ./build/seticore <arguments-for-seticore>
perf script > perf.script
~/FlameGraph/stackcollapse-perf.pl perf.script | ~/FlameGraph/flamegraph.pl > perf.svg
```

Then look at `perf.svg` in a browser.

## Viewer

A Python notebook is provided to show how to analyze the outputs.
Set up Python 3.9 in whichever way you prefer (such as a fresh conda env) and then:

```
pip install -r requirements.txt
jupyter-lab example.ipynb
```

The notebook demonstrates some basic reading of hit files and stamp files.

You can also read the schemas for the [hit file](hit.capnp) and the
[stamp file](stamp.capnp).
