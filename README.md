# seticore
A high-performance implementation of some core SETI algorithms, in particular turboSETI's
dedoppler algorithm that creates .dat files from .h5 files.

## Quickstart

You need g++, cmake 3.10 or above, a CUDA toolkit, and hdf5 with the bitshuffle
plugin installed. `nvcc` should be in your `$PATH`.

Git clone this repo, then update submodules and run the make script:

```
git submodule init
git submodule update
./make.sh
```

then run on an h5 file:

```
./build/seticore /path/to/your.h5
```

It will behave roughly like turboseti. If something doesn't immediately work, try the more
detailed instructions below.

## More Detailed Installation Instructions

The build is using CMake. I've been developing on Ubuntu 18. You need
to install the CUDA dev kit, go follow Nvidia's instructions for that.

You're probably better off deactivating any conda environment before installing. Conda sets
up the environment to use its own copies of `hdf5` and CUDA, and it's simpler to just
use the global ones. If you get errors around multiple copies of the `hdf5` library found, or
notice paths from your conda installation sneaking into the output, this is likely to be
the source.

Then you need to install `hdf5` along with the `bitshuffle` plugin. Two
paths you can take here - what I did, or what I recommend.

## What I did

I manually installed the `hdf5` library:

* Download according to the instructions
[here](https://portal.hdfgroup.org/display/support/HDF5+1.12.1#files)
* Configure and make install according to the instructions
[here](https://github.com/mokus0/hdf5/blob/master/release_docs/INSTALL)
under "Quick installation"
* add `/usr/local/hdf5` to `PATH`

This, however, doesn't have the plugins, and HDF5 doesn't provide
binary plugins for Ubuntu 18 or any easy way to compile them yourself, as far as I
can tell. So, I installed the binary plugins into a `conda`
environment and then copied the contents of the `plugin` directory
into a new directory in `/usr/local/hdf5/lib/plugin/`.

```
conda create --name hdf5
conda activate hdf5
conda install -c conda-forge hdf5plugin
conda deactivate
```

This doesn't seem like a great way to install hdf5 but it's the first
way I got it working.

## What I recommend

If I were doing it again on Ubuntu 18, I'd try installing the `libhdf5-dev` Ubuntu package,
and then finding wherever it put the `.so` files, and then doing the
part where you copy the plugin directory from a `conda`
environment. The plugin directory should be in the same directory that
`libhdf5.so` is in, and you probably have to create the plugin
directory. So for example when I do

```
$ locate libhdf5.so
/usr/local/hdf5/lib/libhdf5.so
/usr/local/hdf5/lib/libhdf5.so.200
/usr/local/hdf5/lib/libhdf5.so.200.1.0
```

that means the plugin directory should be `/usr/local/hdf5/lib/plugin`.

On Ubuntu 20 there's a `bitshuffle` package that may just work if you
install it along with the `libhdf5-dev` package.

## Checking that it works on your machine

To check if your installation is okay:

```
./make.sh
```

You may need to set `CUDACXX` to be the path to your `nvcc` binary.

The plugin installation won't be validated until runtime.

Then grab a Green Bank data file from

```
https://bldata.berkeley.edu/pipeline/AGBT21B_999_31/blc17_blp17/blc17_guppi_59544_62191_HIP99317_0059.rawspec.0000.h5
```

and put it in the `data` directory. It's about a 3 GB file.

You can run on the sample data with:

```
./run.sh
```

To enable regression testing, run

```
cp data/output.txt data/golden.txt
```

and then subsequent calls to `run.sh` will show any variations from the golden output.

## Selecting a Specific GPU Device ID

Reference: https://shawnliu.me/post/nvidia-gpu-id-enumeration-in-linux/

When there is competition for GPU devices, it is often desirable to run CUDA programs like seticore on a specific device. This is a multi-step process to first select a specific device and then run seticore:
* ```export CUDA_DEVICE_ORDER=PCI_BUS_ID``` so that the IDs are consistent with what you will see in ```nvidia-smi```.
* Run ```nvidia-smi``` to view the current loading on the individual GPU devices.  
* ```export CUDA_VISIBLE_DEVICES=<n>``` where ```<n>``` is the ID of the specific GPU device available for use with seticore.
* Run seticore.

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
