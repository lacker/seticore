# seticore
A high-performance implementation of some core SETI algorithms, in particular turboSETI's
dedoppler algorithm that creates .dat files from .h5 files.

## Quickstart

Git clone this repo, run the make script

```
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
