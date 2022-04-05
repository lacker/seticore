# seticore
A high-performance implementation of some core SETI algorithms.

## Installation

Only Ubuntu is supported. You need to install Boost:

```
sudo apt-get install libboost-all-dev
```

## Running things

First grab some Green Bank test data from

```
https://bldata.berkeley.edu/pipeline/AGBT21B_999_31/blc17_blp17/blc17_guppi_59544_62191_HIP99317_0059.rawspec.0000.h5
```

and put it in the `data` directory. It's about a 3 GB file.

Then you can build and run on the sample data with:

```
./make.sh
./run.sh
```
