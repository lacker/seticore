# How to analyze seticore data products

There are currently two data formats that we keep around a long time from Meerkat.

* "hits" are everything found by dedoppler. Each hit is in a single
  beam. They contain beamformed power from a small region in which the
  hit was found, as well as a bunch of metadata about the hit.

* "stamps" are a subset of the raw data that we generate for some
  hits. We FFT the raw data to upchannelize, then select a region
  around the hit, and store the voltage information for all antennas.

Both of these are stored in protocol buffers, specifically the
professionally named [Cap'n Proto](https://capnproto.org/) data
format. I know, you're excited to add yet another data format to your
repertoire.

These data formats are defined by `.capnp` schema files.

* [hit.capnp](hit.capnp)
* [stamp.capnp](stamp.capnp)

The schema ensures that we track precisely what data format was output
by each past version of the software, and ensures backward
compatibility. Code using an older version of the schema can still
read newer files, and vice versa.

# Python instructions

You need the `pycapnp` library.

```
pip install pycapnp
```

Then, the most convenient thing to do for now is just to copy the `.capnp`
files into your working directory. If there's enough demand, I'll wrap
this up into a library that you can `pip install`.

In your Python code, you can then do:

```
import capnp
import hit_capnp
import stamp_capnp
```

`hit_capnp.Hit` and `stamp_capnp.Stamp` then contain methods for
dealing with these objects. For example, to read a file of hits
and print out the beams:

```
def read_hits(filename):
    with open(filename) as f:
        hits = hit_capnp.Hit.read_multiple(f)
        for hit in hits:
            print("the beam is:", hit.filterbank.beam)
```

More docs for pycapnp are [here](http://capnproto.github.io/pycapnp/).

# Languages that aren't Python

There are capnp libraries for many other languages, see [this
list](https://capnproto.org/otherlang.html).

# Julia

capnp doesn't have Julia support out of the box; consider using [David
M's wrapper](https://github.com/david-macmahon/SeticoreCapnp.jl).
