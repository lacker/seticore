"""
A tool to check that the .h5 output of seticore and the .fil output of hpguppi_proc match.

Usage:
  ./crosscheck.py <dir>
or:
  ./crosscheck.py foo.h5 foo.fil

The directory probably looks like:
  ./crosscheck.py /scratch/data/20220915/0013/
"""

from blimpy import Waterfall
import h5py
import numpy as np
import os
import sys

def force_decode(s):
    if hasattr(s, "decode"):
        return s.decode()
    return s

def approx_equal(a, b):
    if type(a) is float:
        return abs(a - b) / abs(a) < 0.01
    return force_decode(a) == force_decode(b)


def check_files(h5_filename, fil_filename):
    print("comparing", h5_filename, "to", fil_filename)
    h5 = h5py.File(h5_filename, mode="r")
    data = h5["data"]
    fil = Waterfall(fil_filename)

    # Check metadata. (tstart is kind of busted still)
    mismatches = 0
    for key in ["fch1", "foff", "source_name", "tsamp"]:
        h5_val = data.attrs[key]
        fil_val = fil.header[key]
        if not approx_equal(h5_val, fil_val):
            print(f"mismatch for {key}:")
            print(f"h5 value: {h5_val}")
            print(f"fil value: {fil_val}")
            return False

    fil_rows = fil.data.shape[0]
    h5_rows = data.shape[0]

    # Allow the implementations to differ on the last row, ie truncation policy
    if abs(fil_rows - h5_rows) > 1:
        print("shape mismatch:")
        print("fil shape:", fil.data.shape)
        print("h5 shape:", data.shape)
        return False

    # Check actual data
    for row in range(min(fil_rows, h5_rows)):
        fil_row = fil.data[row][0]
        h5_row = data[row][0]
        if len(fil_row) != len(h5_row):
            print(f"row {row} length mismatch:")
            print(f"{len(fil_row)=}")
            print(f"{len(h5_row)=}")
            return False
        
        diffs = (fil_row - h5_row) / fil_row.mean()

        # Allow 1% disagreement for numerical instability
        diff_count = np.count_nonzero(diffs > 0.01)
        if diff_count > 0:
            print(f"row {row} diff: {diff_count}")
            return False
        else:
            print(f"row {row} ok")
    
    return True


def guppi_prefix(fname):
    return "_".join(fname.split("_")[:4])


def chop_word(s, word):
    assert s.startswith(word)
    return int(s[len(word):])


def check_dir(dirname):
    # Keyed by (guppi prefix, subband, beam) tuples
    h5map = {}
    filmap = {}

    h5_dir = os.path.join(dirname, "seticore_beamformer")
    for basename in os.listdir(h5_dir):
        if not basename.endswith(".h5"):
            continue
        prefix = guppi_prefix(basename)
        s_band, s_beam, _ = basename.split(".")[-3:]
        band = chop_word(s_band, "band")
        if s_beam == "incoherent":
            beam = 64
        else:
            beam = chop_word(s_beam, "beam")
        h5map[(prefix, band, beam)] = os.path.join(h5_dir, basename)
    
    fil_dir = os.path.join(dirname, "hpguppi_beamformer")
    for basename in os.listdir(fil_dir):
        if not basename.endswith(".fil"):
            continue
        prefix = guppi_prefix(basename)
        s_band, s_beam, _ = basename.split(".")[-3:]
        band = chop_word(s_band, "SB")
        beam = chop_word(s_beam, "B")
        filmap[(prefix, band, beam)] = os.path.join(fil_dir, basename)

    keys = sorted(set(h5map.keys()).intersection(filmap.keys()))
    for key in keys:
        check_files(h5map[key], filmap[key])
        
    
def main():
    if len(sys.argv) == 2:
        check_dir(sys.argv[1])
    elif len(sys.argv) == 3:
        check_files(sys.argv[1], sys.argv[2])
    else:
        print("usage: ./crosscheck.py <dir> or ./crosscheck.py foo.h5 foo.fil")
