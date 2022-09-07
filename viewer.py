#!/usr/bin/env python

import matplotlib
matplotlib.rc("figure", max_open_warning=0)
from matplotlib import pyplot as plt

import capnp
import hit_capnp
import numpy as np
import stamp_capnp

def read_hits(filename):
    with open(filename) as f:
        hits = hit_capnp.Hit.read_multiple(f)
        for hit in hits:
            yield hit

def beam_name(hit):
    n = hit.filterbank.beam
    if n < 0:
        return "incoherent beam"
    return f"beam {n}"

def show_array(arr):
    fig, ax = plt.subplots(figsize=arr.shape)
    ax.imshow(arr, rasterized=True, interpolation="nearest", cmap="viridis")
    display(fig)
    plt.close()
    
def show_hit(hit):
    fb = hit.filterbank
    data = np.array(fb.data).reshape((fb.numTimesteps, fb.numChannels))
    print(f"hit with source {fb.sourceName}, {beam_name(hit)}, " +
          f"{hit.signal.frequency:.5f} MHz, " +
          f"{hit.signal.snr:.1f} SNR, {hit.signal.driftRate:.3f} Hz/s drift:")
    show_array(data)

    
class Stamp(object):
    def __init__(self, filename):
        """
        self.stamp stores the proto data.
        """
        with open(filename) as f:
            stamps = stamp_capnp.Stamp.read_multiple(f)
            self.stamp = list(stamps)[0]

    def real_array(self):
        dimensions = (self.stamp.numTimesteps,
                      self.stamp.numChannels,
                      self.stamp.numPolarities,
                      self.stamp.numAntennas,
                      2)
        return np.array(self.stamp.data).reshape(dimensions)

    def complex_array(self):
        real = self.real_array()
        return real[:, :, :, :, 0] + 1.0j * real[:, :, :, :, 1]

    
def main():
    for hit in read_hits("data/voyager.hits"):
        print(hit.filterbank.numChannels, "x", hit.filterbank.numTimesteps,
              "=", len(hit.filterbank.data))
    
if __name__ == "__main__":
    main()
