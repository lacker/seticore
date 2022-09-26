#!/usr/bin/env python

import matplotlib
matplotlib.rc("figure", max_open_warning=0)
from matplotlib import pyplot as plt

import capnp
import hit_capnp
import math
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
    fig, ax = plt.subplots(figsize=(15, 15))
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
    def __init__(self, stamp):
        """
        self.stamp stores the proto data.
        """
        self.stamp = stamp
        self._real_array = None
        
    def real_array(self):
        if self._real_array is None:
            dimensions = (self.stamp.numTimesteps,
                          self.stamp.numChannels,
                          self.stamp.numPolarities,
                          self.stamp.numAntennas,
                          2)
            self._real_array = np.array(self.stamp.data).reshape(dimensions)
        return self._real_array

    def complex_array(self):
        real = self.real_array()
        return real[:, :, :, :, 0] + 1.0j * real[:, :, :, :, 1]

    def show_incoherent(self):
        incoherent = np.square(self.real_array()).sum(axis=(2, 3, 4))
        show_array(incoherent)

    def show_antenna(self, index):
        voltages = self.real_array()[:, :, :, index, :]
        powers = np.square(voltages).sum(axis=(2, 3))
        show_array(powers)
        
    def show_antennas(self):
        antennas = np.square(self.real_array()).sum(axis=(2, 4))
        cols = 12
        rows = math.ceil(antennas.shape[2] / cols)
        fig, axs = plt.subplots(rows, cols, figsize=(20, 20))
        for i in range(rows * cols):
            row = i // cols
            col = i % cols
            ax = axs[row, col]
            if i < antennas.shape[2]:
                ax.set_title(f"antenna {i}", fontsize=10)
                ax.imshow(antennas[:, :, i], rasterized=True, interpolation="nearest",
                          cmap="viridis")
            else:
                ax.axis("off")
            
        fig.tight_layout()

            
def read_stamps(filename):
    with open(filename) as f:
        stamps = stamp_capnp.Stamp.read_multiple(f)
        for s in stamps:
            yield Stamp(s)
    
def main():
    for hit in read_hits("data/voyager.hits"):
        print(hit.filterbank.numChannels, "x", hit.filterbank.numTimesteps,
              "=", len(hit.filterbank.data))
    
if __name__ == "__main__":
    main()
