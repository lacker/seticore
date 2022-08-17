#!/usr/bin/env python

import matplotlib
matplotlib.rc("figure", max_open_warning=0)
from matplotlib import pyplot as plt

import capnp
import hit_capnp
import numpy as np

def read_hits(filename):
    with open(filename) as f:
        hits = hit_capnp.Hit.read_multiple(f)
        for hit in hits:
            yield hit
            
def show_hit(hit):
    fb = hit.filterbank
    data = np.array(fb.data).reshape((fb.numTimesteps, fb.numChannels))
    print(f"hit with source {fb.sourceName}, beam {fb.beam}, " +
          f"{hit.signal.frequency:.5f} MHz, " +
          f"{hit.signal.snr:.1f} SNR, {hit.signal.driftRate:.3f} Hz/s drift:")
    fig, ax = plt.subplots(figsize=data.shape)
    ax.imshow(data, rasterized=True, interpolation="nearest", cmap="viridis")
    display(fig)
    plt.close()
    
def main():
    for hit in read_hits("data/voyager.hits"):
        print(hit.filterbank.numChannels, "x", hit.filterbank.numTimesteps,
              "=", len(hit.filterbank.data))
    
if __name__ == "__main__":
    main()
