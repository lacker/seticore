#!/usr/bin/env python

import matplotlib
matplotlib.rc("figure", max_open_warning=0)
from matplotlib import pyplot as plt

import capnp

import h5py
import hit_capnp
import math
import numpy as np
import os
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
    if arr.shape[1] > 1000:
        size = 30
    else:
        size = 15
    fig, ax = plt.subplots(figsize=(size, size))
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

    
class Recipe(object):
    def __init__(self, filename):
        self.h5 = h5py.File(filename)
        self.ras = self.h5["/beaminfo/ras"][()]
        self.decs = self.h5["/beaminfo/decs"][()]
        self.obsid = self.h5["/obsinfo/obsid"][()]
        self.src_names = self.h5["/beaminfo/src_names"][()]
        self.delays = self.h5["/delayinfo/delays"][()]
        self.time_array = self.h5["/delayinfo/time_array"][()]
        self.npol = self.h5["/diminfo/npol"][()]
        self.nbeams = self.h5["/diminfo/nbeams"][()]
        self.cal_all = self.h5["/calinfo/cal_all"][()]
        self.nants = self.h5["/diminfo/nants"][()]
        self.nchan = self.h5["/diminfo/nchan"][()]
        
        # Validate shapes of things
        assert self.delays.shape == (len(self.time_array), self.nbeams, self.nants)
        assert self.cal_all.shape == (self.nchan, self.npol, self.nants)
        
        
    def time_array_index(self, time):
        """Return the index in time_array closest to time."""
        dist_tuples = [(i, abs(val - time)) for i, val in enumerate(self.time_array)]
        i, _ = min(dist_tuples)
        return i
    
    
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
        if antennas.shape[1] > antennas.shape[0]:
            cols = 4
        else:
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
        plt.show()
        plt.close()

    def times(self):
        """Returns a list of the times for each timestep."""
        return [self.stamp.tstart + n * self.stamp.tsamp
                for n in range(self.stamp.numTimesteps)]

    def frequencies(self):
        """Returns a list of the frequencies in MHz for each channel."""
        return [self.stamp.fch1 + n * self.stamp.foff
                for n in range(self.stamp.numChannels)]

    def coefficients(self, recipe, beam):
        """Returns a numpy array of beamforming coefficients.

        This does not conjugate, so it should match up with c++ generateCoefficients.
        Output dimensions are [time, chan, pol, ant]
        """
        recipe_channel_index = self.stamp.schan + self.stamp.coarseChannel

        answer = np.zeros((self.stamp.numTimesteps,
                           self.stamp.numChannels,
                           self.stamp.numPolarities,
                           self.stamp.numAntennas),
                          dtype=np.cdouble)

        for timestep, time_value in enumerate(self.times()):
            for chan, freq_value in enumerate(self.frequencies()):
                time_array_index = recipe.time_array_index(time_value)
                ghz = freq_value * 0.001
                tau = recipe.delays[time_array_index, beam, :]
                angles = tau * (ghz * 2 * np.pi * 1.0j)
                for pol in range(self.stamp.numPolarities):
                    cal = recipe.cal_all[recipe_channel_index, pol, :]
                    answer[timestep, chan, pol, :] = cal * np.exp(angles)

        return answer

    def beamform_voltage(self, recipe, beam):
        """Beamforms, leaving the result in complex voltage space.
        
        Output dimensions are [time, chan]
        """
        coeffs = self.coefficients(recipe, beam)
        inputs = self.complex_array()

        # Sum along polarity and antenna
        return (np.conjugate(coeffs) * inputs).sum(axis=(2, 3))

    def beamform_power(self, recipe, beam):
        voltage = self.beamform_voltage(recipe, beam)
        return np.square(np.real(voltage)) + np.square(np.imag(voltage))

    
def read_stamps(filename):
    with open(filename) as f:
        stamps = stamp_capnp.Stamp.read_multiple(f, traversal_limit_in_words=2**30)
        for s in stamps:
            yield Stamp(s)
    
def main():
    for hit in read_hits("data/voyager.hits"):
        print(hit.filterbank.numChannels, "x", hit.filterbank.numTimesteps,
              "=", len(hit.filterbank.data))

def find_stamp_files(directory):
    "Find all stamp files under this directory."
    for root, dirs, files in os.walk(directory, topdown=False):
        for f in files:
            full = os.path.join(root, f)
            if full.endswith(".stamps"):
                yield full
                
def scan_dir(directory):
    """
    Display antenna data for each stamp file under this directory.
    """
    count = 0
    for stamp_filename in find_stamp_files(directory):
        try:
            for (i, stamp) in enumerate(read_stamps(stamp_filename)):
                print(f"stamp {i} from {stamp_filename}")
                stamp.show_antennas()
                count += 1
        except Exception as e:
            print(f"error opening {stamp_filename}: {e}")
    print(count, "stamps shown total")
        
if __name__ == "__main__":
    main()
