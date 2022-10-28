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
import pandas as pd
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

def show_array(arr, cmap="viridis", tick_spacing=None):
    if arr.shape[1] > 1000:
        size = 30
    else:
        size = 15
    fig, ax = plt.subplots(figsize=(size, size))
    ax.imshow(arr, rasterized=True, interpolation="nearest", cmap=cmap)
    if tick_spacing is not None:
        ax.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(tick_spacing))
        ax.yaxis.set_major_locator(matplotlib.ticker.MultipleLocator(tick_spacing))
    display(fig)
    plt.close()
    
def show_hit(hit):
    fb = hit.filterbank
    data = np.array(fb.data).reshape((fb.numTimesteps, fb.numChannels))
    print(f"hit with source {fb.sourceName}, {beam_name(hit)}, " +
          f"{hit.signal.frequency:.5f} MHz, " +
          f"{hit.signal.snr:.1f} SNR, {hit.signal.driftRate:.3f} Hz/s drift:")
    show_array(data)

def show_multiple(named_waterfalls):
    """
    Show multiple waterfalls.

    named_waterfalls is a list of (name, waterfall).
    waterfall is an array indexed like [time, chan]
    """
    if not named_waterfalls:
        return
    first_waterfall = named_waterfalls[0][1]
    if first_waterfall.shape[1] > first_waterfall.shape[0]:
        cols = 4
    else:
        cols = 12
    rows = math.ceil(len(named_waterfalls) / cols)
    fig, axs = plt.subplots(rows, cols, figsize=(20, 20))
    for i in range(rows * cols):
        row = i // cols
        col = i % cols
        ax = axs[row, col]
        if i < len(named_waterfalls):
            name, waterfall = named_waterfalls[i]
            ax.set_title(name, fontsize=10)
            ax.imshow(waterfall, rasterized=True, interpolation="nearest",
                      cmap="viridis")
        else:
            ax.axis("off")

    fig.tight_layout()
    plt.show()
    plt.close()

def round_up_power_of_two(n):
    assert n >= 1
    answer = 1
    while answer < n:
        answer *= 2
    return answer
    
def interpolate_drift(total_drift, timesteps):
    """Returns a list of drifts, one for each timestep, to get to a total of total_drift.
    Must start with 0 and end with total_drift.
    """
    rounded_up = round_up_power_of_two(timesteps)
    if rounded_up > timesteps:
        return interpolate_drift(total_drift, rounded_up)[:timesteps]
        
    assert("{0:b}".format(timesteps).count("1") == 1)
    if timesteps == 2:
        return [0, total_drift]
    assert timesteps > 2
    assert timesteps % 2 == 0

    if total_drift < 0:
        return [-x for x in interpolate_drift(-total_drift, timesteps)]
    
    shift = total_drift // (timesteps - 1)
    if shift > 0:
        post_shift_drift = total_drift % (timesteps - 1)
        post_shift_interpolate = interpolate_drift(post_shift_drift, timesteps)
        return [i * shift + x for (i, x) in enumerate(post_shift_interpolate)]
    
    parity = total_drift % 2
    half_drift = total_drift // 2
    assert half_drift * 2 + parity == total_drift
    half_interpolate = interpolate_drift(half_drift, timesteps // 2)
    second_half_start = half_interpolate[-1] + parity
    return half_interpolate + [x + second_half_start for x in half_interpolate]
    
    
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
    def __init__(self, stamp, recipe=None):
        """
        self.stamp stores the proto data.
        """
        self.stamp = stamp
        self.recipe = recipe
        self._real_array = None
        
    def real_array(self):
        if self._real_array is None:
            dimensions = (self.stamp.numTimesteps,
                          self.stamp.numChannels,
                          self.stamp.numPolarizations,
                          self.stamp.numAntennas,
                          2)
            self._real_array = np.array(self.stamp.data).reshape(dimensions)
        return self._real_array

    def complex_array(self):
        real = self.real_array()
        return real[:, :, :, :, 0] + 1.0j * real[:, :, :, :, 1]

    def show_incoherent(self):
        incoherent = np.square(self.real_array()).sum(axis=(2, 3, 4))
        print("SNR:", self.snr(incoherent))
        show_array(incoherent)

    def show_antenna(self, index):
        voltages = self.real_array()[:, :, :, index, :]
        powers = np.square(voltages).sum(axis=(2, 3))
        show_array(powers)

    def show_antennas(self):
        antennas = np.square(self.real_array()).sum(axis=(2, 4))
        show_multiple([(f"antenna {i}", antennas[:, :, i])
                       for i in range(self.stamp.numAntennas)])

    def times(self):
        """Returns a list of the times for each timestep."""
        return [self.stamp.tstart + n * self.stamp.tsamp
                for n in range(self.stamp.numTimesteps)]

    def frequencies(self):
        """Returns a list of the frequencies in MHz for each channel."""
        return [self.stamp.fch1 + n * self.stamp.foff
                for n in range(self.stamp.numChannels)]

    def coefficients(self, beam):
        """Returns a numpy array of beamforming coefficients.

        This does not conjugate, so it should match up with c++ generateCoefficients.
        Output dimensions are [time, chan, pol, ant]
        """
        recipe_channel_index = self.stamp.schan + self.stamp.coarseChannel

        answer = np.zeros((self.stamp.numTimesteps,
                           self.stamp.numChannels,
                           self.stamp.numPolarizations,
                           self.stamp.numAntennas),
                          dtype=np.cdouble)

        for timestep, time_value in enumerate(self.times()):
            for chan, freq_value in enumerate(self.frequencies()):
                time_array_index = self.recipe.time_array_index(time_value)
                ghz = freq_value * 0.001
                tau = self.recipe.delays[time_array_index, beam, :]
                angles = tau * (ghz * 2 * np.pi * 1.0j)
                for pol in range(self.stamp.numPolarizations):
                    cal = self.recipe.cal_all[recipe_channel_index, pol, :]
                    answer[timestep, chan, pol, :] = cal * np.exp(angles)

        return answer

    def beamform_voltage(self, beam):
        """Beamforms, leaving the result in complex voltage space.
        
        Output dimensions are [time, chan]
        """
        coeffs = self.coefficients(beam)
        inputs = self.complex_array()

        # Sum along polarization and antenna
        return (np.conjugate(coeffs) * inputs).sum(axis=(2, 3))

    def beamform_power(self, beam):
        voltage = self.beamform_voltage(beam)
        return np.square(np.real(voltage)) + np.square(np.imag(voltage))

    def show_beam(self, beam):
        power = self.beamform_power(beam)
        print("SNR:", self.snr(power))
        show_array(power)

    def show_beams(self):
        charts = []
        for beam in range(self.recipe.nbeams):
            power = self.beamform_power(beam)
            snr = self.snr(power)
            charts.append((f"beam {beam}, snr {snr:.1f}", power))
        show_multiple(charts)

    def signal_mask(self):
        """A bool array flagging which spots are the signal we detected"""
        # We currently don't handle STI
        assert self.stamp.signal.numTimesteps * 2 > self.stamp.numTimesteps

        mask = np.zeros((self.stamp.numTimesteps, self.stamp.numChannels),
                        dtype=np.bool)
        drifts = interpolate_drift(self.stamp.signal.driftSteps,
                                   self.stamp.signal.numTimesteps)
        hit_offset = self.stamp.signal.index - self.stamp.startChannel
        for timestep, drift in enumerate(drifts):
            mask[timestep, hit_offset + drift] = True
        return mask

    def show_mask(self):
        show_array(self.signal_mask())

    def snr(self, data):
        # Calculate the noise based on the first and last 20 column sums
        left_column_sums = data[:, :20].sum(axis=1)
        right_column_sums = data[:, -20:].sum(axis=1)
        column_sums = np.concatenate((left_column_sums, right_column_sums))
        mean = column_sums.mean()
        std = column_sums.std()

        # Calculate the signal based on the mask
        signal = (data * self.signal_mask()).sum()

        return (signal - mean) / std

    def masked_antenna_values(self):
        """Returns a data[2 * time, antenna] array of complex values.
        We have two points for each timestep because we have multiple polarizations.
        """
        # time,chan,pol,ant dimension order
        raw = self.complex_array()
        mask = self.signal_mask()
        answer = np.zeros((2 * self.stamp.numTimesteps,
                           self.stamp.numAntennas),
                          dtype=np.cdouble)
        for ant in range(self.stamp.numAntennas):
            pol0 = raw[:, :, 0, ant][mask]
            pol1 = raw[:, :, 1, ant][mask]
            answer[:, ant] = np.concatenate((pol0, pol1))
        return answer

    def correlations(self):
        """Return a matrix where [i, j] is the correlation coefficient between
        antennas i and j.
        """
        vectors = self.masked_antenna_values()
        answer = np.zeros((self.stamp.numAntennas, self.stamp.numAntennas))
        for i in range(self.stamp.numAntennas):
            for j in range(i, self.stamp.numAntennas):
                vi = vectors[:, i]
                vj = vectors[:, j]
                cc = abs(np.vdot(vi, vj)) / (np.linalg.norm(vi) * np.linalg.norm(vj))
                answer[i, j] = cc
                answer[j, i] = cc
        return answer

    def show_correlations(self):
        corr = self.correlations()
        print("median correlation:", np.median(corr))
        show_array(corr, cmap="plasma", tick_spacing=2)

    def show_correlations_text(self):
        corr = self.correlations()
        nant = self.stamp.numAntennas
        print("   " + " ".join(f"{n:4d}" for n in range(nant)))
        for i in range(self.stamp.numAntennas):
            print(f"{i:2d} " + " ".join(f"{corr[i, j]:.2f}" for j in range(nant)))
        
        
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
