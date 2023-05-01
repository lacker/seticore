#!/usr/bin/env python

import matplotlib
matplotlib.rc("figure", max_open_warning=0)
from matplotlib import pyplot as plt

import capnp
capnp.remove_import_hook()

import h5py
import math
import numpy as np
import os
import pandas as pd

SETICORE_DIR = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
hit_capnp = capnp.load(SETICORE_DIR + "/hit.capnp")
stamp_capnp = capnp.load(SETICORE_DIR + "/stamp.capnp")

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

def plot_array(arr, cmap="viridis"):
    # TODO: decide size intelligently
    fig, ax = plt.subplots(figsize=(10, 2), dpi=300)
    ax.imshow(arr, rasterized=True, interpolation="nearest", cmap=cmap, aspect="auto")
    return fig, ax

def show_array(arr, cmap="viridis"):
    fig, ax = plot_array(arr, cmap=cmap)
    display(fig)
    plt.close()
    
def show_hit(hit):
    fb = hit.filterbank
    data = np.array(fb.data).reshape((fb.numTimesteps, fb.numChannels))
    print(f"hit with source {fb.sourceName}, {beam_name(hit)}, " +
          f"{hit.signal.frequency:.5f} MHz, " +
          f"{hit.signal.snr:.1f} SNR, {hit.signal.driftRate:.3f} Hz/s drift:")
    show_array(data)

def plot_multiple(named_waterfalls):
    """
    Plot multiple waterfalls.

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
    fig, axs = plt.subplots(rows, cols, figsize=(cols*4, rows*3), dpi=300)
    for i in range(rows * cols):
        row = i // cols
        col = i % cols
        ax = axs[row, col]
        if i < len(named_waterfalls):
            name, waterfall = named_waterfalls[i]
            ax.set_title(name, fontsize=10)
            ax.imshow(waterfall, rasterized=True, interpolation="nearest",
                      cmap="viridis", aspect="auto")
        else:
            ax.axis("off")

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.8)
    return fig, axs

def show_multiple(named_waterfalls):
    """
    Show multiple waterfalls.

    named_waterfalls is a list of (name, waterfall).
    waterfall is an array indexed like [time, chan]
    """
    fig, axs = plot_multiple(named_waterfalls)
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

        self.antenna_names = [s.decode("utf-8")
                              for s in self.h5["/telinfo/antenna_names"][()]]
        
        # Validate shapes of things
        assert self.delays.shape == (len(self.time_array), self.nbeams, self.nants)
        if self.cal_all.shape != (self.nchan, self.npol, self.nants):
            print("cal_all shape:", self.cal_all.shape)
            print("nchan, npol, nants:", (self.nchan, self.npol, self.nants))
            raise ValueError("unexpected cal_all size")
        
        
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

    def show_classic_incoherent(self):
        incoherent = np.square(self.real_array()).sum(axis=(2, 3, 4))
        snr, sig = self.snr_and_signal(incoherent)
        print(f"recalculated power: {sig:e}")
        print("local SNR:", snr)
        show_array(incoherent)

    def weighted_incoherent(self):
        # Start off like we're beamforming beam 0
        coeffs = self.coefficients(0)
        inputs = self.complex_array()
        presum = np.conjugate(coeffs) * inputs

        # But do the power calculation before summing
        power = np.square(np.real(presum)) + np.square(np.imag(presum))

        # Then sum along polarization and antenna
        return power.sum(axis=(2, 3))

    def show_weighted_incoherent(self):
        incoherent = self.weighted_incoherent()
        snr, sig = self.snr_and_signal(incoherent)
        print(f"recalculated power: {sig:e}")
        print("local SNR:", snr)
        show_array(incoherent)
    
    def show_antenna(self, index):
        voltages = self.real_array()[:, :, :, index, :]
        powers = np.square(voltages).sum(axis=(2, 3))
        fig, ax = plot_array(powers)
        
        yticks = ax.get_yticks()
        ax.set_yticklabels([
            f"{tick*self.stamp.tsamp*1e6}"
            for tick in yticks
        ])
        ax.set_ylabel("us")
        
        xticks = ax.get_xticks()
        ax.set_xticklabels([
            f"{self.stamp.fch1 + tick*self.stamp.foff}"
            for tick in xticks
        ])
        ax.set_xlabel("Frequency (MHz)")
        
        display(fig)
        plt.close()

    def show_antennas(self):
        antennas = np.square(self.real_array()).sum(axis=(2, 4))
        fig, axs = plot_multiple([(f"antenna {i}", antennas[:, :, i])
                       for i in range(self.stamp.numAntennas)])

        for ax_r in range(axs.shape[0]):
            for ax_c in range(axs.shape[1]):
                ax = axs[ax_r, ax_c]
                    
                yticks = ax.get_yticks()
                if ax_c == 0:
                    ax.set_yticklabels([
                        f"{tick*self.stamp.tsamp*1e3:0.3f}"
                        for tick in yticks
                    ])
                    ax.set_ylabel("ms")
                else:
                    ax.set_yticklabels([])
                
                xticks = ax.get_xticks()
                ax.set_xticklabels([
                    f"{(tick*self.stamp.foff*1000):0.1f}"
                    for tick in xticks
                ])
                ax.set_xlabel(f"Frequency (kHz + {self.stamp.fch1:0.6f} MHz)")
        
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
        
        Output dimensions are [time, chan, pol]
        """
        coeffs = self.coefficients(beam)
        inputs = self.complex_array()

        # Sum along polarization and antenna
        return (np.conjugate(coeffs) * inputs).sum(axis=3)

    def beamform_power(self, beam):
        """ Converts voltage to power and combines across polarities.

        Output dimensions are [time, chan]
        """
        voltage = self.beamform_voltage(beam)
        squared = np.square(np.real(voltage)) + np.square(np.imag(voltage))
        return squared.sum(axis=2)

    def show_beam(self, beam):
        power = self.beamform_power(beam)
        snr, sig = self.snr_and_signal(power)
        print(f"recalculated power: {sig:e}")
        print("local SNR:", snr)
        show_array(power)

    def show_best_beam(self):
        beam = self.stamp.signal.beam
        if beam < 0:
            print("best beam is incoherent")
            print(f"original power: {self.stamp.signal.power:e}")
            print(f"original SNR: {self.stamp.signal.snr}")
            self.show_weighted_incoherent()
            return
        print("best beam is", beam)
        print(f"original power: {self.stamp.signal.power:e}")
        print("original SNR:", self.stamp.signal.snr)
        self.show_beam(beam)

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

    def snr_and_signal(self, data):
        """Returns a (snr, signal) tuple."""
        # Calculate the noise based on the first and last 20 column sums
        left_column_sums = data[:, :20].sum(axis=1)
        right_column_sums = data[:, -20:].sum(axis=1)
        column_sums = np.concatenate((left_column_sums, right_column_sums))
        mean = column_sums.mean()
        std = column_sums.std()

        signal = (data * self.signal_mask()).sum()

        return ((signal - mean) / std, signal)

    def snr(self, data):
        snr, _ = self.snr_and_signal()
        return snr

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

    def best_correlations(self):
        corr = self.correlations()
        nants = self.stamp.numAntennas
        possible = []
        for i in range(nants):
            for j in range(i + 1, nants):
                possible.append((corr[i, j], (i, j)))
        possible.sort()
        possible.reverse()
        for score, (i, j) in possible[:5]:
            name_i = self.recipe.antenna_names[i]
            name_j = self.recipe.antenna_names[j]
            print(f"{name_i} and {name_j} have correlation {corr[i, j]:.3f}")
    
    def show_correlations(self):
        corr = self.correlations()
        nants = self.stamp.numAntennas
        print("median correlation:", np.median(corr))
        size = 15
        fig, ax = plt.subplots(figsize=(size, size))
        ax.imshow(corr, rasterized=True, interpolation="nearest", cmap="plasma", vmin=0)
        ax.set_yticks(list(range(nants)))
        ax.set_yticklabels(self.recipe.antenna_names)
        ax.set_xticks([])
        display(fig)
        plt.close()

    def show_correlations_text(self):
        corr = self.correlations()
        nants = self.stamp.numAntennas
        print("   " + " ".join(f"{n:4d}" for n in range(nants)))
        for i in range(nants):
            print(f"{i:2d} " + " ".join(f"{corr[i, j]:.2f}" for j in range(nants)))

            
        
def read_stamps(filename):
    with open(filename) as f:
        stamps = stamp_capnp.Stamp.read_multiple(f, traversal_limit_in_words=2**30)
        for s in stamps:
            yield Stamp(s)

def read_events(filename):
    with open(filename) as f:
        events = hit_capnp.Event.read_multiple(f)
        for e in events:
            yield e
            
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
            elif full.endswith(".stamp"):
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
