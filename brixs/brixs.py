#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python package for analysis of RIXS spectra.

Galdino"""

# standard libraries
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import copy
from matplotlib.transforms import Bbox
import warnings

# backpack
import brixs
import brixs.arraymanip as manip
from brixs.arraymanip import index
import brixs.filemanip as fmanip
import brixs.figmanip as figmanip
from brixs.model_functions import fwhmGauss, fwhmAreaGauss

class photon_events():

    def __init__(self, filepath=None, data=None, x_max=None, y_max=None):
        """ x_max and y_max are important for:
            calculate_spectrum()
            apply_spread()
            bining()
        """
        self.data       = None
        self.x_max      = None
        self.y_max      = None
        self.columns    = None
        self.x_edges    = None
        self.y_edges    = None
        self.x_centers  = None
        self.y_centers  = None



        if filepath is not None:
            self.load(filepath)
        else:
            if data is None:
                pass#raise AttributeError('File and Data cannot be None at the same time.')
            else:
                self.data = copy.deepcopy(data)
                # self.data_original = copy.copy(self.data)

        if self.x_max is None and x_max is None:
            self.x_max = max(self.data[:, 0])
        elif self.x_max is None:
            self.x_max = x_max

        if self.y_max is None and y_max is None:
            self.y_max = max(self.data[:, 1])
        elif self.y_max is None:
            self.y_max = y_max

        # overlap
        self.x_min_between_events = None
        self.y_min_between_events = None
        self.data_overlaped = None
        self.n_overlaps = None

        # curvature
        self.curvature_binx = None
        self.curvature_biny = None
        self.curvature_ref = None
        self.curvature_poly_order = None
        self.curvature_offsets = None
        self.curvature_parameters = None
        self.curvature_min        = None
        self.curvature_max        = None

        # spectrum
        self.spectrum = None


    def load(self, filepath):
        for row in fmanip.getComments(filepath=Path(filepath)):

            if row.startswith('# x_max '):
                self.x_max = float(row.split('# x_max ')[-1])
            elif row.startswith('# y_max '):
                self.y_max = float(row.split('# y_max ')[-1])
        self.data = fmanip.load_data(filepath=Path(filepath), force_array=True)


    def save(self, filepath):
        header  = f'x_max {self.x_max}\n'
        header += f'y_max {self.y_max}\n'
        header += f'x y I'
        fmanip.save_data(self.data, filepath=Path(filepath), header=header)


    # def apply_spread(self, type='gaussian', x_spread=5e-6, y_spread=5e-6):
    #     for photon_event in self.data:
    #         photon_event[0] += np.random.normal(0, x_spread)
    #         photon_event[1] += np.random.normal(0, y_spread)
    #
    #         if photon_event[0] > self.x_max:
    #             photon_event[0] = self.x_max
    #         elif photon_event[0] < 0:
                # photon_event[0] = 0
    #
    #         if photon_event[1] > self.y_max:
    #             photon_event[1] = self.y_max
    #         elif photon_event[1] < 0:
    #             photon_event[1] = 0


    # def calculate_overlaps(self, x_min_between_events=5e-6, y_min_between_events=5e-6):
    #
    #     overlap_x = [abs(self.data[:, 0]-x)<x_min_between_events for x in self.data[:, 0]]
    #     overlap_y = [abs(self.data[:, 1]-y)<y_min_between_events for y in self.data[:, 1]]
    #
    #     overlaps = [sum(x*y)>1 for x, y in zip(overlap_x, overlap_y)]
    #     self.n_overlaps = sum(overlaps)/2
    #
    #     data_overlaped = []
    #     for idx, photon_event in enumerate(data):
    #         if overlaps[idx]:
    #             data_overlaped.append(photon_event)
    #     self.data_overlaped = np.array(data_overlaped)
    #
    #     self.x_min_between_events = x_min_between_events
    #     self.y_min_between_events = y_min_between_events


    def bining(self, binx, biny):
        self.columns, self.x_edges, self.y_edges = np.histogram2d(self.data[:, 0],
                                                                 self.data[:, 1],
                                                                 bins=(binx, biny),
                                                                 weights=self.data[:, 2],
                                                                 range=((0, self.x_max), (0, self.y_max))
                                                                )

        self.x_centers = self._movingaverage(self.x_edges, window_size=2)[1:]
        self.y_centers = self._movingaverage(self.y_edges, window_size=2)[1:]


    def _movingaverage(self, interval, window_size):
        window = np.ones(int(window_size))/float(window_size)
        return np.convolve(interval, window, 'same')


    # def plot_overlaped(self, ax=None, pointsize=5):
    #
    #     if ax is None:
    #         fig = figmanip.figure()
    #         ax = fig.add_subplot(111)
    #         ax.set_facecolor('black')
    #
    #     ax.errorbar(self.data_overlaped[:, 0]*10**3, self.data_overlaped[:, 1]*10**3,
    #                  linewidth=0,
    #                  fmt='o',
    #                  mfc = 'red',
    #                  elinewidth = 1,
    #                  yerr=self.y_min_between_events *10**3,
    #                  xerr=self.x_min_between_events *10**3,
    #                  marker='o',
    #                  ms=pointsize)
    #     return ax


    def plot(self, ax=None, pointsize=5, show_binx=False, show_biny=False, show_curvature=False, shift=0, show_lim=False):

        if ax is None:
            fig = figmanip.figure()
            ax = fig.add_subplot(111)
            ax.set_facecolor('black')

        ax.plot(self.data[:, 0]*10**3, self.data[:, 1]*10**3,
                     linewidth=0,
                     marker='o',
                     mfc = 'white',
                     ms=pointsize)

        if show_binx:
            ax.vlines(self.x_edges*10**3, 0, self.y_max*10**3, color='red', linewidth=0.8, zorder=3)
        if show_biny:
            ax.hlines(self.y_edges*10**3, 0, self.x_max*10**3, color='red', linewidth=0.5, zorder=3)

        try:
            if show_curvature:
                shift += (self.curvature_min+self.curvature_max)/2
                self.plot_curvature(ax=ax, shift=shift, color='green', linewidth=2)
            if show_lim:
                ax.hlines(self.curvature_min*10**3, 0, self.x_max*10**3, color='white', linewidth=2, zorder=4)
                ax.hlines(self.curvature_max*10**3, 0, self.x_max*10**3, color='white', linewidth=2, zorder=4)
        except:
            pass

        return ax


    def calculate_curvature(self, ref=0, type='cross-correlation', poly_order=2, binx=25, biny=700, y_min=None, y_max=None):
        """

        Args:
            binx = 25  # Must be as large as possible. The larger, the more outliers we got
            biny = 700  # Must be as small as possible, but there must not be more than one point per "row"
        """

        self.bining(binx, biny)
        self.curvature_binx = copy.copy(binx)
        self.curvature_biny = copy.copy(biny)

        if y_min == None:
            y_min = 0

        if y_max == None:
            y_max = self.y_max

        offsets = np.array([])
        for column in self.columns:
            y_centers, ref_column = manip.extract(self.y_centers, self.columns[ref], ranges=((y_min, y_max),) )
            y_centers, column     = manip.extract(self.y_centers, column, ranges=((y_min, y_max),) )

            if type == 'cross-correlation':
                cross_correlation = np.correlate(column, ref_column, mode='same')
                offsets           = np.append(offsets, y_centers[np.argmax(cross_correlation)])
            elif type == 'max':
                offsets           = np.append(offsets, y_centers[np.argmax(column)] - y_centers[np.argmax(ref_column)])

        if type == 'cross-correlation': offsets -= offsets[ref]

        self.curvature_ref        = copy.copy(ref)
        self.curvature_poly_order = copy.copy(poly_order)
        self.curvature_offsets    = copy.copy(offsets)
        self.curvature_parameters = np.polyfit(self.x_centers, offsets, deg=poly_order)
        self.curvature_min        = copy.copy(y_min)
        self.curvature_max        = copy.copy(y_max)


    def plot_offsets(self, ax=None):
        if ax is None:
            fig = figmanip.figure()
            ax = fig.add_subplot(111)

        ax.scatter(self.x_centers*10**3, self.curvature_offsets*10**3)
        return ax


    def plot_columns(self, ref='all', ax=None, show_lim=True, **kwargs):
        if ax is None:
            fig = figmanip.figure()
            ax = fig.add_subplot(111)

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'ms' not in kwargs:
            kwargs['ms'] = 5


        if ref == 'all':
            for col in self.columns:
                ax.plot(self.y_centers*10**3, col, **kwargs)
                try:
                    if max(col) > max_intensity:
                        max_intensity = max(col)
                    if min(col) < min_intensity:
                        min_intensity = min(col)
                except:
                    min_intensity = min(col)
                    max_intensity = max(col)
        else:
            ax.plot(self.y_centers*10**3, self.columns[ref], **kwargs)
            max_intensity = max(self.columns[ref])

        if show_lim:
            ax.vlines(self.curvature_min*10**3, min_intensity, max_intensity, color='red', linewidth=0.8, zorder=3)
            ax.vlines(self.curvature_max*10**3, min_intensity, max_intensity, color='red', linewidth=0.8, zorder=3)


    def plot_curvature(self, ax=None, shift=0, **kwargs):

        if 'color' not in kwargs:
            kwargs['color'] = 'red'
        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = 1
        if 'zorder' not in kwargs:
            kwargs['zorder'] = 10

        if ax is None:
            fig = figmanip.figure()
            ax = fig.add_subplot(111)

        x_curvature = np.linspace(0, self.x_max, 1000)
        y_curvature = np.polyval(self.curvature_parameters, x_curvature)

        # shift = max(self.columns[ref_column]) == self.columns[ref_column]
        ax.plot(x_curvature*10**3, (y_curvature+shift)*10**3, **kwargs)

        return ax


    def apply_curvature(self):
        self.data[:, 1] -= np.polyval(self.curvature_parameters, self.data[:, 0])
        self.bining(binx=self.curvature_binx, biny=self.curvature_biny)


    def calculate_spectrum(self, biny=10000):
        h, _, y_edges = np.histogram2d(self.data[:, 0],
                                     self.data[:, 1],
                                     bins=(1, biny),
                                     weights=self.data[:, 2],
                                     range=((0, self.x_max), (0, self.y_max))
                                    )
        y_centers = self._movingaverage(y_edges, window_size=2)[1:]
        self.spectrum = spectrum(data=np.vstack((y_centers, sum(h))).transpose())


class spectrum():


    def __init__(self, filepath=None, data=None):
        self.data = None
        self.filepath = filepath

        if filepath is not None:
            self.load(filepath)
        else:
            if data is not None:
                self.data = copy.deepcopy(data)


    def save(self, filepath=None, header='x, y'):
        """ add extension.
        """
        if filepath is None:
            filepath = self.filepath
        fmanip.save_data(self.data, filepath=Path(filepath), header=header)


    def load(self, filepath=None):
        if filepath is None:
            filepath = self.filepath
        self.data = fmanip.load_data(filepath, force_array=True)
        self.filepath = filepath


    def plot(self, ax=None, normalized=True, vertical_increment=0, shift=0, factor=1, **kwargs):
        if ax is None:
            fig = figmanip.figure()
            ax = fig.add_subplot(111)

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'ms' not in kwargs:
            kwargs['ms'] = 5

        if normalized:
            ax.plot((self.data[:, 0] + shift)*factor, self.data[:, 1]/max(self.data[:, 1]) + vertical_increment, **kwargs)
        else:
            ax.plot((self.data[:, 0] + shift)*factor, self.data[:, 1] + vertical_increment, **kwargs)

        return ax


    def calibrate(self, dispersion, zero=0, normalized=True):

        calibrated_spectrum = copy.copy(self.data)
        calibrated_spectrum[:, 0] = calibrated_spectrum[:, 0]*dispersion + zero
        if normalized:
            calibrated_spectrum[:, 1] = self.data[:, 1]/max(self.data[:, 1])

        return spectrum(data=calibrated_spectrum)


    def interpolate(self, start=None, stop=None, num=10000, step=None):
        if start is None:
            start = min(self.data[:, 0])
        if stop is None:
            stop = max(self.data[:, 0])
        # print(np.all(np.diff(self.data[:, 0]) > 0))

        if step is not None:   # step overwrites num
            temp = np.arange(start, stop, step=step)
            temp = np.zeros((len(temp), self.data.shape[1]))
            temp[:, 0] = np.arange(start, stop, step=step)
        else:
            temp = np.zeros((num, self.data.shape[1]))
            temp[:, 0] = np.linspace(start, stop, num=num)
        temp[:, 1] = np.interp(temp[:, 0], self.data[:, 0], self.data[:, 1])

        return spectrum(data=temp)


    def shift(self, shift, mode='hard'):
        temp = np.zeros(self.data.shape)
        temp[:, 0], temp[:, 1] = manip.shift(self.data[:, 0], self.data[:, 1], shift=shift, mode=mode)
        return spectrum(data=temp)


    def peak_fit(self, start=None, stop=None, **kwargs):
        if start is None:
            start = min(self.data[:, 0])
        if stop is None:
            stop = max(self.data[:, 0])

        # index
        start = manip.index(self.data[:, 0], start)
        stop  = manip.index(self.data[:, 0], stop)

        if len(self.data[start:stop, 0]) == 0:
            warnings.warn('start/stop outside data range', stacklevel=2)
            return None, None

        # guess
        if 'guess_A' not in kwargs:       kwargs['guess_A'] = max(self.data[start:stop, 1])
        if 'guess_c' not in kwargs:       kwargs['guess_c'] = self.data[start:stop, 0][np.where(self.data[start:stop, 1] == kwargs['guess_A'])[0][0]]
        if 'guess_w' not in kwargs:       kwargs['guess_w'] = (self.data[1, 0] - self.data[0, 0])*10
        if 'guess_offset' not in kwargs:  kwargs['guess_offset'] = np.mean(self.data[start:stop, 1])
        if 'asymmetry' not in kwargs:     kwargs['asymmetry'] = True
        if 'fixed_m' not in kwargs:         kwargs['fixed_m'] = 0

        _, arr100, popt_2 = manip.peak_fit(self.data[start:stop, 0], self.data[start:stop, 1], **kwargs)

        return arr100, popt_2


class spectra():

    def __init__(self, folderpath=None, data=None):
        self.spectrum = []
        self.folderpath = folderpath

        if folderpath is not None:
            self.load(folderpath)
        else:
            if data is None:
                pass#raise AttributeError('File and Data cannot be None at the same time.')
            else:
                self.spectrum = copy.deepcopy(data)

        self.shifts = []
        self.shift_type = 'fixed'
        self.sum = None


    def spectra_count(self):
        return len(self.spectrum)


    def append(self, filepath=None, data=None):
        # check if spectrum is type spectrum


        if filepath is not None:
            self.spectrum.append(spectrum(filepath=filepath))
        elif data is not None:
            self.spectrum.append(data)


    def exclude(self, idx):
        del self.spectrum[idx]


    def save(self, folderpath=None, prefix='spectrum'):
        if folderpath is None:
            folderpath = self.folderpath

        n_digits = figmanip.n_digits(len(self.spectrum)-1)[0]
        for idx, s in enumerate(self.spectrum):
            filename = f'{prefix}_' + f'{idx}'.zfill(n_digits) + '.dat'
            s.save(filepath=folderpath/filename)


    def load(self, folderpath=None):
        if folderpath is None:
            folderpath = self.folderpath
        for filepath in fmanip.filelist(folderpath):
            self.append(filepath=filepath)
        self.folderpath = folderpath


    def plot(self, ax=None, n='all', normalized=True, vertical_increment=0, shift=0, **kwargs):
        if ax is None:
            fig = figmanip.figure()
            ax = fig.add_subplot(111)

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'ms' not in kwargs:
            kwargs['ms'] = 5

        if n == 'all':
            for i in range(len(self.spectrum)):
                self.spectrum[i].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment*i, shift=shift, **kwargs)
        else:
            self.spectrum[n].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, **kwargs)

        return ax


    def calculate_shifts(self, ref=0, type='cross-correlation', start=None, stop=None):


        self._check_array()

        if type == 'cross-correlation':
            self.shift_type = 'fixed'
            zero_shift = np.argmax(np.correlate(self.spectrum[ref].data[:, 1], self.spectrum[ref].data[:, 1], mode='Same'))

            self.shifts = []
            for i, spectrum in enumerate(self.spectrum):
                cross_corr = np.correlate(self.spectrum[ref].data[:, 1], spectrum.data[:, 1], mode='Same')
                self.shifts.append(np.argmax(cross_corr) - zero_shift)
                print(f'spectrum {i} shift calculated!')
        elif type == 'fit':
            self.shift_type = 'float'

            spectrum = self.spectrum[ref]

            _, popt = spectrum.peak_fit(start, stop)
            cen = popt[1]
            for i, spectrum in enumerate(self.spectrum):
                _, popt = spectrum.peak_fit(start, stop)
                self.shifts.append(popt[1] - cen)
                print(f'spectrum {i} shift calculated!')


    def apply_shifts(self):
        if self.shift_type == 'fixed':
            for i, spectrum in enumerate(self.spectrum):
                step = spectrum.data[1, 0] - spectrum.data[0, 0]
                spectrum.data = np.vstack((spectrum.data[:, 0] + self.shifts[i]*step, spectrum.data[:, 1])).T

        elif self.shift_type == 'float':
            for i, spectrum in enumerate(self.spectrum):
                spectrum.data[:, 0] -= self.shifts[i]


    def interpolate(self, start=None, stop=None, num=30000, step=None):
        if start is None:
            start = max([spectrum.data[0, 0] for spectrum in self.spectrum])
        if stop is None:
            stop = min([spectrum.data[-1, 0] for spectrum in self.spectrum])


        for i, spectrum in enumerate(self.spectrum):
            self.spectrum[i] = spectrum.interpolate(start, stop, num, step)


    def crop(self, start=None, stop=None):
        if start is None:
            start = max([spectrum.data[0, 0] for spectrum in self.spectrum])
        if stop is None:
            stop = min([spectrum.data[-1, 0] for spectrum in self.spectrum])


        for spectrum in self.spectrum:
            step = spectrum.data[1, 0] - spectrum.data[0, 0]
            data_range = (start + step/2, stop - step/2)
            a, b = manip.extract(spectrum.data[:, 0], spectrum.data[:, 1], ranges=(data_range,) )

            temp = np.zeros((len(a), 2))
            temp[:, 0] = a
            temp[:, 1] = b
            spectrum.data = copy.copy(temp)


    def _check_array(self):

        # check length
        for idx, spectrum in enumerate(self.spectrum):
            try:
                if len(spectrum.data) != len(self.spectrum[idx+1].data):
                    raise ValueError(f"Spectrum {idx} and {idx+1} have the different length.")
            except IndexError:
                pass

        # check x
        for idx, spectrum in enumerate(self.spectrum):
            try:
                step = spectrum.data[1, 0] - spectrum.data[0, 0]
                if max(abs(spectrum.data[:, 0] - self.spectrum[idx+1].data[:, 0]))*100/step > 0.01:
                    raise ValueError(f"Spectrum {idx} and {idx+1} seems to be different.")
                # print(max(abs(spectrum.data[:, 0] - self.spectrum[idx+1].data[:, 0]))*100/step)
            except IndexError:
                pass


    def calculate_sum(self):

        self._check_array()

        temp = copy.copy(self.spectrum[0])
        for i in range(1, self.spectra_count()):
            temp.data[:, 1] += self.spectrum[i].data[:, 1]
        self.sum = temp


    def calibrate(self, dispersion=None, zero=None, normalized=True):
        s = []
        for spectrum in self.spectrum:
            s.append(spectrum.calibrate(dispersion=dispersion, zero=zero, normalized=normalized))
        return spectra(data=s)
