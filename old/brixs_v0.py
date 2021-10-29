#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module for analysis of RIXS spectra.

This module is based on three class, one for dealing with photon events lists,
one for dealing with single spectra, a another one to deal with many spectra at
a time.


.. autosummary::

    photon_events
    spectrum
    spectra

"""

# standard libraries
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import copy
from matplotlib.transforms import Bbox
import warnings
from scipy.optimize import curve_fit

# backpack
from . import arraymanip as am
from .arraymanip import index
from . import filemanip as fm
from . import figmanip as figm





class photon_events():
    """Creates a ``photon_event`` class type object to deal with photon events lists.

    Args:
        filepath (string or pathlib.Path, optional): filepath to file. It overwrites
            the ``data`` argument.
        data (list or array, optional): two (x, y) or three (x, y, intensity) list or array with
            photon events.
        delimiter (str, optional): The string used to separate values. If whitespaces are used,
            consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).
        x_max (float, optional): maximum x value. If ``None``, it will be infered
            by the data.
        y_max (float, optional): maximum y value. If ``None``, it will be infered
            by the data.

    """


    def __init__(self, filepath=None, data=None, delimiter=',', x_max=None, y_max=None):
        # basic attr
        self.data       = None
        self.filepath   = None
        self.x_max      = None
        self.y_max      = None

        # binning attr
        self.hist       = None
        self.bins       = None
        self.bins_size  = None
        self.x_edges    = None
        self.y_edges    = None
        self.x_centers  = None
        self.y_centers  = None

        # offset attr
        self.offsets         = None
        self.offsets_ref     = None
        self.offsets_ranges  = None
        self.offsets_func    = None
        self.offsets_par     = None

        # spectrum attr
        self.spectrum = None

        self.load(filepath=filepath, data=data, delimiter=delimiter, x_max=x_max, y_max=y_max)
        self.set_binning(bins=1)
        self.calculate_offsets()
        self.fit_offsets(deg=0)
        self.calculate_spectrum()


    def load(self, filepath=None, data=None, delimiter=',', x_max=None, y_max=None):
        """Load photon events data from file or assign data directly.

        The file/data must have three columns, x, y, and intensity.

        One can pass the values for x_max and y_max through the file header by
        using the tag ``# x_max <value>`` and ``# y_max <value>``, respectively.

        args:
            data (list or array, optional): three column list (or array) with photon
                events. Column order should be x, y, and intensity.
            filepath (string or pathlib.Path, optional): filepath to file. It overwrites
                the ``data`` argument.
            delimiter (str, optional): The string used to separate values. If whitespaces are used,
                consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).
            x_max (float, optional): maximum x value. If ``None``, it will be infered
                by the data.
            y_max (float, optional): maximum y value. If ``None``, it will be infered
                by the data.

        returns:
            None

        See Also:
            :py:func:`photon_events.save`.
        """

        x_max_flag = False
        y_max_flag = False
        if filepath is not None:
            self.filepath = Path(filepath)

            # check x_max and y_max
            for row in fm.load_Comments(filepath=self.filepath):

                if row.startswith('# x_max') or row.startswith('#x_max'):
                    self.x_max = float(row.split('x_max')[-1])
                    x_max_flag = True
                elif row.startswith('# y_max') or row.startswith('#y_max'):
                    self.y_max = float(row.split('y_max')[-1])
                    y_max_flag = True

            # get data
            self.data = fm.load_data(filepath=Path(filepath), force_array=True, delimiter=delimiter)

        else:
            self.filepath   = None
            if data is None:
                warnings.warn('No filepath or data to load.', stacklevel=2)
                return
            else:
                self.data = copy.deepcopy(data)




        # infer x_max and y_max if necessary
        if x_max_flag is False:
            if x_max is None:
                self.x_max = max(self.data[:, 0])
            else:
                self.x_max = copy.deepcopy(x_max)

        if y_max_flag is False:
            if y_max is None:
                self.y_max = max(self.data[:, 1])
            else:
                self.y_max = copy.deepcopy(y_max)


    def save(self, filepath, delimiter=','):
        """Saves photon events data to a file.

        args:
            filepath (string or pathlib.Path, optional): filepath to file.
            delimiter (str, optional): The string used to separate values. If whitespaces are used,
                consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).

        returns:
            None

        note:
            x_max and y_max values are saved in the header.

        See Also:
            :py:func:`photon_events.load`.
        """
        header  = f'x_max {self.x_max}\n'
        header += f'y_max {self.y_max}\n'
        header += f'x y I'
        fm.save_data(self.data, filepath=Path(filepath), delimiter=delimiter, header=header)


    def set_binning(self, bins=None, bins_size=None):
        """Compute the histogram of the data set (binning of the data).

        args:
            bins (int or tuple, optional): number of bins. If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the x and y directions, respectively.
            bins_size (int or tuple, optional): size of the bins. This overwrites
                the argument ``bins``. If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the x and y directions, respectively.

        return:
            None
        """
        if bins_size is not None:
            try:
                if len(bins_size) == 1:
                    x_bins_size, y_bins_size = bins_size[0], bins_size[0]
                else:
                    x_bins_size, y_bins_size = bins_size[0], bins_size[1]
            except TypeError:
                x_bins_size, y_bins_size = bins_size, bins_size

            x_bins = int(self.x_max/x_bins_size)
            y_bins = int(self.y_max/y_bins_size)

        else:
            if bins is None:
                warnings.warn('Bins not defined.', stacklevel=2)
                return
            else:
                try:
                    if len(bins) == 1:
                        x_bins, y_bins = bins[0], bins[0]
                    else:
                        x_bins, y_bins = bins[0], bins[1]
                except TypeError:
                    x_bins, y_bins = bins, bins

        bins_size = (self.x_max/x_bins, self.y_max/y_bins)

        self.bins = (x_bins, y_bins)
        self.bins_size = copy.deepcopy(bins_size)

        self.hist, self.x_edges, self.y_edges = np.histogram2d(self.data[:, 0],
                                                                 self.data[:, 1],
                                                                 bins=(x_bins, y_bins),
                                                                 weights=self.data[:, 2],
                                                                 range=((0, self.x_max), (0, self.y_max))
                                                                )

        self.x_centers = am.movingaverage(self.x_edges, window_size=2, remove_boundary_effects=True)
        self.y_centers = am.movingaverage(self.y_edges, window_size=2, remove_boundary_effects=True)


    def calculate_offsets(self, ref=0, mode='cross-correlation', ranges=None):
        """Calculate the offset of each column relative to a reference column.

        args:
            ref (int, optional): reference column. The offset of all other columns
                is calculated based on the reference column. Default is 0.
             mode (string, optional): method used to calculate the offsets.
                The current options are: 'cross-correlation', and 'max'.
            ranges (list, optional): a pair of x values or a list of pairs. Each pair represents
                the start and stop of a data range. If None, the whole data set is used.

        returns:
            None
        """

        if self.hist is None:
            warnings.warn('Data not binned yet. Use binning()', stacklevel=2)
            return

        if ranges is None:
            ranges = [[0, self.y_max]]

        # check ranges
        for r in ranges:
            if r[0] > max(self.y_centers) or r[-1] < min(self.y_centers):
                raise ValueError('Selected ranges outside data range (y_centers).')

        self.offsets = np.zeros(self.hist.shape[0])
        y_centers, ref_column = am.extract(self.y_centers, self.hist[ref], ranges=ranges)
        for i in range(self.hist.shape[0]):

            y_centers, column     = am.extract(self.y_centers, self.hist[i], ranges=ranges)

            if mode == 'cross-correlation':
                cross_correlation = np.correlate(column, ref_column, mode='same')
                self.offsets[i]   = y_centers[np.argmax(cross_correlation)]
            elif mode == 'max':
                self.offsets[i]   = y_centers[np.argmax(column)]

        # if mode == 'cross-correlation': self.offsets -= self.offsets[ref]
        # if mode == 'max': self.offsets -= self.offsets[ref]

        self.offsets -= self.offsets[ref]
        self.offsets_ref        = copy.copy(ref)
        self.offsets_ranges     = copy.copy(ranges)


    def fit_offsets(self, deg=1, f=None):
        """Find the curve that fits the offset values.

        args:
            deg (int, optional): degree for the polnomial fit. The default is 1.
            f (function, optional):  a function y = f(x, a, b, ...) that returns the
                value of y as a function of x and other parameters to be fitted.
                This overwrites the polynomal fit based on the argument ``deg``.

        returns:
            None
        """

        x2fit = self.x_centers
        if x2fit is None:
            warnings.warn('Data not binned yet. Use binning()', stacklevel=2)
            return

        y2fit = self.offsets
        if y2fit is None:
            warnings.warn('Offsets not defined. Use get_offsets().', stacklevel=2)
            return

        if f is None:
            if deg < 0:
                warnings.warn('deg must be a positive value or zero.', stacklevel=2)
                return
            popt = np.polyfit(x2fit, y2fit, deg=deg)
            f = np.poly1d(popt)

        else:
            popt, pcov = curve_fit(f, x2fit, y2fit)
            f = lambda x: f(x, *popt)

        self.offsets_func = lambda x: f(x)
        self.offsets_par = popt


    def offsets_correction(self):
        """Uses the offsets fitted curve to adjust the photon events."""
        f = lambda x, y: (x, y-self.offsets_func(x))
        self.apply_correction(f=f)


    def apply_correction(self, f):
        """Changes the values of x, y based on a function.

        args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        returns:
            None
        """
        self.data[:, 0], self.data[:, 1] = f(self.data[:, 0], self.data[:, 1])
        self.x_max, self.y_max  = f(self.x_max, self.y_max)
        self.set_binning(bins=self.bins)


    def calculate_spectrum(self, y_bins=None, y_bins_size=None):
        """Sum the photon events in the x direction.

        args:
            y_bins (int, optional): number of y bins. If None, the current binning is used.
            bins_size (int or tuple, optional): size of the y bins. This overwrites
                the argument ``y_bins``. If None, the current binning is used.

        returns:
            None
        """
        if y_bins is None:
            self.spectrum = spectrum(data=np.vstack((self.y_centers, sum(self.hist))).transpose())
        else:
            temp = photon_events(data=self.data)
            if y_bins_size is not None:
                temp.set_binning(bins_size=(self.x_max+1, y_bins_size))
            elif y_bins is not None:
                temp.set_binning(bins=(1, y_bins))
            self.spectrum = spectrum(data=np.vstack((temp.y_centers, sum(temp.hist))).transpose())


    def plot(self, ax=None, pointsize=1, show_bins=(False, False), show_offsets=False, show_offsets_fit=False, **kwargs):
        """Plot photon events.

        args:
            ax (matplotlib.axes, optional): axes for plotting on.
            pointsize (int, optional): photon events point size. Default is 1.
            show_bins (bool, optional): if True, bins edges are displayed in cyan.
                A tuple of bools can also be used to display only x bins or y bins,
                e. g., ``show_bins = (True, False)`` will display only x bins.
                The default is (False, False).
            show_offsets (bool, optional): if True, offsets are displayed in yellow
                over its respectively bin. The ranges of data used to calculate
                the offsets are marked by green and red lines.
            show_offsets_fit (bool, optional): if True, the fit of the curve
                defined by the offsets values is displayed in yellow.
                The ranges of data used to calculate
                the offsets are marked by green and red lines.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data (photon events).

        returns:
            matplotlib.axes
        """
        if ax is None:
            fig = figm.figure()
            ax = fig.add_subplot(111)
            ax.set_facecolor('black')

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'ms' not in kwargs:
            kwargs['ms'] = pointsize
        if 'mfc' not in kwargs:
            kwargs['mfc'] = 'white'
        if 'markeredgewidth' not in kwargs:
            kwargs['markeredgewidth'] = 0

        ax.plot(self.data[:, 0], self.data[:, 1],
                     linewidth=0,
                     **kwargs)

        try:
            if len(show_bins) == 1:
                show_bins = (show_bins[0], show_bins[0])
        except TypeError:
            show_bins = (show_bins, show_bins)

        if show_bins[0]:
            plt.vlines(self.x_edges, 0, self.y_max, color='cyan', linewidth=0.8, zorder=3)
        if show_bins[1]:
            plt.hlines(self.y_edges, 0, self.x_max, color='cyan', linewidth=0.5, zorder=3)

        if show_offsets:
            if self.offsets is None:
                warnings.warn('Offsets not defined. Use get_offsets().', stacklevel=2)
            else:
                c = self.y_centers[np.argmax(self.hist[0])]
                self.plot_offsets(ax=ax, shift=(-self.offsets[0] + c), color='yellow', zorder=10)

        if show_offsets_fit:
            if self.offsets is None:
                warnings.warn('Offsets not fitted. Use fit_offsets().', stacklevel=2)
            else:
                c = self.y_centers[np.argmax(self.hist[0])]
                self.plot_offsets_fit(ax=ax,  shift=(-self.offsets[0] + c), linewidth=1, color='yellow', zorder=10)
                # x, y = self.offsets_f(np.linspace(0, self.x_max, 1000), np.zeros(1000))
                # plt.plot(x, y-y[0]+c, linewidth=1, color='yellow')

        if show_offsets or show_offsets_fit:
            for r in self.offsets_ranges:
                plt.axhline(r[0], color='green', linewidth=2, zorder=10)
                plt.axhline(r[1], color='red', linewidth=2, zorder=10)

        return ax


    def plot_offsets(self, ax=None, shift=0, **kwargs):
        """Plot offsets as function of x values (center of x bins).

        args:
            ax (matplotlib.axes, optional): axes for plotting on.
            shift (int, optional): vertical shift. Default is 0.
            **kwargs: kwargs are passed to ``plt.scatter()`` that plots the data.

        returns:
            matplotlib.axes
        """
        if ax is None:
            fig = figm.figure()
            ax = fig.add_subplot(111)

        if self.offsets is None:
            warnings.warn('Offsets not defined. Use get_offsets().', stacklevel=2)
        else:
            ax.scatter(self.x_centers, self.offsets+shift, **kwargs)
        return ax


    def plot_offsets_fit(self, ax=None,  shift=0, **kwargs):
        """Plot the offsets fitted curve as function of x values.

        args:
            ax (matplotlib.axes, optional): axes for plotting on.
            shift (int, optional): vertical shift. Default is 0.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        returns:
            matplotlib.axes
        """
        if ax is None:
            fig = figm.figure()
            ax = fig.add_subplot(111)

        if self.offsets_func is None:
            warnings.warn('Offsets not defined. Use get_offsets().', stacklevel=2)
        else:
            x = np.linspace(0, self.x_max, 200)
            y = self.offsets_func(x)
            ax.plot(x, y+shift, **kwargs)

        return ax


    def plot_columns(self, ax=None, columns='all', show_ranges=False, vertical_increment=0, **kwargs):
        """Plot columns (intensity as function of y values (center of y bins).

        args:
            ax (matplotlib.axes, optional): axes for plotting on.
            columns (int, string or list, optional): number of the columns to plot.
                It can be a single int or a list of int's. If
                ``columns = 'all'``, all columns are ploted.
            vertical_increment (int, optional): if one column is plotted, it
                adds a vertical offset to the plotted curve. If many columns are ploted,
                ``vertical_increment`` defines
                the vertical offset between each curve. Default is 0.
            show_ranges (bool, optional): show ranges in which offsets were calculated.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        returns:
            matplotlib.axes
        """
        if ax is None:
            fig = figm.figure()
            ax = fig.add_subplot(111)

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'ms' not in kwargs:
            kwargs['ms'] = 5

        if columns == 'all':
            columns = np.arange(0, self.hist.shape[0])

        i = 0
        try:
            if len(columns) == 1:
                ax.plot(self.y_centers, self.hist[columns[0]]+vertical_increment, label=columns[0], **kwargs)
            else:
                for i in range(len(columns)):
                    ax.plot(self.y_centers, self.hist[columns[i]]-i*vertical_increment, label=columns[i], **kwargs)
        except TypeError:
                ax.plot(self.y_centers, self.hist[columns]+vertical_increment, label=columns,**kwargs)

        plt.legend()

        if show_ranges:
            if self.offsets_func is None:
                warnings.warn('Offsets range not defined. Use get_offsets().', stacklevel=2)
            else:
                for r in self.offsets_ranges:
                    plt.axvline(r[0], color='green', linewidth=1.2, zorder=10)
                    plt.axvline(r[1], color='red', linewidth=1.2, zorder=10)

        return ax

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

    # def plot_overlaped(self, ax=None, pointsize=5):
    #
    #     if ax is None:
    #         fig = figm.figure()
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

    # Attributes:
    #     data (list or array): three column list (or array) with photon
    #         events. Column order should be x, y, and intensity.
    #     filepath (string or pathlib.Path): filepath to file.
    #     x_max (float): maximum x value.
    #     y_max (float): maximum y value.
    #
    #     hist (list): data histogram.
    #     bins (list): number of x, y bins.
    #     bins_size (list): size of x, y bins
    #     x_edges (list): edges of x bins.
    #     y_edges (list): edges of y bins.
    #     x_centers (list): center of x bins.
    #     y_centers (list): center of y bins.
    #
    #     offsets         = None
    #     offsets_ref     = None
    #     offsets_ranges  = None
    #     offsets_func       = None
    #     offsets_par     = None
    #
    #     # spectrum attr
    #     spectrum = None

    #
    # The first thing that photon_events does is to call :py:func:`photon_events.load`,
    # which loads photon events data from file or assign data directly. The
    # file/data must have three columns, x, y, and intensity. One can pass the values
    # for x_max and y_max through the file header by
    # using the tag ``# x_max <value>`` and ``# y_max <value>``, respectively.
    #
    # note:
    #     One can assign the data directly by editing the ``data`` attribute, but
    #     it is advised to used the :py:func:`photon_events.load` method
    #     as it adjusts other related attributes (x_max, y_max, filepath).

class spectrum():
    """Creates a ``spectrum`` class type object to deal with (x, y) data types.

    Args:
        filepath (string or pathlib.Path, optional): filepath to file. It overwrites
            the ``data`` argument.
        data (list or array, optional): three column list (or array) with photon
            events. Column order should be x, y, and intensity.
        delimiter (str, optional): The string used to separate values. If whitespaces are used,
            consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).

    """

    def __init__(self, filepath=None, data=None, delimiter=','):
        self.data     = None
        self.x        = None
        self.y        = None
        self.filepath = None

        self.load(filepath=filepath, data=data, delimiter=delimiter)


    def load(self, filepath=None, data=None, delimiter=','):
        """Load spectrum from file or assign data directly.

        The file/data must have two columns, x (energy or distance) and intensity.

        args:
            filepath (string or pathlib.Path, optional): filepath to file. It overwrites
                the ``data`` argument.
            data (list or array, optional): three column list (or array) with photon
                events. Column order should be x, y, and intensity.
            delimiter (str, optional): The string used to separate values. If whitespaces are used,
                consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).

        returns:
            None

        See Also:
            :py:func:`spectrum.save`.
        """

        if filepath is not None:
            self.filepath = Path(filepath)
            self.data = fm.load_data(self.filepath, delimiter=delimiter, force_array=True)

        else:
            self.filepath   = None
            if data is None:
                warnings.warn('No filepath or data to load.', stacklevel=2)
                return
            else:
                self.data = copy.deepcopy(data)

        self.x = self.data[:, 0]
        self.y = self.data[:, 1]


    def save(self, filepath, delimiter=',', header=None):
        r"""Saves spectrum to a file.

        args:
            filepath (string or pathlib.Path, optional): filepath to file.
            delimiter (str, optional): The string used to separate values.
                If whitespaces are used, consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).
            header (string, optional): text to add at the beginning of the file. Use ``\n`` for new line. Comment flag (#) is added automatically.

        returns:
            None

        See Also:
            :py:func:`spectrum.load`.
        """
        fm.save_data(self.data, filepath=Path(filepath), delimiter=delimiter, header=header)


    def apply_correction(self, f):
        """Changes the values of x, y based on a function.

        args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        returns:
            None
        """
        self.data[:, 0], self.data[:, 1] = f(self.data[:, 0], self.data[:, 1])


    def calib(self, dispersion, position_energy_pair=(0, 0), normalized=False):
        """Calibrate data (from length to energy).

        args:
            dispersion (number): dispersion of the diffraction grating in
                units of [energy/lenght].
            position_energy_pair (tuple, optional): a y position and its energy value
                of the isoenergetic line at that position.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.

        returns:
            None
        """

        const = position_energy_pair[1] - dispersion*position_energy_pair[0]

        if normalized:
            f = lambda x, y: (x*dispersion + const, y/max(y))
        else:
            f = lambda x, y: (x*dispersion + const, y)
        self.apply_correction(f=f)


    def interp(self, x=None, start=None, stop=None, num=1000, step=None):
        """Interpolate data.

        args:
            x (list or array, optional): The x-coordinates at which to
                evaluate the interpolated values. This overwrites all other arguments.
            start (number, optional): The starting value of the sequence. If None,
                the minium x value will be used.
            stop (number, optional): The end value of the sequence. If None,
                the maximum x value will be used.
            num (int, optional): Number of samples to generate.
            step (number, optional): Spacing between values. This overwrites ``num``.

        returns:
            None
        """
        if x is None:
            if start is None:
                start = min(self.data[:, 0])
            if stop is None:
                stop = max(self.data[:, 0])

            if step is not None:   # step overwrites num
                # temp = np.arange(start, stop, step=step)
                # temp = np.zeros((len(temp), self.data.shape[1]))
                x = np.arange(start, stop, step=step)
            else:
                # temp = np.zeros((num, self.data.shape[1]))
                x = np.linspace(start, stop, num=num)

            # temp[:, 1] =

        self.data = np.column_stack((x, np.interp(x, self.data[:, 0], self.data[:, 1])))
        # return spectrum(data=temp)

    def apply_shift(self, shift, mode='hard'):
        """Shift data.

        Args:
            shift (float or int): shift value.
            mode (string, optional): If ``mode='x'`` or ``mode='hard'``, y is fully preserved
                while x is shifted. If ``mode='y'``, ``'interp'``, or ``'soft'``, x is preserved
                while y is interpolated with a shift. If ``mode='roll'``, x is also preserved
                and y elements are rolled along the array (``shift`` value must be an integer).

        Warning:
            It is always better to use ``mode='hard'`` or ``'roll'`` since the form of y is fully
            preserved (no interpolation). After applying a shift using the ``mode='interp'``,
            one can apply a
            'inverse' shift to retrieve the original data. The diference between the retrieved
            y data and the original data will give an ideia of the information loss
            caused by the interpolation.

        Returns:
            None
        """
        # print(
        # x = self.data[:, 0]
        # y = self.data[:, 1]
        #
        # if mode == 'y' or mode == 'interp' or mode=='soft':
        #     y = np.interp(x, x + shift, y)
        #
        # elif mode == 'x' or mode == 'hard':
        #     x = np.array(x) + shift
        #
        # elif mode == 'roll' or mode == 'rotate':
        #     y = np.roll(y, shift)
        #     if shift > 0:
        #         y[:shift] = 0
        #     elif shift < 0:
        #         y[shift:] = 0

            # self.apply_correction(self, f)
        x, y = am.shift(self.data[:, 0], self.data[:, 1], shift=shift, mode=mode)
        self.data = np.column_stack((x, y))


    def peak_fit(self, ranges=None, **kwargs):
        r"""Fit a peak with a pseudo-voigt curve.

        .. math:: y(x) = A \left[ m  \frac{w^2}{w^2 + (x-c)^2}   + (1-m) e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

        Args:
            ranges (list): a pair of x values or a list of pairs. Each pair represents
                the start and stop of a data range.
            **kwargs: kwargs are passed to :py:func:`brixs.arraymanip.peak_fit`.

        Note:
            By default, peak assimetry is taken into account.

        Returns:
            1) 2 column (x,y) array with "Smoothed" fitted peak. This is just the
                fitted peak array with a linear interpolation with 100 times more data points.
            2) An array with the optimized parameters for Amplitude, Center, FWHM and offset.
        """

        if ranges is None:
            ranges = [[0, self.y_max]]

        x, y = am.extract(self.data[:, 0], self.data[:, 1], ranges)

        # guess
        if 'guess_A' not in kwargs:       kwargs['guess_A'] = max(y)
        if 'guess_c' not in kwargs:       kwargs['guess_c'] = x[am.index(y, max(y))]
        if 'guess_w' not in kwargs:       kwargs['guess_w'] = (x[1] - x[0])*10
        if 'guess_offset' not in kwargs:  kwargs['guess_offset'] = np.mean(y)
        if 'asymmetry' not in kwargs:     kwargs['asymmetry'] = True
        if 'fixed_m' not in kwargs:         kwargs['fixed_m'] = 0.5


        _, arr100, popt_2 = am.peak_fit(x, y, **kwargs)

        return arr100, popt_2


    def plot(self, ax=None, normalized=True, vertical_increment=0, shift=0, factor=1, **kwargs):
        """Plot spectrum.

        args:
            ax (matplotlib.axes, optional): axes for plotting on.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.
            vertical_increment (float or int, optional): defines the vertical offset. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            show_ranges (bool, optional): show ranges in which offsets were calculated.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        returns:
            matplotlib.axes
        """

        if ax is None:
            fig = figm.figure()
            ax = fig.add_subplot(111)

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'ms' not in kwargs:
            kwargs['ms'] = 5

        if normalized:
            ax.plot((self.data[:, 0] + shift), (self.data[:, 1]/max(self.data[:, 1]))*factor + vertical_increment, **kwargs)
        else:
            ax.plot((self.data[:, 0] + shift), self.data[:, 1]*factor + vertical_increment, **kwargs)

        return ax


class spectra():
    """Creates a ``spectra`` class type object to deal with many spectrum at a time.

    Args:
        folderpath (string or pathlib.Path, optional): path to folder with
            spectra data files. It overwrites the ``data`` argument.
        data (list or array, optional): list of :py:class:`spectrum` objects.
        delimiter (str, optional): The string used to separate values. If whitespaces are used,
            consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).

    """

    def __init__(self, folderpath=None, data=None, delimiter=','):
        # basic attr
        self.spectrum = None
        self.folderpath = folderpath

        # shift attr
        self.shift_mode = None
        self.shifts = None
        self.shift_ranges = None

        # sum attr
        self.sum = None

        self.load(folderpath=folderpath, data=data, delimiter=delimiter)


    def get_spectra_count(self):
        """Returns the number of spectra."""
        return len(self.spectrum)


    def get_filelist(self):
        """Returns a filepath list of all spectra."""
        return [x.filepath for x in self.spectrum]


    def get_specrum_by_filename(self, filename):
        """Returns a idx list of spectra associated with filename.
        """
        filelist = self.get_filelist()
        filelist = [Path(x) if x is not None else None for x in filelist]
        return [idx for idx, s in enumerate(filelist) if filename == filelist.name]


    def append(self, filepath=None, data=None, delimiter=','):
        """Add spectrum to the spectrum list.

        args:
            filepath (string or pathlib.Path, optional): filepath to file. It overwrites
                the ``data`` argument.
            data (spectrum obj, optional): spectrum object to be added.
            delimiter (str, optional): The string used to separate values. If whitespaces are used,
                consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).

        returns:
            None

        See Also:
            :py:func:`spectra.exclude`.
        """
        if filepath is not None:
            self.spectrum.append(spectrum(filepath=filepath, delimiter=','))
        elif data is not None:
            self.spectrum.append(data)
        else:
            warnings.warn('No filepath or data to load.', stacklevel=2)
            return


    def exclude(self, idx):
        """Exclude spectrum from the spectrum list.

        args:
            idx (int): index of the spectrum.

        returns:
            None

        See Also:
            :py:func:`spectra.append`.
        """
        del self.spectrum[idx]


    def save(self, folderpath=None, prefix='', suffix='_spectrum', delimiter=',', header=None):
        r"""Saves spectra in a folder.

        args:
            folderpath (string or pathlib.Path, optional): path to folder.
            delimiter (str, optional): The string used to separate values.
                If whitespaces are used, consecutive whitespaces act as delimiter.
                Use ``\\t`` for tab. The default is comma (,).
            header (string, optional): text to add at the beginning of each file.
                Use ``\n`` for new line. Comment flag (#) is added automatically.

        returns:
            None

        See Also:
            :py:func:`spectra.load`.
        """
        n_digits = figm.n_digits(self.get_spectra_count-1)[0]
        for idx, s in enumerate(self.spectrum):
            filename = f'{prefix}' + f'{idx}'.zfill(n_digits) + f'{suffix}'
            s.save(filepath=folderpath/filename, delimiter=',', header=None)


    def load(self, folderpath=None, data=None, delimiter=','):
        """Load all spectra from folder or assign data directly.

        Each file/data must have two columns, x (energy or distance) and intensity.

        args:
            folderpath (string or pathlib.Path, optional): path to folder with
                spectra data files. It overwrites the ``data`` argument.
            data (list or array, optional): list of :py:class:`spectrum` objects.
            delimiter (str, optional): The string used to separate values. If whitespaces are used,
                consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).

        returns:
            None

        See Also:
            :py:func:`spectra.save`.
        """
        if folderpath is not None:
            self.folderpath = Path(folderpath)
            for filepath in fm.filelist(self.folderpath):
                self.append(filepath=filepath, delimiter=delimiter)

        else:
            self.filepath   = None
            if data is None:
                warnings.warn('No filepath or data to load.', stacklevel=2)
                return
            else:
                self.spectrum = copy.deepcopy(data)


    def calib(self, dispersion, position_energy_pair=(0, 0), normalized=False):
        """Calibrate data (from length to energy).

        args:
            dispersion (number): dispersion of the diffraction grating in
                units of [energy/lenght].
            position_energy_pair (tuple, optional): a y position and its energy value
                of the isoenergetic line at that position.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.

        returns:
            None
        """
        for spectrum in self.spectrum:
            spectrum.calib(dispersion=dispersion, position_energy_pair=position_energy_pair, normalized=normalized)


    def interp(self, x=None, start=None, stop=None, num=1000, step=None):
        """Interpolate data.

        args:
            x (list or array, optional): The x-coordinates at which to
                evaluate the interpolated values. This overwrites all other arguments.
            start (number, optional): The starting value of the sequence. If None,
                the minium x value will be used.
            stop (number, optional): The end value of the sequence. If None,
                the maximum x value will be used.
            num (int, optional): Number of samples to generate.
            step (number, optional): Spacing between values. This overwrites ``num``.

        returns:
            None
        """
        if x is None:
            if start is None:
                start = max([min(s.data[:, 0]) for s in self.spectrum])
            if stop is None:
                stop = min([max(s.data[:, 0]) for s in self.spectrum])

        for spectrum in self.spectrum:
            spectrum.interp(x=x, start=start, stop=stop, num=num, step=step)


    def check_x(self, max_error=0.001):
        """Compare spectra to see if they have same x-coordinates.

            args:
                max_error (number, optional): percentage value of the max error.

            Three checks are performed:

                1) Checks if all spectra have same lenght.

                2) Checks if the x step between two data points is the same through out all x-axis.

                3) checks if the max diference between x arrays of two spectra is
                    less then ``max_error`` percentage of the step between points.

            raises:
                ValueError: If any x-coodinate of any two spectrum is different.
        """

        # check length
        for idx, spectrum in enumerate(self.spectrum):
            try:
                if len(spectrum.data) != len(self.spectrum[idx+1].data):
                    raise ValueError(f"Spectrum {idx} and {idx+1} have the different length.")
            except IndexError:
                pass

        # check step
        for idx, spectrum in enumerate(self.spectrum):
            d = np.diff(spectrum.data[:, 0])
            if (max(d) - min(d))*100/np.mean(np.diff(spectrum.data[:, 0])) > max_error:
                raise ValueError(f"Step in the x-coordinate of spectrum {idx} seems not to be uniform.")

        # check x
        for idx, spectrum in enumerate(self.spectrum):
            try:
                step = spectrum.data[1, 0] - spectrum.data[0, 0]
                if max(abs(spectrum.data[:, 0] - self.spectrum[idx+1].data[:, 0]))*100/step > max_error:
                    raise ValueError(f"Spectrum {idx} and {idx+1} seems to be different.")
                # print(max(abs(spectrum.data[:, 0] - self.spectrum[idx+1].data[:, 0]))*100/step)
            except IndexError:
                pass


    def calculate_shifts(self, ref=0, mode='cross-correlation', ranges=None, check_x=True, verbose=True):
        """Calculate the shift of each spectrum relative to a reference spectrum.

        args:
            ref (int, optional): index of reference spectrum. The shift of all other spectra
                is calculated based on the reference spectrum. Default is 0.
             mode (string, optional): method used to calculate the offsets.
                The current options are: 'cross-correlation', and 'max'.
            ranges (list, optional): a pair of x values or a list of pairs. Each pair represents
                the start and stop of a data range. If None, the whole data set is used.
            check_x (bool, optional): if True, it will check if the x-coordinate
                of all spectra is the same.
            verbose (bool,optional): turn verbose on/off.

        returns:
            None
        """
        if ranges is None:
            ranges = [[min(self.spectrum[ref].data[:, 0]), max(self.spectrum[ref].data[:, 0])]]


        if mode == 'cross-correlation':
            if check_x:
                self.check_x()

            self.shifts = np.zeros(len(self.spectrum))
            _, y_ref = am.extract(self.spectrum[ref].data[:, 0], self.spectrum[ref].data[:, 1], ranges=ranges)
            for i, spectrum in enumerate(self.spectrum):
                _, y = am.extract(spectrum.data[:, 0], spectrum.data[:, 1], ranges=ranges)
                cross_correlation = np.correlate(y_ref, y, mode='Same')
                self.shifts[i] = np.argmax(cross_correlation)

                if verbose:
                    print(f'spectrum {i} shift calculated!')

        elif mode == 'fit':

            for i, spectrum in enumerate(self.spectrum):
                _, popt = spectrum.peak_fit(ranges)
                self.shifts[i] = popt[1]

                if verbose:
                    print(f'spectrum {i} shift calculated!')

        self.shifts -= self.shifts[ref]
        # print(self.shifts)
        # self.shifts = self.shifts*(self.spectrum[ref].data[1, 0] - self.spectrum[ref].data[0, 0])
        # print(self.shifts)
        # print([int(x) for x in self.shifts/(self.spectrum[ref].data[1, 0] - self.spectrum[ref].data[0, 0])])

        self.shift_mode = mode
        self.shift_ranges = ranges


    def shifts_correction(self, mode=None):
        """Shift data.

        Args:
            shift (float or int): shift value.
            mode (string, optional): If None, the best mode will be selected.
                If ``mode='x'`` or ``mode='hard'``, y is fully preserved
                while x is shifted. If ``mode='y'``, ``'interp'``, or ``'soft'``, x is preserved
                while y is interpolated with a shift. If ``mode='roll'``, x is also preserved
                and y elements are rolled along the array (``shift`` value must be an integer).

        Warning:
            It is always better to use ``mode='hard'`` or ``'roll'`` since the form of y is fully
            preserved (no interpolation). After applying a shift using the ``mode='interp'``,
            one can apply a
            'inverse' shift to retrieve the original data. The diference between the retrieved
            y data and the original data will give an ideia of the information loss
            caused by the interpolation.

        Returns:
            None
        """

        if mode is None:
            if self.shift_mode == 'cross-correlation':
                mode = 'roll'
            elif self.shift_mode == 'fit':
                mode = 'soft'
            else:
                warnings.warn(f'Shift mode ({self.shift_mode}) not recognized.', stacklevel=2)
                return

        for i in range(self.get_spectra_count()):
            self.spectrum[i].apply_shift(shift=self.shifts[i], mode=mode)



    def crop(self, start=None, stop=None):
        """Crop spectra ends.

        args:
            start (number, optional): The starting value. If None,
                the minium x value will be used.
            stop (number, optional): The end value. If None,
                the maximum x value will be used.

        returns:
            None
        """
        if start is None:
            start = max([min(s.data[:, 0]) for s in self.spectrum])
        if stop is None:
            stop = min([max(s.data[:, 0]) for s in self.spectrum])

        for spectrum in self.spectrum:
            step = spectrum.data[1, 0] - spectrum.data[0, 0]
            data_range = (start + step/2, stop - step/2)
            a, b = am.extract(spectrum.data[:, 0], spectrum.data[:, 1], ranges=(data_range,) )

            temp = np.zeros((len(a), 2))
            temp[:, 0] = a
            temp[:, 1] = b
            spectrum.data = copy.deepcopy(temp)


    def calculate_sum(self):
        """Sum all spectra."""

        self.check_x()

        temp = copy.deepcopy(self.spectrum[0])
        for i in range(1, self.get_spectra_count()):
            temp.data[:, 1] += self.spectrum[i].data[:, 1]
        self.sum = spectrum(data=temp.data)


    def plot(self, ax=None, idx='all', normalized=True, vertical_increment=0, shift=0, factor=1, show_ranges=False, **kwargs):
        """Plot spectra.

        args:
            ax (matplotlib.axes, optional): axes for plotting on.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.
            vertical_increment (int, optional): if one spectrum is plotted, it
                adds a vertical offset to the plotted curve. If many spectra are ploted,
                ``vertical_increment`` defines
                the vertical offset between each curve. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            show_ranges (bool, optional): show ranges in which shifts were calculated.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        returns:
            matplotlib.axes
        """
        if ax is None:
            fig = figm.figure()
            ax = fig.add_subplot(111)

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'ms' not in kwargs:
            kwargs['ms'] = 5

        if idx == 'all':
            for i in range(len(self.spectrum)):
                self.spectrum[i].plot(ax=ax, normalized=normalized, vertical_increment=-vertical_increment*i, shift=shift, factor=factor, label=i, **kwargs)
        else:
            try:
                if len(idx) == 1:
                    self.spectrum[idx[0]].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, label=idx[0], **kwargs)
                else:
                    for i in range(len(idx)):
                        self.spectrum[i].plot(ax=ax, normalized=normalized, vertical_increment=-vertical_increment*i, shift=shift, factor=factor, label=i, **kwargs)
            except TypeError:
                self.spectrum[idx].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, label=idx, **kwargs)

        plt.legend()


        if show_ranges:
            if self.shifts is None:
                warnings.warn('Shift range not defined. Use calculate_shifts().', stacklevel=2)
            else:
                for r in self.shift_ranges:
                    plt.axvline(r[0], color='green', linewidth=1.2, zorder=10)
                    plt.axvline(r[1], color='red', linewidth=1.2, zorder=10)

        return ax

    #
    # Example:
    #     Simple usage: Gets a photon event list and remove the rotation of the
    #     detector.
    #
    #     >>> import brixs
    #     >>> import matplotlib.pyplot as plt
    #     >>> import numpy as np
    #     >>> plt.ion()
    #     >>> # simulating a generic spectrum
    #     >>> I = brixs.dummy_spectrum(0, 0.2, excitations=[[0.5, 2, 2], [0.5, 4, 2]])
    #     >>> # simulating the photon_event list
    #     >>> data = brixs.dummy_photon_events(I, noise=0.02, background=0.01, y_zero_energy=-20, angle=2)
    #     >>> # initializing photon_events object
    #     >>> p = brixs.photon_events(data=data)
    #     >>> # set binning
    #     >>> p.set_binning((10, 50))
    #     >>> p.plot(show_bins=True)
    #
    #     .. image:: _figs/bins.png
    #         :target: _static/bins.png
    #         :width: 600
    #         :align: center
    #
    #     >>> # plot columns
    #     >>> p.plot_columns(columns='all', shift=100)
    #
    #     .. image:: _figs/columns.png
    #         :target: _static/columns.png
    #         :width: 600
    #         :align: center
    #
    #     >>> # fitting offsets
    #     >>> p.set_binning((10, 1000))
    #     >>> p.calculate_offsets(ranges=[[0,  0.005]])
    #     >>> p.fit_offsets()
    #     >>> p.plot(show_offsets=True, show_offsets_fit=True)
    #
    #     .. image:: _figs/offsets.png
    #         :target: _static/offsets.png
    #         :width: 600
    #         :align: center
    #
    #     .. image:: _figs/offsets_zoom.png
    #         :target: _static/offsets_zoom.png
    #         :width: 600
    #         :align: center
    #
    #     >>> # remove offsets
    #     >>> p.offsets_correction()
    #     >>> p.plot()
    #
    #     .. image:: _figs/final.png
    #         :target: _static/final.png
    #         :width: 600
    #         :align: center
    #
    #     .. image:: _figs/final_zoom.png
    #         :target: _static/final_zoom.png
    #         :width: 600
    #         :align: center
