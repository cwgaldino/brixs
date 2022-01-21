#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module for analysis of RIXS spectra.

This module is based on three class, one for dealing with photon events lists,
one for dealing with single spectra, a another one to deal with many spectra at
a time.


.. autosummary::

    PhotonEvents
    Spectrum
    Spectra

IDEAS:
- when you crop the data in Spectra, you modify the spectrum. There's no way to
go back (like we do with the shift and offset). Maybe figure a away to go back.

"""

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
from collections.abc import Iterable
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

# backpack
from .backpack.filemanip import save_data, save_obj, load_obj
from .backpack.arraymanip import index, moving_average, extract, shifted, sort
from .backpack.figmanip import n_digits
from .backpack.model_functions import gaussian_fwhm, voigt_fwhm

cc = ['cross-correlation', 'cc']
roll = ['roll', 'rotate', 'r', 'rot']
hard = ['hard', 'x', 'h', 'Hard']
soft = ['soft', 'Soft', 'interp', 'y', 's']


class _Meta(type):
    """Metaclass to facilitate creation of read-only attributes."""
    def __new__(self, class_name, bases, attrs):

        def lazy_read_only(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    raise AttributeError('Attribute is "read only". Cannot set attribute.')
                def deleter(self):
                    raise AttributeError('Attribute is "read only". Cannot delete object.')
            return getter, setter, deleter, 'read only attribute'

        def lazy_non_removable(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    return setattr(self, variable, value)
                def deleter(self):
                    raise AttributeError('Attribute cannot be deleted.')
            return getter, setter, deleter, 'non-removable attribute'

        new_attrs = {}
        for name, value in attrs.items():
            if name == '_read_only':
                for attr in value:
                    _property = property(*lazy_read_only(attr))
                    new_attrs[attr] = _property
            elif name == '_non_removable':
                for attr in value:
                    _property = property(*lazy_non_removable(attr))
                    # print(_property.__doc__)
                    new_attrs[attr] = _property
            else:
                new_attrs[name] = value

        return type(class_name, bases, new_attrs)

class PhotonEvents(metaclass=_Meta):
    """Creates a ``PhotonEvents`` class type object to deal with photon events lists.

    Args:
        events (list or array): two (x, y) or three (x, y, intensity) list or array with
            photon events.
            data (list or array, optional): three column list (or array) with photon
                events. Column order should be x, y, (and intensity if applicable).
        x_max (float, optional): maximum x value. If ``None``, it will be infered
            by the data.
        y_max (float, optional): maximum y value. If ``None``, it will be infered
            by the data.

        bins (list): cannot be negative. It is always converted to int values

    """

    _read_only = ['hist', 'x_edges', 'y_edges', 'x_centers', 'y_centers',
                  'offsets_ranges', 'offsets_mode', 'offsets_ref',
                  'offsets_func', 'offsets', 'offsets_popt', 'spectrum', 'spectrum_bins']

    def __init__(self, events, x_max=None, y_max=None):
        # basic attr
        self.events     = events
        # self.events0    = events
        self._x_max     = None
        # self.x_max      = x_max
        self._y_max     = None
        # self.y_max      = y_max

        # binning attr
        self._bins       = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._hist       = None
        self._x_edges    = None
        self._y_edges    = None
        self._x_centers  = None
        self._y_centers  = None

        # offset attr
        self._offsets         = None
        self._offsets_ranges  = None
        self._offsets_mode    = None
        self._offsets_ref     = None

        # offset fit attr
        self._offsets_func    = None
        self._offsets_popt    = None

        # spectrum attr
        self._spectrum_bins = np.array((-1, -1))
        self._spectrum = None

        self.x_max= x_max
        self.y_max= y_max


    @property
    def events(self):
        return self._events
    @events.setter
    def events(self, value):
        if value.shape[1] == 3:
            self._events = np.array([event for event in np.array(value) if not event[2]<=0])
        elif value.shape[1] == 2:
            events = np.array(value)
            self._events = np.c_[events, np.ones(events.shape[0])]
        else: raise ValueError("data must have 2 or 3 columns.")
        x_max = None
        y_max = None
    @events.deleter
    def events(self):
        raise AttributeError('Cannot delete object.')


    @property
    def x_max(self):
        return self._x_max
    @x_max.setter
    def x_max(self, value):
        if value is None:
            _x_max = max(self.events[:, 0])
        elif value > 0 :
             _x_max = value
        else:
            raise ValueError('invalid x_max value.')
        if self.x_max != _x_max:
            self._x_max = _x_max
            if -1 not in self.bins:
                self._binning()
    @x_max.deleter
    def x_max(self):
        self._x_max = max(self.events[:, 0])


    @property
    def y_max(self):
        return self._y_max
    @y_max.setter
    def y_max(self, value):
        if value is None:
            _y_max = max(self.events[:, 1])
        elif value > 0 :
             _y_max = value
        else:
            raise ValueError('invalid y_max value.')
        if self.y_max != _y_max:
            self._y_max = _y_max
            if -1 not in self.bins:
                self._binning()
    @y_max.deleter
    def y_max(self):
        self._y_max = max(self.events[:, 1])


    @property
    def bins(self):
        return self._bins
    @bins.setter
    def bins(self, value):
        if type(value) != str:
            try:
                if len(value) == 1:
                    if value[0] > 0:
                        x_bins, y_bins = int(value[0]), int(value[0])
                    else:
                        raise ValueError('Cannot set negative binning.')
                else:
                    if value[0] > 0 and value[1] > 0:
                        x_bins, y_bins = int(value[0]), int(value[1])
                    else:
                        raise ValueError('Cannot set negative binning.')
            except TypeError:
                if value > 0:
                    x_bins, y_bins = int(value), int(value)
                else:
                    raise ValueError('Cannot set negative binning.')

            if sum(self.bins != np.array((x_bins, y_bins))) > 0:
                self._bins = np.array((x_bins, y_bins))
                self._binning()

        elif value == 'guess':
            guess_bins = self.guess_bins()
            self.bins = guess_bins

        else:
            raise ValueError("Not a valid option ('guess', number or tuple).")
    @bins.deleter
    def bins(self):
        raise AttributeError('Cannot delete object.')


    @property
    def bins_size(self):
        return self._bins_size
    @bins_size.setter
    def bins_size(self, value):
        try:
            if len(value) == 1:
                if value[0] > 0:
                    x_bins_size, y_bins_size = value[0], value[0]
                else:
                    raise ValueError('Cannot set negative binning.')
            else:
                if value[0] > 0 and value[1] > 0:
                    x_bins_size, y_bins_size = value[0], value[1]
                else:
                    raise ValueError('Cannot set negative binning.')
        except TypeError:
            if value > 0:
                x_bins_size, y_bins_size = value, value
            else:
                raise ValueError('Cannot set negative binning.')
        if sum(self.bins_size != np.array((x_bins_size, y_bins_size))) > 0:
            self._bins_size = np.array((x_bins_size, y_bins_size))
            self._binning(use_bins_size=True)
    @bins_size.deleter
    def bins_size(self):
        raise AttributeError('Cannot delete object.')

    def get_bins(self):
        return self.bins


    def set_bins(self, *args):
        if len(args) == 2:
            self.bins = (args[0], args[1])
        elif len(args) == 1:
            self.bins = args[0]
        else:
            raise ValueError('wrong number of arguments.')


    def get_bins_size(self):
        return self.bins_size


    def set_bins_size(self, *args):
        if len(args) == 2:
            self.bins_size = (args[0], args[1])
        elif len(args) == 1:
            self.bins_size = args[0]
        else:
            raise ValueError('wrong number of arguments.')


    def save(self, filepath, delimiter=', ', fmt='%.4f'):
        """Saves photon events data to a file.

        Args:
            filepath (string or pathlib.Path, optional): filepath to file.
            delimiter (str, optional): The string used to separate values. If whitespaces are used,
                consecutive whitespaces act as delimiter. Use ``\\t`` for tab. The default is comma (,).

        Returns:
            None

        note:
            x_max and y_max values are saved in the header.

        See Also:
            :py:func:`photon_events.load`.
        """
        header  = f'x_max {self.x_max}\n'
        header += f'y_max {self.y_max}\n'
        header += f'x y I'
        save_data(self.events, filepath=Path(filepath), delimiter=delimiter, header=header, fmt=fmt)


    def _binning(self, use_bins_size=False):
        if use_bins_size:
            self._bins = np.array((int(self.x_max/self.bins_size[0]), int(self.y_max/self.bins_size[1])))
            # print('binned from bins_size bitch')

        else:
            self._bins_size = np.array((self.x_max/self.bins[0], self.y_max/self.bins[1]))
            # print('binned from bins bitch')

        self._hist, self._x_edges, self._y_edges = np.histogram2d(self.events[:, 0],
                                                                 self.events[:, 1],
                                                                 bins=self.bins,
                                                                 weights=self.events[:, 2],
                                                                 range=((0, self.x_max), (0, self.y_max))
                                                                )
        self._x_centers = moving_average(self.x_edges, n=2)
        self._y_centers = moving_average(self.y_edges, n=2)


    def binning(self, bins=None, bins_size=None):
        """Compute the histogram of the data (binning of the data).

        Args:
            bins (int or tuple, optional): number of bins. If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the x and y directions, respectively.
            bins_size (int or tuple, optional): size of the bins. This overwrites
                the argument ``bins``. If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the x and y directions, respectively.

        example:
            ``bins = 10``
            ``bins = (10)``
            ``bins = (10, 5)``


        return:
            None
        """
        if bins_size is None:
            if bins is None:
                raise ValueError('Must define bins parameter.')
            else:
                try:
                    if len(bins) == 1:
                        x_bins, y_bins = bins[0], bins[0]
                    else:
                        x_bins, y_bins = bins[0], bins[1]
                except TypeError:
                    x_bins, y_bins = bins, bins
                self._bins = np.array((x_bins, y_bins))
                self._binning()
        else:
            try:
                if len(bins_size) == 1:
                    x_bins_size, y_bins_size = bins_size[0], bins_size[0]
                else:
                    x_bins_size, y_bins_size = bins_size[0], bins_size[1]
            except TypeError:
                x_bins_size, y_bins_size = bins_size, bins_size
            self._bins_size = np.array((x_bins_size, y_bins_size))
            self._binning(use_bins_size=True)


    def guess_bins(self, bins_initial=(9, 100), bins_step=(1, 50), max_x_iter=3, max_iter=1000, mode='cross-correlation', ranges=None):
        """
        Args:
            max_x_iter (int, optional): every iteration the algorithm calculates
                the offsets for different x_bins. Sometimes, subsequant x_bins
                will give the same result (i.e., same number of equal offsets).
                The algorithm will keep trying max_x_iter different x_bins until
                finaly change to a different y_bins.
            max_iter (int, optional): maximum number of different binnigs to try.
        """
        # offset parameters
        ref = int(bins_initial[0]/2)

        # inital iteration
        # temp = PhotonEvents(self.events, x_max=self.x_max, y_max=self.y_max)
        # temp.bins = bins_initial
        # temp.calculate_offsets(mode=mode, ranges=ranges, ref=ref)
        diff = [0, 0]
        diff_sum = sum(diff)

        # bins_initial[0] -= 1
        x_counter = 0
        total_counter = 0
        bins = (bins_initial[0], bins_initial[1])
        while 0 in diff and total_counter < max_iter:
            temp = PhotonEvents(self.events, x_max=self.x_max, y_max=self.y_max)
            temp.bins = bins
            temp.calculate_offsets(mode=mode, ranges=ranges, ref=ref)
            diff = np.diff(temp.offsets)

            if x_counter < max_x_iter-1:
                if diff_sum == sum(diff):
                    bins = (bins[0]+bins_step[0], bins[1])
                    x_counter += 1
            else:
                bins = (bins_initial[0], bins[1]+bins_step[1])
                x_counter = 0

            diff_sum = sum(diff)
            total_counter += 1
        if total_counter >= max_iter:
            raise RuntimeError('cannot find solution.')
        else:
            return temp.bins


    def best_bins(self, y_bins=None, x_bins=None, deg=2, spectrum_bins=6000, **kwargs):

        if y_bins is None:
            y_bins = [100, 150, 300, 400, 600, 1000, 1500, 2500, 5000, 6000, 8000, 12000, 15000]

        if x_bins is None:
            x_bins = [3, 5, 7, 9, 12, 15, 20, 30, 60]

        # bins_old = self.bins

        fwhm = {y_bin: {x_bin:None for x_bin in x_bins} for y_bin in y_bins}
        best_bins = (0, 0)
        best = 10000
        print(f'Starting iteration...')
        for i, y_bin in enumerate(y_bins):
            print(f'({i}/{len(y_bins)-1}) Trying y_bin = {y_bin}.')
            for x_bin in x_bins:
                temp = PhotonEvents(self.events, x_max=self.x_max, y_max=self.y_max)
                temp.bins = (x_bin, y_bin)
                # print(self.bins)
                temp.calculate_offsets(ref=int(temp.bins[0]/2))
                temp.fit_offsets(deg=deg)
                temp.offsets_correction()
                s = temp.calculate_spectrum(bins=spectrum_bins)
                try:
                    s.guess_elastic_peak(**kwargs)
                    fwhm[y_bin][x_bin] = s.elastic_fwhm
                    if s.elastic_fwhm < best:
                        best_bins = (x_bin, y_bin)
                        best = s.elastic_fwhm
                except RuntimeError:
                    pass


        print(f'Done. Best bins is {best_bins}.')

        return best_bins, best, fwhm


    def calculate_offsets(self, mode='cross-correlation', ranges=None, ref=0):
        """Calculate the offset of each column relative to a reference column.

        Args:
            ref (int, optional): reference column. The offset of all other columns
                is calculated based on the reference column. Default is 0.
             mode (string, optional): method used to calculate the offsets.
                The current options are: 'cross-correlation', 'max', 'elastic',
                and 'elastic-cross-correlation'.
            ranges (list, optional): a pair of x-coordinate values or a list of pairs. Each pair represents
                the start and stop of a data range. If `None`, the whole data set is used.

        Returns:
            None
        """
        # check if data is binned
        if self.hist is None:
            raise ValueError('Data not binned yet. Use binning()')

        # check ranges
        if ranges is None:
            ranges = [[0, self.y_max]]
        else:
            try:
                for r in ranges:
                    if r[0] > max(self.y_centers) or r[-1] < min(self.y_centers):
                        raise ValueError('Selected ranges outside data range (y_centers).')
            except TypeError:
                if ranges[0] > max(self.y_centers) or ranges[-1] < min(self.y_centers):
                    raise ValueError('Selected ranges outside data range (y_centers).')
                ranges = (ranges, )

        # pre allocate
        offsets = np.zeros(self.hist.shape[0])

        # ref
        y_centers, columns = extract(self.y_centers, self.hist, ranges=ranges)

        # calculate
        for i, _ in enumerate(offsets):
            if mode == 'cross-correlation' or mode == 'cc':
                cross_correlation = np.correlate(columns[i], columns[ref], mode='same')
                offsets[i]   = y_centers[np.argmax(cross_correlation)]

            elif mode == 'max' or mode == 'm':
                offsets[i]   = y_centers[np.argmax(column)]
            elif mode == 'elastic' or mode == 'e':
                raise NotImplementedError('Not implementeed yet.')
            elif mode == 'elastic-cross-correlation' or mode == 'ecc':
                raise NotImplementedError('Not implementeed yet.')
            else:
                raise ValueError("valid mode options are: 'cross-correlation', 'max', 'elastic', and 'elastic-cross-correlation'.")

        self._offsets        =  offsets - offsets[ref]
        self._offsets_ref    = ref
        self._offsets_ranges = ranges
        self._offsets_mode   = mode



    def fit_offsets(self, deg=2, f=None, **kwargs):
        """Find the curve that fits the offset values.

        Args:
            deg (int, optional): degree for the polnomial fit. The default is 1.
            f (function, optional):  a function y = f(x, a, b, ...) that returns the
                value of y as a function of x and other parameters to be fitted.
                This overwrites the standard polynomal fit and ``deg``
                is ignored.
            **kwargs: kwargs are passed to ``scipy.optimize.curve_fit()`` that
                that optimizes the parameters of f. If f is not None, kwargs are
                passed to ``numpy.polyfit()``.

        Returns:
            None
        """
        x2fit = self.x_centers
        y2fit = self.offsets

        if x2fit is None:
            raise ValueError('Data not binned yet. Use binning()')
        if y2fit is None:
            raise ValueError('Offsets not defined. Use calculate_offsets().')

        if f is None:
            if deg < 0:
                raise ValueError('deg must be a positive value or zero.')

            popt = np.polyfit(x2fit, y2fit, deg=deg, **kwargs)
            f = np.poly1d(popt)

        else:
            popt, pcov = curve_fit(f, x2fit, y2fit, **kwargs)
            f = lambda x: f(x, *popt)

        self._offsets_func = f
        self._offsets_popt = popt

    def save_popt(self, filepath='./popt.dat', comments='', check_overwrite=True):
        d = {'popt': self.offsets_popt.tolist(), 'comments': comments}
        save_obj(d, filepath=Path(filepath), check_overwrite=check_overwrite)

    def load_popt(self, filepath='./popt.dat', f=None):
        d = load_obj(filepath=Path(filepath))

        if f is None:
            self._offsets_func = np.poly1d(d['popt'])
            self._offsets_popt = np.array(d['popt'])
        else:
            self._offsets_func = lambda x: f(x, d['popt'])
            self._offsets_popt = np.array(d['popt'])

    def offsets_correction(self, f=None):
        """Uses the offsets fitted curve to adjust the photon events.

        Args:
            f (function, optional): f = f(x). If None, function at self.offsets_func
            will be used.
        """
        if f is None:
            f2 = lambda x, y: (x, y-self.offsets_func(x))
        else:
            f2 = lambda x, y: (x, y-f(x))

        self.apply_correction(f=f2)


    def apply_correction(self, f):
        """Changes the values of x, y based on a function.

        Example:
            f = lambda x, y: (x, y**2)

        Args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        Returns:
            None
        """
        self._events[:, 0], self._events[:, 1] = f(self.events[:, 0], self.events[:, 1])

        try:
            self.binning(bins=self.bins)
        except ValueError:
            pass


    def best_bins_for_spectrum(self, y_bins=None, **kwargs):

        if y_bins is None:
            y_bins = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000]

        fwhm = {y_bin: None for y_bin in y_bins}
        best = 10000
        print(f'Starting iteration...')
        for i, y_bin in enumerate(y_bins):
            print(f'({i}/{len(y_bins)-1}) Trying y_bin = {y_bin}')
            s = self.calculate_spectrum(bins=y_bin, mode='length')
            try:
                s.guess_elastic_peak(**kwargs)
                fwhm[y_bin] = s.elastic_fwhm
                if s.elastic_fwhm < best:
                    best_bins = (1, y_bin)
                    best = s.elastic_fwhm
            except RuntimeError:
                fwhm[y_bin] = -1


        return best_bins, best, fwhm


    def calculate_spectrum(self, bins=None, bins_size=None, mode='bins'):
        """Sum the photon events in the x direction.

        Args:
            y_bins (int, optional): number of y bins. If None, the current binning is used.
            bins_size (int or tuple, optional): size of the y bins. This overwrites
                the argument ``y_bins``. If None, the current binning is used.
            mode (string, optional): Type of the x axis. If 'bins', the x axis
                is numbered from 0 to the number of y_bins. If 'length', the
                detector x values are used.

        mode = length is usefull when comparing data with different binnings.

        Returns:
            spectrum type
        """

        temp = PhotonEvents(events=self.events, x_max=self.x_max, y_max=self.y_max)
        if bins_size is None and bins is not None:
            try:
                if len(bins) == 1:
                    y_bins = bins[0]
                else:
                    y_bins = bins[1]
            except TypeError:
                y_bins = bins
            temp.binning(bins=(1, y_bins))
        elif bins is None and bins_size is not None:
            try:
                if len(bins_size) == 1:
                    y_bins_size = bins_size[0]
                else:
                    y_bins_size = bins_size[1]
            except TypeError:
                y_bins_size = bins_size
            temp.binning(bins_size=(temp.x_max+1, bins_size))
        elif self.bins[0] != -1:
            temp.binning(bins=(1, self.bins[1]))
        else:
            raise ValueError('bining for spectrum calculation not defined.')

        # x
        if mode == 'bins':
            x = np.arange(0, temp.bins[1])
        elif mode == 'length':
            x = temp.y_centers
        else:
            raise ValueError('mode can only be `bins` or `length`.')
        self._spectrum = Spectrum(x=x, y=temp.hist[0])
        self._spectrum_bins = temp.bins[:]

        return self._spectrum


    def plot(self, ax=None, pointsize=1, cutoff=0, show_bins=(False, False), show_offsets=False, show_fit=False, offsets_kwargs={}, bins_kwargs={}, fit_kwargs={}, **kwargs):
        """Plot photon events.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            pointsize (int, optional): photon events point size. Default is 1.
            show_bins (bool, optional): if True, bins edges are displayed in cyan.
                A tuple of bools can also be used to display only x bins or y bins,
                e. g., ``show_bins = (True, False)`` will display only x bins.
                The default is (False, False).
            show_offsets (bool, optional): if True, offsets are displayed over
                its respectively bin (it is positined based on the max calue of
                the first bin). The ranges of data used to calculate
                the offsets are marked by green and red lines.
            show_fit (bool, optional): if True, the fit of the curve
                defined by the offsets values is displayed.
                The ranges of data used to calculate
                the offsets are marked by green and red lines.
            bins_kwargs (dict): kwargs that are passed to ``plt.plot()`` that plots the binning lines.
            fit_kwargs (dict): kwargs that are passed to ``plt.plot()`` that plots the binning lines.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data (photon events).


        Returns:
            matplotlib.axes
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            # ax.set_facecolor('black')

        # kwargs
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'ms' not in kwargs and 'markersize' not in kwargs:
            kwargs['ms'] = pointsize
        if 'mfc' not in kwargs and 'markerfacecolor' not in kwargs:
            kwargs['mfc'] = 'black'
            # kwargs['mfc'] = 'white'
        if 'markeredgewidth' not in kwargs and 'mew' not in kwargs:
            kwargs['markeredgewidth'] = 0
        if 'linewidth' not in kwargs and 'lw' not in kwargs:
            kwargs['linewidth'] = 0

        # bins_kwargs
        if 'linewidth' not in bins_kwargs and 'lw' not in bins_kwargs:
            bins_kwargs['linewidth'] = 0.5
        if 'color' not in bins_kwarg and 'c' not in bins_kwargs:
            bins_kwargs['color'] = 'cyan'
        if 'zorder' not in bins_kwargs:
            bins_kwargs['zorder'] = 3

        # offsets_kwargs
        if 'color' not in offsets_kwargs and 'c' not in offsets_kwargs:
            offsets_kwargs['color'] = 'green'
        if 'zorder' not in offsets_kwargs:
            offsets_kwargs['zorder'] = 4

        # fit_kwargs
        if 'color' not in fit_kwargs and 'c' not in fit_kwargs:
            fit_kwargs['color'] = 'red'
        if 'linewidth' not in fit_kwargs and 'lw' not in fit_kwargs:
            fit_kwargs['linewidth'] = 2
        if 'zorder' not in fit_kwargs:
            fit_kwargs['zorder'] = 5

        # plot
        if cutoff != 0:
            temp = np.array([event for event in self.events if event[2]>cutoff])
            ax.plot(temp[:, 0], temp[:, 1], **kwargs)
        else:
            ax.plot(self.events[:, 0], self.events[:, 1], **kwargs)

        # show_bins
        try:
            if len(show_bins) == 1:
                show_bins = (show_bins[0], show_bins[0])
        except TypeError:
            show_bins = (show_bins, show_bins)
        if show_bins != (False, False):
            if self.x_edges is not None:
                if show_bins[0]:
                    plt.vlines(self.x_edges, 0, self.y_max, **bins_kwargs)
                if show_bins[1]:
                    plt.hlines(self.y_edges, 0, self.x_max, **bins_kwargs)
            else:
                raise ValueError('Data not binned yet. Cannot show bins. Use binning()')


        # show_offsets
        if show_offsets:
            if self.offsets is None:
                raise ValueError('Offsets not defined. Cannot show offsets. Use calculate_offsets().')
            else:
                x2, y2 = extract(self.y_centers, self.hist[self.offsets_ref], ranges=self.offsets_ranges)
                c = x2[np.argmax(y2)]
                self.plot_offsets(ax=ax, shift=(-self.offsets[self.offsets_ref] + c), **offsets_kwargs)

        # show offsets fit
        if show_fit:
            if self.offsets is None:
                raise ValueError('Offsets not fitted. Use fit_offsets().')
            else:
                x2, y2 = extract(self.y_centers, self.hist[self.offsets_ref], ranges=self.offsets_ranges)
                c = x2[np.argmax(y2)]
                self.plot_fit(ax=ax,  shift=(-self.offsets[self.offsets_ref] + c), **fit_kwargs)

        if show_offsets or show_fit:
            for r in self.offsets_ranges:
                plt.axhline(r[0], color='green', linewidth=2, zorder=10)
                plt.axhline(r[1], color='red', linewidth=2, zorder=10)

        plt.xlim(0, self.x_max)
        plt.ylim(0, self.y_max)

        return ax


    def plot_offsets(self, ax=None, shift=0, **kwargs):
        """Plot offsets as function of x values (center of x bins).

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            shift (int, optional): vertical shift. Default is 0.
            **kwargs: kwargs are passed to ``plt.scatter()`` that plots the data.

        Returns:
            matplotlib.axes
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if self.offsets is None:
            raise ValueError('Offsets not defined. Use get_offsets().')
        else:
            ax.scatter(self.x_centers, self.offsets+shift, **kwargs)
        return ax


    def plot_fit(self, ax=None,  shift=0, **kwargs):
        """Plot the offsets fitted curve as function of x values.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            shift (int, optional): vertical shift. Default is 0.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if self.offsets_func is None:
            raise ValueError('Offsets not defined. Use calculate_offsets().')
        else:
            x = np.linspace(0, self.x_max, 200)
            y = self.offsets_func(x)
            ax.plot(x, y+shift, **kwargs)

        return ax


    def plot_columns(self, ax=None, columns='all', show_ranges=False, vertical_increment=0, **kwargs):
        """Plot columns (intensity as function of y values (center of y bins).

        Args:
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

        Returns:
            matplotlib.axes
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if self.hist is not None:
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
        else:
            raise ValueError('data not binned yet. Please, bin data.')

def peak_function(asymmetry, fixed_m, idx=0):
    if fixed_m == False and type(fixed_m) == bool:  # variable m
        if asymmetry:
            def function2fit(x, amp, c, w1, m1, w2, m2):
                f = np.heaviside(x-c, 0)*voigt_fwhm(x, amp, c, w1, m1) +\
                    np.heaviside(c-x, 0)*voigt_fwhm(x, amp, c, w2, m2)
                return f
            f_str = f'np.heaviside(x-c_{idx}, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w1_{idx}, m1_{idx}) + np.heaviside(c_{idx}-x, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w2_{idx}, m2_{idx})'
            args_str = f'amp_{idx}, c_{idx}, w1_{idx}, m1_{idx}, w2_{idx}, m2_{idx}'
        else:
            def function2fit(x, amp, c, w, m, offset):
                return voigt_fwhm(x, amp, c, w, m)
            f_str = f'voigt_fwhm(x, amp_{idx}, c_{idx}, w_{idx}, m_{idx})'
            args_str = f'amp_{idx}, c_{idx}, w_{idx}, m_{idx}'
    else:
        if fixed_m > 1:
            m = 1
        elif fixed_m < 0:
            m = 0
        else:
            m = fixed_m
        if asymmetry:
            def function2fit(x, A, c, w1, w2):
                f = np.heaviside(x-c, 0)*voigt_fwhm(x, A, c, w1, fixed_m) +\
                    np.heaviside(c-x, 0)*voigt_fwhm(x, A, c, w2, fixed_m)
                return f
            f_str = f'np.heaviside(x-c_{idx}, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w1_{idx}, {m}) + np.heaviside(c_{idx}-x, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w2_{idx}, {m})'
            args_str = f'amp_{idx}, c_{idx}, w1_{idx}, w2_{idx}'
        else:
            def function2fit(x, A, c, w):
                return voigt_fwhm(x, A, c, w, fixed_m)
            f_str = f'voigt_fwhm(x, amp_{idx}, c_{idx}, w_{idx}, {m})'
            args_str = f'amp_{idx}, c_{idx}, w_{idx}'
    return function2fit, f_str, args_str

class Spectrum(metaclass=_Meta):
    """Creates a ``spectrum`` class type object to deal with (x, y) data types.

    Args:
        data (list or array, optional): two column list (or array).
        x (list or array, optional): x values (1D list/array). Overwrites `data`.
        y (list or array, optional): y values (1D list/array). Overwrites `data`.

    Attributes:
        data (array): three column array (x, y, and intensity)
        x (array): vector with x-coordinate values
        y (array): vector with y-coordinate values
        step (number, read only): step size for x-coordinates.
        shift( number): shift value (value will be added to x-coordinates).
        calib (number): calibration value (x-coordinates will be multiplied by this value).
        offset( number): offset value (value will be added to y-coordinates).
        factor (number): ultiplicative factor (y-coordinates will be multiplied by this value).
        peaks (dict): each entry represents a peak.
        fit_data (dict): fit data for each peak.
        fit_popt (list): optimized parameters of the full fit (only for inspection).
        fit_func (function): function of the fitted peaks.
        fit_ranges (list): data range used to fit the peaks.

    """

    _read_only = ['peaks', 'shift_mode', 'step']
    _non_removable = ['calib', 'fit', 'fit_data', 'fit_popt', 'fit_func',
                      'fit_ranges',]

    def __init__(self, data=None, x=None, y=None):
        # self._calib = 1
        # self._shift = 0
        # self._shift_mode = None
        # self._offset = 0
        # self._factor = 1

        if x is None and y is None:
            self.data = data
        else:
            if x is None:
                self._x = np.arange(0, len(y))
                self._step = 1
            else:
                if len(x) == len(y):
                    self._x = x
                    self._step = None
                else:
                    raise ValueError('x and y data are not the same length.')
            self._y = y
            self._data = np.vstack((x, y)).transpose()
            self._restart_attr()

    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        try:
            if value.shape[1] != 2:
                raise ValueError('Data must have two columns (x, y).')
            elif value is None:
                raise ValueError('No data to load.')
            elif value.shape[1] == 1:   # one column data (list)
                self._x = np.arange(0, len(value))
                self._y = value
                self._data = np.vstack((self.x, self.y)).transpose()
                self._step = 1
                self._restart_attr()
            else:
                self._data = value
                self._x = self.data[:, 0]
                self._y = self.data[:, 1]
                self._step = None
                self._restart_attr()
        except IndexError:  # one column data (list)
            self._x = np.arange(0, len(value))
            self._y = value
            self._data = np.vstack((self.x, self.y)).transpose()
            self._step = 1
            self._restart_attr()
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, value):
        if len(value) != len(self.x):
            raise ValueError('length of x is not compatible with data.')
        else:
            self._x = value
            self._data[:, 0] = value
            self._step = None
            self._restart_attr()
            if self.calib != 1:
                print('Calibration not 1. Applaying calibration to new data.')
                f = lambda x, y: (x*self.calib, y)
                self.apply_correction(f=f)
    @x.deleter
    def x(self):
        raise AttributeError('Cannot delete object.')

    @property
    def y(self):
        return self._y
    @y.setter
    def y(self, value):
        if len(value) != len(self.y):
            raise ValueError('length of y is not compatible with data.')
        else:
            self._y = value
            self._data[:, 1] = value
            self._restart_attr()
    @y.deleter
    def y(self):
        raise AttributeError('Cannot delete object.')

    @property
    def calib(self):
        return self._calib
    @calib.setter
    def calib(self, value):
        self.set_calib(value)
    @calib.deleter
    def calib(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift(self):
        return self._shift
    @shift.setter
    def shift(self, value):
        self.set_shift(value, mode='x')
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        return self._offset
    @offset.setter
    def offset(self, value):
        self.set_offset(value)
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        return self._factor
    @factor.setter
    def factor(self, value):
        self.set_factor(value)
    @factor.deleter
    def factor(self):
        raise AttributeError('Cannot delete object.')

    def set_x(self, x):
        """Set x-coordinates."""
        self.x = x

    def get_x(self):
        """Returns x-coordinates."""
        return self.x

    def set_y(self, y):
        """Set y-coordinates."""
        self.y = y

    def get_y(self):
        """Returns y-coordinates."""
        return self.y

    def _restart_attr(self):
        """Set relevant attributes to initial value."""
        self._calib = 1
        self._shift = 0
        self._shift_mode = None
        self._offset = 0
        self._factor = 1

        self.fit = None
        self.fit_data     = None
        self.fit_popt     = None
        self.fit_func     = None
        self.fit_ranges   = None

        self._peaks = None

    def save(self, filepath, delimiter=', ', header='', fmt='%.8f'):
        r"""Saves spectrum to a file.

        Args:
            filepath (string or pathlib.Path, optional): filepath to file.
            delimiter (str, optional): The string used to separate values.
                If whitespaces are used, consecutive whitespaces act as delimiter.
                Use ``\\t`` for tab. The default is comma (', ').
            header (string, optional): text to add at the beginning of each file.
                Use ``\n`` for new line. Comment flag ('# ') is added automatically.
            fmt (string, or list, optional): format for saving data.
                If string, the value is used for x- and y-coordinates. If tuple
                of strings, the first string is used for x-coordinates and the
                second for y-coordinates.

                    fmt = (%[flag]width[.precision]specifier)

                * `flag` can be:
                    1. '-' for left justify
                    2. '+', puts + or - in front of the numbers
                    3. '0' to Left pad the number with zeros instead of space (see width).

                * `width` is the minimum number of characters to be printed.

                * `precision` is the number of significant digits.

                * `specifier` is the type of notation. Tipically, either 'e' for
                scientific notation of 'f' for decimal floating point.

                * a common `fmt` strings is: '%.3f' for 3 decimal places.

                *  for more information see `np.savetxt <https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html?highlight=savetxt#numpy.savetxt>`_ documentation::

        Returns:
            None

        """
        save_data(self.data, filepath=Path(filepath), delimiter=delimiter, header=header, fmt=fmt)

    def check_step_x(self, max_error=0.1):
        """Checks vector uniformity of the x-coordinates.

            If the step between two data points is the same through out the
            x vector, it sets :py:attr:`step` with the value of the step size.

            Args:
                max_error (number, optional): percentage (of the x step) value of the max error.

            Returns:
                None

            Raises:
                ValueError: If x-ccordinates are not uniform.
        """

        # check step uniformity
        d = np.diff(self.x)
        if (max(d) - min(d))*100/np.mean(np.diff(self.x)) > max_error:
            self._step = None
            raise ValueError(f"Step in the x-coordinate seems not to be uniform.")

        self._step = np.mean(d)

    def find_peaks(self, prominence=None, points=4, moving_average_window=8, ranges=None):
        """Find peaks using scipy.signal.find_peaks() function.

        Sets :py:attr:`peaks` attribute (sets it to None if no peak is found).

        Args:
            prominence (number, optional): minimum prominence of peaks. This
                parameter will be passed directly to the ``scipy.signal.find_peaks()``
                function. If ``None``, prominence is set as 5 % of the max intensity
                value (which should be fine for most RIXS spectra).
            points (number, optional): minimum number of peaks defining a peak
                (width).
            moving_average_window (int, optional): window size for smoothing the data for finding peaks.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.

        Returns:
            None

        Raises:
            ValueError: if points or moving_average_window have an improper value.
        """
        # check points and moving_average_window
        if points < 1:
            raise ValueError('points must be 1 or higher.')
        if isinstance(points, int) == False:
            if points.is_integer() == False:
                raise ValueError('points must be an integer.')
        if moving_average_window < 1:
            raise ValueError('moving_average_window must be 1 or higher.')
        if isinstance(moving_average_window, int) == False:
            if moving_average_window.is_integer() == False:
                raise ValueError('moving_average_window must be an integer.')

        # step
        if self.step is None:
            self.check_step_x()

        # extracting data
        if ranges is None:
            ranges = [[min(self.x), max(self.x)]]
            x = self.x
            y = self.y
        else:
            if isinstance(ranges[0], Iterable) == False:
                ranges = (ranges, )
            x, y = extract(self.x, self.y, ranges=ranges)

        # data smoothing
        if moving_average_window > 1:
            y2 = moving_average(y, moving_average_window)
            x2 = moving_average(x, moving_average_window)
        else:
            x2 = x[:]
            y2 = y[:]

        # parameters
        if prominence is None:
            prominence = max(y2)*0.05

        # find peaks
        try:
            peaks, d = find_peaks(y2, prominence=prominence, width=points)
            self._peaks = {i: {'amp': d['prominences'][i], 'c': x2[peaks[i]],  'fwhm': d['widths'][i]*self.step} for i in range(len(peaks))}
        except IndexError:
            peaks = []

        self._check_peak_order()

    def _check_peak_order(self):
        if self.peaks == []:
            return
        else:
            l = max(self.peaks.keys())+1
            c = [0]*l
            for i in self.peaks:
                c[i] = self.peaks[i]['c']
            if sum(1 for test in np.diff(c) if test < 0) > 0:
                i_new = sort(c, [i for i in range(l)])
                print(i_new)
                peaks_old = copy.deepcopy(self.peaks)
                for i in range(l):
                    self._peaks[i] = peaks_old[i_new[i]]

    def append_peak(self, amp, c, fwhm):
        l = max(self.peaks.keys())+1
        self._peaks[l] = {'amp': amp, 'c': c,  'fwhm': fwhm}
        self._check_peak_order()

    def remove_peak(self, idx):
        del self._peaks[idx]

    def fit_peak(self, idx='all', ranges=None, asymmetry=False, fixed_m=0, multiplicity=1, offset=True):
        """

        Args:
            idx (number, list, or 'all'): idx of the peak to fit.
            ranges (list, optional) if None, ranges is automatically set, based
                on the first and last peak to be fitted (given on `idx`).
                If ranges is not None, only data inside ranges is considered in
                the fitting. Also, it raises an error if any peak in idx is
                outside of ranges.
            assymetry (bool or list of bools, optional): if True, fits each half of the
                with a different width.
            fixed_m (False, number, or list, optional): m is the amount of lorentzian
                contribution for a peak. If False, m will be fitted for each peak.
                If a number (from 0 to 1), this will be used as the value of m
                (fixed m).
            multiplicity (int or list, optional): use this if peaks have a shoulder
                or some spliting. The peak will be fitted with `multiplicity` peaks.
            offset (bool, optional): if True, a offset value will be fitted.

        Returns:
            None

        """
        # check if peaks is defined ============================================
        if self.peaks is None:
            raise ValueError('Spectrum.peaks is not defined. Run Spectrum.find_peaks().')
        if self.peaks == []:
            raise ValueError('Spectrum.peaks = []. No peaks to fit.')

        # fix idx ==============================================================
        if idx == 'all':
            idx = [k for k in self.peaks]
        if isinstance(idx, Iterable) == False:
            idx = [idx, ]
        idx = [i if i>=0 else len(self.peaks)+i for i in idx]

        # check multiplicity ======================================================
        if isinstance(multiplicity, dict):
            for i in idx:
                if i not in multiplicity:
                    multiplicity[i] = 1
                else:
                    if multiplicity[i] < 1 and type(multiplicity) != int:
                        raise ValueError('multiplicity must be a positive integer.')
            for i in multiplicity:
                if i not in idx:
                    raise ValueError(f'Peak {i} not in idx (idx = {idx}) (multiplicity = {multiplicity}).')
        else:
            if isinstance(multiplicity, Iterable):
                raise ValueError('multiplicity must be a number or dictionary.')
            multiplicity = {i: multiplicity for i in idx}

        # check asymmetry ======================================================
        if isinstance(asymmetry, dict):
            for i in idx:
                if i not in asymmetry:
                    asymmetry[i] = False
                else:
                    if type(asymmetry[i]) != bool:
                        raise ValueError('asymmetry must be bool (True or False)')
            for i in asymmetry:
                if i not in idx:
                    raise ValueError(f'Peak {i} not in idx (idx = {idx}) (asymmetry = {asymmetry}).')
        else:
            if isinstance(asymmetry, Iterable):
                raise ValueError('asymmetry must be a number or dictionary.')
            asymmetry = {i: asymmetry for i in idx}

        # check fixed_m ========================================================
        if isinstance(fixed_m, dict):
            for i in idx:
                if i not in fixed_m:
                    fixed_m[i] = 0
                else:
                    if type(fixed_m[i]) == False or (fixed_m[i] > 0 and fixed_m[i] < 1):
                        raise ValueError('fixed_m must be a False or a number from 0 to 1.')
            for i in fixed_m:
                if i not in idx:
                    raise ValueError(f'Peak {i} not in idx (idx = {idx}) (fixed_m = {fixed_m}).')
        else:
            if isinstance(fixed_m, Iterable):
                raise ValueError('fixed_m must be a number, False, or dictionary.')
            fixed_m = {i: fixed_m for i in idx}

        # fix ranges ===========================================================
        if ranges is None:    # try fitting assigned peaks (range is automatically set)
            ranges = [min(self.x),  max(self.x)]
            p_max = self.peaks[max(idx)]['c']+self.peaks[max(idx)]['fwhm']
            p_min = self.peaks[min(idx)]['c']-self.peaks[min(idx)]['fwhm']

            if max(idx)+1 < len(self.peaks):
                next_min = self.peaks[max(idx)+1]['c']-self.peaks[max(idx)+1]['fwhm']
                if p_max > next_min:
                    ranges[1] = p_max + self.peaks[max(idx)]['fwhm']/2
                else:
                    x_temp, y_temp = extract(self.x, self.y, (p_max, next_min))
                    ranges[1] = x_temp[np.argmin(y_temp)]
            else:
                next_min = max(self.x)
                ranges[1] = max(self.x)

            if min(idx)-1 >= 0:
                previous_max = self.peaks[min(idx)-1]['c']-self.peaks[min(idx)-1]['fwhm']
                if p_min < previous_max:
                    ranges[0] = p_min - self.peaks[min(idx)]['fwhm']/2
                else:
                    x_temp, y_temp = extract(self.x, self.y, (previous_max, p_min))
                    ranges[0] = x_temp[np.argmin(y_temp)]
            else:
                previous_max = min(self.x)
                ranges[0] = min(self.x)

            ranges = (ranges, )
            x2fit, y2fit = extract(self.x, self.y, ranges=ranges)
        else:
            if isinstance(ranges[0], Iterable) == False:
                ranges = (ranges, )
            x2fit, y2fit = extract(self.x, self.y, ranges=ranges)

            # check if peaks are inside ranges
            for i in idx:
                flag = True
                for r in ranges:
                    if self.peaks[i]['c'] > r[0] and self.peaks[i]['c'] < r[1]:
                        flag = False
                if flag:
                    raise ValueError(f'peak {i} outside of range.')

        # fit function =====================================================
        f_temp = {}

        model_str = ''
        args_str = ''
        bounds_min = []
        bounds_max = []
        p0 = []

        total_i = -1
        for i in idx:
            for extra in range(multiplicity[i]):
                total_i += 1
                f, f_str, a_str = peak_function(asymmetry=asymmetry[i], fixed_m=fixed_m[i], idx=total_i)
                if fixed_m[i] == False and type(fixed_m[i]) == bool:  # variable m
                    if asymmetry[i]:
                        n_args = 6
                        p0 = np.append(p0, [self.peaks[i]['amp'], self.peaks[i]['c']+self.peaks[i]['fwhm']/4*extra, self.peaks[i]['fwhm'], 0.5, self.peaks[i]['fwhm'], 0.5])
                        bounds_min = np.append(bounds_min, [            0,          self.peaks[i]['c']-self.peaks[i]['fwhm'],            0,          0,       0,               0])
                        bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*3, self.peaks[i]['c']+self.peaks[i]['fwhm'], self.peaks[i]['fwhm']*4, 1, self.peaks[i]['fwhm']*4, 1])
                    else:
                        n_args = 4
                        p0 = np.append(p0, [self.peaks[i]['amp'], self.peaks[i]['c']+self.peaks[i]['fwhm']/4*extra, self.peaks[i]['fwhm'], 0.5])
                        bounds_min = np.append(bounds_min, [            0,          self.peaks[i]['c']-self.peaks[i]['fwhm'],            0,          0])
                        bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*3, self.peaks[i]['c']+self.peaks[i]['fwhm'], self.peaks[i]['fwhm']*4, 1])
                else:
                    if asymmetry[i]:
                        n_args = 4
                        p0 = np.append(p0, [self.peaks[i]['amp'], self.peaks[i]['c']+self.peaks[i]['fwhm']/4*extra, self.peaks[i]['fwhm'], self.peaks[i]['fwhm']])
                        bounds_min = np.append(bounds_min, [            0,          self.peaks[i]['c']-self.peaks[i]['fwhm'],            0,                0,             ])
                        bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*3, self.peaks[i]['c']+self.peaks[i]['fwhm'], self.peaks[i]['fwhm']*4, self.peaks[i]['fwhm']*4])
                    else:
                        n_args = 3
                        p0 = np.append(p0, [self.peaks[i]['amp'], self.peaks[i]['c']+self.peaks[i]['fwhm']/4*extra, self.peaks[i]['fwhm']])
                        bounds_min = np.append(bounds_min, [            0,          self.peaks[i]['c']-self.peaks[i]['fwhm'],            0,        ])
                        bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*3, self.peaks[i]['c']+self.peaks[i]['fwhm'], self.peaks[i]['fwhm']*4])
                f_temp[total_i] = {'asymmetry':asymmetry[i], 'fixed_m':fixed_m[i], '_func':f, 'n_args':n_args}
                model_str += f_str + ' + '
                args_str += a_str + ', '

        if offset:
            model_str = f'lambda x, {args_str}offset: {model_str}offset'
            p0 = np.append(p0, min(y2fit))
            bounds_min = np.append(bounds_min, min(y2fit))
            bounds_max = np.append(bounds_max, np.mean(y2fit))
        else:
            model_str = f'lambda x, {args_str[:-2]}: {model_str[:-3]}'
        model = eval(model_str)

        # fit peaks ============================================================
        popt, pcov = curve_fit(model, x2fit, y2fit, p0=p0, bounds = (bounds_min, bounds_max))

        # final ============================================================
        x_temp = np.linspace(self.x[0], self.x[-1], len(x2fit)*3)
        stop = 0
        for i in range(total_i+1):
            start = stop
            stop  = start + f_temp[total_i]['n_args']
            f_temp[i]['popt'] = popt[start:stop]
            f_temp[i]['pcov'] = pcov[start:stop]
            f_temp[i]['amp'] = popt[start]
            f_temp[i]['c'] = popt[start+1]
            if f_temp[i]['fixed_m'] == False and type(f_temp[i]['fixed_m']) == bool:  # variable m
                if f_temp[i]['asymmetry']:
                    f_temp[i]['fwhm'] = popt[start+2]/2 + popt[start+4]/2
                    f_temp[i]['fwhm1'] = popt[start+2]
                    f_temp[i]['fwhm2'] = popt[start+4]
                    f_temp[i]['m1'] = popt[start+3]
                    f_temp[i]['m2'] = popt[start+5]
                else:
                    f_temp[i]['fwhm'] = popt[start+2]
                    f_temp[i]['m'] = popt[start+3]
            else:
                if f_temp[i]['asymmetry']:
                    f_temp[i]['fwhm'] = popt[start+2]/2 + popt[start+3]/2
                else:
                    f_temp[i]['fwhm'] = popt[start+2]
                f_temp[i]['m'] = f_temp[i]['fixed_m']

            if offset:
                def peak(*args):
                    return lambda x: (f_temp[i]['_func'](x/self.calib - self.shift, *args) + popt[-1])*self.factor + self.offset
                f_temp[i]['func'] = peak(*f_temp[i]['popt'])
            else:
                def peak(*args):
                    return lambda x: (f_temp[i]['_func'](x/self.calib - self.shift, *args))*self.factor + self.offset
                f_temp[i]['func'] = peak(*f_temp[i]['popt'])
            f_temp[i]['curve'] = Spectrum(x=x_temp, y=f_temp[i]['func'](x_temp))


        self.fit_data = f_temp
        self.fit_popt = popt
        self.fit_func = lambda x: (model(x/self.calib - self.shift, *popt))*self.factor + self.offset
        self.fit_ranges = ranges
        self.fit = Spectrum(x=x_temp, y=self.fit_func(x_temp))

    def _multiplicative_x_fix(self, value):
        # peaks =================================
        if self.peaks is not None:
            for i in self.peaks:
                self.peaks[i]['c'] *= value
                self.peaks[i]['fwhm'] *= value

        # fit_popt ==============================
        self.fit_popt = None  # cannot fix

        # fit_data =============================
        if self.fit_data is not None:
            for i in self.fit_data:
                self.fit_data[i]['popt'] = None # cannot fix
                self.fit_data[i]['c'] *= value
                self.fit_data[i]['fwhm'] *= value
                if 'fwhm1' in self.fit_data[i]:
                    self.fit_data[i]['fwhm1'] *= value
                if 'fwhm2' in self.fit_data[i]:
                    self.fit_data[i]['fwhm2'] *= value
                self.fit_data[i]['curve'].calib = value

        # fit_ranges =============================
        if self.fit_ranges is not None:
            self.fit_ranges = [(r[0]*value, r[1]*value) for r in self.fit_ranges]

    def _additive_x_fix(self, value, mode):
        if mode == 'roll' or mode == 'r' or mode == 'rotate':
            value = self.step*value

        # peaks =================================
        if self.peaks is not None:
            for i in self.peaks:
                self.peaks[i]['c'] += value

        # fit_popt ==============================
        self.fit_popt = None  # cannot fix

        # fit ==================================
        if self.fit is not None:
            self.fit.shift = value

        # fit_data =============================
        if self.fit_data is not None:
            for i in self.fit_data:
                self.fit_data[i]['popt'] = None # cannot fix
                self.fit_data[i]['c'] += value
                self.fit_data[i]['curve'].shift = value

        # fit_ranges =============================
        if self.fit_ranges is not None:
            self.fit_ranges = [(r[0]+value, r[1]+value) for r in self.fit_ranges]

    def _additive_y_fix(self, value):
        # peaks =================================
        if self.peaks is not None:
            for i in self.peaks:
                self.peaks[i]['amp'] += value

        # fit_popt ==============================
        self.fit_popt = None  # cannot fix

        # fit ==================================
        if self.fit is not None:
            self.fit.offset = value

        # fit_data =============================
        if self.fit_data is not None:
            for i in self.fit_data:
                self.fit_data[i]['popt'] = None # cannot fix
                self.fit_data[i]['amp'] += value
                self.fit_data[i]['curve'].offset = value

    def _multiplicative_y_fix(self, value):
        # peaks =================================
        if self.peaks is not None:
            for i in self.peaks:
                self.peaks[i]['amp'] *= value

        # fit_popt ==============================
        self.fit_popt = None  # cannot fix

        # fit ==================================
        self.fit.factor = value

        # fit_data =============================
        if self.fit_data is not None:
            for i in self.fit_data:
                self.fit_data[i]['popt'] = None # cannot fix
                self.fit_data[i]['amp'] *= value
                self.fit_data[i]['curve'].factor = value

    def set_calib(self, value):
        """Calibrate data.

        Args:
            value (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].

        Returns:
            None
        """
        if self.calib != value:
            if self.calib != 1:
                self._x = self.x*self.calib**-1
                self.data[:, 0] = self._x
                self._multiplicative_x_fix(self.calib**-1)
            if value != 1:
                self._x = self.x*value
                self.data[:, 0] = self._x
                self._multiplicative_x_fix(value)
            self._calib = value

    def set_shift(self, value, mode):
        """Shift data.

        Args:
            value (float or int): shift value.
            mode (string, optional): If ``mode='x'`` or ``mode='hard'``, y is fully preserved
                while x is shifted. If ``mode='y'``, ``'interp'``, or ``'soft'``, x is preserved
                while y is interpolated with a shift. If ``mode='roll'`` (or rotate or r), x is also preserved
                and y elements are rolled along the array (``shift`` value must be an integer).
                The form of y-coordinates is fully preserved, but the edge of y-coordinates are lost.

        Returns:
            None

        Warning:
            It is always better to use ``mode='hard'`` or ``'roll'`` since the form of y is fully
            preserved (no interpolation). After applying a shift using the ``mode='interp'``,
            one can apply a
            'inverse' shift to retrieve the original data. The diference between the retrieved
            y data and the original data will give an ideia of the information loss
            caused by the interpolation. Also, multiple shifts using interp may cause
            built up of data loss. It's not drastic if data is suficiently well
            resolved, but it should be avoided.
        """
        if mode in ['roll', 'r', 'rotate', 'rot'] and self.step is None:
            try:
                self.check_step_x()
            except ValueError:
                raise ValueError(f'Cannot shift data using mode = {mode}, because x-coordinates are not uniform.')

        if self.shift != value:
            if self.shift != 0:
                self._x, self._y = shifted(self.x, self.y, value=-self.shift, mode=self.shift_mode)
                self._additive_x_fix(self.shift, mode=self.shift_mode)
            if value != 0:
                self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                self._shift_mode = mode
                self._additive_x_fix(value, mode=mode)
            else:
                self._shift_mode = None

            self._shift = value
            self._data[:, 0] = self._x
            self._data[:, 1] = self._y

    def set_offset(self, value):
        if self.offset != value:
            if self.offset != 0:
                self._y = self.y - self.offset
                self.data[:, 1] = self._y
                self._additive_y_fix(-self.offset)
            if value != 0:
                self._y = self.y + value
                self.data[:, 1] = self._y
                self._additive_y_fix(value)
            self._offset = value

    def set_factor(self, value):
        if self.factor != value:
            if self.factor != 1:
                self._y = self.y*self.factor**-1
                self.data[:, 1] = self._y
                self._multiplicative_y_fix(self.factor**-1)
            if value != 1:
                self._y = self.y*value
                self.data[:, 1] = self._y
                self._multiplicative_y_fix(value)
            self._factor = value

    def apply_correction(self, f, verbose=True):
        """Changes the values of x, y based on a function.

        Example:
            f = lambda x, y: (x, y**2)

        Args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        Returns:
            None

        Warning:
            Spectra.calib is set to 1.
            Spectra.step is recalculated if possible.
            Fit information is lost.
        """
        # core ================================
        self._x, self._y = f(self.x, self.y)
        self._data[:, 0] = self._x
        self._data[:, 1] = self._y

        self._restart_attr()
        if self.step is not None:
            try:
                self.check_step_x()
            except ValueError:
                pass

    def interp(self, x=None, start=None, stop=None, num=1000, step=None):
        """Interpolate data.

        Args:
            x (list or array, optional): The x-coordinates at which to
                evaluate the interpolated values. This overwrites all other arguments.
            start (number, optional): The starting value of the sequence. If None,
                the minium x value will be used.
            stop (number, optional): The end value of the sequence. If None,
                the maximum x value will be used.
            num (int, optional): Number of samples to generate.
            step (number, optional): Spacing between values. This overwrites ``num``.

        Returns:
            None
        """
        if x is None:
            if start is None:
                start = min(self.x)
            if stop is None:
                stop = max(self.x)

            if step is not None:   # step overwrites num
                x = np.arange(start, stop, step=step)
            else:
                x = np.linspace(start, stop, num=num)

        self._y = np.interp(x, self.x, self.y)
        self._x = x
        self._data = np.vstack((self._x, self._y)).transpose()

        self._restart_attr()
        if self.step is not None:
            try:
                self.check_step_x()
            except ValueError:
                pass

    def compress(self, ranges):
        """Extract data range from full data."""

        self._x, self._y  = extract(self.x, self.y, ranges)
        self._data = np.vstack((self._x, self._y)).transpose()

        self._restart_attr()
        if self.step is not None:
            try:
                self.check_step_x()
            except ValueError:
                pass

    def crop(self, start=None, stop=None):
        """crop data."""
        if start is None and stop is None:
            return
        if start is None:
            start = min(self.x)
        if stop is None:
            stop = max(self.x)
        self._x, self._y  = extract(self.x, self.y, (start, stop))
        self._data = np.vstack((self._x, self._y)).transpose()

    def floor(self, value=None, n=8, ranges=None):
        """Sets zero value for y-coordinates (shifts data verticaly).

        Args:
            value (number, optional): x-coordinate to set y-coordinate value as
                zero. If `None`, the minium value of x is used.
            n (int, optional): number of data points to average and setting to
                zero.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                ranges overwrites value and n.

        Returns:
            None
        """
        if ranges is None:
            if value is None:
                i = 0
            else:
                i = index(self.x, value)
            start = i-int(n/2)
            end   = i+int(n/2)
            if start < 0:
                end   = i+int(n/2)-start
                start = 0
            if end > len(self.y)-1:
                start = i-int(n/2) - (end - len(self.y)-1)
                end   = len(self.y)-1
            if start == end:
                end = start + 1
            self.offset = self.offset-np.mean(self.y[start:end])
        else:
            _, y= extract(self.x, self.y, ranges)
            self.offset = self.offset-np.mean(y)

    def plot(self, ax=None, offset=0, shift=0, factor=1, ranges=None, **kwargs):
        """Plot spectrum.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            vertical_increment (float or int, optional): defines the vertical offset. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            ranges (list, optional): ranges.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            None
        """
        if ax is None:
            ax = plt

        if ranges is None:
            ax.plot(self.x + shift, self.y*factor + offset, **kwargs)
        else:
            x, y = extract(self.x, self.y, ranges=ranges)
            ax.plot(x + shift, y*factor + offset, **kwargs)

    def plot_detected_peaks(self, ax=None, offset=0, shift=0, factor=1, **kwargs):
        """Plut a marker where peaks were detected.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            **kwargs: kwargs are passed to ``plt.errorbar()`` or ``plt.scatter()``.

        Returns:
            None

        Raises:
            ValueError: if :py:attr:`Spectrum.peaks` is not defined.
        """
        if self.peaks is None:
            raise ValueError('No peaks defined.')

        if ax is None:
            ax = plt

        # data
        c = np.array([self.peaks[i]['c'] for i in self.peaks])
        amp = np.array([self.peaks[i]['amp'] for i in self.peaks])
        xerr = np.array([self.peaks[i]['fwhm']/2 for i in self.peaks])

        if 'lw' not in kwargs and 'linewidth' not in kwargs:
            kwargs['lw'] = 0
        if 'elinewidth' not in kwargs :
            kwargs['elinewidth'] = 2
        if 'marker' not in kwargs :
            kwargs['marker'] = 'o'
        if 'markersize' not in kwargs and 'ms' not in kwargs:
            kwargs['markersize'] = 5

        # for i in self.peaks:
        ax.errorbar(c+shift, amp*factor+offset, xerr=xerr, **kwargs)

    def plot_fitted_ranges(self, ax=None, **kwargs):
        if self.fit_ranges is None:
            raise ValueError('data was not fitted. Use fit_peak().')


        if ax is None:
            ax = plt


        if 'lw' not in kwargs and 'linewidth' not in kwargs:
            kwargs['lw'] = 2
        if 'zorder' not in kwargs:
            kwargs['zorder'] = 10

        kwargs1 = copy.deepcopy(kwargs)
        kwargs2 = copy.deepcopy(kwargs)
        if 'color' not in kwargs1 and 'c' not in kwargs1:
            kwargs1['color'] = 'green'
        if 'color' not in kwargs2 and 'c' not in kwargs2:
            kwargs2['color'] = 'red'

        for r in self.fit_ranges:
            ax.vlines(r[0], min(self.y), max(self.y), **kwargs1)
            ax.vlines(r[1], min(self.y), max(self.y), **kwargs2)

    # OBSOLETE =================================================================
    def _plot_old(self, ax=None, normalized=False, vertical_increment=0, shift=0, factor=1, show_fit_ranges=False, show_fit=False, show_fit_partial=False, mark_peaks=False, kwargs_fit={}, **kwargs):
        """Plot spectrum.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.
            vertical_increment (float or int, optional): defines the vertical offset. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            show_fit_ranges (bool, optional): show ranges used to fit peaks.
            show_fit (bool, optional): show fitted curve.
            show_fit_partial (bool, optional): same as `show_fit`, but it shows all
                peaks used to fit.
            mark_peaks (book, optional): put a marker where peaks were detected.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes
        """

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if normalized:
            factor = factor/max(self.y)


        ax.plot((self.x + shift), self.y*factor + vertical_increment, **kwargs)
        if show_fit:
            self.plot_fit(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, **kwargs_fit)
        if show_fit_partial:
            self.plot_fit_partial(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, color='black', lw=0.5, ls='--')

        if mark_peaks:
            self.mark_peaks(ax=ax)

        if show_fit_ranges:
            for r in self.fit_ranges:
                ax.vlines(r[0], min(self.y), max(self.y), color='green', linewidth=2, zorder=10)
                ax.vlines(r[1], min(self.y), max(self.y), color='red', linewidth=2, zorder=10)

        return ax

    def _plot_fit_old(self, ax=None, normalized=False, vertical_increment=0, shift=0, factor=1, **kwargs):
        """Plot spectrum.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.
            vertical_increment (float or int, optional): defines the vertical offset. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            show_ranges (bool, optional): show ranges in which offsets were calculated.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes

        Raises:
            ValueError: if :py:attr:`Spectrum.fit_func` is not defined.
        """
        if self.fit_func is None:
            raise ValueError('No fit data to plot.')

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if normalized:
            factor = factor/max(self.y)

        # x,  y = extract(self.x, self.y, self.elastic.ranges)
        x = np.linspace(min(self.x), max(self.x), len(self.x)*10)
        ax.plot(x + shift, self.fit_func(x)*factor + vertical_increment, **kwargs)

        return ax

    def _plot_fit_partial_old(self, ax=None, normalized=False, vertical_increment=0, shift=0, factor=1, **kwargs):
        """Plot spectrum.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.
            vertical_increment (float or int, optional): defines the vertical offset. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            show_ranges (bool, optional): show ranges in which offsets were calculated.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes

        Raises:
            ValueError: if :py:attr:`Spectrum.fit_data` is not defined.
        """
        if self.fit_data is None:
            raise ValueError('No fit data to plot.')

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if normalized:
            factor = factor/max(self.y)

        # x,  y = extract(self.x, self.y, self.elastic.ranges)
        x = np.linspace(min(self.x), max(self.x), len(self.x)*10)
        for i in self.fit_data:
            ax.plot(x + shift, self.fit_data[i]['func'](x)*factor + vertical_increment, **kwargs)

        return ax

    def _mark_peaks_old(self, ax=None, show_width=True, **kwargs):
            """Plut a marker where peaks were detected.

            Args:
                ax (matplotlib.axes, optional): axes for plotting on.
                show_width (bool, optional): if True, it will use ``plt.errorbar()`` to
                    plot a horizontal line with the marker the size of the estimated
                    width. If False, it uses ``plt.scatter()`` to draw the markers.
                **kwargs: kwargs are passed to ``plt.errorbar()`` or ``plt.scatter()``.

            Returns:
                matplotlib.axes

            Raises:
                ValueError: if :py:attr:`Spectrum.peaks` is not defined.
            """
            if self.peaks is None:
                raise ValueError('No peaks defined.')

            if ax is None:
                fig = plt.figure()
                ax = fig.add_subplot(111)

            # data
            c = [self.peaks[i]['c'] for i in self.peaks]
            amp = [self.peaks[i]['amp'] for i in self.peaks]
            xerr = [self.peaks[i]['fwhm']/2 for i in self.peaks]



            if show_width:
                if 'lw' not in kwargs and 'linewidth' not in kwargs:
                    kwargs['lw'] = 0
                if 'elinewidth' not in kwargs :
                    kwargs['elinewidth'] = 2
                if 'marker' not in kwargs :
                    kwargs['marker'] = 'o'
                if 'markersize' not in kwargs and 'ms' not in kwargs:
                    kwargs['markersize'] = 5

                for i in self.peaks:
                    ax.errorbar(c, amp, xerr=xerr, **kwargs)
            else:
                if 'lw' not in kwargs and 'linewidth' not in kwargs:
                    kwargs['lw'] = 0
                if 'marker' not in kwargs :
                    kwargs['marker'] = 'o'
                if 's' not in kwargs:
                    kwargs['s'] = 100

                for i in self.peaks:
                    ax.scatter(c, amp, **kwargs)

            return ax

class Spectra(metaclass=_Meta):
    """Returns a ``spectra`` class type object for dealing with many spectrum at a time.

    Args:
        n (int, optional): when the number of spectra `n` is known, one can preallocate
            the data array.
        data (list or array, optional): list of :py:class:`spectrum` objects.
            Overwrites ``n``.

    Attributes:
        data (list of :py:attr:`Spectrum`): list with :py:attr:`Spectrum` objects.
        step (number, read only): step size for x-coordinates.
        length (number, read only): length of x-coordinates vector.
        x (number, read only): x-coordinates (only applicable if x-coordinates
            are the same for all spectra.
        shifts (list of int): Calculated shifts as int. Must always have the same lenght as data.
        shift_ref (number, read only): index of the spectrum used as reference
            for shift calculation.
        shift_ref_value (number, read only): value of the shift for the reference
            spectrum, e.g.:
                * if shift_calc_mode is 'peak'
                    shift_ref_value returns the position of the elastic peak of the reference spectrum.
                * if shift_calc_mode is 'cross-correlation'
                    shift_ref_value returns the max value of the cross-correlation of the reference spectrum with itself.
                * if shift_calc_mode is 'max'
                    shift_ref_value returns the position of the maximum value of the reference spectrum.
        shift_calc_mode (string, read only): mode used to calculate shifts.
        shift_ranges (list, read only): data range used to calculate shifts.
        disp_values (list, read only): value of the 'relative spectrum position'
            for each spectrum. See :py:attr:`shift_ref_value`
        disp_energies (number, read only): energy value associated with each spectrum.
        disp_popt (number, read only): optimal parameters of the fit of
            disp_values vs disp_energies.
        disp_func (number, read only): function of the fitted curve.
        disp (number, read only): final calculated dispersion.
        sum (number, read only): spectrum resulting from the sumation of all spectra.
    """

    _read_only     = ['step', 'length', 'x',
                      'shift_ref', 'shift_ref_value', 'shift_calc_mode', 'shift_ranges',
                      'disp_values', 'disp_shifts', 'disp_popt', 'disp', 'disp_func',
                      'sum']


    def __init__(self, n=None, data=None):

        if n is not None and data is None:
            self._data = [-1]*n
            self._restart_sum_attr()
            self._restart_check_attr()
            self._restart_shift_attr()
            self._restart_disp_attr()
        else:
            self.data = data

    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        if value is None:
            self._data = []
        else:
            if isinstance(value, Iterable) == False:
                # for s in value:
                #     if isinstance(s, Spectrum) == False:
                #         raise ValueError('all entries must be of type brixs.spectrum.')
                self._data = copy.deepcopy(value)
            else:
                raise ValueError('data must be a list.')
        self._restart_sum_attr()
        self._restart_check_attr()
        self._restart_shift_attr()
        self._restart_disp_attr()
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shifts(self):
        return self._shifts
    @shifts.setter
    def shifts(self, value):
        try:
            if len(value) == self.get_spectra_count():
                self._shifts = np.array(value)
                self._shift_calc_mode = 'user-defined'
            else:
                raise ValueError(f'Shifts must have the same length as number of spetra. Number of spectra = {self.get_spectra_count()}, number of shifts = {len(value)}.')
        except TypeError:
            if self.get_spectra_count() == 1:
                self._shifts = np.array([value])
                self._shift_calc_mode = 'user-defined'
            else:
                raise ValueError(f'Shifts must have the same length as number of spetra. Number of spectra = {self.get_spectra_count()}, number of shifts = {len(value)}.')
    @shifts.deleter
    def shifts(self):
        raise AttributeError('Cannot delete object.')

    def __next__(self):
        result = self.data[self._index]
        self._index +=1
        return result
        # End of Iteration
        raise StopIteration

    def __getitem__(self, item):
         return self.data[item]

    def __setitem__(self, item, value):
        # if isinstance(value, Spectrum) == False:
        #     raise ValueError('value must be brixs.spectrum.')
        self.data[item] = value
        self._restart_attr()

    def _restart_sum_attr(self):
        # sum attr
        self._sum = None

    def _restart_check_attr(self):
        # check x attrs
        self._length = None
        self._step = None
        self._x = None

    def _restart_shift_attr(self):
        # shifts attrs
        if self.data is None:
            self._shifts        = np.array([])
        elif self.get_spectra_count() > 0:
            self._shifts = np.array([np.NaN]*self.get_spectra_count())
        else:
            self._shifts        = np.array([])
        self._shift_ref = None
        self._shift_ref_value = None
        self._shift_calc_mode   = None
        self._shift_ranges = None

    def _restart_disp_attr(self):
        """Set relevant attributes to initial value."""
        # dispersion attrs
        self._disp_values = None
        self._disp_energies = None
        self._disp_popt  = None
        self._disp       = None
        self._disp_func  = None

    def get_spectra_count(self):
        """Returns the number of spectra."""
        return len(self.data)

    def append(self, s=None, x=None, y=None, data=None):
        """Add spectrum to the spectrum list.

        Args:
            s (spectrum obj): :py:class:`Spectrum` object to be added.
            data (list or array, optional): two column list (or array).
            x (list or array, optional): x values (1D list/array). Overwrites `data`.
            y (list or array, optional): y values (1D list/array). Overwrites `data`.

        Returns:
            None

        See Also:
            :py:func:`Spectra.remove`.
        """
        if s is None:
            self.append(s=Spectrum(x=x, y=y, data=data))
            self._restart_check_attr()
            self._restart_sum_attr()
            self._restart_disp_attr()
        elif s is not None:
            if isinstance(s, Iterable):
                # if isinstance(s, Spectrum) == False:
                #     raise ValueError('spectrum must be of type brixs.Spectrum.')
                self._data += s
                self._shifts        = np.append(self.shifts, [np.NaN]*len(s))
                # self._shifts_length = np.append(self.shifts_length, [np.NaN]*len(s))
            else:
                # for temp in s:
                #     if isinstance(temp, Spectrum) == False:
                #         raise ValueError('all entries must be of type brixs.spectrum.')
                self._data += [s]
                self._shifts        = np.append(self.shifts, np.NaN)
                # self._shifts_length = np.append(self.shifts_length, np.NaN)
            self._restart_check_attr()
            self._restart_sum_attr()
            self._restart_disp_attr()
        else:
            raise ValueError('No data to append.')

    def remove(self, idx):
        """Exclude spectrum from the spectrum list.

        Args:
            idx (int): index of the spectrum.

        Returns:
            None

        See Also:
            :py:func:`Spectra.append`.
        """
        del self.data[idx]
        self._restart_attr()

    def check_length(self):
        """Checks if all spectra has the same length.

        If all spectra have the same length, it sets :py:attr:`Spectra.length` = lenght.
        Otherwise, it raises an error.

        Returns:
            None

        Raises:
            ValueError: spectra does not have the same length.

        See Also:
            :py:func:`Spectra.check_step_x`, :py:func:`Spectra.check_same_x`.
        """
        for idx, s in enumerate(self.data):
            try:
                if len(s.data) != len(self.data[idx+1].data):
                    self._length = None
                    raise ValueError(f"Spectrum {idx} and {idx+1} have the different length.")
            except IndexError:
                pass
        self._length = len(self.data[0].x)

    def check_step_x(self, max_error=0.1):
        """Check step between data points in the x-coordinates.

            If data has a well defined step size, it sets :py:attr:`Spectra.step` = step.
            Otherwise, it raises an error.

            Args:
                max_error (number, optional): percentage value (in terms of the
                    average x step) of the maximum allowed error.

            Three checks are performed:

            1. If not yet, checks if all spectra have same length by using
            :py:func:`Spectra.check_length`.

            2. Checks if the step between two data points is the same
            through out the x vector for each spectrum (vector uniformity).

            3. Checks if the step is the same between different spectra.

            Returns:
                None

            Raises:
                ValueError: If condition 1, 2, or 3 are not satisfied.

            See Also:
                :py:func:`Spectra.check_length`, :py:func:`Spectra.check_same_x`.
        """

        # 1) check length
        if self.length is None:
            self.check_length()

        # 2) check step uniformity
        steps = np.zeros(self.get_spectra_count())
        for idx, s in enumerate(self.data):
            if s.step is None:
                try:
                    s.check_step_x()
                except ValueError:
                    raise ValueError(f"Step in the x-coordinate of spectrum {idx} seems not to be uniform.")
            steps[idx] = s.step

        # 3) check step between spectra
        avg_step = np.mean(steps)
        if sum([steps[i]-steps[i+1] > avg_step*max_error/100 for i in range(self.get_spectra_count()-1)]) > 0:
            self._step = None
            raise ValueError(f"Spectra seems to have different step size. Calculated step sizes = {steps}")
        self._step = avg_step

    def check_same_x(self, max_error=0.1):
        """Check if spectra have same x-coordinates.

        If data has same x-coordinates, it sets :py:attr:`Spectra.x` = x.
        Otherwise, it raises an error.

        Args:
            max_error (number, optional): percentage value (in terms of the
                x step) of the maximum allowed error.

        Returns:
            None

        Raises:
            ValueError: If any x-coodinate of any two spectrum is different.

        See Also:
            :py:func:`Spectra.check_length`, :py:func:`Spectra.check_step_x`.
        """
        if self.step is None:  # calculate step if it hasn't yet
            self.check_step_x(max_error=max_error)

        # check x between spectra
        for idx, s in enumerate(self.data):
            try:
                if max(abs(s.x - self.data[idx+1].x))*100/self.step > max_error:
                    self._x = None
                    raise ValueError(f"x axis of spectrum {idx} and {idx+1} seem to be different.")
            except IndexError:
                pass
        self._x = self.data[0].x

    def _gather_ys(self, ranges=None):
        """Return two lists, x and y's"""
        # check x
        if self.x is None:
            self.check_same_x()

        ys = np.zeros((self.length, self.get_spectra_count()))
        for i in range(self.get_spectra_count()):
            ys[:, i] = self.data[i].y
        if ranges is None:
            x = self.data[0].x
        else:
            x, ys = extract(self.data[0].x, ys, ranges=ranges)
        return x, ys

    def find_peaks(self, prominence=None, points=4, moving_average_window=8, ranges=None):
        for s in self.data:
            s.find_peaks(prominence=prominence, points=points, moving_average_window=moving_average_window, ranges=ranges)

    def calculate_shifts(self, ref=0, mode='cross-correlation', ranges=None, verbose=False, **kwargs):
        """Calculate how much shifted spectra are relatively to a reference spectrum.

        If :py:attr:`Spectra.step` is defined, it also sets :py:attr:`Spectra.shifts`.

        Args:
            ref (int, optional): index of reference spectrum. The shift of all
                other spectra is calculated based on the reference spectrum.
                Default is 0.
            mode (string, optional): method used to calculate the offsets.
                The current options are:

                1. 'cross-correlation' ('cc')
                2. 'max' ('m')
                3. 'elastic' ('e')

            ranges (list, optional): a pair of x values or a list of pairs. Each pair represents
                the start and stop of a data range. If None, the whole data set is used.
                If mode is peak, ranges is passed to Spectrum.fit_peak()
            verbose (bool,optional): turn verbose on/off.
            **kwargs: if ``mode = 'elastic'``, kwargs are passed to
                :py:func:`Spectra.guess_elastic_peak`.

        Returns:
            None

        Raises:
            ValueError: if ref is not an integer or a non-valid mode is selected.

        Note:
            For ``mode = 'cc'``, spectra must have the same x-coordinates (this
            is checked before execution).

        See Also:
            :py:func:`Spectra.apply_shifts`
        """
        if type(ref) != int:
            raise ValueError('ref must be an integer.')

        # try checking if x is the same (execution is much faster if it is)
        if self.x is None:
            try:
                self.check_same_x()
            except ValueError:
                pass
        # else: # if x is uniform, shifts whithin each spectrum can be saved
        #     for s in self.data:
        #         if s.step is None:
        #             try:
        #                 s.check_step_x()
        #             except ValueError:
        #                 pass

        # CROSS-CORRELATION ====================================================
        if mode in cc:
            # x must be the same for cc
            x, ys = self._gather_ys(ranges=ranges)
            # calculate cross-correlation
            for i in range(self.get_spectra_count()):
                if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                cross_correlation = np.correlate(ys[:, ref], ys[:, i], mode='full')
                self.shifts[i] = np.argmax(cross_correlation)
                if self.data[i].shift != 0: # fix shift in case there was a previous shift set
                    if self.data[i].shift_mode in roll:
                        self.shifts[i] += self.data[i].shift
                    else:
                        self.shifts[i] += round(self.data[i].shift/self.data[i])
            # save values
            self._shift_ref_value = self.shifts[ref]
            self.shifts        -= self.shift_ref_value
            self.shifts        *= -1
        # MAX ==================================================================
        elif mode == 'max':
            if self.x is not None: # if all data has the same x, execution is much faster
                x, ys = self._gather_ys(ranges=ranges)
                j_ref = np.argmax(ys[:, ref])
                self._shift_ref_value = x[j_ref]
                for i in range(self.get_spectra_count()):
                    if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                    j = np.argmax(ys[:, i])
                    self.shifts[i] = x[j] - x[j_ref]
                    if self.data[i].shift != 0: # fix shift in case there was a previous shift set
                        if self.data[i].shift_mode in roll:
                            self.shifts[i] += self.data[i].shift*self.data[i].step
                        else:
                            self.shifts[i] += self.data[i].shift
            else:  # if x is not the same, execution is slower
                if ranges is None:   # if ranges is not defined, execution is faster
                    j_ref = np.argmax(self.data[ref].y)
                    self._shift_ref_value = self.data[ref].x[j_ref]
                    for i in range(self.get_spectra_count()):
                        if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                        j = np.argmax(self.data[i].y)
                        self.shifts[i] = self.data[i].x[j] - self.data[ref].x[j_ref]
                        if self.data[i].shift != 0: # fix shift in case there was a previous shift set
                            if self.data[i].shift_mode in roll:
                                self.shifts[i] += self.data[i].shift*self.data[i].step
                            else:
                                self.shifts[i] += self.data[i].shift
                else: # if ranges is defined, execution is slower
                    x_ref, y_ref = extract(self.data[ref].x, self.data[ref].y, ranges=ranges)
                    j_ref = np.argmax(y_ref)
                    self._shift_ref_value = x_ref[j_ref]
                    for i in range(self.get_spectra_count()):
                        x, y = extract(self.data[i].x, self.data[i].y, ranges=ranges)
                        if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                        j = np.argmax(y)
                        self.shifts[i] = x[j] - x_ref[j_ref]
                        if self.data[i].shift != 0: # fix shift in case there was a previous shift set
                            if self.data[i].shift_mode in roll:
                                self.shifts[i] += self.data[i].shift*self.data[i].step
                            else:
                                self.shifts[i] += self.data[i].shift
        # PEAK =================================================================
        elif mode == 'peak':

            if 'idx' not in kwargs:
                idx = 0
            else:
                if isinstance(kwargs['idx'], Iterable):
                    raise ValueError('idx needs to be a number (not iterable).')
                else:
                    idx = kwargs['idx']

            if ranges is not None:
                kwargs['ranges'] = ranges

            self.data[ref].fit_peak(**kwargs)
            self._shift_ref_value = self.data[ref].fit_data[idx]['c']
            for i in range(self.get_spectra_count()):
                if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                self.data[i].fit_peak(**kwargs)
                self.shifts[i] = self.data[i].fit_data[idx]['c'] - self.data[ref].fit_data[idx]['c']
                if self.data[i].shift != 0:  # fix shift in case there was a previous shift set
                    if self.data[i].shift_mode in roll:
                        self.shifts[i] += self.data[i].shift*self.data[i].step
                    else:
                        self.shifts[i] += self.data[i].shift
        else:
            raise ValueError('mode not valid (cross-correlation, max, elastic).')

        # save ranges
        if ranges is None:
            self._shift_ranges = ((self.get_min_x(), self.get_max_x()), )
        elif isinstance(ranges[0], Iterable):
            self._shift_ranges = ranges
        else:
            self._shift_ranges = (ranges, )

        # save stuff
        self._shift_ref  = ref
        self._shift_calc_mode = mode

        # finish
        if verbose:
            print('done!')

    def set_shifts(self, mode=None):
        """Shift data according to calculated shifts.

        Args:
            mode (string, optional): Three modes can be selected (default is 'roll'):

                1. 'x' or 'hard': y coordinates are fully preserved, while
                    x-coordinates are shifted.
                2. 'y', 'interp', or 'soft': x-coordinates are preserved, while
                    y-coordinates are interpolated with a shift.
                3. 'roll', 'rotate', or 'r': both x- and y-coordinates are preserved
                    as y-coordinates are simply rolled along the array. In this case,
                    shift values must be an integer and x-coordinates must be
                    evenly spaced.

        Returns:
            None

        Raises:
            ValueError: if a non-valid mode is selected, if shifts have not been
                calculated, or if x-coordinates are not evenly space for mode =
                'roll'.

        Warning:
            It is always better to use mode='hard' or 'roll' since y-coordinates are fully
            preserved (no interpolation). After applying a shift using the 'interp',
            one can apply a 'inverse' shift to retrieve the original data.
            The diference between the retrieved y-coordinates and the original
            data will give an ideia of the information loss caused by the interpolation.
            If statiscs are good enough, interpolating the data should not be a problem.

        See Also:
            :py:func:`Spectra.calculate_shifts`
        """
        if  np.all(self.shifts == 0):
            return

        if mode is None:
            if self.shift_calc_mode in cc:
                mode = 'roll'
            elif self.shift_calc_mode == 'max':
                mode = 'x'
            elif self.shift_calc_mode == 'peak':
                mode = 'x'
            elif self.shift_calc_mode == 'user-defined':
                text = 'shift mode not define (please, select "roll", "soft", "hard").' +\
                        'Note that shifts were calculated elsewhere by the user, ' +\
                        'the script can has no way of knowing the most appropriate shift mode.'
                raise ValueError(text)

        if self.shift_calc_mode in cc and mode not in roll:
            raise ValueError('shifts were calculated using cross-correlation. Only shift mode possible is "roll".')

        for i in range(self.get_spectra_count()):
            self.data[i].set_shift(value=-self.shifts[i], mode=mode)

        self._restart_shift_attr()
        self._restart_check_attr()
        self._restart_sum_attr()

    def calculate_sum(self):
        """Returns Spestrum object with the sum of all spectra.

        It also save it at :py:attr:`Spectra.sum`.

        Returns:
            :py:class:`Spectra` object.

        Note:
            All spectra have to have the same x-coordinates. This is verified
            before summing up the spectra.
        """

        if self.x is None:
            self.check_same_x()

        y = np.zeros(len(self.x))
        for i in range(self.get_spectra_count()):
            y += self.data[i].y

        self._sum = Spectrum(x=self.x, y=y)
        return self._sum

    def save(self, folderpath='./', prefix='spectrum_', suffix='.dat', zfill=None, delimiter=', ', header='', fmt='%.8f', verbose=False):
        r"""Saves all spectra in a folder.

        Args:
            folderpath (string or pathlib.Path, optional): path to directory.
                Default is the current directory.
            prefix (string, optional): prefix used for naming the files.
            suffix (string, optional): suffix used for naming the files.
            zfill (int, optional): number of digits for file numbering. If `None`,
                zfill will be determined.
            delimiter (str, optional): The string used to separate values.
                If whitespaces are used, consecutive whitespaces act as delimiter.
                Use ``\\t`` for tab. The default is comma (', ').
            header (string, optional): text to add at the beginning of each file.
                Use ``\n`` for new line. Comment flag ('# ') is added automatically.
            fmt (string, or list, optional): format for saving data.
                If string, the value is used for x- and y-coordinates. If tuple
                of strings, the first string is used for x-coordinates and the
                second for y-coordinates.

                    fmt = (%[flag]width[.precision]specifier)

                * `flag` can be:
                    1. '-' for left justify
                    2. '+', puts + or - in front of the numbers
                    3. '0' to Left pad the number with zeros instead of space (see width).

                * `width` is the minimum number of characters to be printed.

                * `precision` is the number of significant digits.

                * `specifier` is the type of notation. Tipically, either 'e' for
                scientific notation of 'f' for decimal floating point.

                * a common `fmt` strings is: '%.3f' for 3 decimal places.

                *  for more information see `np.savetxt <https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html?highlight=savetxt#numpy.savetxt>`_ documentation::

            verbose (bool, optional): turn verbose on and off. Default is `False`.

        Returns:
            None

        """
        folderpath = Path(folderpath)

        if zfill is None:
            zfill = n_digits(self.get_spectra_count()-1)[0]

        for idx, s in enumerate(self.data):
            filename = f'{prefix}' + f'{idx}'.zfill(zfill) + f'{suffix}'
            if verbose:  print(f'({i}/{self.get_spectra_count()-1}) saving {filename}')
            s.save(filepath=folderpath/filename, delimiter=delimiter, header=header, fmt='%.4f')
        if verbose: print('Done!')

    def calib(self, value):
        """Calibrate data (from x-coordinates to energy).

        Args:
            value (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].

        Returns:
            None
        """
        for s in self.data:
            s.set_calib(value=value)

        self._restart_sum_attr()
        self._restart_check_attr()
        self._restart_shift_attr()
        self._restart_disp_attr()

    def interp(self, x=None, start=None, stop=None, num=None, step=None):
        """Interpolate data.

        Args:
            x (list or array, optional): The x-coordinates at which to
                evaluate the interpolated values. This overwrites all other arguments.
            start (number, optional): The starting value of the sequence. If `None`,
                the minium x value will be used.
            stop (number, optional): The end value of the sequence. If `None`,
                the maximum x value will be used.
            num (int, optional): Number of samples to generate.
            step (number, optional): Spacing between values. This overwrites ``num``.

        Returns:
            None
        """
        if x is None:
            if self.x is None:
                if start is None:
                    start = max([min(s.x) for s in self.data])
                if stop is None:
                    stop = min([max(s.x) for s in self.data])
                if num is None:
                    num = max([len(s.x) for s in self.data])
            else:
                if start is not None or stop is not None or num is not None or step is not None:
                    if start is None:
                        start = max(self.x)
                    if stop is None:
                        stop = min(self.x)
                    if num is None:
                        num = len(self.x)
                else:
                    return

        for s in self.data:
            s.interp(x=x, start=start, stop=stop, num=num, step=step)

        self._restart_sum_attr()
        self.check_same_x()

    def compress(self, ranges):
        """Slices data.

        Args:
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.

        Returns:
            None
        """
        if isinstance(ranges, Iterable):
            if isinstance(ranges[0], Iterable):
                for s in self.data:
                    s.compress(ranges=ranges)
                self._restart_sum_attr()
                self._restart_check_attr()
            else:
                self.crop(start=ranges[0], stop=ranges[0])

    def get_min_x(self):
        if self.x is None:
            return min([min(s.x) for s in self.data])
        else:
            return min(self.x)

    def get_max_x(self):
        if self.x is None:
            return max([max(s.x) for s in self.data])
        else:
            return max(self.x)

    def crop(self, start=None, stop=None):
        if self.x is None:
            print('here')
            if start is None:
                start = min(self.data[0].x)
                for s in self.data:
                    temp = min(s.x)
                    if temp > start:
                        start = temp
            if stop is None:
                stop = max(self.data[0].x)
                for s in self.data:
                    temp = max(s.x)
                    if temp < stop:
                        stop = temp
            for s in self.data:
                s.crop(start, stop)
            self._restart_sum_attr()
            self._restart_check_attr()
        elif start is not None or stop is not None:
            for s in self.data:
                s.crop(start, stop)
            self._restart_sum_attr()
            self._restart_check_attr()

    def floor(self, value=None, n=8, ranges=None):
        """Sets zero value for y-coordinates (shifts data verticaly).

        Args:
            value (number, optional): x-coordinate to set y-coordinate value as
                zero. If `None`, the minium value of x is used.
            n (int, optional): number of data points to average and setting to
                zero.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                ranges overwrites value and n.

        Returns:
            None
        """
        for s in self.data:
            s.floor(value=value, n=n, ranges=ranges)

    def plot(self, ax=None, vertical_increment=0, shift=0, factor=1, ranges=None, **kwargs):
        """Plot spectra.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            vertical_increment (int, optional): if one spectrum is plotted, it
                adds a vertical offset to the plotted curve. If many spectra are ploted,
                ``vertical_increment`` defines
                the vertical offset between each curve. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            None
        """
        if ax is None:
            ax = plt


        for i in range(self.get_spectra_count()):
            self.data[i].plot(ax=ax, offset=vertical_increment*-i, shift=shift, factor=factor, ranges=ranges, label=i, **kwargs)

        plt.legend()

    def plot_shift_ranges(self, ax=None, **kwargs):

        if self.shift_ranges is None:
            raise ValueError('Shift range not defined. Use calculate_shifts().')

        if ax is None:
            ax = plt

        if 'lw' not in kwargs and 'linewidth' not in kwargs:
            kwargs['lw'] = 2
        if 'zorder' not in kwargs:
            kwargs['zorder'] = 10
        kwargs1 = copy.deepcopy(kwargs)
        kwargs2 = copy.deepcopy(kwargs)
        if 'color' not in kwargs1 and 'c' not in kwargs1:
            kwargs1['color'] = 'green'
        if 'color' not in kwargs2 and 'c' not in kwargs2:
            kwargs2['color'] = 'red'

        for r in self.shift_ranges:
            ax.axvline(r[0], **kwargs1)
            ax.axvline(r[1], **kwargs2)

    def plot_detected_peaks(self, ax=None, vertical_increment=0, shift=0, factor=1, **kwargs):

        if ax is None:
            ax = plt

        for i in range(len(idx)):
            if self.data[i].peaks is None:
                raise ValueError('find_peaks needs to be ran in Spectrum {i}.')

        for i in range(len(idx)):
            self.data[i].plot_detected_peaks(ax=ax, offset=vertical_increment*-i, shift=shift, factor=factor,ranges=ranges, label=f'peaks for spectrun {i}', **kwargs)
        plt.legend()

    def calculate_calib(self, start=None, stop=None, values=None, **kwargs):
        """Calculate calibration factor (dispersion).

        It calculates the shift between spectra and infer the calibration factor.

        Args:
            start (number): value used for measuring the first spectrum.
            stop (number): value used for measuring the last spectrum.
            values (list): value list. It overwrites start and
                stop. Must be the same lenght as data.
            **kwargs: parameters to be passed to :py:func:`Spectra.calculate_shifts()`

        Returns:
            calibration factor (number)
        """
        # check number of values matches the numbre of spectra
        if values is None:
            values = np.linspace(start, stop, self.get_spectra_count())
        if self.get_spectra_count() != len(values):
            raise ValueError(f'number of values ({len(values)}) do not match the number of spectra ({self.get_spectra_count()})')

        # self.calculate_shifts() parameters
        if 'asymmetry' not in kwargs:
            kwargs['asymmetry'] = False
        if 'ref' not in kwargs:
            kwargs['ref'] = 0
        if 'mode' not in kwargs:
            kwargs['mode'] = 'cc'

        # calculate
        if self.shift_calc_mode != kwargs['mode']:
            self.calculate_shifts(**kwargs)
        if self.shift_calc_mode in cc:
            self._disp_shifts = (self.shifts + self.shift_ref_value)*self.step
        else:
            self._disp_shifts = self.shifts + self.shift_ref_value
        self._disp_values = values

        # fit elastic positions
        popt = np.polyfit(self.disp_values, self.disp_shifts, deg=1)
        self._disp_popt = popt
        self._disp   = 1/popt[0]
        self._disp_func   = np.poly1d(popt)

        return self.disp

    def plot_disp(self, ax=None, show_fit=True, **kwargs):
        """Plot position of elastic peak (or data shift) as function of energy.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            show_fit (bool, optional): show fit.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes
        """
        if ax is None:
            ax = plt

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'linewidth' not in kwargs and 'lw' not in kwargs:
            kwargs['linewidth'] = 0


        ax.plot(self.disp_values, self.disp_shifts, **kwargs)
        if show_fit:
            x = np.linspace(self.disp_values[0], self.disp_values[-1], self.get_spectra_count()*10)
            y = self.disp_func(x)
            ax.plot(x, y, color='black')

    def plot_map(self, start=None, stop=None, values=None, ranges=None, vlines=False, **kwargs):
        """Plot position of elastic peak (or data shift) as function of energy.

        x-coordinates must be the same for all spectra.

        Args:
            start (number): value (energy or momentum) used for measuring the first spectrum.
            stop (number): value (energy or momentum) used for measuring the last spectrum.
                Energies are them extrapolated from start to stop in
                steps of one.
            values (list): values list. It overwrites start and
                stop. Must be the same lenght as data.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                ranges overwrites value and n.
            vlines (bool, optional): show separation lines between spectra.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes
        """

        if values is None:
            values = np.arange(start, stop+1)
        if self.get_spectra_count() != len(values):
            raise ValueError(f'number of values ({len(values)}) do not match the number of spectra ({self.get_spectra_count()})')
        if self.x is None:
            self.check_same_x()

        x, ys = self._gather_ys(ranges=ranges)

        if 'vmin' not in kwargs:
            kwargs['vmin'] = min([min(k) for k in ys])
        if 'vmax' not in kwargs:
            kwargs['vmax'] = max([max(k) for k in ys])
        if 'shading' not in kwargs:
            kwargs['shading'] = 'nearest'

        fig = plt.figure()
        ax = fig.add_subplot(111)

        ax.pcolormesh(values, self.x, ys, **kwargs)

        if vlines:
            d = (np.diff(values)/2)[0]
            ax.vlines(values-d, min(self.x), max(self.x),)
            for i, x in enumerate(values):
                ax.text(x-d, max(self.x)-(max(self.x)-min(self.x))*0.1, i,)


        return ax



# %%
