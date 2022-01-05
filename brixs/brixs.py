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

"""

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
from collections.abc import Iterable
from scipy.signal import find_peaks

# backpack
from .backpack.filemanip import save_data, save_obj, load_obj
from .backpack.arraymanip import index, moving_average, extract, shifted, peak_fit
from .backpack.figmanip import n_digits




class Meta(type):
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
            return getter, setter, deleter, 'read only variable'

        def lazy_non_removable(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    return setattr(self, variable, value)
                def deleter(self):
                    raise AttributeError('Attribute cannot be deleted.')
            return getter, setter, deleter, 'non removable variable'

        new_attrs = {}
        for name, value in attrs.items():
            if name == '_read_only':
                for attr in value:
                    _property = property(*lazy_read_only(attr))
                    new_attrs[attr] = _property
            elif name == '_non_removable':
                for attr in value:
                    _property = property(*lazy_non_removable(attr))
                    new_attrs[attr] = _property
            else:
                new_attrs[name] = value

        return type(class_name, bases, new_attrs)


class MetaIterator(type):
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
            return getter, setter, deleter, 'read only variable'

        def lazy_non_removable(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    return setattr(self, variable, value)
                def deleter(self):
                    raise AttributeError('Attribute cannot be deleted.')
            return getter, setter, deleter, 'non removable variable'

        new_attrs = {}
        for name, value in attrs.items():
            if name == '_read_only':
                for attr in value:
                    _property = property(*lazy_read_only(attr))
                    new_attrs[attr] = _property
            elif name == '_non_removable':
                for attr in value:
                    _property = property(*lazy_non_removable(attr))
                    new_attrs[attr] = _property
            else:
                new_attrs[name] = value

        return type(class_name, bases, new_attrs)

    def __next__(self):
        result = self.data[self._index]
        self._index +=1
        return result
        # End of Iteration
        raise StopIteration


class PhotonEvents(metaclass=Meta):
    """Creates a ``PhotonEvents`` class type object to deal with photon events lists.

    Args:
        events (list or array): two (x, y) or three (x, y, intensity) list or array with
            photon events.
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


    def save(self, filepath, delimiter=', ', data_format='%.4f'):
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
        save_data(self.events, filepath=Path(filepath), delimiter=delimiter, header=header, data_format=data_format)


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

        args:
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
        args:
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
                    fwhm[y_bin][x_bin] = s.elastic_w
                    if s.elastic_w < best:
                        best_bins = (x_bin, y_bin)
                        best = s.elastic_w
                except RuntimeError:
                    pass


        print(f'Done. Best bins is {best_bins}.')

        return best_bins, best, fwhm


    def calculate_offsets(self, mode='cross-correlation', ranges=None, ref=0):
        """Calculate the offset of each column relative to a reference column.

        args:
            ref (int, optional): reference column. The offset of all other columns
                is calculated based on the reference column. Default is 0.
             mode (string, optional): method used to calculate the offsets.
                The current options are: 'cross-correlation', 'max', 'elastic',
                and 'elastic-cross-correlation'.
            ranges (list, optional): a pair of x values or a list of pairs. Each pair represents
                the start and stop of a data range. If None, the whole data set is used.

        returns:
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

        args:
            deg (int, optional): degree for the polnomial fit. The default is 1.
            f (function, optional):  a function y = f(x, a, b, ...) that returns the
                value of y as a function of x and other parameters to be fitted.
                This overwrites the standard polynomal fit and ``deg``
                is ignored.
            **kwargs: kwargs are passed to ``scipy.optimize.curve_fit()`` that
                that optimizes the parameters of f. If f is not None, kwargs are
                passed to ``numpy.polyfit()``.

        returns:
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

        args:
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

        args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        returns:
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
            s = self.calculate_spectrum(bins=y_bin, mode='lenght')
            try:
                s.guess_elastic_peak(**kwargs)
                fwhm[y_bin] = s.elastic_w
                if s.elastic_w < best:
                    best_bins = (1, y_bin)
                    best = s.elastic_w
            except RuntimeError:
                fwhm[y_bin] = -1


        return best_bins, best, fwhm


    def calculate_spectrum(self, bins=None, bins_size=None, mode='bins'):
        """Sum the photon events in the x direction.

        args:
            y_bins (int, optional): number of y bins. If None, the current binning is used.
            bins_size (int or tuple, optional): size of the y bins. This overwrites
                the argument ``y_bins``. If None, the current binning is used.
            mode (string, optional): Type of the x axis. If 'bins', the x axis
                is numbered from 0 to the number of y_bins. If 'lenght', the
                detector x values are used.

        mode = lenght is usefull when comparing data with different binnings.

        returns:
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
        elif mode == 'lenght':
            x = temp.y_centers
        else:
            raise ValueError('mode can only be `bins` or `lenght`.')
        self._spectrum = Spectrum(x=x, y=temp.hist[0])
        self._spectrum_bins = temp.bins[:]

        return self._spectrum


    def plot(self, ax=None, pointsize=1, cutoff=0, show_bins=(False, False), show_offsets=False, show_fit=False, offsets_kwargs={}, bins_kwargs={}, fit_kwargs={}, **kwargs):
        """Plot photon events.

        args:
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


        returns:
            matplotlib.axes
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            # ax.set_facecolor('black')

        # kwargs
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'ms' not in kwargs:
            kwargs['ms'] = pointsize
        if 'mfc' not in kwargs:
            kwargs['mfc'] = 'black'
            # kwargs['mfc'] = 'white'
        if 'markeredgewidth' not in kwargs:
            kwargs['markeredgewidth'] = 0
        if 'linewidth' not in kwargs:
            kwargs['linewidth'] = 0

        # bins_kwargs
        if 'linewidth' not in bins_kwargs:
            bins_kwargs['linewidth'] = 0.5
        if 'color' not in bins_kwargs:
            bins_kwargs['color'] = 'cyan'
        if 'zorder' not in bins_kwargs:
            bins_kwargs['zorder'] = 3

        # offsets_kwargs
        if 'color' not in offsets_kwargs:
            offsets_kwargs['color'] = 'green'
        if 'zorder' not in offsets_kwargs:
            offsets_kwargs['zorder'] = 4

        # fit_kwargs
        if 'color' not in fit_kwargs:
            fit_kwargs['color'] = 'red'
        if 'linewidth' not in fit_kwargs:
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

        args:
            ax (matplotlib.axes, optional): axes for plotting on.
            shift (int, optional): vertical shift. Default is 0.
            **kwargs: kwargs are passed to ``plt.scatter()`` that plots the data.

        returns:
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

        args:
            ax (matplotlib.axes, optional): axes for plotting on.
            shift (int, optional): vertical shift. Default is 0.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        returns:
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


class Spectrum(metaclass=Meta):
    """Creates a ``spectrum`` class type object to deal with (x, y) data types.

    Args:
        data (list or array, optional): three column list (or array) with photon
            events. Column order should be x, y, and intensity.
        x (list or array, optional): x values (1D list/array). Overwrites `data`.
        y (list or array, optional): y values (1D list/array). Overwrites `data`.

    """

    _non_removable = ['elastic_amp', 'elastic_c', 'elastic_w', 'elastic_func', 'elastic_ranges']

    def __init__(self, data=None, x=None, y=None):
        if x is None and y is None:
            self.data = data
        else:
            if x is None:
                self._x = np.arange(0, len(y))
            else:
                if len(x) == len(y):
                    self._x = x
                else:
                    raise ValueError('x and y data are not the same lenght.')
            self._y = y
            self._data = np.vstack((x, y)).transpose()

        self.elastic_amp    = None
        self.elastic_c      = None
        self.elastic_w      = None
        self.elastic_func   = None
        self.elastic_ranges = None

    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        try:
            if value.shape[1] > 1:
                raise ValueError('Data must have two columns (x, y).')
            elif value is None:
                raise ValueError('No data to load.')
            elif value.shape[1] == 1:   # one column data (list)
                self._x = np.arange(0, len(value))
                self._y = value
                self._data = np.vstack((self.x, self.y)).transpose()
                self.elastic_amp    = None
                self.elastic_c      = None
                self.elastic_w      = None
                self.elastic_func   = None
                self.elastic_ranges = None
            else:
                self._data = value
                self._x = self.data[:, 0]
                self._y = self.data[:, 1]
                self.elastic_amp    = None
                self.elastic_c      = None
                self.elastic_w      = None
                self.elastic_func   = None
                self.elastic_ranges = None
        except IndexError:  # one column data (list)
            self._x = np.arange(0, len(value))
            self._y = value
            self._data = np.vstack((self.x, self.y)).transpose()
            self.elastic_amp    = None
            self.elastic_c      = None
            self.elastic_w      = None
            self.elastic_func   = None
            self.elastic_ranges = None
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, value):
        if len(value) != len(self.x):
            raise ValueError('Lenght of x is not compatible with data.')
        else:
            self._x = value
            self._data[:, 0] = value
            self.elastic_amp    = None
            self.elastic_c      = None
            self.elastic_w      = None
            self.elastic_func   = None
            self.elastic_ranges = None
    @x.deleter
    def x(self):
        raise AttributeError('Cannot delete object.')

    @property
    def y(self):
        return self._y
    @y.setter
    def y(self, value):
        if len(value) != len(self.y):
            raise ValueError('Lenght of y is not compatible with data.')
        else:
            self._y = value
            self._data[:, 1] = value
            self.elastic_amp    = None
            self.elastic_c      = None
            self.elastic_w      = None
            self.elastic_func   = None
            self.elastic_ranges = None
    @y.deleter
    def y(self):
        raise AttributeError('Cannot delete object.')

    def set_x(self, x):
        self.x = x

    def get_x(self):
        return self.x

    def set_y(self, y):
        self.y = y

    def get_y(self):
        return self.y

    def save(self, filepath, delimiter=', ', header='', data_format='%.4f'):
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
        save_data(self.data, filepath=Path(filepath), delimiter=delimiter, header=header, data_format='%.4f')

    def guess_elastic_peak(self, prominence=None, points=5, moving_average_winddow=8, ranges=None, asymmetry=True, side='right', mode='simple'):
        """

        if prominence is none, prominence is 5% of the max.

        1) data is smoothed until max of the smoothed data is 20 smaller than max
        of data.

        2) use find_peaks

            3) if only one peak is found, fit peak

            4) if more than one peak is found, get the two peak to the far right
            (or left) sides and try to find a minimum between than.

                5) if a minimum can be found. Fit peak until this minimum.

                6) if it can't be found, fit until center of peak + width/2

        """
        right = ['right', 'r', 'RIGHT', 'R', 'max']
        left = ['left', 'l', 'LEFT', 'L', 'min']

        if ranges is None:
            ranges = [[min(self.x), max(self.x)]]
            x = self.x
            y = self.y
        else:
            x, y = extract(self.x, self.y, ranges=ranges)

        # data smoothing
        y2 = copy.deepcopy(y)
        # n = 8
        # while max(y)*0.8 < max(y2) and n < 20:
        #     n *= 1.2
        #     n = round(n)
        y2 = moving_average(y, moving_average_winddow)
        # print(n)
        x2 = moving_average(x, moving_average_winddow)

        if prominence is None:
            prominence = max(y2)*0.05

        # find peaks
        try:
            if side in right:
                i = -1
            elif side in left:
                i = 0
            else:
                raise ValueError('side not valid. Choose "left" or "right".')

            peaks, d = find_peaks(y2, prominence=prominence, width=points)
            peaks = [index(x, x2[p]) for p in peaks]

            # initial parameters for fitting
            guess_c = x[peaks[i]]
            guess_w = d['widths'][i]*(x[1]-x[0])
            guess_A = d['prominences'][i]
        except IndexError:
            raise RuntimeError('Cannot find peaks. Try changing ``points`` and ``prominence``.')

        # find peak next to elastic line
        if len(peaks) > 1:
            if side in right:
                if peaks[-2]+d['widths'][-2]/2 > peaks[-1]-d['widths'][-1]/2:
                    ranges = ((peaks[-1]-d['widths'][-1]/2, x[-1]), )
                else:
                    ranges_temp = (x[peaks[-2]], x[peaks[-1]])
                    x_temp, y_temp = extract(x, y, ranges_temp)
                    ranges = ((x_temp[np.argmin(y_temp)], x[-1]), )
            else:
                if peaks[1]+d['widths'][1]/2 > peaks[2]-d['widths'][2]/2:
                    ranges = ((x[0], peaks[1]-d['widths'][1]/2), )
                else:
                    ranges_temp = (x[peaks[1]], x[peaks[2]])
                    x_temp, y_temp = extract(x, y, ranges_temp)
                    ranges = ((x[0], x_temp[np.argmin(y_temp)]), )

            x2fit, y2fit = extract(x, y, ranges)
            _, popt, _, f = peak_fit(x2fit, y2fit, guess_c=guess_c, guess_A=guess_A, guess_w=guess_w, guess_offset=0, fixed_m=False, asymmetry=asymmetry)

        else:
            _, popt, _, f = peak_fit(x, y, guess_c=guess_c, guess_A=guess_A, guess_w=guess_w, guess_offset=0, fixed_m=False, asymmetry=asymmetry)

        self.elastic_amp    = popt[0]
        self.elastic_c      = popt[1]
        self.elastic_w      = popt[2]
        self.elastic_ranges = ranges
        self.elastic_func   = f

    def apply_shift(self, value, mode='roll'):
        """Shift data.

        Args:
            value (float or int): shift value.
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
        self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
        self._data[:, 0] = self._x
        self._data[:, 1] = self._y
        if self.elastic_c is not None:
            self.elastic_c      += value
        if self.elastic_func is not None:
            elastic = copy.deepcopy(self.elastic_func)
            def elastic2(x, value):
                return elastic(x-value)
            self.elastic_func = lambda x: elastic2(x, value)
        if self.elastic_ranges is not None:
            self.elastic_ranges = [(r[0]+value, r[1]+value) for r in self.elastic_ranges]

    def calib(self, value):
        """Calibrate data.

        args:
            value (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].

        returns:
            None
        """
        f = lambda x, y: (x*value, y)
        self.apply_correction(f=f)

    def apply_correction(self, f):
        """Changes the values of x, y based on a function.

        Example:
            f = lambda x, y: (x, y**2)

        args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        returns:
            None
        """
        self._x, self._y = f(self.x, self.y)
        self._data[:, 0] = self._x
        self._data[:, 1] = self._y
        if self.elastic_c is not None:
            self.elastic_c, self.elastic_amp = f(self.elastic_c, self.elastic_amp)
        if self.elastic_func is not None:
            elastic = copy.deepcopy(self.elastic_func)
            self.elastic_func = lambda x: f(x, elastic(x))
        if self.elastic_ranges is not None:
            def temp2(x):
                x2, _ = f(x, np.interp(x, self.x, self.y))
                return x2
            self.elastic_ranges = [(temp2(r[0]), temp2(r[1])) for r in self.elastic_ranges]

    def plot(self, ax=None, normalized=False, vertical_increment=0, shift=0, factor=1, show_elastic_range=False, show_elastic=False, kwargs_elastic={}, **kwargs):
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
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if normalized:
            factor = factor/max(self.y)


        ax.plot((self.x + shift), self.y*factor + vertical_increment, **kwargs)
        if show_elastic:
            self.plot_elastic(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, **kwargs_elastic)
        if show_elastic_range:
            for r in self.elastic_ranges:
                ax.vlines(r[0], min(self.y), max(self.y), color='green', linewidth=2, zorder=10)
                ax.vlines(r[1], min(self.y), max(self.y), color='red', linewidth=2, zorder=10)

        return ax

    def plot_elastic(self, ax=None, normalized=False, vertical_increment=0, shift=0, factor=1, **kwargs):
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
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if normalized:
            factor = factor/max(self.y)

        # x,  y = extract(self.x, self.y, self.elastic.ranges)
        x = np.linspace(min(self.x), max(self.x), len(self.x)*10)
        ax.plot(x + shift, self.elastic_func(x)*factor + vertical_increment, **kwargs)

        return ax

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
                start = min(self.x)
            if stop is None:
                stop = max(self.x)

            if step is not None:   # step overwrites num
                x = np.arange(start, stop, step=step)
            else:
                x = np.linspace(start, stop, num=num)

        self._y = np.interp(x, self.x, self.y)
        self._x = x

    def compress(self, ranges):
        self._x, self._y  = extract(self.x, self.y, ranges)

    def floor(self, value=None, n=8):
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
        f = lambda x, y: (x, y - np.mean(y[start:end]))
        self.apply_correction(f)


class Spectra(metaclass=MetaIterator):
    """Creates a ``spectra`` class type object to deal with many spectrum at a time.

    Args:
        data (list or array, optional): list of :py:class:`spectrum` objects.
    """

    _non_removable = ['step',]
    _read_only     = ['shift_mode', 'shift_ranges', 'x', 'sum']

    def __init__(self, data=None):
        # basic attr
        if data is None:
            self._data = []
            self._shifts         = np.array([])
            self._shifts_lenght  = np.array([])
        else:
            self.data = data
        self._step = None

        # shift attr
        self._shift_mode   = None
        self._shift_ranges = None

        self._x = None

        # sum attr
        self._sum = None


    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        if value is None:
            self._data = []
            self._shifts        = np.array([])
            self._shifts_lenght  = np.array([])
            self._step = None
            self._x = None
            self._shift_mode = None
            self._shift_ranges = None
            self._sum = None
        else:
            try:
                # for s in value:
                #     if isinstance(s, Spectrum) == False:
                #         raise ValueError('all entries must be of type brixs.spectrum.')
                self._data = copy.deepcopy(value)
                self._shifts = np.array([np.NaN]*len(value))
                self._shifts_lenght = np.array([np.NaN]*len(value))
                self._step = None
                self._x = None
                self._shift_mode = None
                self._shift_ranges = None

            except TypeError:
                raise ValueError('data must be a list.')
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
            else:
                raise ValueError(f'Shifts must have the same lenght as number of spetra. Number of spectra = {self.get_spectra_count()}, number of shifts = {len(value)}.')
        except TypeError:
            if self.get_spectra_count() == 1:
                self._shifts = np.array([value])
            else:
                raise ValueError(f'Shifts must have the same lenght as number of spetra. Number of spectra = {self.get_spectra_count()}, number of shifts = {len(value)}.')
    @shifts.deleter
    def shifts(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shifts_lenght(self):
        return self._shifts_lenght
    @shifts_lenght.setter
    def shifts_lenght(self, value):
        try:
            if len(value) == self.get_spectra_count():
                self._shifts_lenght = np.array(value)
            else:
                raise ValueError(f'Shifts must have the same lenght as number of spetra. Number of spectra = {self.get_spectra_count()}, number of shifts = {len(value)}.')
        except TypeError:
            if self.get_spectra_count() == 1:
                self._shifts_lenght = np.array([value])
            else:
                raise ValueError(f'Shifts must have the same lenght as number of spetra. Number of spectra = {self.get_spectra_count()}, number of shifts = {len(value)}.')
    @shifts_lenght.deleter
    def shifts_lenght(self):
        raise AttributeError('Cannot delete object.')

    def __getitem__(self, item):
         return self.data[item]

    def get_spectra_count(self):
        """Returns the number of spectra."""
        return len(self.data)

    # def get_filelist(self):
    #     """Returns a filepath list of all spectra."""
    #     return [x.filepath for x in self.spectrum]
    #
    # def get_specrum_by_filename(self, filename):
        """Returns a idx list of spectra associated with filename.
        """
        filelist = self.get_filelist()
        filelist = [Path(x) if x is not None else None for x in filelist]
        return [idx for idx, s in enumerate(filelist) if filename == filelist.name]

    def append(self, s):
        """Add spectrum to the spectrum list.

        args:
            s (spectrum obj): spectrum object to be added.

        returns:
            None

        See Also:
            :py:func:`spectra.exclude`.
        """
        if s is not None:
            if isinstance(s, Iterable):
                # if isinstance(s, Spectrum) == False:
                #     raise ValueError('spectrum must be of type brixs.Spectrum.')
                self._data += s
                self._shifts        = np.append(self.shifts, [np.NaN]*len(s))
                self._shifts_lenght = np.append(self.shifts_lenght, [np.NaN]*len(s))
            else:
                # for temp in s:
                #     if isinstance(temp, Spectrum) == False:
                #         raise ValueError('all entries must be of type brixs.spectrum.')
                self._data += [s]
                self._shifts        = np.append(self.shifts, np.NaN)
                self._shifts_lenght = np.append(self.shifts_lenght, np.NaN)
            self._step = None
            self._x    = None
            self._sum  = None
        else:
            raise ValueError('No data to append.')

    def remove(self, idx):
        """Exclude spectrum from the spectrum list.

        args:
            idx (int): index of the spectrum.

        returns:
            None

        See Also:
            :py:func:`Spectra.append`.
        """
        del self.data[idx]
        np.delete(self.shifts, idx)
        np.delete(self.shifts_lenght, idx)
        self._sum = None

    def check_step_x(self, max_error=0.001):
        """Check the step between data point in the x-coordinates

            args:
                max_error (number, optional): percentage (of the x step) value of the max error.

            Three checks are performed:

                1) Checks if all spectra have same lenght.

                2) Checks if the x step between two data points is the same
                    through out all x-axis.

                3) checks if the step between two data points is the same between
                    different spectra.

            raises:
                ValueError: If condition 1, 2, and 3 are not satisfied.
        """

        # check length
        for idx, s in enumerate(self.data):
            try:
                if len(s.data) != len(self.data[idx+1].data):
                    raise ValueError(f"Spectrum {idx} and {idx+1} have the different length.")
            except IndexError:
                pass

        # check step uniformity
        for idx, s in enumerate(self.data):
            d = np.diff(s.x)
            if (max(d) - min(d))*100/np.mean(np.diff(s.x)) > max_error:
                raise ValueError(f"Step in the x-coordinate of spectrum {idx} seems not to be uniform.")

        # check step between spectra
        for idx, s in enumerate(self.data):
            try:
                avg_step = (s.x[1] - s.x[0])
                avg_step_2 = (self.data[idx+1].x[1] - self.data[idx+1].x[0])
                if  (avg_step - avg_step_2)*100/avg_step > max_error:
                    raise ValueError(f"Spectrum {idx} ({avg_step}) and {idx+1} ({avg_step_2}) seems to have different step size.")
            except IndexError:
                pass

        self._step = avg_step
        if np.isnan(sum(self.shifts)) and np.isnan(sum(self.shifts_lenght)) == False:
            self.shifts = np.round(self.shifts_lenght/self.step)

    def calculate_shifts(self, ref=0, mode='cross-correlation', ranges=None, verbose=False, **kwargs):
        """Calculate the shift of each spectrum relative to a reference spectrum.

        Data needs evenly spaced data points. The scripts checks that.

        args:
            ref (int, optional): index of reference spectrum. The shift of all other spectra
                is calculated based on the reference spectrum. Default is 0.
             mode (string, optional): method used to calculate the offsets.
                The current options are:

                    1. 'cross-correlation'
                    2. 'max'
                    3. 'elastic_peak'

            ranges (list, optional): a pair of x values or a list of pairs. Each pair represents
                the start and stop of a data range. If None, the whole data set is used.
            verbose (bool,optional): turn verbose on/off.

        returns:
            None
        """
        # gather y's
        ys = np.zeros((len(self.data[ref].x), self.get_spectra_count()))
        for i in range(self.get_spectra_count()):
            ys[:, i] = self.data[i].y
        if ranges is None:
            ranges = [[min(self.data[ref].x), max(self.data[ref].x)]]
            x = self.data[ref].x
        else:
            x, ys = extract(self.data[ref].x, ys, ranges=ranges)

        if mode == 'cross-correlation' or mode == 'cc':
            if self.step is None:  # calculate step if it hasn't yet
                self.check_step_x()

            # calculate cross-correlation
            for i in range(self.get_spectra_count()):
                cross_correlation = np.correlate(ys[:, ref], ys[:, i], mode='Same')
                self.shifts[i]        = np.argmax(cross_correlation)
                self.shifts_lenght[i] = self.shifts[i]*self.step

                if verbose:
                    print(f'spectrum {i} shift calculated!')
            self.shifts        -= self.shifts[ref]
            self.shifts        *= -1
            self.shifts_lenght -= self.shifts_lenght[ref]
            self.shifts_lenght *= -1


        elif mode == 'max':
              # if self.step is None:
                # print("Warning: Step check hasn't been run yet. ``Spectra.shifts`` is not trustworth. Use ``Spectra.shifts_lenght`` instead or run  Spectra.check_step_x()")
            j_ref = np.argmax(ys[:, ref])
            for i in range(self.get_spectra_count()):
                j = np.argmax(ys[:, i])
                if self.step is None:
                    self.shifts[i] = np.NaN
                else:
                    self.shifts[i]    = j - j_ref
                self.shifts_lenght[i] = x[j] - x[j_ref]

                if verbose:
                    print(f'spectrum {i} shift calculated!')

        elif mode == 'elastic':
            # if self.step is None:
            #     print("Warning: Step check hasn't been run yet. ``Spectra.shifts`` is not trustworth. Use ``Spectra.shifts_lenght`` instead or run  Spectra.check_step_x()")
            s_ref = Spectrum(x=x, y=ys[:, ref])
            s_ref.guess_elastic_peak(**kwargs)
            for i in range(self.get_spectra_count()):
                s = Spectrum(x=x, y=ys[:, i])
                s.guess_elastic_peak(**kwargs)

                self.shifts_lenght[i] = s.elastic_c - s_ref.elastic_c
                if self.step is None:
                    self.shifts[i] = np.NaN
                else:
                    self.shifts[i]  = int(round(self.shifts_lenght[i]/self.step))

                if verbose:
                    print(f'spectrum {i} shift calculated!')
        else:
            raise ValueError('mode not valid (cross-correlation, max, elastic).')

        self._shift_mode   = mode

        try: ranges[0][0]
        except TypeError: (ranges, )
        self._shift_ranges = ranges

    def apply_shifts(self, mode='roll'):
        """Shift data.

        Args:
            shift (float or int): shift value.
            mode (string, optional): If None, the best mode will be selected.
                If ``mode='x'`` or ``mode='hard'``, y is fully preserved
                while x is shifted. If ``mode='y'``, ``'interp'``, or ``'soft'``, x is preserved
                while y is interpolated with a shift. If ``mode='roll'``, x is also preserved
                and y elements are rolled along the array (``shift`` value must be an integer.
                x must be evenly spaced for roll to make sense. The script does not check that.

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
        if mode == 'y' or mode == 'interp' or mode=='soft':
            shifts_lenght = True
        elif mode == 'x' or mode == 'hard':
            shifts_lenght = True
        elif mode == 'roll' or mode == 'rotate':
            shifts_lenght = False
        else:
            raise ValueError("mode not recognized (valid: 'y', 'x', 'roll').")

        if shifts_lenght:
            if np.NaN in self.shifts_lenght:
                raise ValueError('shifts not fully defined (or not defined at all). Use Spectra.calculate_shifts().')
            for i in range(self.get_spectra_count()):
                self.data[i].apply_shift(value=-self.shifts_lenght[i], mode=mode)
        else:
            if np.NaN in self.shifts_lenght and np.NaN in self.shifts:
                raise ValueError('shifts not fully defined (or not defined at all). Use Spectra.calculate_shifts().')
            elif np.NaN in self.shifts:
                raise ValueError('Integer shifts not fully defined (or not defined at all). However, float shifts seems defined. Maybe try using Spectra.check_step_x(). If possible, integer shifts will be calculated.')
            for i in range(self.get_spectra_count()):
                self.data[i].apply_shift(value=-self.shifts[i], mode=mode)

    def check_same_x(self, max_error=0.001):
        """Compare spectra to see if they have same x-coordinates.

            args:
                max_error (number, optional): percentage value of the max error.

            raises:
                ValueError: If any x-coodinate of any two spectrum is different.
        """
        if self.step is None:  # calculate step if it hasn't yet
            self.check_step_x(max_error=max_error)

        # check x between spectra
        for idx, s in enumerate(self.data):
            try:
                step = s.x[1] - s.x[0]
                if max(abs(s.x - self.data[idx+1].x))*100/step > max_error:
                    self._x = None
                    raise ValueError(f"x axis of spectrum {idx} and {idx+1} seem to be different.")
            except IndexError:
                pass

        self._x = self.data[0].x

    def calculate_sum(self):
        """Sum all spectra."""

        if self.x is None:
            self.check_same_x()

        y = np.zeros(len(self.x))
        for i in range(self.get_spectra_count()):
            y += self.data[i].data[:, 1]
        # print(type(y))
        self._sum = Spectrum(x=self.x, y=y)
        return self._sum

    def save(self, folderpath='./', prefix='spectrum_', suffix='.dat', zfill=None, delimiter=', ', header='', data_format='%.4f', verbose=False):
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

        """
        folderpath = Path(folderpath)

        if zfill is None:
            zfill = n_digits(self.get_spectra_count()-1)[0]

        for idx, s in enumerate(self.data):
            filename = f'{prefix}' + f'{idx}'.zfill(zfill) + f'{suffix}'
            if verbose:
                print(f'Saving {filename}')
            s.save(filepath=folderpath/filename, delimiter=delimiter, header=header, data_format='%.4f')
        if verbose:
            print('Done')

    def calib(self, value):
        """Calibrate data (from length to energy).

        args:
            value (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].

        returns:
            None
        """
        for s in self.data:
            s.calib(value=value)

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
                start = max([min(s.x) for s in self.data])
            if stop is None:
                stop = min([max(s.x) for s in self.data])

        for s in self.data:
            s.interp(x=x, start=start, stop=stop, num=num, step=step)

    def compress(self, ranges):
        for s in self.data:
            s._x, s._y  = extract(s.x, s.y, ranges)

    def floor(self, value=None, n=8):
        for s in self.data:
            s.floor(value=value, n=n)

    def plot(self, ax=None, idx='all', normalized=False, vertical_increment=0, shift=0, factor=1, show_ranges=False, **kwargs):
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
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if idx == 'all':
            for i in range(self.get_spectra_count()):
                self.data[i].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment*-i, shift=shift, factor=factor, label=i, **kwargs)
        else:
            try:
                if len(idx) == 1:  # (1)
                    self.data[idx[0]].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment*0, shift=shift, factor=factor, label=idx[0], **kwargs)
                else:  # (1, 2, ...)
                    for i in range(len(idx)):
                        self.data[i].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment*-i, shift=shift, factor=factor, label=i, **kwargs)
            except TypeError:  # 1
                self.data[idx].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, label=idx, **kwargs)

        plt.legend()

        if show_ranges:
            if self.shift_ranges is None:
                raise ValueError('Shift range not defined. Use calculate_shifts().')
            else:
                for r in self.shift_ranges:
                    plt.axvline(r[0], color='green', linewidth=1.2, zorder=10)
                    plt.axvline(r[1], color='red', linewidth=1.2, zorder=10)

        return ax
