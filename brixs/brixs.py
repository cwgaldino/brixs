#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module for analysis of RIXS spectra.

.. autosummary::

    Image
    PhotonEvents
    Spectrum
    Spectra
"""

# IDEAS:
# - when you crop the data in Spectra, you modify the spectrum. There's no way to
# go back (like we do with the shift and offset). Maybe figure a away to go back.
#
#
# Todo:
#     * flip() or transpose() method
#     * noise_filter() method (remove cosmic rays)
#     * save image as tiff
#
#

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
import numbers
from collections.abc import Iterable, MutableMapping
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

# backpack
from .backpack.filemanip import load_Comments, filelist
from .backpack.arraymanip import index, moving_average, extract, shifted, sort, get_attributes, derivative
from .backpack.arraymanip import is_integer, all_equal, factors, flatten, check_monotonicity, fix_monotonicity
from .backpack.figmanip import n_digits, n_decimal_places, figure, set_window_position
from .backpack.interact import query
from .peaks import Peak, Peaks, Collection
from .backpack.model_functions import voigt_fwhm, dirac_delta

# common definitions ===========================================================
from .config import settings

cc   = ['cross-correlation', 'cc']
roll = ['roll', 'rotate', 'r', 'rot']
hard = ['hard', 'x', 'h', 'Hard']
soft = ['soft', 'Soft', 'interp', 'y', 's']
relative = ['relative', 'r', 'rel']
absolute = ['a', 'abs', 'absolute']
increasing = ['inc', 'i', 'up', 'increasing', 'increasingly']
decreasing = ['dec', 'd', 'down', 'decreasing', 'decreasingly']

# support class ================================================================
class _Meta(type):
    """Metaclass to facilitate creation of read-only and non-removable attributes."""
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

# Interpreters =================================================================
def _axis_interpreter(axis):
    """Allows for flexibility when signaling axis direction.

    It follows numpy's convention. ``0`` is the vertical axis, while ``1`` is the horizontal one.

    Args:
        axis (int or str) = r, row, rows, hor, horizontal, x, or any
            string that starts with 'c' or 'h', will return 1. On the other hand,
            c, col, column, columns, v, ver, vertical, y, or any string that starts with
            'r' or 'v', will return 0.

    Returns:
        0 or 1
    """
    if type(axis)!=str:
        if (axis != 0 and axis != 1):
            raise ValueError("axis must be 0 ('y') or 1 ('x')")
        elif axis == 1:
            return 1
        elif axis == 0:
            return 0
    elif type(axis)==str and (axis=='x' or axis.startswith('r') or axis.startswith('h')):
        return 1
    elif axis == 0:
        return 0
    elif type(axis)==str and (axis=='y' or axis.startswith('c') or axis.startswith('v')):
        return 0
    else:
        raise ValueError("axis must be 0 ('y') or 1 ('x')")

def _bins_interpreter(*args, **kwargs):
    """Allows for flexible binning parameters.

    Keyword arguments (kwargs) cannot be mixed with positional arguments. If
    positional arguments are used, they are assumend to be the number of bins.
    Size of bins cannot be assigned via positional arguments.

    Args:
        shape (tuple): shape of data. MUST BE PASSED AS KEYWORD argument.
        factor_check (bool, optional): if True, it will assert that nbins or
            bins_size are a factor of shape (divisible). MUST BE PASSED AS
            KEYWORD argument.

    kwargs:
        nbins (int or tuple, optional): number of bins. If one value is given,
            this is used for both vertical (number of rows) and horizontal
            (number of columns) directions. If two values are given,
            they are used separetely for the number of rows and number of
            columns, respectively. If the number of pixels in the image
            cannot be divided by the selected number of bins, it will raise an error.
            ``bins_size`` is calculated. nbins overwrites all other arguments.
        nrows (int, optional): number of rows.
        ncolumns (int, optional): number of columns.
        bins_size (int or tuple, optional): size of the bins. If one value is given,
            this is used for both vertical (number of rows) and horizontal
            (number of columns) directions. If two values are given,
            they are used separetely for rows and columns size, respectively.
            If number of pixels cannot be divided by ``bins_size``, the value of
            ``bins_size`` will be recalculated to the closest possible value.
            ``nbins`` is calculated.
        rows_size (int or float, optional): size of the rows.
        columns_size (int or float, optional): size of the columns.

    Returns:
        bins, bins_size
    """
    shape         = kwargs['shape']
    factor_check  = kwargs['factor_check']
    del kwargs['factor_check']
    del kwargs['shape']

    # initial check
    if kwargs != {} and args != ():
        raise AttributeError('Cannot mix keyword arguments with positional arguments.')

    # initialization
    nbins        = None
    nrows        = None
    ncolumns     = None
    bins_size    = None
    rows_size    = None
    columns_size = None

    # keyword arguments
    if 'nbins' in kwargs:
        nbins = kwargs['nbins']
        if isinstance(nbins, Iterable):
            if len(nbins) == 1:
                nrows, ncolumns = nbins[0], nbins[0]
            else:
                nrows, ncolumns = nbins[0], nbins[1]
        else:
            nrows, ncolumns = nbins, nbins
    else:
        if 'nrows' in kwargs:
            nrows = kwargs['nrows']
        if 'ncolumns' in kwargs:
            ncolumns = kwargs['ncolumns']

    if 'bins_size' in kwargs:
        bins_size = kwargs['bins_size']
        if isinstance(bins_size, Iterable):
            if len(bins_size) == 1:
                rows_size, columns_size = bins_size[0], bins_size[0]
            else:
                rows_size, columns_size = bins_size[0], bins_size[1]
        else:
            rows_size, columns_size = bins_size, bins_size
    else:
        if 'rows_size' in kwargs:
            rows_size = kwargs['rows_size']
        if 'columns_size' in kwargs:
            columns_size = kwargs['columns_size']

    # positional arguments
    if len(args) == 1:
        if isinstance(args[0], Iterable):
            if len(args[0]) == 1:
                nrows, ncolumns = args[0][0], args[0][0]
            elif len(args[0]) > 1:
                nrows, ncolumns = args[0][0], args[0][1]
            else:
                raise AttributeError(f'Cannot understand the arguments passed to the function.\nArguments passed: {args[0]}.\nExamples:\n10\n(10)\n10, 5\n(10, 5)\n')
        else:
            nrows, ncolumns = args[0], args[0]
        nbins = np.array((nrows, ncolumns))
    elif len(args) == 2:
        if isinstance(args[0], Iterable) or isinstance(args[1], Iterable):
            raise AttributeError(f'Cannot understand the arguments passed to the function.\nArguments passed: {args[0]}.\nExamples:\n10\n(10)\n10, 5\n(10, 5)\n')
        else:
            nrows, ncolumns = args[0], args[1]
        nbins = np.array((nrows, ncolumns))
    elif len(args) > 2:
        raise AttributeError(f'Cannot understand the arguments passed to the function.\nArguments passed: {args[0]}.\nExamples:\n10\n(10)\n10, 5\n(10, 5)\n')


    # after this part, either nbins, nrows, and ncolumns are defined or bins_size, row_size, and column_size are defined

    # final
    if bins_size is None and nbins is None:
        raise AttributeError('Number of bins or bin size must be passed as an argument.')
    elif nbins is not None:
        if nrows is None:    nrows    = shape[0]
        if ncolumns is None: ncolumns = shape[1]
        if nrows <= 0 or ncolumns <= 0 or is_integer(nrows)==False or is_integer(ncolumns)==False:
            raise ValueError("Number of bins must be a positive integer.")
        else:
            if factor_check:
                assert shape[1] % ncolumns == 0, f"The {shape[1]} pixels in a row is not evenly divisible by {ncolumns}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[1])))}"
                assert shape[0] % nrows    == 0, f"The {shape[0]} pixels in a column is not evenly divisible by {nrows}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[0])))}"

            nbins = np.array((nrows, ncolumns))
            # bins_size = np.array((shape[0]//nbins[0], shape[1]//nbins[1]))
            bins_size = np.array((shape[0]/nbins[0], shape[1]/nbins[1]))

    else:
        if rows_size is None:    rows_size =    shape[1]
        if columns_size is None: columns_size = shape[0]
        if rows_size <= 0 or columns_size <= 0 or is_integer(rows_size)==False or is_integer(columns_size)==False:
            raise ValueError("Size of bins must be a positive integer.")
        else:
            if factor_check:
                assert shape[1] % columns_size == 0, f"The {shape[1]} pixels in a row is not evenly divisible by {columns_size}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[1])))}"
                assert shape[0] % rows_size == 0,    f"The {shape[0]} pixels in a column is not evenly divisible by {rows_size}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[0])))}"

            nbins = np.array((round(shape[0]/rows_size), round(shape[1]/columns_size)))
            # bins_size = np.array((shape[0]//_nbins[0], shape[1]//_nbins[1]))
            bins_size = np.array((shape[0]/nbins[0], shape[1]/nbins[1]))

    return nbins, bins_size

def _check_ranges(ranges, vmin, vmax):
    """check if ranges is the right format.

    If any item of ranges is None, this item is replaced by the min or max
        value of the data.
    """
    text = 'Ranges should be a pair (x_init1, x_final1) or a list of pairs like this: ((x_init1, x_final1), (x_init2, x_final2), ...)\nUse None to indicate the minimum or maximum x value of the data.'

    # check format
    if ranges is None:
        ranges = ((vmin, vmax),)
    elif isinstance(ranges, Iterable) == True:
        if isinstance(ranges[0], Iterable) == False:
                ranges = (ranges, )
    else:
        raise ValueError(text)

    # check pairs
    ranges = list(ranges)
    for i in range(len(ranges)):
        if len(ranges[i]) == 2:
            if None in ranges[i]:
                ranges[i] = list(ranges[i])
                if ranges[i][0] is None:
                     ranges[i][0] = vmin
                if ranges[i][1] is None:
                     ranges[i][1] = vmax
        else:
            raise ValueError(f'Ranges pair {r} is not a valid pair.\n'+text)
    return tuple(ranges)

# BRIXS ========================================================================
class Image(metaclass=_Meta):
    """Image object.

    Args:
        data (2D array, optional): Image.
        filepath (string or path object, optional): filename or file handle.
            If the filename ends in .gz, the file is automatically saved in
            compressed gzip format. This is overwriten by data.

    Attributes:
        data (2D np.array): This is where we store the Image.
        vmin, vmax (float, read only): Minimum value in data.
        shape (tuple, read only): Shape of data (vertical size, horizontal size).
        x_centers, y_centers (np.array): pixel center labels.
        x_edges, y_edges (np.array): pixel edges labels.
        histogram (brixs.Spectrum): Data intensity histogram.

        nbins (tuple): Number of bins (number of rows, number of columns).
        bins_size (tuple): Bins size (size of rows, size of columns).
        reduced (brixs.Image, read only): Binned image.

        shifts_v, shifts_h (np.array): Shift values in the vertical and
            horizontal direction.
        p (np.array, read only): polynomial values returned by calculated_shifts().
        f (function, read only): funcion returned by calculated_shifts().
        calculated_shifts (brixs.Spectrum): Calculated shifts.

        spectrum_v, spectrum_h (brixs.Spectrum): Spectrum obtained by integrating
            pixels in the vertical and horizontal direction.
        spectrum (brixs.Spectrum): same as spectrum_h.
        columns, rows (brixs.Spectra): Spectra obtained from each pixel columns
            or row.

    Methods:
        save()
        load()

        floor()
        crop()

        pcolormesh()
        imshow()
        plot()

        binning()
        calculate_histogram()
        calculate_spectrum()
        calculate_shifts()
        set_shifts()
        fix_curvature()
    """
    # read only and non-removable arguments
    _read_only = ['shape', 'vmin', 'vmax', 'reduced', 'calculated_shifts']
    _non_removable = []

    def __init__(self, *args, **kwargs):
        # besic attr
        self._data      = None
        self._vmin      = None
        self._vmax      = None
        self._shape     = None
        self._x_centers = None
        self._y_centers = None
        self._x_edges   = None
        self._y_edges   = None

        # binning attr
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None

        # shift attr
        self._calculated_shifts = None
        self._shifts_v = None
        self._shifts_h = None

        # set data
        data, filepath = self._sort_args(args, kwargs)
        if data is not None:
            self.data = copy.deepcopy(data)
        elif filepath is not None:
            self.load(filepath)

    def _sort_args(self, args, kwargs):
        """checks initial arguments.

        Keyword arguments (kwargs) cannot be mixed with positional arguments.

        For positional arguments, if one data set is passed, it assumes it is
            `data`. If this one argument is of type string or Pathlib.Path, it
            assumes it is a filepath.

        Raises:
            AttributeError: if kwargs and args cannot be read.

        Returns:
            data, filepath
        """
        # initial check
        if kwargs != {} and args != ():
            raise AttributeError('Cannot mix keyword arguments with positional arguments. Keyword arguents are `data`, and `filepath`.')

        # initialization
        data = None
        filepath = None

        # keyword arguments
        if 'data' in kwargs:
            data = kwargs['data']
        if 'filepath' in kwargs:
            filepath = kwargs['filepath']

        # positional arguments
        if len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepath = args[0]
            elif isinstance(args[0], Iterable):
                data = args[0]
        elif len(args) > 2:
            raise AttributeError('brixs.Image() cannot figure out the data out of the arguments passed. Maybe use keyword arguments (data, filepath).')

        return data, filepath

    def _transfer_attributes(self, object):
        """Transfer user defined attributes to output objects."""

        # # list of attributes NOT to transfer
        # do_not_transfer = ['_data', '_vmin', '_vmax', '_shape', 
        #                    '_nbins', '_bins_size', 
        #                    '_reduced', '_shifts_v', '_shifts_h', 
        #                    '_p', '_f', '_calculated_shifts', 
        #                    '_x_centers', '_y_centers', 
        #                    '_x_edges', '_y_edges']
        
        # dict = get_attributes(self)
        # for n in dict:
        #     if n not in do_not_transfer:
        #         object.__setattr__(n, self.__getattribute__(n))
        # return object
    
        # copy user defined attributes
        for attr in self.__dict__:
            if attr.startswith('_') == False:
                object.__setattr__(attr, self.__dict__[attr])

        return object

    def __len__(self):
        if self._data is None:
            return 0
        else:
            return len(self._data)

    def __add__(self, object):
        if isinstance(object, Image):
            if self.shape == object.shape:
                if np.issubdtype(self.data.dtype, np.integer) and np.issubdtype(object.data.dtype, np.integer):
                    dtype = 'int'
                else:
                    dtype = 'float'
                final = Image(data=np.add(self.data, object.data, dtype=dtype))
                return self._transfer_attributes(final)
            else:
                raise ValueError(f'Shape is different.\nShape 1: {self.shape}\nShape 2: {im.shape}')
        elif isinstance(object, (np.floating, float, int)):
            if np.issubdtype(self.data.dtype, np.integer) and isinstance(object, (np.floating, float))==False:
                dtype = 'int'
            else:
                dtype = 'float'
            final = Image(data=np.add(self.data, object, dtype=dtype))
            return self._transfer_attributes(final)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with Image')

    def __sub__(self, object):
        if isinstance(object, Image):
            if self.shape == object.shape:
                if np.issubdtype(self.data.dtype, np.integer) and np.issubdtype(object.data.dtype, np.integer):
                    dtype = 'int'
                else:
                    dtype = 'float'
                final = Image(data=np.subtract(self.data, object.data, dtype=dtype))
                return self._transfer_attributes(final)
            else:
                raise ValueError(f'Shape is different.\nShape 1: {self.shape}\nShape 2: {im.shape}')
        elif isinstance(object, (np.floating, float, int)):
            if np.issubdtype(self.data.dtype, np.integer) and isinstance(object, (np.floating, float))==False:
                dtype = 'int'
            else:
                dtype = 'float'
            final = Image(data=np.subtract(self.data, object, dtype=dtype))
            return self._transfer_attributes(final)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with Image')

    def __mul__(self, object):
        if isinstance(object, Image):
            if self.shape == object.shape:
                if np.issubdtype(self.data.dtype, np.integer) and np.issubdtype(object.data.dtype, np.integer):
                    dtype = 'int'
                else:
                    dtype = 'float'
                final = Image(data=np.multiply(self.data, object.data, dtype=dtype))
                return self._transfer_attributes(final)
            else:
                raise ValueError(f'Shape is different.\nShape 1: {self.shape}\nShape 2: {im.shape}')
        elif isinstance(object, (np.floating, float, int)):
            if np.issubdtype(self.data.dtype, np.integer) and isinstance(object, (np.floating, float))==False:
                dtype = 'int'
            else:
                dtype = 'float'
            final = Image(data=np.multiply(self.data, object, dtype=dtype))
            return self._transfer_attributes(final)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with Image')

    def __div__(self, object):
        if isinstance(object, Image):
            if self.shape == object.shape:
                if 0 in object.data:
                    raise ZeroDivisionError(f'Image contain zeros. Cannot divide by zero.')
                else:
                    dtype = 'float'
                    final = Image(data=np.divide(self.data, object.data, dtype=dtype))
                    # final = Image(data = self.data / object.data)
                    return self._transfer_attributes(final)
            else:
                raise ValueError(f'Shape is different.\nShape 1: {self.shape}\nShape 2: {im.shape}')
        elif isinstance(object, (np.floating, float, int)):
            if object == 0:
                raise ZeroDivisionError(f'Cannot divide by zero.')
            else:
                dtype = 'float'
                final = Image(data=np.divide(self.data, object, dtype=dtype))
                return self._transfer_attributes(final)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Image')

    def __truediv__(self, object):
        return self.__div__(object)


    @property
    def data(self):
        return copy.deepcopy(self._data)
    @data.setter
    def data(self, value):
        # basic attr
        self._data  = np.array(value, dtype='float')
        self._vmin  = min([min(x) for x in self.data])
        self._vmax  = max([max(x) for x in self.data])
        self._shape = (self.data.shape[0], self.data.shape[1])
        self.x_centers = None
        self.y_centers = None

        # binning attr
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None

        # shift attr
        self._calculated_shifts = None
        self._shifts_v = np.zeros(self.data.shape[1])
        self._shifts_h = np.zeros(self.data.shape[0])
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x_centers(self):
        return self._x_centers
    @x_centers.setter
    def x_centers(self, value):
        if value is None:
            value = np.arange(0, self.data.shape[1])
        elif isinstance(value, Iterable):
            if len(value) != self.shape[1]:
                raise ValueError(f"Length of x must be the same as the number of pixels in the horizontal direction.\nNumber of pixels = {self.data.shape[1]}\nLength of the array = {len(value)}")
        else:
            raise ValueError(f"x must be None or an iterable (list, tuple, or 1D array)")
        # self._x_centers = Spectrum(x=value, y=np.zeros(len(value)))
        self._x_centers = np.array(value, dtype='float')
        self._x_edges = None
    @x_centers.deleter
    def x_centers(self):
        self._x_centers = np.arange(0, self.data.shape[1])

    @property
    def y_centers(self):
        return self._y_centers
    @y_centers.setter
    def y_centers(self, value):
        if value is None:
            value = np.arange(0, self.data.shape[0])
        elif isinstance(value, Iterable):
            if len(value) != self.shape[0]:
                raise ValueError(f"Length of y must be the same as the number of pixels in the vertical direction.\nNumber of pixels = {self.data.shape[0]}\nLength of the array = {len(value)}")
        else:
            raise ValueError(f"y must be None or an iterable (list, tuple, or 1D array)")
        # self._y_centers = Spectrum(x=value, y=np.zeros(len(value)))
        self._y_centers = np.array(value, dtype='float')
        self._y_edges = None
    @y_centers.deleter
    def y_centers(self):
        self._y_centers  = np.arange(0, self.data.shape[0])

    @property
    def x_edges(self):
        return self._x_edges
    @x_edges.setter
    def x_edges(self, value):
        if value is None:
            self._x_edges = None
            return
        elif isinstance(value, Iterable):
            if len(value) != self.shape[1]+1:
                raise ValueError(f"Length of x_edges must be the same as the number of edges (number of pixels + 1) in the horizontal direction.\nNumber of edges = {self.data.shape[1]+1}\nLength of the array = {len(value)}")
        else:
            raise ValueError(f"x_edges must be None or an iterable (list, tuple, or 1D array)")
        # self._x_edges   = Spectrum(x=value, y=np.zeros(len(value)))
        self._x_edges   = np.array(value, dtype='float')
        centers = moving_average(value, n=2)
        # self._x_centers = Spectrum(x=centers, y=np.zeros(len(centers)))
        self._x_centers = np.array(centers, dtype='float')
    @x_edges.deleter
    def x_edges(self):
        self._x_edges = None

    @property
    def y_edges(self):
        return self._y_edges
    @y_edges.setter
    def y_edges(self, value):
        if value is None:
            self._y_edges = None
            return
        elif isinstance(value, Iterable):
            if len(value) != self.shape[0]+1:
                raise ValueError(f"Length of y_edges must be the same as the number of edges (number of pixles + 1) in the vertical direction.\nNumber of edges = {self.data.shape[0]+1}\nLength of the array = {len(value)}")
        else:
            raise ValueError(f"y must be None or an iterable (list, tuple, or 1D array)")
        # self._y_edges   = Spectrum(x=value, y=np.zeros(len(value)))
        self._y_edges   = np.array(value, dtype='float')
        centers = moving_average(value, n=2)
        # self._y_centers = Spectrum(x=centers, y=np.zeros(len(centers)))
        self._y_centers = np.array(centers, dtype='float')
    @y_centers.deleter
    def y_centers(self):
        self._y_edges  = None

    @property
    def shifts_v(self):
        return copy.deepcopy(self._shifts_v)
    @shifts_v.setter
    def shifts_v(self, value):
        self.set_shifts(value=value, axis=0)
    @shifts_v.deleter
    def shifts_v(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shifts_h(self):
        return copy.deepcopy(self._shifts_h)
    @shifts_h.setter
    def shifts_h(self, value):
        self.set_shifts(value=value, axis=1)
    @shifts_h.deleter
    def shifts_h(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nbins(self):
        return self._nbins
    @nbins.setter
    def nbins(self, value):
        if type(value) != str:
            self.binning(nbins=value)
        elif value == 'guess':
            raise NotImplementedError('not implemented yet')
        elif value == 'possible':
            print(f"The {self.shape[1]} pixels in a row is evenly divisible by: {np.sort(list(factors(self.shape[1])))}\nThe {self.shape[0]} pixels in a column is evenly divisible by: {np.sort(list(factors(self.shape[0])))}")
        else:
            raise ValueError("Not a valid option of `nbins`.\nValid options are: 'guess', a number, a tuple, or a list.")
    @nbins.deleter
    def nbins(self):
        raise AttributeError('Cannot delete object.')

    @property
    def bins_size(self):
        return self._bins_size
    @bins_size.setter
    def bins_size(self, value):
        self.binning(bins_size=value)
    @bins_size.deleter
    def bins_size(self):
        raise AttributeError('Cannot delete object.')

    @property
    def histogram(self):
        return self.calculate_histogram()
    @histogram.setter
    def histogram(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @histogram.deleter
    def histogram(self):
        raise AttributeError('Cannot delete object.')

    @property
    def columns(self):
        ss = Spectra(n=self.shape[1])
        for i in range(self.shape[1]):
            ss[i] = Spectrum(x=self.y_centers, y=self.data[:, i])
        return ss
    @columns.setter
    def columns(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @columns.deleter
    def columns(self):
        raise AttributeError('Cannot delete object.')

    @property
    def rows(self):
        ss = Spectra(n=self.shape[0])
        for i in range(self.shape[0]):
            ss[i] = Spectrum(x=self.x_centers, y=self.data[i, :])
        return ss
    @rows.setter
    def rows(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @rows.deleter
    def rows(self):
        raise AttributeError('Cannot delete object.')

    @property
    def spectrum(self):
        return self.spectrum_h
    @spectrum.setter
    def spectrum(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @spectrum.deleter
    def spectrum(self):
        raise AttributeError('Cannot delete object.')

    @property
    def spectrum_v(self):
        return self.calculate_spectrum(axis=0)
    @spectrum_v.setter
    def spectrum_v(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @spectrum_v.deleter
    def spectrum_v(self):
        raise AttributeError('Cannot delete object.')

    @property
    def spectrum_h(self):
        return self.calculate_spectrum(axis=1)
    @spectrum_h.setter
    def spectrum_h(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @spectrum_h.deleter
    def spectrum_h(self):
        raise AttributeError('Cannot delete object.')


    def save(self, filepath=None, only_data=False,  check_overwrite=True, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format. If None, filepath can be inferred from a
                user defined attr 'im.filepath'. Default is None.  Last used 
                filepath is saved to im.filepath.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user wants to overwrite file.
            **kwargs: kwargs are passed to ``np.savetxt()`` that saves the data.

        If not specified, the following parameters are passed to `numpy.savetxt()`_:

        Args:
            fmt (string, or list, optional): A single format (like ``%10.5f``), or a
                sequence of formats. See numpy's documentation for more information.
                If not specified, best fmt is calculated based on the
                number of decimal places of the data. Specifying fmt makes the code
                runs a little faster (not much tough, but it might make a difference
                if saving a lot of files).
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (", ").
            newline (str, optional): String or character separating lines.
                Default is ``\n``.
            header (bool, optional): String that will be written at the beginning of the file.
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        # filepath
        if filepath is None:
            try: 
                filepath = Path(self.filepath)
            except AttributeError:
                raise TypeError("Missing 1 required argument: 'filepath'")
        
        # check overwrite
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')
                
        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([n_decimal_places(x) for x in flatten(self._data)])
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # save
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
            np.savetxt(Path(filepath), self._data, **kwargs)
        else:
            if 'header' not in kwargs:
                kwargs['header'] = ''
            else:
                kwargs['header'] += '\n'
            dict = get_attributes(self)
            kwargs['header'] += '==== brixs Image ===='  + '\n'
            for n in dict:
                if n not in ['_data', '_reduced', '_calculated_shifts', '_histogram', '_spectrum_h', '_spectrum_v',]:
                    if isinstance(dict[n], Iterable):
                        kwargs['header'] += f'{n}: {list(dict[n])}'  + '\n'
                    else:
                        kwargs['header'] += f'{n}: {dict[n]}'  + '\n'
            np.savetxt(Path(filepath), self._data, **kwargs)

        # save filepath
        self.filepath = filepath

    def load(self, filepath, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first 
                decompressed. Last used filepath is saved to im.filepath.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``# ``. Attributes picked up
                from the header will be loaded too.

        Returns:
            None

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # read data
        data = np.genfromtxt(Path(filepath), **kwargs)
        self.data = data

        # read header
        header = load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
        attr_start = 0
        binning_flag = False  # check if data needs binning

        # find where attributes listing starts
        for i, line in enumerate(header):
            if '==== brixs Image ====' in line:
                attr_start = i
                break
            attr_start = -1

        # read attributes
        if attr_start != -1:
            for i, line in enumerate(header[attr_start+1:-1]):

                # extract name and value
                name = line[1:-1].split(':')[0].strip()
                value = eval(line[1:-1].split(':')[1].strip())


                ### DEALING WITH ATTRIBUTES THAT NEED TO RUN SOMETHING ###
                if name in ['_nbins']:
                    if value[0] > 0 and value[1] > 0:
                        binning_flag = True
                elif name in ['_bins_size']:
                    if value[0] > 0 and value[1] > 0:
                        binning_flag = True
                ### DEALING WITH OTHER ATTRIBUTES ###
                elif name not in ['_vmin', '_vmax', '_shape']:  # except these attrs
                    try:
                        setattr(self, name, value)
                    except Exception as e:
                        print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')

        ### RUN ###
        if binning_flag:
            try:
                self.binning(nbins=self.nbins, bins_size=self.bins_size)
            except:
                pass
        
        # save filepath
        self.filepath = str(filepath)

    def get_user_defined_attrs(self):
        """return attrs that are user defined."""
        default_attrs =  ['_data', '_vmin', '_vmax', '_shape', 
                           '_nbins', '_bins_size', 
                           '_reduced', '_shifts_v', '_shifts_h', 
                           '_calculated_shifts', 
                           '_x_centers', '_y_centers', 
                           '_x_edges', '_y_edges']
        return [key for key in self.__dict__.keys() if key not in default_attrs]

    # plot and visualization
    def pcolormesh(self, ax=None, colorbar=False, **kwargs):
        """Display data as a mesh. Wrapper for `matplotlib.pyplot.pcolormesh()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
            **kwargs: kwargs are passed to `matplotlib.pyplot.pcolormesh()`_.

        If not specified, the following parameters are passed to `matplotlib.pyplot.pcolormesh()`_:

        Args:
            cmap: The Colormap instance. Default is 'jet'.
            vmin: Minimum intensity that the colormap covers. The intensity histogram is
                calculated and vmin is set on the position of the maximum.
            vmax: Maximmum intensity that the colormap covers.  The intensity histogram is
                calculated and vmax is set to the value where the :del:` integral of the
                intensity histogram is 0.998 of the total integral after vmin.`
                intensity drops below 0.01 % of the maximum.
                **THIS MIGHT CHANGE IN THE FUTURE**.

        Returns:
            `matplotlib.collections.QuadMesh`_

        .. _matplotlib.pyplot.pcolormesh(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.pcolormesh.html
        .. _matplotlib.collections.QuadMesh: https://matplotlib.org/3.5.0/api/collections_api.html#matplotlib.collections.QuadMesh
        """
        # initialization
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figure()
                if settings.FIGURE_POSITION is not None:
                    try:
                        set_window_position(settings.FIGURE_POSITION)
                    except:
                        pass
            elif plt.get_fignums() == [] and settings.FIGURE_POSITION is not None:
                try:
                    set_window_position(settings.FIGURE_POSITION)
                except:
                    pass

        # kwargs
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'
        if 'vmin' not in kwargs:
            kwargs['vmin'] = self.histogram.x[np.argmax(self.histogram.y)]
        if 'vmax' not in kwargs:
            # vamx is set when histogram drops below 0.01% of the max
            x2fit = self.histogram.x[np.argmax(self.histogram.y)+1:]
            y2fit = self.histogram.y[np.argmax(self.histogram.y)+1:]
            data2fit = np.array([[i, j] for i, j in zip(x2fit, y2fit) if j > max(y2fit)*0.0001])  # clean zeros
            kwargs['vmax'] = data2fit[-1, 0]

        # fix monotonicity of labels x
        if check_monotonicity(self.x_centers) != 1:
            x, ordering = fix_monotonicity(self.x_centers, np.arange(len(self.x_centers)), mode='increasing')
            assert len(x)==len(self.x_centers), f'Cannot plot when Image.x have repeated elements.\nEither fix Image.x or set it to None.\nx: {self.x_centers}'
            ordering = [int(i) for i in ordering]
            data = copy.deepcopy(self.data)
            for i in ordering:
                if i != ordering[i]:
                    data[:, i] = self.data[:, ordering[i]]
        else:
                # check if x_edges is defined
            if self.x_edges is not None and self.y_edges is not None:
                x = self.x_edges
            else:
                x = self.x_centers
            data = self.data

        # fix monotonicity of labels y
        if check_monotonicity(self.y_centers) != 1:
            y, ordering = fix_monotonicity(self.y_centers, np.arange(len(self.y_centers)), mode='increasing')
            assert len(y)==len(self.y_centers), f'Cannot plot when Image.y have repeated elements.\nEither fix Image.x or set it to None.\ny: {self.y_centers}'
            ordering = [int(i) for i in ordering]
            data2 = copy.deepcopy(data)
            for i in ordering:
                if i != ordering[i]:
                    data2[i, :] = data[ordering[i], :]
        else:
            # check if y_edges is defined
            if self.x_edges is not None and self.y_edges is not None:
                y = self.y_edges
            else:
                y = self.y_centers
            data2 = data

        # plot
        X, Y = np.meshgrid(x, y)
        pos = ax.pcolormesh(X, Y, data2, **kwargs)

        # colorbar
        if colorbar:
            plt.colorbar(pos, aspect=50)

        return pos

    def imshow(self, ax=None, colorbar=False, hindexes=False, vindexes=False, verbose=True, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.imshow()`_.

        Warning:
            Pixels are always square. For irregular pixel row/columns, see Image.pcolormesh()

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
            verbose (bool, optional): if True, a warning will show up if data
                has iregular pixel sizes. Default is true.
            hindexes, vindexes (bool, optional): If True, the index number of a
                pixel column/row will be ploted. Default is False.
            **kwargs: kwargs are passed to `matplotlib.pyplot.imshow()`_.

        If not specified, the following parameters are passed to `matplotlib.pyplot.imshow()`_:

        Args:
            cmap: The Colormap instance. Default is 'jet'.
            aspect: The aspect ratio of the Axes. Default is 'auto'. If 'equal',
                an aspect ratio of 1 will be used (pixels will be square).
            origin: Location of the [0, 0] index. default is 'lower'.
            interpolation: The interpolation method used. Default is 'none'.
                Supported values are 'none', 'antialiased', 'nearest', 'bilinear',
                'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite',
                'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell',
                'sinc', 'lanczos', 'blackman'.
            extent: minimun and maximum x and y values. Default will be given by
                the Image.x and Image.y attributes.
            vmin: Minimum intensity that the colormap covers. The intensity histogram is
                calculated and vmin is set on the position of the maximum.
            vmax: Maximmum intensity that the colormap covers.  The intensity histogram is
                calculated and vmax is set to the value where the :del:` integral of the
                intensity histogram is 0.998 of the total integral after vmin.`
                intensity drops below 0.01 % of the maximum.
                **THIS MIGHT CHANGE IN THE FUTURE**.

        Returns:
            `matplotlib.image.AxesImage`_

        .. _matplotlib.pyplot.imshow(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.imshow.html
        .. _matplotlib.image.AxesImage: https://matplotlib.org/3.5.0/api/image_api.html#matplotlib.image.AxesImage
        """
        # initialization
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figure()

        # default arguments
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'
        if 'aspect' not in kwargs:
            kwargs['aspect'] = 'auto'
        if 'origin' not in kwargs:
            kwargs['origin'] = 'lower'
        if 'interpolation' not in kwargs:
            kwargs['interpolation'] = 'none'
        if 'vmin' not in kwargs:
            kwargs['vmin'] = self.histogram.x[np.argmax(self.histogram.y)]
        if 'vmax' not in kwargs:
            # vamx is set when histogram drops below 0.01% of the max
            x2fit = self.histogram.x[np.argmax(self.histogram.y)+1:]
            y2fit = self.histogram.y[np.argmax(self.histogram.y)+1:]
            data2fit = np.array([[i, j] for i, j in zip(x2fit, y2fit) if j > abs(max(y2fit)*0.0001)])  # clean zeros

            try:
                kwargs['vmax'] = data2fit[-1, 0]
            except IndexError:  # in case the max of y is too high
                kwargs['vmin'] = self.vmin
                kwargs['vmax'] = self.vmax


        # # x and y
        # assert check_monotonicity(self.x) == 1, f'x axis (Image.x) must be increasingly monotonic.\nData cannot be plotted by Image.imshow().\nPlease, use Image.plot() or set Image.x = None\nx: {self.x}'
        # assert check_monotonicity(self.y) == 1, f'y axis (Image.y) must be increasingly monotonic.\nData cannot be plotted by Image.imshow().\nPlease, use Image.plot() or set Image.y = None\ny: {self.y}'
        # x = np.linspace(self.x[0], self.x[-1], len(self.x))
        # y = np.linspace(self.y[0], self.y[-1], len(self.y))
        # if 'extent' not in kwargs:
        #     kwargs['extent'] = [x[0], x[-1], y[0], y[-1]]

        # fix monotonic of labels x
        if check_monotonicity(self.x_centers) != 1:
            x, ordering = fix_monotonicity(self.x_centers, np.arange(len(self.x_centers)), mode='increasing')
            assert len(x)==len(self.x_centers), f'Cannot plot when Image.x have repeated elements.\nEither fix Image.x or set it to None.\nx: {self.x_centers}'
            ordering = [int(i) for i in ordering]
            data = copy.deepcopy(self.data)
            for i in ordering:
                if i != ordering[i]:
                    data[:, i] = self.data[:, ordering[i]]
                extent_x = [min(x), max(x)]
        else:
            # check if x_edges is defined
            if self.x_edges is not None:
                x = self.x_edges
                extent_x = [min(x), max(x)]
            else:
                x = self.x_centers
                x = np.linspace(x[0], x[-1], len(x))
                dx  = np.mean(np.diff(x))
                extent_x = [x[0]-dx/2, x[-1]+dx/2]
            data = self.data

        # fix monotonic of labels y
        if check_monotonicity(self.y_centers) != 1:
            y, ordering = fix_monotonicity(self.y_centers, np.arange(len(self.y_centers)), mode='increasing')
            assert len(y)==len(self.y_centers), f'Cannot plot when Image.y have repeated elements.\nEither fix Image.x or set it to None.\ny: {self.y_centers}'
            ordering = [int(i) for i in ordering]
            data2 = copy.deepcopy(data)
            for i in ordering:
                if i != ordering[i]:
                    data2[i, :] = data[ordering[i], :]
            extent_y = [min(y), max(y)]
        else:
            # check if y_edges is defined
            if self.y_edges is not None:
                y = self.y_edges
                extent_y = [min(y), max(y)]
            else:
                y = self.y_centers
                y = np.linspace(y[0], y[-1], len(y))
                dy  = np.mean(np.diff(y))
                extent_y = [y[0]-dy/2, y[-1]+dy/2]

            data2 = data

        if 'extent' not in kwargs:
            kwargs['extent'] = np.append(extent_x, extent_y)

        # check irregular spacing (issues a warning)
        if verbose:
            sx = Spectrum(x=self.x_centers, y=self.x_centers)
            sy = Spectrum(x=self.y_centers, y=self.y_centers)
            try:
                sx.check_step()
                sy.check_step()
            except ValueError:
                print('Data seems to have irregular pixel size. Maybe plot it using Image.pcolormesh().\nTo turn off this warning set verbose to False.')

        # plot
        pos = ax.imshow(data2, **kwargs)

        # column label
        if hindexes:
            for i in range(self.data.shape[1]):
                plt.text(self.x_centers[i], self.y_centers[-1], i, ha='center')
        if vindexes:
            for i in range(self.data.shape[0]):
                plt.text(self.y_centers[i], self.x_centers[-1], i, va='center')

        # colorbar
        if colorbar:
            # delta = self.vmax - self.vmin
            # delta_digits =n_digits(self.vmax)[0] - n_digits(delta)[0]
            # # scientific notation
            # if self.vmax > 1e4:
            #     letter = 'e'
            #     n = delta_digits + 1
            # elif self.vmax < 1e-4:
            #     letter = 'e'
            #     d = str(self.vmax).split('.')[1]
            #     d1 = str(self.vmax).split('.')[1].lstrip('0')
            #     power = d - d1
            #
            #     d = str(delta).split('.')[1]
            #     d1 = str(delta).split('.')[1].lstrip('0')
            #     power1 = d - d1
            #
            #     n = power - power1 + 1
            # else:
            #     letter = 'f'
            #     if delta > 10:
            #         n = 0
            #     else:
            #         n = 4
            # print(delta)
            # print(delta_digits)
            #
            # fmt = f'%.{n}{letter}'

            plt.colorbar(pos, aspect=50)

        return pos

    def plot(self, *args, **kwargs):
        """Same as Image.imshow.

        See:
            :py:func:`Image.imshow` """
        return self.imshow(*args, **kwargs)

    def vlines(self, ax=None, mode='edges', *args, **kwargs):
        """edges or centers"""
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figure()

        assert self.y_edges is not None, 'y_edges cannot be None.'
        if mode == 'edges':
            plt.vlines(self.x_edges, self.y_edges[0], self.y_edges[-1], *args, **kwargs)
        elif mode == 'centers':
            plt.vlines(self.x_centers, self.y_edges[0], self.y_edges[-1], *args, **kwargs)
        else:
            raise AttributeError('mode must be either `edges` or `centers`.')
    
    def hlines(self, ax=None, mode='edges', *args, **kwargs):
        """edges or centers"""
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figure()

        assert self.x_edges is not None, 'x_edges cannot be None.'
        if mode == 'edges':
            plt.hlines(self.y_edges, self.x_edges[0], self.x_edges[-1], *args, **kwargs)
        elif mode == 'centers':
            plt.hlines(self.y_centers, self.x_edges[0], self.x_edges[-1], *args, **kwargs)
        else:
            raise AttributeError('mode must be either `edges` or `centers`.')

    def possible_nbins(self):
        """return possible values for nbins in the y (nrows) and x (ncols) directions."""
        return np.sort(list(factors(self.shape[0]))), np.sort(list(factors(self.shape[1])))

    def binning(self, *args, **kwargs):
        """Compute the 2D histogram of the data (binning of the data).

        Args:
            nbins (int or tuple, optional): number of bins. If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the number of rows and number of
                columns, respectively. If the number of pixels in the image
                cannot be divided by the selected number of bins, it will raise an error.
                ``bins_size`` is calculated.
            bins_size (int or tuple, optional): size of the bins. This overwrites
                the argument ``bins``. If one value is given,
                it is used for both x and y directions. If two values are given,
                they are used separetely for rows and columns size, respectively.
                If number of pixels cannot be divided by ``bins_size``, the value of
                ``bins_size`` will be recalculated to the closest possible value.
                ``nbins`` is calculated.

        Example:
            .. code-block:: python

                nbins = 10       # (10 rows, 10 columns)
                nbins = (10)     # (10 rows, 10 columns)
                nbins = (10, 5)  # (10 rows, 5 columns)

        Returns:
            binned image
        """
        kwargs['shape']        = self.shape
        kwargs['factor_check'] = True
        _nbins, _bins_size = _bins_interpreter(*args, **kwargs)
        # reduced = Image(np.add.reduceat(np.add.reduceat(self._data,            np.arange(0, self.shape[0], _bins_size[0]), axis=0),             np.arange(0, self.shape[1], _bins_size[1]), axis=1))
        # reduced = Image(np.add.reduceat(np.add.reduceat(self._data, map(float, np.arange(0, self.shape[0], _bins_size[0])), axis=0), map(float, np.arange(0, self.shape[1], _bins_size[1]), axis=1)))
        reduced = Image(np.add.reduceat(np.add.reduceat(self._data, list(map(float, np.arange(0, self.shape[0], _bins_size[0]))), axis=0), list(map(float, np.arange(0, self.shape[1], _bins_size[1]))), axis=1))

        _x_edges   = np.arange(0, self.shape[1]+_bins_size[1]/2, _bins_size[1])
        _y_edges   = np.arange(0, self.shape[0]+_bins_size[0]/2, _bins_size[0])
        # _x_centers = moving_average(_x_edges, n=2)
        # _y_centers = moving_average(_y_edges, n=2)

        # saving
        self._nbins     = _nbins
        self._bins_size = _bins_size
        self._reduced   = reduced
        # self.reduced._x_centers = _x_centers
        # self.reduced._y_centers = _y_centers
        self.reduced.x_edges   = _x_edges
        self.reduced.y_edges   = _y_edges

        return self.reduced

    def calculate_histogram(self, **kwargs):
        """Compute the histogram of data. Wrapper for `numpy.histogram()`_.

        If not specified, the following parameters are passed to `numpy.histogram()`_:

        Args:
            bins (int or tuple, optional): number of bins. If not specified, it will
                be set to 1000. If 1000 is too high (maximum value of the histogram
                is less than 5 % of the total integrated intensity) bins will be
                reduced 10 by 10 until the criteria is satisfied.

        Returns
            brixs.Spectrum

        .. _numpy.histogram(): https://numpy.org/doc/stable/reference/generated/numpy.histogram.html
        """
        if 'bins' not in kwargs:
            kwargs['bins'] = 1000
            hist, bin_edges = np.histogram(flatten(self._data), **kwargs)
            while max(hist) < self.shape[0]*self.shape[1]*0.05:
                kwargs['bins'] -= 10
                hist, bin_edges = np.histogram(flatten(self._data), **kwargs)
                if kwargs['bins'] < 50:
                    break
        else:
            hist, bin_edges = np.histogram(flatten(self._data), **kwargs)
        x = moving_average(bin_edges, 2)
        s = Spectrum(x=x, y=hist)
        return s

    def calculate_spectrum(self, axis=1):
        """Integrate data in one direction (sum columns or rows).

        Args:
            axis (int or string, optional): Axis along which elements are integrated.
                By default, data is integrated in the horizontal direction.

        Returns:
            :py:class:`Spectrum`.
        """
        axis = _axis_interpreter(axis)

        if axis == 0:
            return Spectrum(x=self.x_centers, y=np.sum(self._data, axis=0))
        elif axis == 1:
            return Spectrum(x=self.y_centers, y=np.sum(self._data, axis=1))

    def calculate_shifts(self, axis=0, mode='cc', limit_size=1000):
        """Calculate intensity misalignments via cross-correlation.

        Args:
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.
            mode (string, optional): method used to calculate the shifts.
                The current options are: 'cross-correlation' ('cc').
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images.
                Default is 1000. Set to False to bypass this limit.

        Returns:
            None
        """
        axis = _axis_interpreter(axis)

        # assert self.reduced is None, 'Image was not binned yet.\nPlease, use Image.binning()'

        # print('ff')
        # print(self.reduced)
        # mode = 'cross-correlation'
        peak = 0
        bkg_check = True

        # select axis
        if axis == 0:
            if limit_size:
                if len(self.x_centers) > limit_size:
                    raise ValueError(f'Number of columns is bigger than limit_size.\nImage is seems to be too big.\nAre you sure you want to calculate shifts for such a big image.\nIf so, either set limit_size to False or a higher value.\nNumber of columns: {len(self.x_centers)}\nlimit size: {limit_size}')
            ss = self.columns
            centers = self.x_centers
        elif axis == 1:
            if limit_size:
                if len(self.rows) > limit_size:
                    raise ValueError(f'Number of rows is bigger than limit_size.\nImage is seems to be too big.\nAre you sure you want to calculate shifts for such a big image.\nIf so, either set limit_size to False or a higher value.\nNumber of columns: {len(self.y_centers)}\nlimit size: {limit_size}')
            ss = self.rows
            centers = self.y_centers

        # peaks
        if mode == 'fitted peaks' or mode == 'peaks':
            ss.fit_peak()

        # calculate
        ss.calculate_shifts(mode=mode)
        if mode in cc:
            self._calculated_shifts = ss.calculated_shifts
            self.calculated_shifts.factor = ss.step
        else:
            self._calculated_shifts = ss.calculated_shifts
            self._calculated_shifts.y = [int(round(y)) for y in self._calculated_shifts.y]
        self._calculated_shifts.x = centers

    def set_shifts(self, value=None, p=None, f=None, axis=0, type_='absolute'):
        """Roll array of pixels along a given axis.

        Elements that roll beyond the edge are NOT re-introduced at the first.
        They are lost. This might change in the future.

        Args:
            value (int or tuple): The number of pixels by which the data are shifted.
                If a tuple, then it must be of the same size as the number of
                columns (``for axis=0``) or rows (``for axis=1``). If elements are
                not int, it will rounded to an integer value.
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.

        Returns:
            None
        """
        axis = _axis_interpreter(axis)

        # calculate value if necessary =========================================
        if value is not None:
            pass
        elif p is not None:
            if axis == 0:
                value = np.polyval(p, self.x_centers)
            elif axis == 1:
                value = np.polyval(p, self.y_centers)
        elif p is not None:
            if axis == 0:
                value = f(self.x_centers)
            elif axis == 1:
                value = f(self.y_centers)

        if type_ in relative:
            if axis == 0:
                value = self.shifts_v + value
            elif axis == 1:
                value = self.shifts_h + value

        # shift is absolute!!!!!!!!! Works ==============================
        if axis == 0:
            if isinstance(value, Iterable):
                assert len(value) == self.shape[1], f'Number of values must be the same as the number of columns ({self.shape[1]})'
                value = [round(k) for k in value]
            else:
                value = [int(round(value))]*self.shape[1]

            # undo last shift
            for i, v in enumerate(self._shifts_v):
                if v != 0:
                    temp = Spectrum(self._data[:, i])
                    temp.shift_roll = -v
                    self._data[:, i] = temp.y
                    self._shifts_v[i] = 0

            # apply shift
            for i, v in enumerate(value):
                if v != 0:
                    temp = Spectrum(self._data[:, i])
                    temp.shift_roll = v
                    self._data[:, i] = temp.y
                    self._shifts_v[i] = v
        elif axis == 1:
            if isinstance(value, Iterable):
                assert len(value) == self.shape[0], f'Number of values must be the same as the number of rows ({self.shape[0]})'
                value = [round(v) for v in value]
            else:
                value = [int(round(value))]*self.shape[0]

            # undo last shift
            for i, h in enumerate(self._shifts_h):
                if h != 0:
                    temp = Spectrum(self._data[i, :])
                    temp.shift_roll = -h
                    self._data[i, :] = temp.y
                    self._shifts_h[i] = 0

            # apply shift
            for i, h in enumerate(value):
                if h != 0:
                    temp = Spectrum(self._data[i, :])
                    temp.shift_roll = h
                    self._data[i, :] = temp.y
                    self._shifts_h[i] = h

        # if self.reduced is not None:
        #     self.binning(nbins=self.nbins)
        self._reduced = None

    def fix_curvature(self, deg=2, axis=0, mode='cc'):
        """Fix curvature.

        Args:
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical direction.

        Returns:
            None
        """
        axis = _axis_interpreter(axis)

        assert self.reduced is not None, 'Image was not binned yet.\nPlease, use Image.binning()'

        # calculate shifts
        self.reduced.floor()
        self.calculate_shifts(axis=axis, mode=mode)

        p, f = self.calculated_shifts.polyfit(deg=deg)
        self._p = p
        self._f = f

        self.set_shifts(p=p, axis=axis)


    def floor(self, x=0, y=0, n=30, nx=None, ny=None):
        """Set background intensity to zero.

        Args:
            x, y (int, optional): x and y position to sample background intensity.
            n, nx, ny (int, optional): size of the pixel window around x, y.

        Returns:
            None
        """
        assert x >= 0 and is_integer(x) and x<self.shape[1], f'x must be a positive integer smaller than {self.shape[1]}.'
        assert y >= 0 and is_integer(y) and y<self.shape[0], f'y must be a positive integer smaller than {self.shape[0]}.'

        # sorting n
        if nx is None: nx = n
        if ny is None: ny = n

        # check if x falls inside the image
        if x-nx/2 >= 0 and x+nx/2 < self.shape[1]:
            x_start = x-nx/2
            x_stop  = x+nx/2
        elif x-nx/2 < 0:
            x_start = 0
            x_stop = x+nx/2-(x-nx/2)
        elif x+n/2 >= self.shape[1]:
            x_start = x-nx/2-(x+nx/2-self.shape[1])
            x_stop  = self.shape[1]-1
        else:
            raise ValueError('Averaging range falls outside of the image. Please, change x or n.')

        # check if y falls inside the image
        if y-ny/2 >= 0 and y+ny/2 < self.shape[0]:
            y_start = y-ny/2
            y_stop  = y+ny/2
        elif y-nx/2 < 0:
            y_start = 0
            y_stop = y+ny/2-(y-ny/2)
        elif y+n/2 >= self.shape[0]:
            y_start = y-ny/2-(y+ny/2-self.shape[0])
            y_stop  = self.shape[0]-1
        else:
            raise ValueError('Averaging range falls outside of the image. Please, change y or n.')

        self._data -= np.mean(self._data[int(y_start):int(y_stop), int(x_start):int(x_stop)]).astype(self.data.dtype)
        self._vmin = min([min(x) for x in self.data])
        self._vmax = max([max(x) for x in self.data])

    def crop(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Crop Image.

        Args:
            x_start, x_stop, y_start, y_stop (int): pixel range. start is
                inclusive and stop is exclusive. Use None to to indicate the
                edge of the image.

        Returns:
            croped image
        """
        # check if None
        if x_start is None: x_start = 0
        if x_stop is None:  x_stop = self.shape[1]
        if y_start is None: y_start = 0
        if y_stop is None:  y_stop = self.shape[0]

        # verification
        assert x_start >= 0 and x_start<=self.shape[1], f'x_start must be a positive integer smaller than {self.shape[1]}.'
        assert is_integer(x_start), f'x_start must be an integer.'
        assert x_stop  >= 0 and x_stop<=self.shape[1],  f'x_stop must be a positive integer smaller than {self.shape[1]}.'
        assert is_integer(x_stop),  f'x_stop must be an integer.'
        assert x_stop > x_start, f'x_start must be smaller than x_stop.'

        assert y_start >= 0 and y_start<=self.shape[0], f'y_start must be a positive integer smaller than {self.shape[0]}.'
        assert is_integer(y_start), f'y_start must be an integer.'
        assert y_stop  >= 0 and y_stop<=self.shape[0],  f'y_stop must be a positive integer smaller than {self.shape[0]}.'
        assert is_integer(y_stop), f'y_stop must be an integer.'
        assert y_stop > y_start, f'y_start must be smaller than y_stop.'

        # crop
        final = Image(data=self.data[int(y_start):int(y_stop), int(x_start):int(x_stop)])
        return self._transfer_attributes(final)


class PhotonEvents(metaclass=_Meta):
    """Photon events object.

    Args:
        data (list or array): two (x, y) or three (x, y, intensity) list or array with
            photon events.
        filepath (string or path object, optional): filename or file handle.
            If the filename ends in .gz, the file is automatically saved in
            compressed gzip format. This is overwriten by data.
        shape (tuple, optional): Shape of data (maximum vertical size, maximum horizontal size).
        x, y, I (list or array): data columns.

    Attributes:
        data (2D array): This is where we store the Image.
        shape (tuple): Shape of data (vertical size, horizontal size).
        x, y, I (1D array): x, y, and intensity arrays.

        nbins (tuple): Number of bins (number of rows, number of columns).
        bins_size (tuple): Bins size (size of rows, size of columns).
        reduced (brixs.Image): Binned image.

        calculated_shifts (brixs.Spectrum): Calculated shifts.
        p (list): Polynomial coefficients of the fitted shift values.
        f (function): Function f(x) of the fitted shift values.

    Methods:
        save()
        load()
        plot()
        binning()
        calculate_histogram()
        calculate_spectrum()
        floor()
        calculate_shifts()
        set_shifts()
        fix_curvature()

    """

    _read_only = ['reduced', 'calculated_shifts']

    def __init__(self, *args, **kwargs):
        # argument parsing
        data, filepath, shape = self._sort_args(args, kwargs)

        # basic attr
        self._data  = None
        self._shape = [None, None]

        # binning attr
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None

        # shifts
        self._calculated_shifts = None
        self._shifts = None

        # set data
        if data is not None:
            self.data = copy.deepcopy(data)
        elif filepath is not None:
            self.load(filepath)
        if shape != (None, None):
            self.shape = shape

    @property
    def data(self):
        return copy.deepcopy(self._data)
    @data.setter
    def data(self, value):
        # basic attr
        if value.shape[1] == 3:
            self._data = np.array([event for event in np.array(value, dtype='float') if not event[2]<=0])
        elif value.shape[1] == 2:
            data = np.array(value, dtype='float')
            self._data = np.c_[data, np.ones(data.shape[0])]
        else:
            raise ValueError("Data must have 2 or 3 columns.")

        # set shape
        if self.shape[0] is None:
            self._shape[0] = max(self.data[:, 1])
        elif max(self.data[:, 1]) > self._shape[0]:
            self._shape[0] = max(self.data[:, 1])
        if self.shape[1] is None:
            self._shape[0] = max(self.data[:, 0])
        elif max(self.data[:, 0]) > self._shape[1]:
            self._shape[1] = max(self.data[:, 0])

        # binning attr
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None

        # shift attr
        self._calculated_shifts = None
        self._shifts = None
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shape(self):
        return copy.deepcopy(self._shape)
    @shape.setter
    def shape(self, value):
        if isinstance(value, Iterable):
            if len(value) == 2:
                if value[0] is None:
                    value[0] = max(self.data[:, 1])
                if value[1] is None:
                    value[1] = max(self.data[:, 0])
                if value[0] < 0 or value[1] < 0:
                    raise ValueError('Shape cannot be negative.')
                self._shape = value
            else:
                raise ValueError(f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}')
        else:
            raise ValueError(f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}')
    @shape.deleter
    def shape(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x(self):
        return copy.deepcopy(self._data[:, 0])
    @x.setter
    def x(self, value):
        raise AttributeError('Cannot set object. This might change in the future.')
    @x.deleter
    def x(self):
        raise AttributeError('Cannot delete object.')

    @property
    def y(self):
        return copy.deepcopy(self._data[:, 1])
    @y.setter
    def y(self, value):
        raise AttributeError('Cannot set object. This might change in the future.')
    @y.deleter
    def y(self):
        raise AttributeError('Cannot delete object.')

    @property
    def I(self):
        return copy.deepcopy(self._data[:, 2])
    @I.setter
    def I(self, value):
        raise AttributeError('Cannot set object. This might change in the future.')
    @I.deleter
    def I(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nbins(self):
        return self._nbins
    @nbins.setter
    def nbins(self, value):
        if type(value) != str:
            self.binning(nbins=value)
        # elif value == 'guess':
        #     guess_bins = self.guess_bins()
        #     self.bins = guess_bins
        else:
            raise ValueError("Not a valid option of `nbins`.\nValid options are: a number, a tuple, or a list.")
    @nbins.deleter
    def nbins(self):
        raise AttributeError('Cannot delete object.')

    @property
    def bins_size(self):
        return self._bins_size
    @bins_size.setter
    def bins_size(self, value):
        self.binning(bins_size=value)
    @bins_size.deleter
    def bins_size(self):
        raise AttributeError('Cannot delete object.')

    def _sort_args(self, args, kwargs):
        """checks initial arguments.

        Keyword arguments (kwargs) cannot be mixed with positional arguments.

        For positional arguments, if one data set is passed, it assumes it is
            `data`. If this one argument is of type string or Pathlib.Path, it
            assumes it is a filepath.

        Raises:
            AttributeError: if kwargs and args cannot be read.

        Returns:
            data, filepath, shape
        """
        # initial check
        if kwargs != {} and args != ():
            raise AttributeError('Cannot mix keyword arguments with positional arguments. Keyword arguents are `x`, `y`, `I`, `shape`, `data`, and `filepath`.')

        # initialization
        data     = None
        filepath = None
        x        = None
        y        = None
        I        = None
        shape    = (None, None)
        error = 'brixs.PhotonEvents() cannot figure out the data out of the arguments passed.\nMaybe use keyword arguments.\nValid arguments: x, y, I, shape, data, filepath.'

        # keyword arguments
        if 'data' in kwargs:
            data = kwargs['data']
        if 'x' in kwargs:
            x = kwargs['x']
        if 'y' in kwargs:
            y = kwargs['y']
        if 'I' in kwargs:
            I = kwargs['I']
        if 'shape' in kwargs:
            shape = kwargs['shape']
        if 'data' in kwargs:
            data = kwargs['data']
        if 'filepath' in kwargs:
            filepath = kwargs['filepath']

        if x is not None and y is not None and data is None:
            assert len(x) == len(y), f'x and y must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}.'
            data = np.c_[np.array(x), np.array(y)]
            if I is not None:
                assert len(x) == len(I), f'x, y and I must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}\nLength of I = {len(I)}.'
                data = np.c_[data, np.array(I)]


        # positional arguments
        if len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepath = Path(args[0])
            else:
                data = args[0]
        elif len(args) == 2:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepath = Path(args[0])
                shape = list(args[1])
                assert len(shape) == 2, f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}'
            else:
                if len(args[0]) == len(args[1]):
                    x = args[0]
                    y = args[1]
                    assert len(x) == len(y), f'x and y must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}.'
                    data = np.c_[np.array(x), np.array(y)]
                else:
                    data = args[0]
                    shape = list(args[1])
                    assert len(shape) == 2, f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}'
        elif len(args) == 3:
            if len(args[0]) == len(args[1]):
                x = args[0]
                y = args[1]
                assert len(x) == len(y), f'x and y must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}.'
                data = np.c_[np.array(x), np.array(y)]
                if len(args[0]) == len(args[2]):
                    I = args[2]
                    data = np.c_[data, np.array(I)]
                else:
                    shape = args[2]
                    assert len(shape) == 2, f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}'
        elif len(args) == 4:
            x = args[0]
            y = args[1]
            assert len(x) == len(y), f'x and y must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}.'
            data = np.c_[np.array(x), np.array(y)]
            I = args[2]
            assert len(x) == len(I), f'x, y and I must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}\nLength of I = {len(I)}.'
            data = np.c_[data, np.array(I)]
            shape = args[3]
            assert len(shape) == 2, f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}'
        elif len(args) > 4:
            raise AttributeError(error)

        return data, filepath, shape

    def __len__(self):
        if self._data is None:
            return 0
        else:
            return len(self._data[:, 0])

    def get_user_defined_attrs(self):
        """return attrs that are user defined."""
        default_attrs =  ['_data', '_shape', 
                           '_nbins', '_bins_size', 
                           '_reduced', '_shifts',
                           '_calculated_shifts']
        return [key for key in self.__dict__.keys() if key not in default_attrs]
    
    def save(self, filepath=None, only_data=False,  check_overwrite=True, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format. If None, filepath can be inferred from a
                user defined attr 'im.filepath'. Default is None.  Last used 
                filepath is saved to pe.filepath.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user wants to overwrite file.
            **kwargs: kwargs are passed to ``np.savetxt()`` that saves the data.

        If not specified, the following parameters are passed to `numpy.savetxt()`_:

        Args:
            fmt (string, or list, optional): A single format (like ``%10.5f``), or a
                sequence of formats. See numpy's documentation for more information.
                If not specified, best fmt is calculated based on the
                number of decimal places of the data. Specifing fmt makes the code
                runs a little faster (not much tough, but it might make a difference
                if saving a lot of files).
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (", ").
            newline (str, optional): String or character separating lines.
                Default is ``\n``.
            header (bool, optional): String that will be written at the beginning of the file.
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        # filepath
        if filepath is None:
            try: 
                filepath = Path(self.filepath)
            except AttributeError:
                raise TypeError("Missing 1 required argument: 'filepath'")

        # check overwrite
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')
                
        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([n_decimal_places(x) for x in flatten(self._data)])
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # save
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
            np.savetxt(Path(filepath), self._data, **kwargs)
        else:
            if 'header' not in kwargs:
                kwargs['header'] = ''
            else:
                kwargs['header'] += '\n'
            dict = get_attributes(self)
            kwargs['header'] += '==== brixs PhotonEvents ===='  + '\n'
            for n in dict:
                if n not in ['_data', '_x', '_y', '_I', '_reduced', '_calculated_shifts', '_f', '_shifts']:
                    if isinstance(dict[n], Iterable):
                        kwargs['header'] += f'{n}: {list(dict[n])}'  + '\n'
                    else:
                        kwargs['header'] += f'{n}: {dict[n]}'  + '\n'
            np.savetxt(Path(filepath), self._data, **kwargs)

        # save filepath
        self.filepath = filepath

    def load(self, filepath, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first 
                decompressed. Last used filepath is saved to pe.filepath.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``# ``. Attributes picked up
                from the header will be loaded too.

        Returns:
            None

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # read data
        data = np.genfromtxt(Path(filepath), **kwargs)
        self.data = data

        # read header
        header = load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
        attr_start = 0
        binning_flag = False  # check if data needs binning

        # find where attributes listing starts
        for i, line in enumerate(header):
            if '==== brixs PhotonEvents ====' in line:
                attr_start = i
                break
            attr_start = -1

        # read attributes
        if attr_start != -1:
            for i, line in enumerate(header[attr_start+1:-1]):

                # extract name and value
                name = line[1:-1].split(':')[0].strip()
                value = eval(line[1:-1].split(':')[1].strip())


                ### DEALING WITH ATTRIBUTES THAT NEED TO RUN SOMETHING ###
                if name in ['_nbins']:
                    if value[0] > 0 and value[1] > 0:
                        binning_flag = True
                elif name in ['_bins_size']:
                    if value[0] > 0 and value[1] > 0:
                        binning_flag = True
                ### DEALING WITH OTHER ATTRIBUTES ###
                elif name not in []:  # except these attrs
                    try:
                        setattr(self, name, value)
                    except Exception as e:
                        print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')

        ### RUN ###
        if binning_flag:
            try:
                self.binning(nbins=self.nbins, bins_size=self.bins_size)
            except:
                pass

        # save filepath
        self.filepath = str(filepath)

    def plot(self, ax=None, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.scatter()`_.

        The limits of the plot is also set to the size of the detector (given in
            the attribute ``shape``.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.

        If not specified, the following parameters are passed to `matplotlib.pyplot.scatter()`_:

        Args:
            s: The marker size in points**2. Default is 0.1.

        Returns:
            `matplotlib.image.AxesImage`_

        .. _matplotlib.pyplot.scatter(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.scatter.html
        .. _matplotlib.image.AxesImage: https://matplotlib.org/3.5.0/api/image_api.html#matplotlib.image.AxesImage
        """
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figure()
                if settings.FIGURE_POSITION is not None:
                    try:
                        set_window_position(settings.FIGURE_POSITION)
                    except:
                        pass
            elif plt.get_fignums() == [] and settings.FIGURE_POSITION is not None:
                try:
                    set_window_position(settings.FIGURE_POSITION)
                except:
                    pass

        # kwargs
        if 's' not in kwargs:
            kwargs['s'] = 0.1

        # plot
        pos = ax.scatter(self.data[:, 0], self.data[:, 1], **kwargs)

        # set limits
        ax.xlim(0, self.shape[1])
        ax.ylim(0, self.shape[0])

        return pos

    def binning(self, *args, **kwargs):
        """Compute the 2D histogram of the data (binning of the data).

        Args:
            nbins (int or tuple, optional): number of bins. If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the number of rows and number of
                columns, respectively. ``bins_size`` is calculated.
            bins_size (int or tuple, optional): size of the bins. This overwrites
                the argument ``bins``. If one value is given,
                it is used for both x and y directions. If two values are given,
                they are used separetely for rows and columns size, respectively.
                ``nbins`` is calculated.

        Example:
            .. code-block:: python

                nbins = 10       # (10 rows, 10 columns)
                nbins = (10)     # (10 rows, 10 columns)
                nbins = (10, 5)  # (10 rows, 5 columns)

        Returns:
            None
        """
        kwargs['shape']        = self.shape
        kwargs['factor_check'] = False
        _nbins, _bins_size = _bins_interpreter(*args, **kwargs)

        temp, _x_edges, _y_edges = np.histogram2d(self.data[:, 0],
                                                            self.data[:, 1],
                                                            bins=_nbins[::-1],
                                                             weights=self.data[:, 2],
                                                             range=((0, self.shape[1]), (0, self.shape[0]))
                                                            )
        self._reduced   = Image(temp.transpose())
        self.reduced._x_centers = moving_average(_x_edges, n=2)
        self.reduced._y_centers = moving_average(_y_edges, n=2)
        self.reduced._x_edges = _x_edges
        self.reduced._y_edges = _y_edges
        self._nbins     = _nbins
        self._bins_size = _bins_size
        return

    def calculate_spectrum(self, nbins=None, bins_size=None, axis=1, xaxis='bins'):
        """Integrate data in one direction (sum columns or rows).

        Args:
            nbins (int, optional): number of bins. nbins overwrites bins_size. The
                nbins of the opposite axis is automatically set to 1.
            bins_size (int or tuple, optional): size of the bins.
            axis (int or string, optional): Axis along which elements are integrated.
                By default, data is integrated in the horizontal direction.

        Returns:
            :py:class:`Spectrum`.
        """
        axis = _axis_interpreter(axis)

        temp = PhotonEvents(data=self.data, shape=self.shape)
        temp.binning(nbins=nbins, bins_size=bins_size)
        s = temp.reduced.calculate_spectrum(axis=axis)
        return s

    def calculate_shifts(self, axis=0):
        """For now, the only mode tested is cc"""
        axis = _axis_interpreter(axis)
        mode = 'cc'
        assert self.reduced is not None, 'Image was not binned yet.\nPlease, use PhotonEvents.binning()'

        self.reduced.calculate_shifts(axis=axis, mode=mode)

        if mode in cc:
            if axis == 0:
                self.reduced.calculated_shifts._y = self.reduced.calculated_shifts.y#*self.bins_size[0]
            elif axis == 1:
                self.reduced.calculated_shifts._y = self.reduced.calculated_shifts.y#*self.bins_size[1]
        self._calculated_shifts = self.reduced.calculated_shifts

    def transform(self, f):
        """Changes the values of x, y based on a function.

        Args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        Example:
            f = lambda x, y: (x, y**2)

        Returns:
            None
        """
        self._data[:, 0], self._data[:, 1] = f(self.data[:, 0], self.data[:, 1])

    def set_shifts(self, value=None, p=None, f=None, axis=0):
        """assuming calculation mode = 'cc'
        value
        p
        f
        axis=0 or 1
        """
        axis = _axis_interpreter(axis)

        if value is not None:
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)
            assert len(value) == len(self), f'value must have the same length as the data\nLength of the data = {len(self)}.'

            if axis == 0:
                self._data[:, 1] = self.data[:, 1] + value
                self._shifts = Spectrum(x=self.data[:, 0], y=value)
            elif axis == 1:
                self._data[:, 0] = self.data[:, 0] + value
                self._shifts = Spectrum(x=self.data[:, 1], y=value)
        elif p is not None:
            if axis == 0:
                self._data[:, 1] = self.data[:, 1] + np.polyval(p, self.data[:, 0])
            elif axis == 1:
                self._data[:, 0] = self.data[:, 0] + np.polyval(p, self.data[:, 1])
            self._p = p
        elif f is not None:
            if axis == 0:
                self._data[:, 1] = self.data[:, 1] + f(self.data[:, 0])
            elif axis == 1:
                self._data[:, 0] = self.data[:, 0] + f(self.data[:, 1])
            self._f = f

    def fix_curvature(self, deg=2, axis=0):
        """mode = 'cc'
        
        return popt
        """
        axis = _axis_interpreter(axis)

        assert self.reduced is not None, 'Image was not binned yet.\nPlease, use Image.binning()'

        # self.reduced.floor()
        self.calculate_shifts(axis=axis)

        popt, model, r2 = self.calculated_shifts.polyfit(deg=deg)

        # set shifts
        self.set_shifts(p=popt)

        return popt, model, r2

    #
    # def guess_bins(self, bins_initial=(9, 100), bins_step=(1, 50), max_x_iter=3, max_iter=1000, mode='cross-correlation', ranges=None):
    #     """
    #     Args:
    #         max_x_iter (int, optional): every iteration the algorithm calculates
    #             the offsets for different x_bins. Sometimes, subsequant x_bins
    #             will give the same result (i.e., same number of equal offsets).
    #             The algorithm will keep trying max_x_iter different x_bins until
    #             finaly change to a different y_bins.
    #         max_iter (int, optional): maximum number of different binnigs to try.
    #     """
    #     # offset parameters
    #     ref = int(bins_initial[0]/2)
    #
    #     # inital iteration
    #     # temp = PhotonEvents(self.data, x_max=self.x_max, y_max=self.y_max)
    #     # temp.bins = bins_initial
    #     # temp.calculate_offsets(mode=mode, ranges=ranges, ref=ref)
    #     diff = [0, 0]
    #     diff_sum = sum(diff)
    #
    #     # bins_initial[0] -= 1
    #     x_counter = 0
    #     total_counter = 0
    #     bins = (bins_initial[0], bins_initial[1])
    #     while 0 in diff and total_counter < max_iter:
    #         temp = PhotonEvents(self.data, x_max=self.x_max, y_max=self.y_max)
    #         temp.bins = bins
    #         temp.calculate_offsets(mode=mode, ranges=ranges, ref=ref)
    #         diff = np.diff(temp.offsets)
    #
    #         if x_counter < max_x_iter-1:
    #             if diff_sum == sum(diff):
    #                 bins = (bins[0]+bins_step[0], bins[1])
    #                 x_counter += 1
    #         else:
    #             bins = (bins_initial[0], bins[1]+bins_step[1])
    #             x_counter = 0
    #
    #         diff_sum = sum(diff)
    #         total_counter += 1
    #     if total_counter >= max_iter:
    #         raise RuntimeError('cannot find solution.')
    #     else:
    #         return temp.bins
    #
    # def best_bins(self, y_bins=None, x_bins=None, deg=2, spectrum_bins=6000, **kwargs):
    #
    #     if y_bins is None:
    #         y_bins = [100, 150, 300, 400, 600, 1000, 1500, 2500, 5000, 6000, 8000, 12000, 15000]
    #
    #     if x_bins is None:
    #         x_bins = [3, 5, 7, 9, 12, 15, 20, 30, 60]
    #
    #     # bins_old = self.bins
    #
    #     fwhm = {y_bin: {x_bin:None for x_bin in x_bins} for y_bin in y_bins}
    #     best_bins = (0, 0)
    #     best = 10000
    #     print(f'Starting iteration...')
    #     for i, y_bin in enumerate(y_bins):
    #         print(f'({i}/{len(y_bins)-1}) Trying y_bin = {y_bin}.')
    #         for x_bin in x_bins:
    #             temp = PhotonEvents(self.data, x_max=self.x_max, y_max=self.y_max)
    #             temp.bins = (x_bin, y_bin)
    #             # print(self.bins)
    #             temp.calculate_offsets(ref=int(temp.bins[0]/2))
    #             temp.fit_offsets(deg=deg)
    #             temp.offsets_correction()
    #             s = temp.calculate_spectrum(bins=spectrum_bins)
    #             try:
    #                 s.guess_elastic_peak(**kwargs)
    #                 fwhm[y_bin][x_bin] = s.elastic_fwhm
    #                 if s.elastic_fwhm < best:
    #                     best_bins = (x_bin, y_bin)
    #                     best = s.elastic_fwhm
    #             except RuntimeError:
    #                 pass
    #
    #
    #     print(f'Done. Best bins is {best_bins}.')
    #
    #     return best_bins, best, fwhm
    #
    # def best_bins_for_spectrum(self, y_bins=None, **kwargs):
    #
    #     if y_bins is None:
    #         y_bins = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000]
    #
    #     fwhm = {y_bin: None for y_bin in y_bins}
    #     best = 10000
    #     print(f'Starting iteration...')
    #     for i, y_bin in enumerate(y_bins):
    #         print(f'({i}/{len(y_bins)-1}) Trying y_bin = {y_bin}')
    #         s = self.calculate_spectrum(bins=y_bin, mode='length')
    #         try:
    #             s.guess_elastic_peak(**kwargs)
    #             fwhm[y_bin] = s.elastic_fwhm
    #             if s.elastic_fwhm < best:
    #                 best_bins = (1, y_bin)
    #                 best = s.elastic_fwhm
    #         except RuntimeError:
    #             fwhm[y_bin] = -1
    #
    #
    #     return best_bins, best, fwhm


class Spectrum(metaclass=_Meta):
    """Creates a ``spectrum`` class type object to deal with (x, y) data types.

    The hierarchy for Keyword arguments is: 1) data, 2) y (and x), and finaly
    3) filepath. For example, if `data` and `filepath` are passed as arguments,
    `filepath` is ignored.

    Keyword arguments (kwargs) cannot be mixed with positional arguments.

    For positional arguments, if one data set is passed, it assumes it is
    `data`. If this one argument is of type string or Pathlib.Path, it
    assumes it is a filepath. If two data sets are passed, it will
    assume one is the x coordinates and the next one is the y coordinates.

    Args:
        data (list or array, optional): two column list (or array).
        x (list or array, optional): x values (1D list/array). Overwrites `data`.
        y (list or array, optional): y values (1D list/array). Overwrites `data`.

    Attributes:
        data (array): three column array (x, y, and intensity)
        x (array): vector with x-coordinate values
        y (array): vector with y-coordinate values
        step (number, read only): step size for x-coordinates.
        shift (number): shift value (value will be added to x-coordinates).
        calib (number): calibration value (x-coordinates will be multiplied by this value).
        offset( number): offset value (value will be added to y-coordinates).
        factor (number): multiplicative factor (y-coordinates will be multiplied by this value).
        peaks (dict): each entry represents a peak.



    Note, peaks can be copyed around without need of deepcopy.
    """
    # read only and non-removable arguments
    # _read_only = ['step', 'monotonicity',  'R2', 'fit', 'residue', 'guess', 'pcov']
    _read_only = []
    _non_removable = []

    def __init__(self, *args, **kwargs):
        # basic
        self._x    = None
        self._y    = None
        self._xlabel   = ''
        self._ylabel   = ''
        self._filepath = ''

        # modifiers
        self._reset_modifiers()

        # check
        self._step         = None
        self._monotonicity = None

        # # fit
        # self._fit     = None
        # self._residue = None
        # self._guess   = None
        # self._R2      = None
        # self._pcov    = None

        # peaks
        self._peaks = Peaks()

        # sorting data
        data, x, y, filepath = self._sort_args(args, kwargs)
        if data is not None:
            self.data = data
        elif filepath is not None:
            self.load(filepath)
        elif y is not None:
            if x is None:
                x = np.arange(0, len(y))
                self.data = np.vstack((x, y)).transpose()
            else:
                if len(x) == len(y):
                    self.data = np.vstack((x, y)).transpose()
                else:
                    raise ValueError('x and y data are not the same length.')
        else:
            self._x = None
            self._y = None

    # basic
    def __len__(self):
        if self.x is None:
            return 0
        else:
            return len(self.x)

    def __add__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot operate on spectra. x axis is different.\nMaybe try interpolating the x axis.')
            final = Spectrum(x=self.x, y=self.y + object.y)
            return self._transfer_attributes(final)
        elif isinstance(object, (np.floating, float, int)):
            final = Spectrum(x=self.x, y=self.y + object)
            return self._transfer_attributes(final)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')

    def __sub__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot operate on spectra. x axis is different.\nMaybe try interpolating the x axis.')
            final = Spectrum(x=self.x, y=self.y - object.y)
            return self._transfer_attributes(final)
        elif isinstance(object, (np.floating, float, int)):
            final = Spectrum(x=self.x, y=self.y - object)
            return self._transfer_attributes(final)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')

    def __mul__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot operate on spectra. x axis is different.\nMaybe try interpolating the x axis.')
            final = Spectrum(x=self.x, y=self.y * object.y)
            return self._transfer_attributes(final)
        elif isinstance(object, (np.floating, float, int)):
            final = Spectrum(x=self.x, y=self.y * object)
            return self._transfer_attributes(final)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')

    def __div__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot divide spectra. x axis is different.\nMaybe try interpolating the x axis.')
            if 0 in object.y:
                raise ZeroDivisionError(f'y axis contain zeros. Cannot divide by zero.')
            else:
                final = Spectrum(x=self.x, y=self.y/object.y)
                return self._transfer_attributes(final)
        elif isinstance(object, (np.floating, float, int)):
            if object == 0:
                raise ZeroDivisionError(f'Cannot divide by zero.')
            else:
                final = Spectrum(x=self.x, y=self.y / object)
                return self._transfer_attributes(final)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')

    def __truediv__(self, object):
        return self.__div__(object)

    # support
    def _sort_args(self, args, kwargs):
        """checks initial arguments.

        Keyword arguments (kwargs) cannot be mixed with positional arguments.

        For positional arguments, if one data set is passed, it assumes it is
            `data`. If this one argument is of type string or Pathlib.Path, it
            assumes it is a filepath. If two data sets are passed, it will
            assume one is the x coordinates and the next one is the y coordinates.

        Raises:
            AttributeError: if kwargs and args cannot be read.

        Returns:
            data, x, y, filepath
        """
        # initial check
        if kwargs != {} and args != ():
            raise AttributeError('cannot mix key word arguments with positional arguments. Key word arguents are `x`, `y`, `data`, and `filepath`.')
        if any([item not in ['data', 'x', 'y', 'filepath'] for item in kwargs.keys()]):
            raise AttributeError(f'invalid attributes.\nValid atributes are `data`, `x`, `y`, and `filepath`\nInput attributes: {kwargs.keys()}')

        # initialization
        data = None
        x    = None
        y    = None
        filepath = None

        # keyword arguments
        if 'data' in kwargs:
            data = kwargs['data']
        if 'y' in kwargs:
            if 'x' in kwargs:
                x = kwargs['x']
            else:
                x = None
            y = kwargs['y']
        elif 'x' in kwargs:
            raise AttributeError('cannot seem to find y data.')
        if 'filepath' in kwargs:
            filepath = kwargs['filepath']

        # positional arguments
        if len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepath = args[0]
            elif isinstance(args[0], Iterable):
                data = args[0]
        elif len(args) == 2:
            x = args[0]
            y = args[1]
        elif len(args) > 2:
            raise AttributeError('brixs.Spectrum() cannot figure out the data out of the arguments passed. Maybe use key word arguments (x, y, data, filepath).')

        return data, x, y, filepath

    def _reset_modifiers(self):
        self._factor       = 1
        self._offset       = 0
        self._calib        = 1
        self._shift        = 0
        self._shift_roll   = 0
        self._shift_interp = 0

    def _transfer_attributes(self, object):
        """Transfer user defined attributes to output objects."""
        
        # copy user defined attributes
        for attr in self.__dict__:
            if attr.startswith('_') == False:
                object.__setattr__(attr, self.__dict__[attr])

        return object

    def get_user_defined_attrs(self):
        """return attrs that are user defined."""
        default_attrs =  ['_x', '_y', '_xlabel', '_ylabel', '_filepath', '_factor', '_offset', '_calib', '_shift', '_shift_roll', 
'_shift_interp', '_step', '_monotonicity', '_peaks']
        return [key for key in self.__dict__.keys() if key not in default_attrs]


    # attributes and properties
    @property
    def data(self):
        return np.vstack((self.x, self.y)).transpose()
    @data.setter
    def data(self, value):
        try:
            value = np.array(value, dtype='float')
            if value is None:
                raise ValueError('No data to load.')
            elif value.shape[1] == 1:   # one column data (list)
                self._x = np.arange(0, len(value), dtype='float')
                self._y = np.array(value, dtype='float')
            elif value.shape[1] != 2:
                value = value.transpose()
                if value.shape[1] != 2:
                    raise ValueError('Data must have two columns (x, y).')
                else:
                    self._x = value[:, 0]
                    self._y = value[:, 1]
            else:
                self._x = value[:, 0]
                self._y = value[:, 1]
        except IndexError:  # one column data (list)
            self._x = np.arange(0, len(value), dtype='float')
            self._y = np.array(value, dtype='float')
        # check
        self.step         = None
        self.monotonicity = None
        self._reset_modifiers()
        #special
        self.peaks.clear()
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, value):
        if self.y is not None:
            assert len(value) == len(self.y), f'Length of x array (len={len(value)}) you are trying to set is not compatible with current length of the y array (len={len(self.y)}).'
        self._x = np.array(value, dtype='float')
        # check
        self.step         = None
        self.monotonicity = None
        self._reset_modifiers()
        #special
        self.peaks.clear()
    @x.deleter
    def x(self):
        raise AttributeError('Cannot delete object.')

    @property
    def y(self):
        return self._y
    @y.setter
    def y(self, value):
        if self.x is not None:
            assert len(value) == len(self.x), f'Length of y array (len={len(value)}) you are trying to set is not compatible with current length of the x array (len={len(self.x)}).'
        else:
            self._x = np.arange(0, len(value))
            self._y = np.array(value, dtype='float')
        # check
        # self.step         = None
        # self.monotonicity = None
        # self._reset_modifiers()
        #special
        self.peaks.clear()
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
    def shifts(self):
        return {'hard': self.shift, 'interp': self.shift_interp, 'roll':self.shift_roll}
    @shifts.setter
    def shifts(self, value):
        raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shift(value, mode).")
    @shifts.deleter
    def shifts(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_roll(self):
        return self._shift_roll
    @shift_roll.setter
    def shift_roll(self, value):
        # raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shifts(value, mode='roll').")
        self.set_shift(value, mode='roll')
    @shift_roll.deleter
    def shift_roll(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_interp(self):
        return self._shift_interp
    @shift_interp.setter
    def shift_interp(self, value):
        # raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shifts(value, mode='interp').")
        self.set_shift(value, mode='interp')
    @shift_interp.deleter
    def shift_interp(self):
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

    @property
    def peaks(self):
        return self._peaks
    @peaks.setter
    def peaks(self, value):
        raise AttributeError('Attribute is "read only".')
    @peaks.deleter
    def peaks(self):
        raise AttributeError('Cannot delete object.')

    @property
    def area(self):
        return self.calculate_area()
    @area.setter
    def area(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @area.deleter
    def area(self):
        raise AttributeError('Cannot delete object.')

    @property
    def xlabel(self):
        return self._xlabel
    @xlabel.setter
    def xlabel(self, value):
        if value is None:
            value = ''
        elif isinstance(value, str):
            forbidden = [':', '#', '\n', '\r']
            if np.sum([x in value for x in forbidden]):
                raise ValueError(f'xlabel cannot contain the follwing characters: {forbidden}')
            self._xlabel = str(value)
        else:
            raise TypeError('Invalid type for xlabel\nxlabel must be a string.')
    @xlabel.deleter
    def xlabel(self):
        raise AttributeError('Cannot delete object.')

    @property
    def ylabel(self):
        return self._ylabel
    @ylabel.setter
    def ylabel(self, value):
        if value is None:
            value = ''
        elif isinstance(value, str):
            forbidden = [':', '#', '\n', '\r']
            if np.sum([x in value for x in forbidden]):
                raise ValueError(f'xlabel cannot contain the follwing characters: {forbidden}')
            self._ylabel = str(value)
        else:
            raise TypeError('Invalid type for ylabel\nylabel must be a string.')
    @ylabel.deleter
    def ylabel(self):
        raise AttributeError('Cannot delete object.')   
    
    @property
    def filepath(self):
        return self._filepath
    @filepath.setter
    def filepath(self, value):
        if value is None:
            value = ''
        elif isinstance(value, str) or isinstance(value, Path):
            self._filepath = value
        else:
            raise TypeError(r'Invalid type ' + str(type(value)) + 'for filepath\nFilepath can only be str or pathlib.Path type.')
    @filepath.deleter
    def filepath(self):
        raise AttributeError('Cannot delete object.')   

    @property
    def step(self):
        return self._step
    @step.setter
    def step(self, value):
        # raise AttributeError('Attribute is "read only". Cannot set attribute.')
        if value is None:
            self._step       = None
            self._shift_roll = 0

            # special
            self.peaks.step  = None
        else:
            raise AttributeError('Attribute is "read only". Cannot set attribute unless set to None.')
    @step.deleter
    def step(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def monotonicity(self):
        return self._monotonicity
    @monotonicity.setter
    def monotonicity(self, value):
        # raise AttributeError('Attribute is "read only". Cannot set attribute.')
        if value is None:
            self._monotonicity = None
            # self._shift_roll = 0
            # special
            # self.peaks.monotonicity  = None
        else:
            raise AttributeError('Attribute is "read only". Cannot set attribute unless set to None.')
    @monotonicity.deleter
    def monotonicity(self):
        raise AttributeError('Cannot delete object.')
    

    # save and load
    def _header(self, verbose):
        """Gather attrs to be saved to a file."""
        header = ''
        temp = get_attributes(self)
        for name in temp:
            if name not in ['_x', '_y', '_peaks', '_xlabel', '_ylabel', '_filepath']:
                if temp[name] is None:
                    header += f'{name}: None'  + '\n'
                elif isinstance(temp[name], str):
                    temp2 = str(temp[name]).replace('\n','\\n')
                    header += f'{name}: \"{temp2}\"'  + '\n'
                elif isinstance(temp[name], Iterable):
                    header += f'{name}: {list(temp[name])}'  + '\n'
                elif isinstance(temp[name], dict):
                    if verbose:
                        type_ = str(type(temp[name]))
                        print(r'Warning: Cannot save attr of type: ' + type_ + r'.\attr: '+ name + r'.\nTo turn off this warning, set verbose to False.')
                elif isinstance(temp[name], numbers.Number):
                    tosave = str(temp[name])
                    if tosave[-1] == '\n':
                        tosave = tosave[:-1]
                    header += f'{name}: {tosave}'  + '\n'
                else:
                    temp2 = str(temp[name]).replace('\n','\\n')
                    header += f'{name}: \"{temp2}\"'  + '\n'
        return header[:-1]

    def save(self, filepath=None, only_data=False, check_overwrite=True, verbose=True, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        User defined attr are saved in the header (except for dictionaries).

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format. If None, filepath can be inferred from a
                user defined attr 's.filepath'. Default is None. Last used 
                dirpath is saved to an attr ss.filepath.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user wants to overwrite file.
            **kwargs: kwargs are passed to ``np.savetxt()`` that saves the data.

        If not specified, the following parameters are passed to `numpy.savetxt()`_:

        Args:
            fmt (string, or list, optional): A single format (like ``%10.5f``), or a
                sequence of formats. See numpy's documentation for more information.
                If not specified, best fmt is calculated based on the
                number of decimal places of the data. Specifing fmt makes the code
                runs a little faster (not much tough, but it might make a difference
                if saving a lot of files).
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (", ").
            newline (str, optional): String or character separating lines.
                Default is ``\n``.
            header (bool, optional): String that will be written at the beginning of the file.
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        # filepath
        if filepath is None:
            if filepath == '':
                raise TypeError("Missing 1 required argument: 'filepath'")
            else:
                filepath = self.filepath               
        filepath = Path(filepath)
        
        # check overwrite
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([n_decimal_places(x) for x in flatten(self.data)])
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # save
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
        else:
            if 'header' not in kwargs:
                kwargs['header'] = self._header(verbose=verbose)
            else:
                if kwargs['header'] == '':
                    kwargs['header'] = self._header(verbose=verbose)
                elif kwargs['header'][-1] != '\n':
                    kwargs['header'] += '\n'
            
            if self.xlabel != '' and self.ylabel != '':
                kwargs['header'] += '\n' + self.xlabel + kwargs['delimiter'] + self.ylabel
        np.savetxt(Path(filepath), self.data, **kwargs)

        # save filepath
        self.filepath = str(filepath)

    def load(self, filepath, only_data=False, verbose=True, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        This a very simple loading function that works well with two column text
        files. If file has more columns than two columns, the first two columns
        will be loades. Use `usecols` to select columns, for example:
        
        usecols = (1, 4) 

        If file was saved by br.Spectrum.save(), then the metadata (comments) can be 
        recovered. If not, only_data must be set to True.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first 
                decompressed. Last used filepath is saved to an attr s.filepath.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is loaded.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``# ``. Attributes picked up
                from the header will be loaded too.

        Returns:
            None

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '#'
        if 'usecols' not in kwargs:
            kwargs['usecols'] = (0, 1)

        # read data
        data = np.genfromtxt(Path(filepath), **kwargs)

        # check data
        x = data[:, 0]
        y = data[:, 1]
        assert len(x) == len(y), f'Length of x array (len={len(x)}) is not compatible with y array (len={len(y)}).'
        self._x = x
        self._y = y

        # read header
        if only_data is False:
            header = load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            if header:
                for line in header:
                    if ':' not in line:
                        temp = line.split(kwargs['delimiter'])
                        self.xlabel = temp[kwargs['usecols'][0]].replace(kwargs['comments'], '').replace('\n', '').replace('\r', '').strip()
                        self.ylabel = temp[kwargs['usecols'][1]].replace(kwargs['comments'], '').replace('\n', '').replace('\r', '').strip()
                    else:
                        # extract name and value
                        name = line[1:-1].split(':')[0].strip()
                        value = eval(str(':'.join(line[1:-1].split(':')[1:])).strip())
                        try:
                            setattr(self, name, value)
                        except Exception as e:
                            if verbose:
                                print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')
        # save filepath
        self.filepath = str(filepath)

    # check
    def check_step(self, max_error=None):
        """Checks vector uniformity of the x-coordinates.

            If the step between two data points is the same through out the
            x vector, it sets :py:attr:`step` with the value of the step size.

            Args:
                max_error (number, optional): percentage (of the x step) value 
                of the max error.
                
                (max(steps) - min(steps))/np.mean(steps) * 100 < max_error

            Returns:
                None

            Raises:
                ValueError: If x-ccordinates are not uniform.
        """
        if max_error is None:
            max_error = settings.MAX_ERROR_STEP_X

        # if data is not monotonic, than it is not uniform
        if self.monotonicity is None:
            try:
                self.check_monotonicity()
            except ValueError:
                raise ValueError(f"Step in the x-coordinate seems not to be uniform. In fact, it is not even monotonic. Use Spectrum.fix_monotonicity().")

        # check step uniformity
        d = np.diff(self.x)
        if abs((max(d) - min(d))*100/np.mean(np.diff(self.x))) > max_error:
            self._step = None
            raise ValueError(f"Step in the x-coordinate seems not to be uniform.")

        # set step
        self._step = np.mean(d)

        # special
        self.peaks._step = self._step

    def check_monotonicity(self):
        """Sets monotonicity attribute to 'increasing' or 'decreasing'.

        Raises:
            ValueError if data is not monotonic.

        Returns:
            None
        """
        if np.all(np.diff(self.x) > 0) == True:
            self._monotonicity = 'increasing'
        elif np.all(np.diff(self.x) < 0) == True:
            self._monotonicity = 'decreasing'
        else:
            self._monotonicity = None
            raise ValueError('x array is not monotonic. Use Spectrum.fix_monotonicity()')

    def fix_monotonicity(self, mode='increasing'):
        """Rearrange (x, y) such as x array is monotonically increasing or decreasing.

        Args:
            mode (str, optional): increasing or decreasing.

        Returns:
            None
        """
        # check mode
        if mode not in increasing and mode not in decreasing:
            raise ValueError('mode should be "decreasing" or "increasing".')
        
        # turn array into monotonic
        if self.monotonicity is None:
            try:
                self.check_monotonicity()
            except ValueError:
                unqa, ID, counts = np.unique(self.x, return_inverse=True, return_counts=True)
                self.data        = np.column_stack(( unqa , np.bincount(ID,self.y)/counts ))
                # check
                self.check_monotonicity()
            
        # make it decreasing or increasing
        if self.monotonicity != mode:
            self._x = self.x[::-1]
            self._y = self.y[::-1]
        self.check_monotonicity()

        # check
        self.step         = None
        # self.monotonicity = None
        # self._reset_modifiers()
        #special
        # self.peaks.clear()

    # modifiers
    def set_calib(self, value, type_='absolute'):
        """Set calibration value.

        Args:
            value (number): calibration value (x-coordinates will be multiplied
                by this value).
            type_ (string, optional): either 'absolute' (default) or 'relative'.

        Returns:
            None
        """
        # calib cannot be zero
        if value == 0:
            raise ValueError('cannot set calib = 0.0')

        # is relative?
        if type_ in relative:
            value = self.calib * value

        # apply calibration
        if self.calib != value:
            if self.calib != 1:
                self._x = self.x*self.calib**-1
            if value != 1:
                self._x = self.x*value
            self._calib = value

        # check
        self.step         = None
        self.monotonicity = None
        # self._reset_modifiers()
        #special
        # self.peaks.clear()
        self.peaks.set_calib(value=value, type_=type_)

    def set_shift(self, value, mode, type_='absolute'):
        """Set shift value.

        Args:
            value (float or int): shift value (value will be added to x-coordinates).
            mode (string, optional): If ``mode='x'`` or ``mode='hard'``, y is fully preserved
                while x is shifted. If ``mode='y'``, ``'interp'``, or ``'soft'``, x is preserved
                while y is interpolated with a shift. If ``mode='roll'`` (or rotate or r), x is also preserved
                and y elements are rolled along the array (``shift`` value must be an integer).
                The form of y-coordinates is fully preserved, but the edge of y-coordinates are lost.
            type_ (string, optional): either 'absolute' (default) or 'relative'.

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
        # ROLL        
        if mode in roll:
            # x axis must be uniform if mode is roll
            if self.step is None:
                try:
                    self.check_step()
                except ValueError:
                    raise ValueError(f'Cannot shift data using mode = {mode}, because x-coordinates are not uniform. Use s.interp() to make data uniform.')

            # if mode is roll, value must be an integer
            if is_integer(value) == False:
                raise ValueError('shift must be an integer for mode = `roll`.')

            # is relative?
            if type_ in relative:
                value = self.shift_roll + value

            # apply shift
            if self.shift_roll != value:
                if self.shift_roll != 0:  # undo last shift
                    self._x, self._y = shifted(self.x, self.y, value=-self.shift_roll, mode='roll')
                if value != 0:
                    self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                self._shift_roll = value
        # HARD (x)
        elif mode in hard:
            # is relative?
            if type_ in relative:
                value = self.shift + value

            # apply shift
            if self.shift != value:
                if self.shift != 0:
                    self._x, self._y = shifted(self.x, self.y, value=-self._shift, mode='x')
                if value != 0:
                    self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                self._shift = value
        # SOFT (y)
        elif mode in soft:
            # data must be monotonic and increasing
            if self.shift_interp != value:
                if self.monotonicity is None:
                    self.check_monotonicity()
                if self.monotonicity != 'increasing':
                    raise ValueError('x array must be monotonicaly increasing.\nTip: use Spectrum.fix_monotonicity()')

                # is relative?
                if type_ in relative:
                    value = self.shift_interp + value
                
                # apply shift
                if self.shift_interp != 0:
                    self._x, self._y = shifted(self.x, self.y, value=-self.shift_interp, mode='interp')
                if value != 0:
                    self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                self._shift_interp = value
        else:
            raise ValueError(f'Invalid mode. Valid options are `roll`, `x`, `interp`.')

        # special
        self.peaks.set_shifts(value=value, mode=mode, type_=type_)

    def set_offset(self, value, type_='absolute'):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).
            type_ (string, optional): either 'absolute' (default) or 'relative'.

        Returns:
            None
        """
        # is relative
        if type_ in relative:
            value = self.offset + value

        # apply offset
        if self.offset != value:
            if self.offset != 0:
                self._y = self.y - self.offset
            if value != 0:
                self._y = self.y + value
            self._offset = value

        # special
        self.peaks.set_offsets(value=value, type_=type_)

    def set_factor(self, value, type_='absolute'):
        """Set y multiplicative factor.

        Args:
            value (number): multiplicative factor (y-coordinates will be
                multiplied by this value).
            type_ (string, optional): either 'absolute' (default) or 'relative'.

        Returns:
            None
        """
        # check
        if value == 0:
            raise ValueError('cannot set factor = 0.')

        # is relative
        if type_ in relative:
            value = self.factor * value

        # apply factor
        if self.factor != value:
            if self.factor != 1:
                self._y = self.y*self.factor**-1
            if value != 1:
                self._y = self.y*value
            self._factor = value

        # special
        self.peaks.set_factors(value=value, type_=type_)


    def floor(self, value=None, n=20, ranges=None):
        """Sets zero value for y-coordinates (shifts data verticaly).

        use ranges = (None, None) to bring the avg of all data point to zero.

        If data is monotonicaly, n will bring the avg of n point around value to
        zero. If data is not monotonicaly, the datapoint at value is brought to
        zero.

        Args:
            value (number, optional): x-coordinate to set y-coordinate value as
                zero. If `None`, the minium value of x is used.
            n (int, optional): number of data points to average and setting to
                zero. This only applies if data is monotonicaly.
            ranges (list, optional): Pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                Use None to indicate the minimum or maximum x value of the data.
                The mean value of the y-coordinates inside the defined range
                will be set to zero. This overwrites value and n.

        Returns:
            None
        """
        if ranges is None:
            if self.monotonicity is None:
                try:
                    self.check_monotonicity()
                except ValueError:
                    pass
            
            # value to index
            if value is None:
                i = index(self.x, min(self.x))
            else:
                i = index(self.x, value)

            # calculate offset
            if self.monotonicity in increasing or self.monotonicity in decreasing:
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
                value2 = -np.mean(self.y[start:end])
            else:
                value2 = -self.y[i]
        else:
            ranges = _check_ranges(ranges, vmin=min(self.x), vmax=max(self.x))
            _, y   = extract(self.x, self.y, ranges)
            value2 = -np.mean(y)
        self.set_offset(self.offset + value2)

    def flip(self):
        """Flip sign of x axis (via Spectrum.calib).

        Returns:
            None
        """
        self.set_calib(value=-self.calib)

    def normalize(self, value=1, ranges=None):
        """Set a factor such as the average y between ranges is equal to value.

        It uses Spectra.set_factors().

        Args:
            value (number): value. Default is 1.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data. If None, the maximum
                y value of the data is set to 1.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_factors`
        """
        if ranges is None:
            self.set_factor(self.factor*value/np.max(s.y))
        else:
            s = self.extract(ranges=ranges)
            self.set_factor(self.factor*value/np.mean(s.y))

    def zero(self, mode='max', ref_spectrum=0, ref_peak=0, shift_mode='x'):
        """[EXPERIMENTAL] Uses Spectrum.set_shift() to shift a specified position to zero.

        It shifts data using shift mode: 'hard' ('x').

        Args:
            mode (string, optional): method for identifing the point which the x-
                coordinate is zero. Default is 'max'. Options are: 
                'second derivative minimum'
                'cross-correlation reference'
                'max', 
                'min',
                'fitted peaks', 
                'peak'.
            ref (int, optional): if mode='peak' or mode='fitted peaks', this
                variable selects which peak to align data. If mode='cc', ref 
                must be a reference spectrum of type brixs.Spectrum.
            bkg_check (bool, optional): If mode = 'cc', data needs to be 
                normalized and shifts must be calculated and if bkg_check is 
                True, data with big bkg will raise an error. Later, I can
                describe here how this bkg check is done.

        Returns:
            None
        """
        # second derivative minimum
        if mode == 'second derivative minimum':
            s = self.derivative(2)
            value = -self.x[np.argmin(s.y)]
        # cross-correlation reference
        elif mode in cc:
            assert isinstance(ref_spectrum, Spectrum), 'ref must be of type brixs.Spectrum.'
            try:
                ref_spectrum.check_step()
            except ValueError:
                raise ValueError(f'Cannot shift data using mode = {mode}, because x-coordinates of reference spectrum are not uniform.')

            # max
            ss = Spectra([copy.deepcopy(ref_spectrum), copy.deepcopy(self)])
            ss.calculate_shifts('max')
            value1 = ss.calculated_shifts.y[1]
            ss.set_shifts()

            # cc
            ss.interp(x=ss[0].x)
            ss.calculate_factors(mode='area')
            ss.set_factors()
            ss.calculate_shifts(mode='cc')
            value2 = ss.calculated_shifts.y[1]*ref_spectrum.step
            value = value1+value2
        # max
        elif mode == 'max':
            value = -self.x[np.argmax(self.y)]
        # min
        elif mode == 'min':
            value = -self.x[np.argmin(self.y)]        
        # peak
        elif mode == 'peak':
            assert len(self.peaks) > 0, 'no peaks defined for this spectrum.'
            value = -self.peaks[ref_peak]['c'].value
        else:
            raise ValueError('mode not valid.\nValid modes: max, min, fitted peaks, peak.')
        
        # if roll, shift value must be an integer
        if shift_mode in roll and self.step is None:
            try:
                self.check_step()
            except ValueError:
                raise ValueError(f'Cannot shift data using shift_mode = {shift_mode}, because x-coordinates are not uniform.')
            value = int(value/self.step)
        
        # apply shift
        self.set_shift(value, shift_mode=shift_mode, type_='relative')

    # persistent modifiers
    def interp(self, x=None, start=None, stop=None, num=None, step=None):
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

        None:
            If none arguments is given, this function will adjust x to have regular spacing.

        Returns:
            None

        Warning:
            Spectrum.interp() will change data for good. There is no way to
                recover the previous state of the array.
        """
        self.check_monotonicity()
        if self.monotonicity != 'increasing':
            raise ValueError('monotonicity of data must be strictly increasing. Use Spectrum.fix_monotonicity().')

        if x is None:
            if start is None:
                start = min(self.x)
            if stop is None:
                stop = max(self.x)

            if step is not None:   # step overwrites num
                x = np.arange(start, stop, step=step)
            if num is not None:
                x = np.linspace(start, stop, num=num)
            else:
                temp, _ = extract(self.x, self.y, (start, stop))
                x = np.linspace(start, stop, num=len(temp))

        self._y = np.interp(x, self.x, self.y)
        self._x = x

        # check
        self.step         = None
        self.monotonicity = None
        # self._reset_modifiers()
        #special
        # self.peaks.clear()

    def remove(self, ranges):
        """Remove data points inside range."""
        ranges = _check_ranges(ranges, vmin=min(self.x), vmax=max(self.x))
        self._x, self._y   = extract(self.x, self.y, ranges, invert=True)

        # check
        self.step         = None
        # self.monotonicity = None
        # self._reset_modifiers()
        #special
        # self.peaks.clear()

    def makeover(self, ranges, mode='poly', deg=1, n=3):
        """[EXPERIMENTAL] Replace y data points inside range with averaged out data points nearby.
        
        mode can be average (avg) or polinomial (poly). If poly, the deg may also 
        be selected. Default is mode='poly'.
        
        n is the number of points to be used for calculating the average before 
        (or after) the data point to be substituted.
        
        it replaces data points from the beggining of the dataset by calculating
        the average of n previous points and by replacing data points from the
        end of the dataset by calculation the average of n subsequent data 
        points. If mode is poly, a polynomial curve is extrapolated and n points
        are used in the extrapolation.
        """
        ranges = _check_ranges(ranges, vmin=min(self.x), vmax=max(self.x))

        # check
        if self.monotonicity is None:
            self.check_monotonicity()
        if self.monotonicity != 'increasing' and self.monotonicity != 'decreasing':
            raise RuntimeError('Data must be monotonic (either increasing or decreasing).')

        if mode == 'poly':
            for r in ranges:
                x, _   = extract(self.x, self.y, r)
                # print(x)
                # # plot data that will be replaced
                # plt.plot(x, y)

                start = index(self.x, x[0])
                stop  = index(self.x, x[-1])+1

                # print(x)
                for i in range(len(x)): 
                    if stop-i<(start+(stop-start)/2):
                        # print(start)
                        # print(stop)
                        # print((start+stop-start)/2)
                        break

                    # # plot step by step by changing i
                    # if i == 1:
                    #     # plot point to be replaced
                    #     plt.plot(self.x[start+i], self.y[start+i],    marker='X', ms=10)
                    #     plt.plot(self.x[stop-i-1],  self.y[stop-i-1], marker='X', ms=10)

                    #     # plot points to be used as reference
                    #     plt.plot(self.x[start+i-n:start+i], self.y[start+i-n:start+i], color='g', marker='o')
                    #     plt.plot(self.x[stop-i:stop-i+n], self.y[stop-i:stop-i+n], color='g', marker='o', lw=3)

                    # calculate polyfit
                    p1 = np.polyfit(self.x[start+i-n:start+i], self.y[start+i-n:start+i], deg)
                    p2 = np.polyfit(self.x[stop-i:stop-i+n], self.y[stop-i:stop-i+n], deg)

                    self._y[start+i]  = np.polyval(p1, self.x[start+i])
                    self._y[stop-i-1]  = np.polyval(p2, self.x[stop-i-1])

                    # # plot step by step by changing i
                    # if i == 1:
                    #     # plot new point
                    #     plt.plot(self.x[start+i], np.polyval(p1, self.x[start+i]), marker='v', ms=10)
                    #     plt.plot(self.x[stop-i-1], np.polyval(p2, self.x[stop-i-1]), marker='v', ms=10)

                    # self._y[start+i] = np.mean(self.y[start+i-n:start+i])
                    # self._y[stop-i]  = np.mean(self.y[stop+i:stop+i+n])
                    pass
        elif mode == 'avg':
            for r in ranges:
                x, _   = extract(self.x, self.y, r)
                # print(x)
                # # plot data that will be replaced
                # plt.plot(x, y)

                start = index(self.x, x[0])
                stop  = index(self.x, x[-1])+1

                # print(x)
                for i in range(len(x)): 
                    if stop-i<(start+(stop-start)/2):
                        # print(start)
                        # print(stop)
                        # print((start+stop-start)/2)
                        break

                    # # plot step by step by changing i
                    # if i == 1:
                    #     # plot point to be replaced
                    #     plt.plot(self.x[start+i], self.y[start+i],    marker='X', ms=10)
                    #     plt.plot(self.x[stop-i-1],  self.y[stop-i-1], marker='X', ms=10)

                    #     # plot points to be used as reference
                    #     plt.plot(self.x[start+i-n:start+i], self.y[start+i-n:start+i], color='g', marker='o')
                    #     plt.plot(self.x[stop-i:stop-i+n], self.y[stop-i:stop-i+n], color='g', marker='o', lw=3)

                    self._y[start+i]  = np.mean(self.y[start+i-n:start+i])
                    self._y[stop-i-1] = np.mean(self.y[stop-i:stop-i+n])

    def crop(self, start=None, stop=None):
        """Crop data.

        Args:
            start (number, optional): start x value. If None, the minimum value of
                x will be used.
            stop (number, optional): final x value. If None, the maximum value of
                x will be used.

        Returns:
            None
        """
        if start is None and stop is None:
            return
        if start is None:
            start = min(self.x)
        if stop is None:
            stop = max(self.x)
        self._x, self._y  = extract(self.x, self.y, (start, stop))

    # extractors
    def extract(self, ranges):
        """Return the extracted data range from full data.

        Args:
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            :py:attr:`Spectrum`
        """
        ranges = _check_ranges(ranges, vmin=min(self.x), vmax=max(self.x))
        x, y   = extract(self.x, self.y, ranges)
        s      = Spectrum(x=x, y=y)

        # modifiers
        s._offset       = self.offset
        s._factor       = self.factor
        s._calib        = self.calib
        s._shift        = self.shift
        s._shift_roll   = self.shift_roll
        s._shift_interp = self.shift_interp

        # special
        # s._peaks        = self.peaks

        return self._transfer_attributes(s)

    def derivative(self, order=1):
        """Returns the derivative of y-coordinates as a function of x-coordinates.

        Args:
            order (number, optional): derivative order. Defaut is 1.

        Returns:
            Derivative spectrum
        """
        x, y = derivative(self.x, self.y, order=order)
        s = Spectrum(x=x, y=y)

        # modifiers
        # s._offset       = self.offset
        # s._factor       = self.factor
        # s._calib        = self.calib
        # s._shift        = self.shift
        # s._shift_roll   = self.shift_roll
        # s._shift_interp = self.shift_interp

        # special
        # s._peaks        = self.peaks

        return self._transfer_attributes(s)
    
    def smooth(self, n=8):
        """Returns the moving average of the data.

        Args:
            n (int, optional): number of points to average. Default is 8.

        Returns:
            spectrum of length given by (len(x)-n+1).
        """
        x = moving_average(self.x, n=n)
        y = moving_average(self.y, n=n)
        
        s = Spectrum(x=x, y=y)

        # modifiers
        s._offset       = self.offset
        s._factor       = self.factor
        s._calib        = self.calib
        s._shift        = self.shift
        s._shift_roll   = self.shift_roll
        s._shift_interp = self.shift_interp

        # special
        s._peaks        = self.peaks

        return self._transfer_attributes(s)

    def polyfit(self, deg=2, ranges=None):
        """Fit data with a polynomial. Wrapper for `numpy.polyfit()`_.

            popt, model, R2 = polyfit()
            opt = Polynomial coefficients, highest power first.
            R2
            model(x) = Model function f(x).

        Args:
            deg (int, optional): degree of the fitting polynomial. Default is 2.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            popt, f(x) model, R2

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        # ranges
        if ranges is None:
            x = self.x
            y = self.y
        else:
            ranges = _check_ranges(ranges, vmin=min(self.x), vmax=max(self.x))
            x, y   = extract(self.x, self.y, ranges)

        # fit
        popt  = np.polyfit(x, y, deg=deg)
        model = lambda x: np.polyval(popt, x)
        R2 =  1- (sum((self.y-model(self.x))**2)/sum((self.y-np.mean(self.y))**2))
    
        return popt, model, R2

    # calculation and info
    def calculate_area(self):
        """Returns the calculated area under the curve. Wrapper for `numpy.trapz()`_.

        Returns:
            Number

        .. _numpy.trapz(): https://numpy.org/doc/stable/reference/generated/numpy.trapz.html
        """
        return np.trapz(y=self.y, x=self.x)

    # plot and visualization
    def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, smooth=1, **kwargs):
        """Plot spectrum. Wrapper for `matplotlib.pyplot.plot()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number, optional): defines a vertical offset. Default is 0.
            shift (number, optional): horizontal shift value. Default is 0.
            factor (number, optional): multiplicative factor on the y axis.
                Default is 1.
            calib (number, optional): multiplicative factor on the x axis.
                Default is 1.
            smooth (int, optional): number of points to average data (moving average).
                Default is 1.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            `Line2D`_

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figure()
        x = (self.x*calib) + shift
        y = self.y*factor + offset

        if 'label' not in kwargs and hasattr(self, 'label'):
            kwargs['label'] = self.label

        if smooth > 1:
            x = moving_average(x, int(smooth))
            y = moving_average(y, int(smooth))
        return ax.plot(x, y, **kwargs)

    # Peaks
    def find_peaks(self, prominence=None, width=4, moving_average_window=8, ranges=None):
        """Find peaks. Wrapper for `scipy.signal.find_peaks()`_.

        Sets :py:attr:`peaks` attribute.

        Args:
            prominence (number, optional): minimum prominence of peaks in percentage
                of the maximum prominence [max(y) - min(y)]. Default is 5.
            width (number, optional): minimum number of data points defining a peak.
            moving_average_window (int, optional): window size for smoothing the
                data for finding peaks. Default is 4.

        Returns:
            None

        .. _scipy.signal.find_peaks(): https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
        """
        if ranges is None:
            s = self
        else:
            s = self.extract(self.x, self.y, ranges)

        # check monotonicity
        if s.monotonicity is None:
            s.check_monotonicity()

        # step
        if s.step is None:
            s.check_step()

        self.peaks.find(x=s.x, y=s.y, prominence=prominence, width=width, moving_average_window=moving_average_window)

    def fit_peak(self, asymmetry=None, moving_average_window=8, ranges=None):
        """Fits one peak. Initial guess is based on the maximum y value.

        Args:
            asymmetry (bool or dict, optional): if True, fits each half of the
                with a different width.
            fixed_m (False, number, or dict, optional): m is the amount of lorentzian
                contribution for a peak. If False, m will be fitted for each peak.
                If a number (from 0 to 1), this will be used as the value of m
                (fixed m).
            offset (bool, optional): if True, a offset value will be fitted.
            amp_bounds, c_bounds, fwhm_bounds (tuple, optional): minimum and
                maximum multiplication factor for boundary values. For amp, the
                bounds are set between amp*amp_bounds[0] and amp*amp_bounds[1].
                For c, bounds are set c+fwhm*c_bound[0] and c+fwhm*c_bound[1].
                Finaly, for fwhm, the bounds are set between fwhm1*fwhm_bounds[0]
                and fwhm1*fwhm_bounds[0]. Note that if fwhm1*fwhm_bounds[0] is
                less than zero, the minimum fwhm boundary is set to zero as it
                cannot be negative. If None, the data max and min limits will be
                used.

        Returns:
            None
        """
        if ranges is None:
            x0 = self.x
            y0 = self.y
        else:
            ranges = _check_ranges(ranges, vmin=min(self.x), vmax=max(self.x))
            x0, y0 = extract(self.x, self.y, ranges)

        if asymmetry is None:
            asymmetry = False

        # smoothing
        x = moving_average(x0, moving_average_window)
        y = moving_average(y0, moving_average_window)
        # plt.plot(x, y)

        # guess amp and c
        amp = max(y)
        c = x[np.argmax(y)]

        # guess fwhm
        try:
            w1 = x[np.argmax(y)] - x[:np.argmax(y)][::-1][index(y[:np.argmax(y)][::-1], max(y)/2)]
        except ValueError:
            w1 = x[np.argmax(y):][index(y[np.argmax(y):], max(y)/2)] - x[np.argmax(y)]
        try:
            w2 = x[np.argmax(y):][index(y[np.argmax(y):], max(y)/2)] - x[np.argmax(y)]
        except ValueError:
            w2 = x[np.argmax(y)] - x[:np.argmax(y)][::-1][index(y[:np.argmax(y)][::-1], max(y)/2)]
        w = w1 + w2

        # x, y = derivative(moving_average(self.x, 10), moving_average(self.y, 10))
        # w = np.abs(x[np.argmin(y)] - x[np.argmax(y)])
        if w == 0:
            w = 0.1*(max(self.x)-min(self.x))

        # peaks
        self.peaks.clear()
        if asymmetry:
            self.peaks.append(Peak(amp=amp, c=c, w1=w1, w2=w2, asymmetry=asymmetry))
        else:
            self.peaks.append(Peak(amp=amp, c=c, w=w, asymmetry=asymmetry))

        # fit
        self.fit_peaks(ranges=ranges)

    def fit_peaks(self, method='least_squares', ranges=None):
        """Fit peaks. Wrapper for `scipy.optimize.curve_fit()`_.

        Args:
            asymmetry (bool or dict, optional): if True, fits each half of the
                with a different width.
            fixed_m (False, number, or dict, optional): m is the amount of lorentzian
                contribution for a peak. If False, m will be fitted for each peak.
                If a number (from 0 to 1), this will be used as the value of m
                (fixed m).
            offset (bool, optional): if True, a offset value will be fitted.

        Returns:
            None

        .. _scipy.optimize.curve_fit(): https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
        """
        if ranges is None:
            self.check_monotonicity()
            x = self.x
            y = self.y
        else:
            ranges = _check_ranges(ranges, vmin=min(self.x), vmax=max(self.x))
            s = self.extract(ranges)
            s.check_monotonicity()
            x = s.x
            y = s.y

        # check if peaks is defined
        if len(self.peaks) == 0:
            raise ValueError('No peaks to fit.\nRun Spectrum.find_peaks().')

        # self.guess = 
        return self.peaks.fit(x, y, method=method, return_fitted_peaks=False)

        # # build model ==========================================================
        # p0, bounds_min, bounds_max, decode = self.peaks.build_guess()
        # model = self.peaks.build_model()

        # # guess ================================================================
        # self._guess = Spectrum(x=self.x, y=model(self.x, *p0))
        # self.guess._offset       = self.offset
        # self.guess._factor       = self.factor
        # self.guess._calib        = self.calib
        # self.guess._shift        = self.shift
        # self.guess._shift_roll   = self.shift_roll
        # self.guess._shift_interp = self.shift_interp
        # self.guess.peaks         = copy.deepcopy(self.peaks)

        # # fitting ==============================================================
        # # print(model)
        # # print(p0)
        # if verbose: print(f'Fitting data...')
        # popt, pcov = curve_fit(model, x, y, p0=p0, bounds=(bounds_min, bounds_max))
        # if verbose: print(f'Done!')

        # # save fitted peaks parameters =========================================
        # psigma = np.sqrt(np.diag(pcov))
        # peaks  = decode(popt, psigma)

        # # modifiers
        # for peak in peaks:
        #     total = self.shift_interp + self.shift
        #     if self.step is not None:
        #         total += self.shift_roll*self.step
        #     peak._shift  = total
        #     peak._offset = self.offset
        #     peak._calib  = self.calib
        #     peak._factor = self.factor

        # # fit ==================================================================
        # x_temp = peaks._find_suitable_x()
        # self._fit = Spectrum(x=x_temp, y=model(x_temp, *popt))
        # self.fit._offset       = self.offset
        # self.fit._factor       = self.factor
        # self.fit._calib        = self.calib
        # self.fit._shift        = self.shift
        # self.fit._shift_roll   = self.shift_roll
        # self.fit._shift_interp = self.shift_interp
        # self.fit._peaks = peaks
        # self.fit._pcov = pcov
        # self.fit._R2 =  1- (sum((self.y-model(self.x, *popt))**2)/sum((self.y-np.mean(self.y))**2))
        # # residue ==============================================================
        # if offset:
        #     self._residue = Spectrum(x=self.x, y=self.y-model(self.x, *popt)+popt[-1])
        # else:
        #     self._residue = Spectrum(x=self.x, y=self.y-model(self.x, *popt))
        # self.residue._offset       = self.offset
        # self.residue._factor       = self.factor
        # self.residue._calib        = self.calib
        # self.residue._shift        = self.shift
        # self.residue._shift_roll   = self.shift_roll
        # self.residue._shift_interp = self.shift_interp

    # OBSOLET
    def apply_correction(self, f):
        """Changes the values of x, y based on a function.

        Warning:
            All attributes are set to the default values.

        Example:
            f = lambda x, y: (x, y**2)

        Args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        Returns:
            None
        """
        # core ================================
        self._x, self._y = f(self.x, self.y)
        self._data[:, 0] = self._x
        self._data[:, 1] = self._y

        # modifiers
        self._factor       = 1
        self._offset       = 0
        self._calib        = 1
        self._shift        = 0
        self._shift_roll   = 0
        self._shift_interp = 0

        # check
        self._step         = None
        self._monotonicity = None

        # fit
        self._fit     = None
        self._residue = None
        self._guess   = None
        self._R2      = None

        # peaks
        self._peaks = Peaks()


class Spectra(metaclass=_Meta):
    """Returns a ``spectra`` class type object for dealing with many spectrum at a time.

    Args:
        n (int, optional): when the number of spectra `n` is known, one can preallocate
            the data array.
        data (list or array, optional): list of :py:class:`spectrum` objects.
            Overwrites ``n``.
        dirpath (string, path object, or list): folderpath, folder handle,
                or a list of filepaths. The list can be comprised of file and
                folderpaths. All files inside a folder where string can be found
                in the filename are imported. If the filename extension is .gz or .bz2,
                the file is first decompressed.

    Attributes:
        data (list of :py:attr:`Spectrum`): list with :py:attr:`Spectrum` objects.
        step (number, read only): step size for x-coordinates.
        length (number, read only): length of x-coordinates vector.
        x (number, read only): x-coordinates (only applicable if x-coordinates
            are the same for all spectra.
        shifts (list of int): Calculated shifts as int. Must always have the same length as data.
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

    _read_only     = ['step', 'length', 'x', 'monotonicity',
                      'calculated_calib', 'calculated_factors',
                      'calculated_offsets', 'calculated_shifts']
    _non_removable = []


    def __init__(self, *args, **kwargs):
        # basic
        self._data  = None
        self._dirpath  = ''
        self._filepath  = ''
        # self._peaks = Collection()

        # argument parsing
        data, dirpath, n = self._sort_args(args, kwargs)

        # sorting data
        if data is not None:
            self.data = data
            # _restart_attr is builtin inside data
        elif dirpath is not None:
            self.load(dirpath)
            # _restart_attr is builtin inside load
        elif n is not None:
            self._data     = [-1]*n
            # check
            self._restart_check_attr()
            self._restart_calculated_modifiers()
        else:
            self._data = []
            # check
            self._restart_check_attr()
            self._restart_calculated_modifiers()

    # Basic
    def __str__(self):
        return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

    def __repr__(self):
        return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

    def __next__(self):
        result = self.data[self._index]
        self._index +=1
        return result
        # End of Iteration
        raise StopIteration

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.data[item]
        elif isinstance(item, slice):
            temp = Spectra(self.data[item])

            dict = get_attributes(self)
            for n in dict:
                if n not in ['_x', '_y', '_data', '_fit', '_guess', '_residue', '_pcov', '_peaks', '_step','_monotonicity','_calculated_shifts','_calculated_factors','_calculated_offsets','_calculated_calib',]:
                    temp.__setattr__(n, dict[n])
            return temp        
        else:
            raise TypeError('Index must be int, not {}'.format(type(item).__name__))

    def __setitem__(self, item, value):
        if isinstance(value, Spectrum) == False:
            raise ValueError(f'value must be of type brixs.spectrum, not {type(value)}')
        self.data[item] = value
        self._restart_check_attr()
        self._calculated_shifts  = None
        self._calculated_calib  = None
        self._calculated_factors = None
        self._calculated_offsets = None

    def __len__(self):
        return len(self.data)

    def __delitem__(self, item):
        self.remove(item)

    # Support
    def _sort_args(self, args, kwargs):
        """checks initial arguments.

         Keyword arguments (kwargs) cannot be mixed with positional arguments.

        The hierarchy for Keyword arguments is: 1) data, 2) y (and x), and finaly
            3) filepath. For example, if `data` and `filepath` are passed as
            arguments, `filepath` is ignored.

        For positional arguments, if one data set is passed, it assumes it is
            `data`. If this one argument is of type string or Pathlib.Path, it
            assumes it is a filepath. If two data sets are passed, it will
            assume one is the x coordinates and the next one is the y coordinates.

        Raises:
            AttributeError: if kwargs and args cannot be read.

        Returns:
            data, x, y, filepath
        """
        # print(kwargs)
        # print(args)
        if kwargs != {} and args != ():
            raise AttributeError('cannot mix key word arguments with positional arguments. Key word arguents are `x`, `y`, `data`, and `filepath`.')
        if any([item not in ['data', 'n', 'dirpath'] for item in kwargs.keys()]):
            raise AttributeError(f'invalid attributes.\nValid atributes are `data`, `n`, and `dirpath`\nInput attributes: {kwargs.keys()}')

        data = None
        n    = None
        dirpath = None
        if 'data' in kwargs:
            data = kwargs['data']
        elif 'dirpath' in kwargs:
            dirpath = kwargs['dirpath']
        elif 'n' in kwargs:
            n = kwargs['n']
        elif len(args) == 1:
            if isinstance(args[0], int):
                n = args[0]
            elif isinstance(args[0], str) or isinstance(args[0], Path):
                dirpath = [args[0], ]
            elif isinstance(args[0], Iterable):
                if isinstance(args[0][0], str) or isinstance(args[0][0], Path):
                    dirpath = args[0]
                else:
                    data = args[0]
        elif len(args) > 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                dirpath = args
            elif isinstance(args[0], Iterable):
                data = args
        return data, dirpath, n

    def _restart_check_attr(self):
        """restart check attributes."""
        # check x attrs
        self._length       = None
        self._step         = None
        self._x            = None
        self._monotonicity = None
    
    def _restart_calculated_modifiers(self):
        """restart calculated modifiers."""
        self._calculated_shifts  = None
        self._calculated_calib   = None
        self._calculated_factors = None
        self._calculated_offsets = None

    def _gather_ys(self, ranges=None):
        """Return two lists, x and y's. This structure speeds up some operations."""
        # check x
        if self.x is None:
            self.check_same_x()

        ys = np.zeros((self.length, len(self)))
        for i in range(len(self)):
            ys[:, i] = self.data[i].y
        if ranges is None:
            x = self.data[0].x
        else:
            x, ys = extract(self.data[0].x, ys, ranges=ranges)
        return x, ys

    def _transfer_attributes(self, object):
        """Transfer user defined attributes to output objects."""

        # copy user defined attributes
        for attr in self.__dict__:
            if attr.startswith('_') == False:
                object.__setattr__(attr, self.__dict__[attr])

        return object
    
    def _transfer_attributes_from_each_spectrum(self, object):
        """Transfer user defined attributes from each spectrum to output objects.

        User defined attr are saved to the summed spectrum. A number is 
            added in front of the attr name indicating the index of the 
            spectrum in which the attr was copied from.
        """
        # copy user defined attributes
        for j in range(len(self)):
            for attr in self[j].get_user_defined_attrs():
                object.__setattr__(attr+f'_{j}', self[j].__dict__[attr])
        return object

    # attributes and properties
    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        if value is None:
            self._data = []
        else:
            if isinstance(value, Iterable):
                for i, s in enumerate(value):
                    if isinstance(s, Spectrum) == False:
                        raise ValueError(f'All entries must be of type brixs.spectrum.\nEntry {i} is of type {type(s)}')
                # self._data = copy.deepcopy(value)
                self._data = value
            else:
                raise ValueError('data must be a list.')

        # check
        self._restart_check_attr()
        self._restart_calculated_modifiers()
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def area(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].area
        return temp
    @area.setter
    def area(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @area.deleter
    def area(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].shift
        return temp
    @shift.setter
    def shift(self, value):
        self.set_shifts(value=value, mode='x')
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_roll(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].shift_roll
        return temp
    @shift_roll.setter
    def shift_roll(self, value):
        self.set_shifts(value=value, mode='roll')
    @shift_roll.deleter
    def shift_roll(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_interp(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].shift_interp
        return temp
    @shift_interp.setter
    def shift_interp(self, value):
        self.set_shifts(value=value, mode='interp')
    @shift_interp.deleter
    def shift_interp(self):
        raise AttributeError('Cannot delete object.')

    @property
    def calib(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].calib
        return temp
    @calib.setter
    def calib(self, value):
        self.set_calib(value)
    @calib.deleter
    def calib(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].factor
        return temp
    @factor.setter
    def factor(self, value):
        self.set_factors(value)
    @factor.deleter
    def factor(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].offset
        return temp
    @offset.setter
    def offset(self, value):
        self.set_offsets(value)
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    @property
    def sum(self):
        return self.calculate_sum()
    @sum.setter
    def sum(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @sum.deleter
    def sum(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def dirpath(self):
        return self._dirpath
    @dirpath.setter
    def dirpath(self, value):
        if value is None:
            value = ''
        elif isinstance(value, str) or isinstance(value, Path):
            self._dirpath = value
        else:
            raise TypeError(r'Invalid type ' + str(type(value)) + 'for dirpath\ndirpath can only be str or pathlib.Path type.')
    @dirpath.deleter
    def dirpath(self):
        raise AttributeError('Cannot delete object.')  
    
    @property
    def filepath(self):
        return self._filepath
    @filepath.setter
    def filepath(self, value):
        if value is None:
            value = ''
        elif isinstance(value, str) or isinstance(value, Path):
            self._filepath = value
        else:
            raise TypeError(r'Invalid type ' + str(type(value)) + 'for filepath\nfilepath can only be str or pathlib.Path type.')
    @filepath.deleter
    def filepath(self):
        raise AttributeError('Cannot delete object.')  
    
    @property
    def xlabels(self):
        return [s.xlabel for s in self]
    @xlabels.setter
    def xlabels(self, value):
        if isinstance(value, str):
            for s in self:
                s.xlabel = value
        elif isinstance(value, Iterable):
            assert len(value) == len(self), f'Lenght of xlabels ({len(value)}) does not match with number of spectra ({len(self)})\nxlabels must be set to a string or a list with the same lenght as the number of spectra.' 
            for i, s in enumerate(self):
                s.xlabel = value[i]
        else:
            raise TypeError('Invalid type for xlabels\nxlabels must be a string of list of strings with the same length as the number of spectra.') 
    @xlabels.deleter
    def xlabels(self):
        raise AttributeError('Cannot delete object.')

    @property
    def ylabels(self):
        return [s.ylabel for s in self]
    @ylabels.setter
    def ylabels(self, value):
        if isinstance(value, str):
            for s in self:
                s.ylabel = value
        elif isinstance(value, Iterable):
            assert len(value) == len(self), f'Lenght of ylabels ({len(value)}) does not match with number of spectra ({len(self)})\nylabels must be set to a string or a list with the same lenght as the number of spectra.' 
            for i, s in enumerate(self):
                s.ylabel = value[i]
        else:
            raise TypeError('Invalid type for ylabel\nylabel must be a string of list of strings with the same length as the number of spectra.')
    @ylabels.deleter
    def ylabels(self):
        raise AttributeError('Cannot delete object.')

    # @property
    # def R2(self):
    #     temp = [0]*len(self)
    #     for i in range(len(self)):
    #         temp[i] = self[i].R2
    #     return temp
    # @R2.setter
    # def R2(self, value):
    #     raise AttributeError('Cannot be set.')
    # @R2.deleter
    # def R2(self):
    #     raise AttributeError('Cannot delete object.')

    # @property
    # def fit(self):

    #     self._check_fit()
    #     ss = Spectra(n=len(self))
    #     for i in range(len(self)):
    #         ss[i] = self[i].fit
    #     return ss
    # @fit.setter
    # def fit(self, value):
    #     raise AttributeError('Attribute is "read only". Cannot set attribute.')
    # @fit.deleter
    # def fit(self):
    #     raise AttributeError('Attribute cannot be deleted.')

    # @property
    # def guess(self):
    #     self._check_fit()
    #     ss = Spectra(n=len(self))
    #     for i in range(len(self)):
    #         ss[i] = self[i].guess
    #     return ss
    # @guess.setter
    # def guess(self, value):
    #     raise AttributeError('Attribute is "read only". Cannot set attribute.')
    # @guess.deleter
    # def guess(self):
    #     raise AttributeError('Attribute cannot be deleted.')

    # @property
    # def residue(self):
    #     self._check_fit()
    #     ss = Spectra(n=len(self))
    #     for i in range(len(self)):
    #         ss[i] = self[i].residue
    #     # if None in temp:
    #     #     print(f'Some spectra do not have fit defined.\nfit={temp}')
    #     return ss
    # @residue.setter
    # def residue(self, value):
    #     raise AttributeError('Attribute is "read only". Cannot set attribute.')
    # @residue.deleter
    # def residue(self):
    #     raise AttributeError('Attribute cannot be deleted.')

    @property
    def peaks(self):
        return Collection([s.peaks for s in self])
    @peaks.setter
    def peaks(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @peaks.deleter
    def peaks(self):
        raise AttributeError('Attribute cannot be deleted.')

    # @property
    # def errors(self):
    #     return self.get_errors()
    # @errors.setter
    # def errors(self, value):
    #     raise AttributeError('Attribute is "read only". Cannot set attribute.')
    # @errors.deleter
    # def errors(self):
    #     raise AttributeError('Attribute cannot be deleted.')

    @property
    def map(self):
        return self.calculate_map()
    @map.setter
    def map(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @map.deleter
    def map(self):
        raise AttributeError('Attribute cannot be deleted.')

    # basic
    def append(self, *args, **kwargs):
        """[Not fully tested] Append spectrum to the spectrum list.

        Args:
            s (Spectrum obj or list): Spectrum object to be added or
                list of Spectrum.
            data (list or array, optional): two column list (or array).
            x (list or array, optional): x values (1D list/array). Overwrites `data`.
            y (list or array, optional): y values (1D list/array). Overwrites `data`.
            filepath (str, optional)

        Returns:
            None

        See Also:
            :py:func:`Spectra.remove`.
        """
        # initial check
        if kwargs != {} and args != ():
            raise AttributeError('cannot mix key word arguments with positional arguments. Key word arguents are `x`, `y`, `data`, and `filepath`.')
        if any([item not in ['data', 'x', 'y', 'filepath'] for item in kwargs.keys()]):
            raise AttributeError(f'invalid attributes.\nValid atributes are `s`, `data`, `x`, `y`, and `filepath`\nInput attributes: {kwargs.keys()}')

        # is spectrum?
        if 's' in kwargs:
            assert isinstance(s, Spectrum), 'Spectrum must be of type brixs.Spectrum.'
            self._data += [s]
        else:
            if len(args) == 1 and isinstance(args[0], Spectrum):
                self._data += [args[0]]
                return
            
        # if not spectrum
        if len(args)>0 and kwargs != {}:
            s = Spectrum(*args, **kwargs)
            self.append(s=s)
        else:
            raise ValueError('No data to append.')
        
        # check
        self._restart_check_attr()
        self._restart_calculated_modifiers()

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

        # check
        self._restart_check_attr()
        self._restart_calculated_modifiers()

    def reorder(self, i1, i2):
        temp     = copy.deepcopy(self[i1])
        self[i1] = copy.deepcopy(self[i2])
        self[i2] = temp

        # check
        # self._restart_check_attr()
        self._restart_calculated_modifiers()
    
    def flip_order(self):
        self._data = self._data[::-1]

        # check
        # self._restart_check_attr()
        self._restart_calculated_modifiers()

    def get_by_attr(self, attr, value, approx=True, verbose=True):
        """Return spectrum with attr closest to value.

        Only works for numerical attr. If multiple spectrum have same value, it
        returns the first one.

        attr: name of the attr
        value: value to look for
        approx: if False, returns spectrum only if value is exact.
        """
        # check all spectra has attr
        has = [hasattr(s, attr) for s in self]
        if False in has:
            raise ValueError(f'Some spectra do not have attr: {attr}.\nhas attr:{has}')

        # get attr
        temp = [getattr(s, attr) for s in self]

        if approx:
            i = index(temp, value)
        else:
            i = temp.index(value)

        if verbose:
            print(f'Spectrum number {i} have {attr}={temp[i]}')
                
        return self[i]
    
    # save and load
    def save(self, dirpath=None, prefix='spectrum_', suffix='.dat', filenames=None, zfill=None, only_data=False, check_overwrite=True, verbose=True, **kwargs):
        r"""Save spectra. Wrapper for `numpy.savetxt()`_.

        Args:
            dirpath (string or pathlib.Path, optional): folderpath, folder handle. 
                If None, filepath can be inferred from a
                user defined attr 'ss.dirpath'. Default is None. Last used 
                dirpath is saved to ss.dirpath.
            prefix (string, optional): prefix used for naming the files.
            suffix (string, optional): suffix used for naming the files. If the
                filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            filenames (string, optional): overwrites prefix and suffix. Use `{<attr>}` to 
                include spectrum attr in the name. Use `{i}` to include the
                spectrum index. Example: 'spectrum_{i}_T{T}K.dat'
            zfill (int, optional): number of digits for file numbering. If `None`,
                zfill will be determined.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            check_overwrite (bool, optional): if True, it will check if files exist
                and ask if user wants to overwrite file.
            verbose (bool, optional): turn verbose on and off. Default is `False`.
            **kwargs: kwargs are passed to ``np.savetxt()`` that saves the data.
            
        If not specified, the following parameters are passed to `numpy.savetxt()`_:

        Args:
            fmt (string, or list, optional): A single format (like ``%10.5f``), or a
                sequence of formats. See numpy's documentation for more information.
                If not specified, best fmt is calculated based on the
                number of decimal places of the data. Specifing fmt makes the code
                runs a little faster (not much tough, but it might make a difference
                if saving a lot of files).
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (", ").
            newline (str, optional): String or character separating lines.
                Default is ``\n``.
            header (bool, optional): String that will be written at the beginning of the file.
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        # dirpath
        if dirpath is None:
            if dirpath == '':
                raise TypeError("Missing 1 required argument: 'dirpath'")
            else:
                dirpath = self.dirpath               
        dirpath = Path(dirpath)

        assert dirpath.exists(), f'dirpath does not exists.\ndirpath: {dirpath}'
        assert dirpath.is_dir(), f'dirpath is not a directory.\ndirpath: {dirpath}'

        # set filenames
        if zfill is None:
            zfill = n_digits(len(self)-1)[0]

        # saving
        # if verbose: print('saving {} files...')
        for i, s in enumerate(self.data):
            if filenames is not None:
                filename = filenames
                for slot in filenames.split('{'):
                    if '}' in slot:
                        name = slot.split('}')[0]
                        if '{' + name + '}' == '{i}':
                            filename = filename.replace('{i}', str(i))
                        else:
                            filename = filename.replace('{' + name + '}', str(getattr(self[i], name)))
            else:
                filename = f'{prefix}' + f'{i}'.zfill(zfill) + f'{suffix}'
            if verbose:  print(f'{i}/{len(self)-1}: {filename}')
            s.save(filepath=dirpath/filename, only_data=only_data, check_overwrite=check_overwrite, verbose=verbose, **kwargs)
        
        # save dirpath
        self.dirpath = dirpath
        if verbose: print('Done!')
    
    def load(self, dirpath, string='*', verbose=False, **kwargs):
        """Load data from text files. Wrapper for `numpy.genfromtxt()`_.

        Args:
            dirpath (string, path object, or list): folderpath, folder handle,
                or a list of filepaths. All files inside a folder where string 
                can be found in the filename are imported. If the filename 
                extension is .gz or .bz2, the file is first decompressed. 
                Last used dirpath is saved to ss.dirpath.
            string (str, optional): file names without this string will be ignored.
                Use '*' for matching anything. Default is '*'.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``# ``. Attributes picked up
                from the header will be loaded too.

        Returns:
            None

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        # reset data
        self._data = []
        self._restart_check_attr()
        self._restart_calculated_modifiers()

        fl = filelist(dirpath=dirpath, string=string)
        for i, filepath in enumerate(fl):
            if verbose: print(f'Loading: {dirpath}')
            if verbose: print(f'    {i+1}/{len(fl)}: {filepath.name}')
            self.append(Spectrum(filepath=filepath))

        # # if dirpath is str
        # if isinstance(dirpath, str) or isinstance(dirpath, Path):
        #     dirpath = Path(dirpath)



        #     if dirpath.is_file():
        #         if verbose: print('dirpath is a file')
        #         if verbose: print('Loading...')
        #         self.append(Spectrum(filepath=dirpath))
        #         if verbose: print('Done!')
        #         return
        #     elif dirpath.is_dir():
        #         dirpath = [dirpath, ]
        #     else:
        #         raise ValueError(f'cannot read dirpath.\ndirpath: {dirpath}')

        # # if dirpath is iterable
        # if isinstance(dirpath, Iterable):
        #     if verbose: print('Loading...')
        #     for j, filepath in enumerate(dirpath):
        #         if verbose: print(f'{j+1}/{len(dirpath)}: {filepath}')

        #         # if Path(filepath).is_dir():
        #         #     fl = filelist(dirpath=filepath, string=string)
        #         #     for i, f in enumerate(fl):
        #         #         if verbose: print(f'    {i+1}/{len(fl)}: {f}')
        #         #         self.append(Spectrum(filepath=f))

        #         elif Path(filepath).is_file():
        #             self.append(Spectrum(filepath=filepath))
        #         else:
        #             raise ValueError(f'cannot read dirpath.\dirpath: {dirpath}')
        
        # save dirpath
        self.dirpath = dirpath

        if verbose: print('Done!')

    def _header(self, verbose):
        """Gather attrs to be saved to a file."""
        header = ''
        temp = get_attributes(self)
        for name in temp:
            if name not in ['_data', '_x', '_filepath', '_dirpath', '_calculated_shifts', '_calculated_calib', '_calculated_factors', '_calculated_offsets', 'calculated_shifts']:
                if temp[name] is None:
                    header += f'{name}: None'  + '\n'
                elif isinstance(temp[name], str):
                    temp2 = str(temp[name]).replace('\n','\\n')
                    header += f'{name}: \"{temp2}\"'  + '\n'
                elif isinstance(temp[name], Iterable):
                    header += f'{name}: {list(temp[name])}'  + '\n'
                elif isinstance(temp[name], dict):
                    if verbose:
                        type_ = str(type(temp[name]))
                        print(r'Warning: Cannot save attr of type: ' + type_ + r'.\attr: '+ name + r'.\nTo turn off this warning, set verbose to False.')
                elif isinstance(temp[name], numbers.Number):
                    tosave = str(temp[name])
                    if tosave[-1] == '\n':
                        tosave = tosave[:-1]
                    header += f'{name}: {tosave}'  + '\n'
                else:
                    temp2 = str(temp[name]).replace('\n','\\n')
                    header += f'{name}: \"{temp2}\"'  + '\n'
        return header[:-1]

    def save_all_single_file(self, filepath='./spectra.dat', only_data=False, ranges=None, verbose=True, **kwargs):
        r"""Save all Spectra in one single file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                The additive factor will be calculated for spectra such the average value of
                data inside ranges is the same as for the first spectrum. Spectra
                outside this range will not be saved.

        If not specified, the following parameters are passed to `numpy.savetxt()`_:

        Args:
            fmt (string, or list, optional): A single format (like ``%10.5f``), or a
                sequence of formats. See numpy's documentation for more information.
                If not specified, best fmt is calculated based on the
                number of decimal places of the data. Specifing fmt makes the code
                runs a little faster (not much tough, but it might make a difference
                if saving a lot of files).
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (", ").
            newline (str, optional): String or character separating lines.
                Default is ``\n``.
            header (bool, optional): String that will be written at the beginning of the file.
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        # dirpath
        if filepath is None:
            if filepath == '':
                raise TypeError("Missing 1 required argument: 'filepath'")
            else:
                filepath = self.filepath               
        filepath = Path(filepath)

        assert filepath.parent.exists(), f'filepath folder does not exists.\ndirpath: {filepath.parent}'

        # check x is the same
        if self.x is None:
            try:
                self.check_same_x()
            except ValueError:
                raise ValueError('Cannot save spectra in one file. x axis are different.\nMaybe try interpolating the x axis (Spectra.interp()) or use Spectra.save() to save spectra in multiple files.')
        
        # gather ys
        x, ys = self._gather_ys(ranges=ranges)

        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([n_decimal_places(x) for x in self.x])
            temp_decimal = max([n_decimal_places(y) for y in flatten(ys)])
            if temp_decimal > decimal:
                decimal = copy.deepcopy(temp_decimal)
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        
        # dataset
        final = np.zeros((self.length, len(self)+1))
        final[:, 0]  = x
        final[:, 1:] = ys
        

        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
            np.savetxt(Path(filepath), final, **kwargs)
        else:
            if 'header' not in kwargs:
                kwargs['header'] = self._header(verbose=verbose)
            else:
                if kwargs['header'] == '':
                    kwargs['header'] = self._header(verbose=verbose)
                elif kwargs['header'][-1] != '\n':
                    kwargs['header'] += '\n'

            # gather ylabels
            if self[0].xlabel != '':
                if '' in self.ylabels:
                    if verbose:
                        print('Warning: Some spectra have no ylabel.\n ylabels: {self.ylabels}.\nFile is being saved without column names.')

                kwargs['header'] += '\n' + self[0].xlabel + kwargs['delimiter'] + kwargs['delimiter'].join(self.ylabels)

        np.savetxt(Path(filepath), final, **kwargs)

    def load_from_single_file(self, filepath, only_data=False, verbose=True, **kwargs):
        """
        if usecols
        """
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '#'
        if 'usecols' not in kwargs:
            kwargs['usecols'] = None

        # read data
        data = np.genfromtxt(Path(filepath), **kwargs)
        x  = data[:, 0]
        ys = data[:, 1:]

        # artificial usecols
        if kwargs['usecols'] is None:
            kwargs['usecols'] = [i for i in range(data.shape[1])]

        # reset data
        self._data     = []
        self._restart_check_attr()
        self._restart_calculated_modifiers()

        # allocate data
        for i in range(ys.shape[1]):
            s = Spectrum(x=x, y=ys[:, i])
            self.append(s)

        # read header
        if only_data is False:
            header = load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            if header:
                for line in header:
                    if ':' not in line:
                        temp = line.split(kwargs['delimiter'])
                        self.xlabels = temp[kwargs['usecols'][0]].replace(kwargs['comments'], '').replace('\n', '').replace('\r', '').strip()
                        self.ylabels = [temp[kwargs['usecols'][i]].replace(kwargs['comments'], '').replace('\n', '').replace('\r', '').strip() for i in kwargs['usecols'][1:]]
                    else:
                        # extract name and value
                        name = line[1:-1].split(':')[0].strip()
                        value = eval(str(':'.join(line[1:-1].split(':')[1:])).strip())
                        try:
                            setattr(self, name, value)
                        except Exception as e:
                            if verbose:
                                print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')
        # save filepath
        self.filepath = str(filepath)

    # check
    def check_monotonicity(self):
        """Sets monotonicity attribute to 'increasing' or 'decreasing'.

        Raises:
            ValueError if data is not monotonic.

        Returns:
            None
        """
        monotonicity = [None]*len(self)
        for i in range(len(self)):
            try:
                self[i].check_monotonicity()
                monotonicity[i] = self.data[i].monotonicity
            except ValueError:
                pass
        if all(x == 'increasing' for x in monotonicity):
            self._monotonicity = 'increasing'
        elif all(x == 'decreasing' for x in monotonicity):
            self._monotonicity = 'decreasing'
        else:
            text = ''
            for i in range(len(self)):
                text += f'spectrum: {i}, motonicity: {monotonicity[i]}\n'
            raise ValueError(f'some spectra have different monotonicity (increasing, decreasing) or no monotonicity at all (None): \n{text}')

    def fix_monotonicity(self, mode='increasing'):
        """Rearrange x, y such as x array is monotonically increasing or decreasing.

        Args:
            mode (str, optional): increasing or decreasing.

        Returns:
            None
        """
        for s in self.data:
            s.fix_monotonicity(mode=mode)
        self.check_monotonicity()

    def check_length(self):
        """Checks if all spectra has the same length.

        If all spectra have the same length, it sets :py:attr:`Spectra.length` = length.
        Otherwise, it raises an error.

        Returns:
            None

        Raises:
            ValueError: spectra does not have the same length.

        See Also:
            :py:func:`Spectra.check_step`, :py:func:`Spectra.check_same_x`.
        """
        # collect
        length = [None]*len(self)
        for i in range(len(self)):
            try:
                length[i] = len(self.data[i])
            except ValueError:
                pass

        # apply or raise error
        if all(x == length[0] for x in length):
            self._length = length[0]
        else:
            text = ''
            for i in range(len(self)):
                text += f'spectrum: {i}, length: {length[i]}\n'
            raise ValueError(f'some spectra have different length: \n{text}')

    def check_step(self, max_error=None):
        """Check step between data points in the x-coordinates.

            If data has a well defined step size, it sets :py:attr:`Spectra.step` = step.
            Otherwise, it raises an error.

            Args:
                max_error (number, optional): percentage value (in terms of the
                    average x step) of the maximum allowed error. If None, the 
                    Default value from settings will be used.

            Two checks are performed:

            1. Checks if the step between two data points is the same
            through out the x vector for each spectrum (vector uniformity).
            See Spectrum.check_step().

            (max(steps) - min(steps))/np.mean(steps) * 100 < max_error

            2. Checks if the step is the same between different spectra.

            (step[i]-step[i+1])/avg(step) * 100 < max_error

            Returns:
                None

            Raises:
                ValueError: If condition 1, or 2 are not satisfied.

            See Also:
                :py:func:`Spectra.check_length`, :py:func:`Spectra.check_same_x`.
        """
        if max_error is None:
            max_error = settings.MAX_ERROR_STEP_X

        if self.x is None:
            # try and see if spectra have the same x
            try:
                self.check_same_x(max_error=max_error)

                # check step uniformity
                temp = Spectrum(x=self.x, y=self.x)
                try:
                    temp.check_step(max_error=max_error)
                except ValueError:
                    raise ValueError(f"Spectra have the same x-coordinates, but it is not uniform.")
                self._step = temp.step
                return

            # if spectra have different x-coordinates
            except ValueError:

                # 1) check step uniformity
                steps = ['not uniform']*len(self)
                for idx, s in enumerate(self.data):
                    if s.step is None:
                        try:
                            s.check_step(max_error=max_error)
                            steps[idx] = s.step
                        except ValueError:
                            pass
                    else:
                        steps[idx] = s.step
                # raise error
                if all(x == steps[0] for x in steps):
                    pass
                else:
                    text = ''
                    for i in range(len(self)):
                        text += f'spectrum: {i}, step: {step[i]}\n'
                    raise ValueError(f'some spectra have different step: \n{text}')

                # 2) check step between spectra
                avg_step = np.mean(steps)
                if sum([abs(steps[i]-steps[i+1]) > abs(avg_step*max_error/100) for i in range(len(self)-1)]) > 0:
                    self._step = None
                    raise ValueError(f"Spectra seems to have different step size. Calculated step sizes = {steps}")
                self._step = avg_step
                return
        
        # if all spectra have the same x, check step becames easier
        else:
            # check step uniformity
            temp = Spectrum(x=self.x, y=self.x)
            try:
                temp.check_step(max_error=max_error)
            except ValueError:
                raise ValueError(f"Spectra have the same x-coordinates, but it is not uniform.")
            self._step = temp.step
            return

    def check_same_x(self, max_error=None):
        """Check if spectra have same x-coordinates.

        If data has same x-coordinates, it sets :py:attr:`Spectra.x` = x.
        Otherwise, it raises an error.

        max(s[i].x - s[i+1].x)/avg(step) * 100 < max_error

        Args:
            max_error (number, optional): percentage value (in terms of the
                average x step) of the maximum allowed error. If None, the 
                Default value from settings will be used.

        Returns:
            None

        Raises:
            ValueError: If any x-coodinate of any two spectrum is different.

        See Also:
            :py:func:`Spectra.check_length`, :py:func:`Spectra.check_step`.
        """
        if max_error is None:
            max_error = settings.MAX_ERROR_STEP_X

        # check length
        self.check_length()

        # average step
        if self.step is None:
            step = 0
            for s in self:
                step += np.mean(np.diff(s.x))
            step = step/len(self)
        else:
            step = self.step

        # check x between spectra
        x = ['same as the previous']*len(self)
        for idx in range(len(self)-1):
            if max(abs(self[idx].x - self[idx+1].x))*100/abs(step) > max_error:
                x[idx] = 'different'

        # apply or raise error
        if 'different' not in x:
            self._x = self[0].x
        else:
            text = ''
            for i in range(len(self)):
                text += f'spectrum: {i}: {x[i]}\n'
            raise ValueError(f'some spectra have different x: \n{text}\n\nUse brixs.Spectra.interp() to interpolate the data and make the x axis for different spectra match.')

        # update length
        self._length = len(self.x)

    # persistent modifiers
    def interp(self, crop=True, x=None, start=None, stop=None, num=None, step=None):
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

        None:
            If none arguments is given, this function will adjust all spectra to
            have the same x.

        Returns:
            None

        Warning:
            Spectra.interp() will change data for good. There is no way to
                recover the previous state of the array.
        """
        if x is None:
            if self.x is None:
                if crop:
                    if start is None:
                        start = max([min(s.x) for s in self.data])
                    if stop is None:
                        stop = min([max(s.x) for s in self.data])
                else:
                    if start is None:
                        start = min([min(s.x) for s in self.data])
                    if stop is None:
                        stop = max([max(s.x) for s in self.data])
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

        for s in self:
            s.interp(x=x, start=start, stop=stop, num=num, step=step)
        # check
        self._restart_check_attr()
        # self._restart_calculated_modifiers()

    def crop(self, start=None, stop=None):
        """Crop data.

        Args:
            start (number, optional): start x value. If None, the minimum value of
                x will be used.
            stop (number, optional): final x value. If None, the maximum value of
                x will be used.

        Returns:
            None
        """
        if self.x is None:
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
        elif start is not None or stop is not None:
            for s in self.data:
                s.crop(start, stop)
        # check
        self._restart_check_attr()
        # self._restart_calculated_modifiers()

    # plot and visualization
    def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, smooth=1, **kwargs):
        """Plot spectra. Wrapper for `matplotlib.pyplot.plot()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number or list, optional): defines a vertical offset. Default is 0.
            shift (number or list, optional): horizontal shift value. Default is 0.
            factor (number or list, optional): multiplicative factor on the y axis.
                Default is 1.
            calib (number or list, optional): multiplicative factor on the x axis.
                Default is 1.
            smooth (int, optional): number of points to average data (moving average).
                Default is 1.
            vertical_increment or vi (number): vertical increment between curves
                in terms of percentage of the maximum delta y.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.


        Returns:
            `Line2D`_ list, offsets list, shifts list

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figure()


        # percentage wise vertical increment ====================
        if 'vi' in kwargs and 'vertical_increment' in kwargs:
            raise SyntaxError('keyword argument repeated: vertical increment and vi')
        elif 'vi' in kwargs or 'vertical_increment' in kwargs:
            if 'vi' in kwargs:
                vi = kwargs['vi']
                del kwargs['vi']
            if 'vertical_increment' in kwargs:
                vi = kwargs['vertical_increment']
                del kwargs['vertical_increment']
            temp = [0]*len(self)
            for i in range(len(self)):
                temp[i] = -(max(self.data[i].y) - min(self.data[i].y))
            vi = max(temp)*factor*vi/100
        else:
            vi = 0

        # percentage wise horizontal increment ====================
        if 'hi' in kwargs and 'horizontal_increment' in kwargs:
            raise SyntaxError('keyword argument repeated: horizontal increment and hi')
        elif 'hi' in kwargs or 'horizontal_increment' in kwargs:
            if 'hi' in kwargs:
                hi = kwargs['hi']
                del kwargs['hi']
            if 'horizontal_increment' in kwargs:
                hi = kwargs['horizontal_increment']
                del kwargs['horizontal_increment']
            temp = [0]*len(self)
            for i in range(len(self)):
                temp[i] = max(self.data[i].x) - min(self.data[i].x)
            hi = max(temp)*factor*hi/100
        else:
            hi = 0

        # offset
        if isinstance(offset, Iterable):
            if len(offset) != len(self):
                raise ValueError(f'offset must be a number of a list with length compatible with the number of spectra.\nnumber of offsets: {len(offset)}\nnumber of spectra: {len(self)}')
        else:
            offset = [offset]*len(self)
            for i in range(len(self)):
                offset[i] = offset[i]+(vi*i)

        # shift
        if isinstance(shift, Iterable):
            if len(shift) != len(self):
                raise ValueError(f'shift must be a number of a list with length compatible with the number of spectra.\nnumber of shift: {len(shift)}\nnumber of spectra: {len(self)}')
        else:
            shift = [shift]*len(self)
            for i in range(len(self)):
                shift[i] = shift[i]+(hi*i)

        # calib
        if isinstance(calib, Iterable):
            if len(calib) != len(self):
                raise ValueError(f'calib must be a number of a list with length compatible with the number of spectra.\nnumber of calib: {len(calib)}\nnumber of spectra: {len(self)}')
        else:
            calib = [calib]*len(self)

        # factor
        if isinstance(factor, Iterable):
            if len(factor) != len(self):
                raise ValueError(f'factor must be a number of a list with length compatible with the number of spectra.\nnumber of factor: {len(factor)}\nnumber of spectra: {len(self)}')
        else:
            factor = [factor]*len(self)

        # factor
        if isinstance(smooth, Iterable):
            if len(smooth) != len(self):
                raise ValueError(f'smooth must be a number of a list with length compatible with the number of spectra.\nnumber of smooth: {len(factor)}\nnumber of spectra: {len(self)}')
        else:
            smooth = [smooth]*len(self)

        # plot
        # temp = [0]*len(self)
        # for i in range(len(self)-1, -1, -1):
        #     temp[i] = self.data[i].plot(ax=ax, offset=offset[i], shift=shift[i], factor=factor[i], calib=calib[i], **kwargs)

        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self.data[i].plot(ax=ax, offset=offset[i], shift=shift[i], factor=factor[i], calib=calib[i], smooth=smooth[i],**kwargs)

        return temp, offset, shift

    def sequential_plot(self, xlim=None, ylim=None, update_function=None):
        """plot where you can use arrows to flip trhough spectra.

        update_function must be a function of type:

            def update_function(ss):
                plt.title(ss.__i)
                ss[ss.__i].plot(color='black', marker='o')
            
        note that the attr __i can be used to index the spectra.

        """
        self.__i = 0

        # core update function
        if update_function is None:
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                plt.title(self.__i)
                ss[self.__i].plot(color='black', marker='o')
                
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)
        else:
            # add counter and xlim/ylim to update function
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                update_function(ss, __i=self.__i)
                
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'right':
                self.__i = self.__i + 1

                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'left':
                self.__i = self.__i - 1
                
                plt.cla()
                _update(ss)
                plt.draw()

        # mouse events
        def mouse(event):
            """Mouse can be used with a keyboard key"""
            # print('mouse')
            # print(event.key)
            # print(event.button)
            pass

        # plotting
        fig = figure()
        fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
        fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))
        _update(self)
        return

    # modifiers
    def set_shifts(self, value=None, mode=None, type_='relative'):
        """Shift data recursively.

        if value is none, calculated values will be used and type_ will be set to relative.

        Args:
            value (number or list, optional): value will be added to x-coordinates.
                If None, it will look for calculated values from Spectra.calculate_shifts().
            mode (string, optional): If ``mode='x'`` or ``mode='hard'``, y is fully preserved
                while x is shifted. If ``mode='y'``, ``'interp'``, or ``'soft'``, x is preserved
                while y is interpolated with a shift. If ``mode='roll'`` (or rotate or r), x is also preserved
                and y elements are rolled along the array (``shift`` value must be an integer).
                The form of y-coordinates is fully preserved, but the edge of y-coordinates are lost.
            type_ (str, optional): set values 'relative' or in 'absolute' units.
                If 'relative', modifications will be done on top of previous
                modifications. Default is 'relative'.

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
            :py:func:`Spectra.calculate_shifts()`
        """
        # if verbose: print(f'Setting shifts...')

        # if value is None, check calculated_shifts
        if value is None:
            if self.calculated_shifts is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_shifts().')

            value = self.calculated_shifts.y
            type_ = 'relative'

            # if mode is not defined, auto pick the right one
            if mode is None:
                if self.calculated_shifts.mode in cc:
                    mode = 'roll'
                elif self.calculated_shifts.mode == 'max':
                    mode = 'x'
                elif self.calculated_shifts.mode == 'peak' or self.calculated_shifts.mode == 'peaks':
                    mode = 'x'
                elif self.calculated_shifts.mode == 'user-defined':
                    text = 'shift mode not define (please, select "roll", "soft", "hard").' +\
                            'Note that shifts were calculated elsewhere by the user, ' +\
                            'the script can has no way of knowing the most appropriate shift mode.'
                    raise ValueError(text)
                elif self.calculated_shifts.mode in None:
                    raise ValueError('calculated shifts not defined. Use Spectra.calculate_shifts().')

            # if shift were calculated by cc, one should use roll
            if self.calculated_shifts.mode in cc and mode not in roll:
                # one could in principle do a roll even though the calculation mode
                # is cc by multiplying the shift value by the step, however,
                # the whole point of using cc is to be able to do a roll shift.
                # Therefore, when cc is used, roll shift is enforced.
                raise ValueError('shifts were calculated using cross-correlation. Only shift mode possible is "roll".')

            # if shift were calculated by peak or max, one can opt to shift via roll
            # as long as self.step is defined
            if self.calculated_shifts.mode in ['peak', 'peaks', 'max'] and mode in roll:
                value = np.array([int(round(x)) for x in value/self.step])
        else:
            # check if value is a number
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)

            # if mode is not defined, raise an error
            if mode is None:
                raise ValueError(f'Invalid mode. Valid options are `roll`, `x`, `interp`.')

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        # if all shifts are zero, does nothing
        if  all(x==0 for x in value):
            return

        # if mode interp, check monotonicity
        if mode in soft:
            if self.monotonicity is None:
                self.check_monotonicity()
            if self.monotonicity != 'increasing':
                raise ValueError('Arrays must be monotonicaly increasing.\nUse Spectra.fix_monotonicity()')

        # if mode is roll, steps must be defined for all spectra
        if mode in roll:
            if self.step is None:
                try:
                    self.check_step()
                except ValueError:
                    temp = [self[i].step for i in  range(len(self))]
                    raise ValueError(f'Cannot apply roll because some spectra have no step defined.\nData must be homogeneous.\nSteps: {temp}\nUse Spectra.check_step()')

        # set values ===========================================================
        for i in range(len(self)):
            self.data[i].set_shift(value=value[i], mode=mode, type_=type_)

        # if mode not interp and shifts are different, reset x
        if mode not in soft:
            if all_equal(value) == False:
                self._x = None

        # check
        # self._restart_check_attr()
        # self._restart_calculated_modifiers()
        self._calculated_shifts = None
        self._calculated_calib  = None
        # self._calculated_factors = None
        # self._calculated_offsets = None

    def set_factors(self, value=None, type_='relative'):
        """Apply multiplicative y factor recursively.

        if value is none, calculated values will be used and type_ will be set to relative.

        Args:
            value (number or list, optional): value will be multiplied to y-coordinates.
                If None, it will look for calculated values from Spectra.calculate_factors().
            type_ (str, optional): set values 'relative' or in 'absolute' units.
                If 'relative', modifications will be done on top of previous
                modifications. Default is 'relative'.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_factors()`
        """
        if value is None:
            if self.calculated_factors is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_factors().')
            value = self.calculated_factors.y
            type_ = 'relative'
        else:
            # check if value is a number
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        # set values ===========================================================
        for i in range(len(self)):
            self[i].set_factor(value=value[i], type_=type_)

        # check
        # self._restart_check_attr()
        # self._restart_calculated_modifiers()
        # self._calculated_shifts  = None
        # self._calculated_calib  = None
        self._calculated_factors = None
        # self._calculated_offsets = None

    def set_calib(self, value, type_='relative'):
        """Apply multiplicative x factor recursively.

        Args:
            value (number or list): value will be multiplied to x-coordinates.
            type (str, optional): set values 'relative' or in 'absolute' units.
                If 'relative', modifications will be done on top of previous
                modifications. Default is 'relative'.

        Returns:
            None
        """
        # check if value is a number
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        # set values ===========================================================
        for i in range(len(self)):
            self[i].set_calib(value=value[i], type_=type_)

        # check
        # self._restart_check_attr()
        # self._restart_calculated_modifiers()
        self._calculated_shifts  = None
        self._calculated_calib  = None
        # self._calculated_factors = None
        # self._calculated_offsets = None

    def set_offsets(self, value=None, type_='relative'):
        """Apply additive y factor recursively.

        if value is none, calculated values will be used and type will be set to relative.

        Args:
            value (number or list, optional): value will be added to y-coordinates.
                If None, it will look for calculated values from Spectra.calculate_offsets().
            type_ (str, optional): set values 'relative' or in 'absolute' units.
                If 'relative', modifications will be done on top of previous
                modifications. Default is 'relative'.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_offsets()`, :py:func:`Spectra.floor()`
        """
        if value is None:
            if self.calculated_offsets is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_offsets().')
            value = self.calculated_offsets.y
            type_ = 'relative'
        else:
            # check if value is a number
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'


        # set values ===========================================================
        for i in range(len(self)):
            self[i].set_offset(value=value[i], type_=type_)

        # check
        # self._restart_check_attr()
        # self._restart_calculated_modifiers()
        # self._calculated_shifts  = None
        # self._calculated_calib   = None
        # self._calculated_factors = None
        self._calculated_offsets = None


    def align(self, mode='cc', shift_mode=None, ref_spectrum=0, ref_value=None, ref_peak=0, ranges=None):
        """Uses calculate_shifts via cross-correlation and align spectra.

        Args:
            bkg_check (bool, optional): if mode is 'cc', bkg_check=true prevents
                from cross correlating spectra with significant background (10 % of
                the maximum y) which can lead to wrong results.
                See :py:func:`Spectra.calculate_shifts()` for more info.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_shifts`
        """
        self.calculate_shifts(mode=mode, ref_spectrum=ref_spectrum, ref_value=ref_value, ref_peak=ref_peak, ranges=ranges)
        self.set_shifts(mode=shift_mode)

    def normalize(self, mode='cross-correlation', ref_spectrum=None, ref_value=None, ref_peak=0, ranges=None):
        """Uses Spectra.calculate_factors() and Spectra.set_factors() to normalize spectra.

        Args:
            mode (string, optional): method used to calculate the shifts. If
                None, it will try to use 'delta'.
            bkg_check (bool, optional): if mode is 'area', bkg_check=true prevents
                from comparing spectra with significant background (10 % of
                the maximum y) which can lead to wrong results.
                The algorithm uses find_peaks to find
                approximately whats is data and what is background, Then the
                background is said to be big if average(bkg) > max(y)*0.1. If
                background region cannot be found clearly, the criteria
                becames average(y) > max(y)*0.1.
            peak (int, optional): if mode='peak' or mode='fitted peaks', this variable selects which peak
                to align data.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_factors`
        """
        if mode is None:
            mode = 'delta'
        self.calculate_factors(mode=mode, ref_spectrum=ref_spectrum, ref_value=ref_value, ref_peak=ref_peak, ranges=ranges)
        self.set_factors()

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
        # reset attributes =====================================================
        # self._restart_check_attr()
        # self._calculated_shifts  = None
        # self._calculated_calib  = None
        self._calculated_factors = None
        self._calculated_offsets = None

    def flip(self):
        """Flip x axis.

        Returns:
            None
        """
        for s in self:
            s.flip()
        self._restart_check_attr()
        self._calculated_shifts  = None
        self._calculated_calib  = None

    # Calculate modifiers
    def calculate_shifts(self, mode='cross-correlation', ref_spectrum=0, ref_value=None, ref_peak=0, ranges=None):
        """Calculate how much spectra must be shifted to align with the first spectrum.

        Args:
            mode (string, optional): method used to calculate the shifts.
                The current options are: 'cross-correlation' ('cc'), 'max',
                'fitted peaks', or 'peak'. For mode='peak' or 
                mode='fitted peaks', data must be fitted via
                Spectrum.fit_peak(), i. e.,
                spectra must have a Spectrum.fit.peaks object defined.
            ref (int, optional): index of the spectrum to be taken as reference. 
                Default is 0.
            ref_value (float, optional): if not None, spectra will be adjusted 
                to ref_value. Mode 'cc' does not work with ref_value and only 
                    ref must be used). ref_value overwrites ref, unless 
                    mode='cc'. Default is None.
            ref_peak (int, optional): if mode='peaks' (or any peak related mode)
                this variable selects the reference peak to be used. Default is
                0.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                The additive factor will be calculated for spectra such the average value of
                data inside ranges is the same as for the first spectrum.

        Returns:
            None

        Note:
            For ``mode = 'cc'``, spectra must have the same x-coordinates (this
            is checked before execution).

        See Also:
            :py:func:`Spectra.set_shifts`
        """
        # Initialization
        values = np.array([0.0]*len(self))

        # CALCULATION ==========================================================
        # cc
        if mode in cc:
            # ranges
            if ranges is not None:
                ss = self.extract(ranges)
            else:
                ss = self

            # x must be the same for cc
            self.check_same_x()

            # x must be uniform (same step between each data point)
            self.check_step()

            # calculate cross-correlation
            for i in range(len(self)):
                cross_correlation = np.correlate(self[ref_spectrum].y, self[i].y, mode='full')
                values[i] = np.argmax(cross_correlation)
            ref_value = values[ref_spectrum]
            values -= ref_value
        # max
        elif mode == 'max':
            # ranges
            if ranges is not None:
                ss = self.extract(ranges)
            else:
                ss = self
            # calculation
            if ref_value is None:
                j_ref = np.argmax(self[ref_spectrum].y)
                ref_value = self[ref_spectrum].x[j_ref]
            for i in range(len(self)):
                j = np.argmax(self[i].y)
                values[i] = -(self[i].x[j] - ref_value)
        # fitted peaks
        # elif mode == 'fitted peaks':
        #     # check if all data has fit
        #     self._check_fit()
        #     # check if number of peaks for each spectrum is the same
        #     self._check_number_of_peaks_per_spectrum(fitted_peaks=True)

        #     # check if peaks are defined
        #     assert len(self[0].fit.peaks) > 0, 'Spectra does not have defined fitted peaks.\nMaybe use Spectra.fit_peaks(), Spectra.find_peaks().'

        #     # get peaks center
        #     peaks     = self.fit.peaks
        #     if ref_value is None:
        #         ref_value = peaks['c'][ref_peak][ref_spectrum]
        #     values    = peaks['c'][ref_peak]
        #     if None in values:
        #         raise ValueError(f'cannot calculate shifts.\npeak {ref_peak} is not defined for all spectra.\ncenter of peak {ref_peak}: {values}')
        #     values    = ref_value-np.array(values)
        # peaks
        elif mode in ['peaks', 'peak']:
            # check if number of peaks for each spectrum is the same
            # self._check_number_of_peaks_per_spectrum()
            # check if peaks are defined
            # assert len(self[0].peaks) > 0, 'Spectra does not have defined peaks.\nMaybe use Spectra.find_peaks().'

            # check peaks
            n_peaks = [0]*len(self)
            for i in range(len(self)):
                n_peaks[i] = len(self[i].peaks)
            
            assert (0 not in n_peaks), f'Some spectra do not have peaks defined.\nn_peaks={n_peaks}'


            # calculate
            if ref_value is None:
                ref_value = self[ref_spectrum].peaks[ref_peak]['c']
            values = self.peaks['c'][ref_peak]
            if None in values:
                raise ValueError(f'cannot calculate shifts.\npeak {ref_peak} is not defined for all spectra.\ncenter of peak {ref_peak}: {values}')
            values = ref_value - np.array(values)
        else:
            raise ValueError('mode not valid.\nValid modes: cross-correlation, max, peaks.')

        # save calculated values ===============================================
        self._calculated_shifts              = Spectrum(y=values)
        # self._calculated_shifts              = values
        self._calculated_shifts.ref_spectrum = ref_spectrum
        self._calculated_shifts.mode         = mode
        self._calculated_shifts.ref_value    = ref_value

    def calculate_factors(self, mode='fitted peaks', ref_spectrum=0, ref_value=None, ref_peak=0, ranges=None):
        """Calculate mult. factor for spectra to be same height as the first spectrum.

        ranges only work for max and delta (for now).

        Args:
            mode (string, optional): method used to calculate the multiplicative factors.
                The current options are: 'max', 'delta', 'area', 'fitted peaks',
                'fitted peaks area', 'peaks', or 'peaks area'.
                For peak related modes, data must be fitted via
                Spectrum.fit_peak(), i. e., spectra must have a
                Spectrum.fit.peaks object defined.
            ref_spectrum (int, optional): index of the spectrum to be taken as reference. 
                Default is 0.
            ref_value (float, optional): if not None, spectra will be adjusted 
                to ref_value. ref_value overwrites ref. Default is None.
            ref_peak (int, optional): if mode='peaks' (or any peak related mode)
                this variable selects the reference peak to be used. Default is
                0.

        Returns:
            None

        See Also:
            :py:func:`Spectra.set_factor`
        """
        # Initialization
        values = np.array([0.0]*len(self))
        
        # max
        if mode == 'max':
            # ranges
            if ranges is not None:
                ss = self.extract(ranges)
            else:
                ss = self
            # calculation
            if ref_value is None:
                ref_value = max(ss.data[ref_spectrum].y)
            for i in range(len(ss)):
                values[i] = ref_value/max(ss.data[i].y)
        # delta
        elif mode == 'delta':
            # ranges
            if ranges is not None:
                ss = self.extract(ranges)
            else:
                ss = self
            # calculation
            if ref_value is None:
                ref_value = max(ss.data[ref_spectrum].y) - min(ss.data[ref_spectrum].y)
            for i in range(len(ss)):
                values[i] = ref_value/(max(ss.data[i].y) - min(ss.data[i].y))
        # area
        elif mode == 'area':
            # ranges
            if ranges is not None:
                ss = self.extract(ranges)
            else:
                ss = self
            # calculation
            if ref_value is None:
                ref_value = self.data[ref_spectrum].area
            for i in range(len(self)):
                values[i] = ref_value/self.data[i].area
        # fitted peaks
        elif mode == 'fitted peaks' or mode == 'fitted peak':
            # check if all data has fit
            self._check_fit()
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum(fitted_peaks=True)
            # check if peaks are defined
            assert len(self[0].fit.peaks) > 0, 'Spectra does not have defined fitted peaks.\nMaybe use Spectra.fit_peaks(), Spectra.find_peaks().'

            # calculate peak
            peaks     = self.fit.peaks
            values    = peaks['amp'][ref_peak]
            if ref_value is None:
                ref_value = peaks['amp'][ref_peak][ref_spectrum]
            if None in values:
                raise ValueError(f'cannot calculate multiplicative factors.\npeak {ref_peak} is not defined for all spectra.\ncenter of peak {ref_peak}: {values}')
            values    = ref_value/np.array(values)
        # fitted peaks area
        elif mode == 'fitted peaks area':
            # check if all data has fit
            self._check_fit()
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum(fitted_peaks=True)
            # check if peaks are defined
            assert len(self[0].fit.peaks) > 0, 'Spectra does not have defined fitted peaks.\nMaybe use Spectra.fit_peaks(), Spectra.find_peaks().'

            # calculate peak
            peaks     = self.fit.peaks
            values    = peaks['area'][ref_peak]
            if ref_value is None:
                ref_value = peaks['area'][ref_peak][ref_spectrum]
            if None in values:
                raise ValueError(f'cannot calculate multiplicative factors.\npeak {ref_peak} is not defined for all spectra.\ncenter of peak {ref_peak}: {values}')
            values    = ref_value/np.array(values)
        # peaks
        elif mode == 'peaks' or mode == 'peak':
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum()
            # check if peaks are defined
            assert len(self[0].peaks) > 0, 'Spectra does not have defined peaks.\nMaybe use Spectra.find_peaks().'

            # calculate peak
            peaks     = self.peaks
            values    = peaks['amp'][ref_peak]
            if ref_value is None:
                ref_value = peaks['amp'][ref_peak][ref_spectrum]
            if None in values:
                raise ValueError(f'cannot calculate multiplicative factors.\npeak {ref_peak} is not defined for all spectra.\ncenter of peak {ref_peak}: {values}')
            values    = ref_value/np.array(values)
        # peaks area
        elif mode == 'peaks area' or mode == 'peak area':
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum()
            # check if peaks are defined
            assert len(self[0].peaks) > 0, 'Spectra does not have defined peaks.\nMaybe use Spectra.find_peaks().'

            # calculate peak
            peaks     = self.peaks
            values    = peaks['area'][ref_peak]
            if ref_value is None:
                ref_value = peaks['area'][ref_peak][ref_spectrum]
            if None in values:
                raise ValueError(f'cannot calculate multiplicative factors.\npeak {ref_peak} is not defined for all spectra.\ncenter of peak {ref_peak}: {values}')
            values    = ref_value/np.array(values)
        else:
            raise ValueError('mode not valid.\nValid modes: max, delta, area, fitted peaks (or fitted peak), fitted peaks area, peak, peak area.')

        # save calculated values ===============================================
        self._calculated_factors           = Spectrum(y=values)
        self._calculated_factors.mode      = mode
        self._calculated_factors.ref_spectrum = ref_spectrum
        self._calculated_factors.ref_value    = ref_value

    def calculate_offsets(self, ref_spectrum=0, ref_value=None, ranges=None):
        """Average data inside range and calculate additive factor to verticaly align data.

        Args:
            ref_spectrum (int, optional): index of the spectrum to be taken as reference. 
                Default is 0.
            ref_value (float, optional): if not None, spectra will be adjusted 
                to ref_value. ref_value overwrites ref. Default is None.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                The additive factor will be calculated for spectra such the average value of
                data inside ranges is the same as for the first spectrum.

        Returns:
            None

        See Also:
            :py:func:`Spectra.set_offset`, :py:func:`Spectra.floor()`,
        """
        # Initialization
        values = np.array([0.0]*len(self))

        # CALCULATION 
        if ref_value is None:
            s_ref = self[ref_spectrum].extract(ranges)
            ref_value = np.mean(s_ref.y)

        for i in range(len(self)):
            s = self[i].extract(ranges)
            values[i] = ref_value-np.mean(s.y)

        # save calculated values ===============================================
        self._calculated_offsets              = Spectrum(y=values)
        self._calculated_offsets.ref_spectrum = ref_spectrum
        self._calculated_offsets.ref_value    = ref_value

    def calculate_calib(self, start_value=None, stop_value=None, values=None, mode='cross-correlation', ref_peak=0, deg=1):
        """Calculate calibration factor via :py:func:`Spectra.calculate_shifts()`.

        The calibration factor is the shift value (x array) as a function of the
            values array.

        Args:
            start_value (number): value used for measuring the first spectrum.
            stop_value (number): value used for measuring the last spectrum.
            values (list): value list. It overwrites start and stop. Must be the
                same length as number of spectra.
            mode (string, optional): method used to calculate the shifts.
                Default is 'cross-correlation'. This
                parameter is passed to :py:func:`Spectra.calculate_shifts()`.
                see :py:func:`Spectra.calculate_shifts()` for more.
            bkg_check (bool, optional): see :py:func:`Spectra.calculate_shifts()` for more.
            ref_peak (int, optional): see :py:func:`Spectra.calculate_shifts()` for more.
            deg (int, ooptional): degree of the fitting polynomial. Default is 1.

        Returns:
            calibration factor (number)

            ss.calculated_calib becames type spectrum with fit results
            ss.calculated_calib.mode
            ss.calculated_calib.popt
            ss.calculated_calib.model
        """
        # check number of values matches the numbre of spectra
        if values is None:
            values = np.linspace(start_value, stop_value, len(self))
        if len(self) != len(values):
            raise ValueError(f'number of values ({len(values)}) do not match the number of spectra ({len(self)})')

        # CALCULATION ==========================================================
        self.calculate_shifts(mode=mode, ref_peak=ref_peak)
        centers = -(self.calculated_shifts.y + self.calculated_shifts.ref_value)

        if mode in cc:
            centers = centers*self.step

        # save calculated values ===============================================
        self._calculated_calib      = Spectrum(x=values, y=centers)
        self._calculated_calib.mode = mode
        # self._calculated_calib.ref_value = ref_value
        popt, model, R2 = self.calculated_calib.polyfit(deg=deg)
        self._calculated_calib.popt = popt
        self._calculated_calib.model = model
        self._calculated_calib.R2 = R2

        return 1/popt[0]

    # Extractors
    def _extract(self, ranges):
        """Same as extract(), but attributes are not copied to the new object."""

        if isinstance(ranges, Iterable):
            temp = Spectra(n=len(self))
            for i, s in enumerate(self.data):
                temp[i] = s.extract(ranges=ranges)
        else:
            raise ValueError(f'Ranges is not a valid value.\nRanges should be a pair (x_init1, x_final1) or a list of pairs like this: ((x_init1, x_final1), (x_init2, x_final2), ...)\nUse None to indicate the minimum or maximum x value of the data.')
        return temp

    def extract(self, ranges):
        """Return the extracted data range from full data.

        Args:
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            :py:attr:`Spectra`
        """
        temp = self._extract(ranges)
        return self._transfer_attributes(temp)

    def concatenate(self):
        """Return spectrum of concatenate spectra.
        
        User defined attr are saved to the summed spectrum. A number is 
            added in front of the attr name indicating the index of the 
            spectrum in which the attr was copied from.
        """
        x = np.concatenate([s.x for s in self.data])
        y = np.concatenate([s.y for s in self.data])
        s = Spectrum(x=x, y=y)
        s = self._transfer_attributes(s)
        return self._transfer_attributes_from_each_spectrum(s)

    def calculate_sum(self, ranges=None):
        """Returns Spectrum object with the sum of all spectra.

        Returns:
            :py:class:`Spectra` object.

        Note:
            All spectra must have the same x-coordinates. This is verified
            before summing up the spectra.

            User defined attr are saved to the summed spectrum. A number is 
            added in front of the attr name indicating the index of the 
            spectrum in which the attr was copied from.
        """
        # gather ys
        x, ys = self._gather_ys(ranges=ranges)

        # calculate sum
        y = np.zeros(len(x))
        for i in range(len(self)):
            y += ys[:, i]

        s = Spectrum(x=x, y=y)
        s = self._transfer_attributes(s)
        return self._transfer_attributes_from_each_spectrum(s)

    def calculate_map(self, axis=0, ranges=None):
        """Return image.

        Args:
            axis (int or string, optional): Axis along which x axis will be laid
                down. By default, x axis will run across the vertical (0) direction.

        Returns:
            :py:class:`Image`.
        """
        axis = _axis_interpreter(axis)

        # gather ys
        y, ys = self._gather_ys(ranges=ranges)
        x = None

        if axis == 1:
            ys = ys.transpose()
            x = y
            y = None

        im = Image(data=ys)
        im.x_centers = x
        im.y_centers = y

        im = self._transfer_attributes(im)
        return self._transfer_attributes_from_each_spectrum(im)

    def polyfit(self, deg=2, ranges=None):
        """Fit data recursively with a polynomial. Wrapper for `numpy.polyfit()`_.

        Args:
            deg (int, optional): degree of the fitting polynomial. Default is 2.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.
                
        Returns:
            list with polynomial coefficients, highest power first.
            list with Model function f(x).
            list with R2.

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        popt = [0]*len(self)
        model = [0]*len(self)
        R2 = [0]*len(self)
        for i in range(len(self)):
            popt[i], model[i], R2[i] = self[i].polyfit(deg=deg, ranges=ranges)

        return popt, model, R2

    # peaks
    def find_peaks(self, prominence=None, width=4, moving_average_window=8, ranges=None):
        """Find peaks recursively. Wrapper for `scipy.signal.find_peaks()`_.

        Args:
            prominence (number, optional): minimum prominence of peaks in percentage
                of the maximum prominence [max(y) - min(y)]. Default is 5.
            width (number, optional): minimum number of data points defining a peak.
            moving_average_window (int, optional): window size for smoothing the
                data for finding peaks. Default is 4.

        Returns:
            None

        .. _scipy.signal.find_peaks(): https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
        """
        for s in self.data:
            s.find_peaks(prominence=prominence, width=width, moving_average_window=moving_average_window, ranges=ranges)

    def fit_peak(self, asymmetry=None, ranges=None, verbose=False):
        if verbose: print('Fitting\n')
        for i in range(len(self)):
            # if verbose: print(f'spectrum {i}')
            self[i].fit_peak(asymmetry=asymmetry, ranges=ranges)
            if verbose: print(f'spectrum {i}: done')
            # if verbose: print('='*20)
            # if verbose: print('\n')

    def fit_peaks(self, method='least_squares', ranges=None, verbose=False):
        """Fit peaks recursively. Wrapper for `scipy.optimize.curve_fit()`_.

        Args:
            asymmetry (bool or dict, optional): if True, fits each half of the
                with a different width.
            fixed_m (False, number, or dict, optional): m is the amount of lorentzian
                contribution for a peak. If False, m will be fitted for each peak.
                If a number (from 0 to 1), this will be used as the value of m
                (fixed m).
            offset (bool, optional): if True, a offset value will be fitted.
            amp_bounds, c_bounds, fwhm_bounds (tuple, optional): minimum and
                maximum multiplication factor for boundary values. For amp, the
                bounds are set between amp*amp_bounds[0] and amp*amp_bounds[1].
                For c, bounds are set c+fwhm*c_bound[0] and c+fwhm*c_bound[1].
                Finaly, for fwhm, the bounds are set between fwhm1*fwhm_bounds[0]
                and fwhm1*fwhm_bounds[0]. Note that if fwhm1*fwhm_bounds[0] is
                less than zero, the minimum fwhm boundary is set to zero as it
                cannot be negative.


        Returns:
            None

        .. _scipy.optimize.curve_fit(): https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
        """
        if verbose: print('Fitting.\n')
        for i in range(len(self)):
            # if verbose: print(f'spectrum {i} ===============================')
            self[i].fit_peaks(method=method, ranges=ranges)
            if verbose: print(f'spectrum {i}: done')
            # if verbose: print('='*20)
            # if verbose: print('\n')




    # OBSOLETE
    def _check_fit(self):
        """Check if all data has fit object. Raises an error."""
        temp = [None for i in range(len(self))]
        for i in range(len(self)):
                temp[i] = self[i].fit
        if None in temp:
            raise RuntimeError(f'Fit is not defined for some spectra.\nList of fits: {temp}')

    def _check_number_of_peaks_per_spectrum(self, fitted_peaks=False):
        """Check if number of peaks for each spectrum is the same."""
        if fitted_peaks:
            peaks = self.fit.peaks
        else:
            peaks = self.peaks

        temp = [len(peaks[i]) for i in range(len(peaks))]
        if np.all(np.diff(temp) != 0) == True:
            temp = str({i: temp[i] for i in range(len(temp))}).replace(', ', '\n').replace('{', '').replace('}', '')
            print(f'number of peaks is different for different spectrum.\nCalculated values based on peaks might be wrong.\nNumber of peaks per spectrum:\n{temp}')

    def get_errors(self):
        """Return a list with peak errors data from all spectra.
    
        Returns:
            list
        """
        peaks_attr = ['amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2', 'area']
        n_peaks = max([len(s.peaks) for s in self])
        if n_peaks == 0:
            return {}
        return [{k: [s.peaks[peak].error[k] if peak in s.peaks else None for s in self] for k in peaks_attr} for peak in range(n_peaks)]

    def get_peaks(self, key='peak'):
        """Return a list with peak data from all spectra.
    
        Args:
            key (string, optional): if 'peaks', list is organized by peak
                indices. If 'spectrum', list is organized by spectrum
                indices.
    
        Returns:
            list
        """
        if key.startswith('p'):
            peaks_attr = ['amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2', 'area']
            n_peaks = max([len(s.peaks) for s in self])
            if n_peaks == 0:
                return {}
            return [{k: [s.peaks[peak][k] if peak in s.peaks else None for s in self] for k in peaks_attr} for peak in range(n_peaks)]
        elif key.startswith('s'):
            return [self[i].peaks for i in range(len(self))]
        else:
            raise ValueError(f"Ivalid key.\nPossible keys are: 'peak' (p) and 'spectrum' (s).")

    def plot_peaks(self, ax=None, offset=0, shift=0, factor=1, calib=1, **kwargs):
        """Place a marker at the maximum of every peak positionfor all spectra.
    
        Wrapper for `matplotlib.pyplot.errorbar()`_.
    
        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number or list, optional): defines a vertical offset. Default is 0.
            shift (number or list, optional): horizontal shift value. Default is 0.
            factor (number or list, optional): multiplicative factor on the y axis.
                Default is 1.
            calib (number or list, optional): multiplicative factor on the x axis.
                Default is 1.
            **kwargs: kwargs are passed to `matplotlib.pyplot.errorbar()`_ that plots the data.
    
        Returns:
            list with dicts with `ErrorbarContainer`_ for each spectrum.
    
        .. matplotlib.pyplot.errorbar(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
        .. ErrorbarContainer: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
        """
        if ax is None:
            ax = plt
    
        # percentage wise increment ====================
        if 'vi' in kwargs and 'vertical_increment' in kwargs:
            raise SyntaxError('keyword argument repeated: vertical increment/vi')
        elif 'vi' in kwargs or 'vertical_increment' in kwargs:
            if 'vi' in kwargs:
                vi = kwargs['vi']
                del kwargs['vi']
            if 'vertical_increment' in kwargs:
                vi = kwargs['vertical_increment']
                del kwargs['vertical_increment']
            temp = [0]*len(self)
            for i in range(len(self)):
                temp[i] = max(self.data[i].y) - min(self.data[i].y)
            vi = max(temp)*factor*vi/100
        else:
            vi = 0
    
        # percentage wise horizontal increment ====================
        if 'hi' in kwargs and 'horizontal_increment' in kwargs:
            raise SyntaxError('keyword argument repeated: horizontal increment and hi')
        elif 'hi' in kwargs or 'horizontal_increment' in kwargs:
            if 'hi' in kwargs:
                hi = kwargs['hi']
                del kwargs['hi']
            if 'horizontal_increment' in kwargs:
                hi = kwargs['horizontal_increment']
                del kwargs['horizontal_increment']
            temp = [0]*len(self)
            for i in range(len(self)):
                temp[i] = max(self.data[i].x) - min(self.data[i].x)
            hi = max(temp)*factor*hi/100
        else:
            hi = 0
    
        # offset
        if isinstance(offset, Iterable):
            if len(offset) != len(self):
                raise ValueError(f'offset must be a number of a list with length compatible with the number of spectra.\nnumber of offsets: {len(offset)}\nnumber of spectra: {len(self)}')
        else:
            offset = [offset]*len(self)
            for i in range(len(self)):
                offset[i] = offset[i]+(vi*i)
    
        # shift
        if isinstance(shift, Iterable):
            if len(shift) != len(self):
                raise ValueError(f'shift must be a number of a list with length compatible with the number of spectra.\nnumber of shift: {len(shift)}\nnumber of spectra: {len(self)}')
        else:
            shift = [shift]*len(self)
            for i in range(len(self)):
                shift[i] = shift[i]+(hi*i)
    
        # calib
        if isinstance(calib, Iterable):
            if len(calib) != len(self):
                raise ValueError(f'calib must be a number of a list with length compatible with the number of spectra.\nnumber of calib: {len(calib)}\nnumber of spectra: {len(self)}')
        else:
            calib = [calib]*len(self)
    
        # factor
        if isinstance(factor, Iterable):
            if len(factor) != len(self):
                raise ValueError(f'factor must be a number of a list with length compatible with the number of spectra.\nnumber of factor: {len(factor)}\nnumber of spectra: {len(self)}')
        else:
            factor = [factor]*len(self)
    
        # get peaks ======================================
        peaks = self.get_peaks(key='spectrum')
    
        # plot
        r = [0]*len(self)
        for i in range(len(self)):
            r[i] = peaks[i].plot(ax=ax, offset=offset[i], shift=shift[i], factor=factor[i], calib=calib[i], **kwargs)
        return r

    def load2(self, dirpath, string='*', verbose=False, **kwargs):
        """Load data from text files. Wrapper for `numpy.genfromtxt()`_.

        Args:
            dirpath (string, path object, or list): folderpath, folder handle,
                or a list of filepaths. All files inside a folder where string 
                can be found in the filename are imported. If the filename 
                extension is .gz or .bz2, the file is first decompressed. 
                Last used dirpath is saved to ss.dirpath.
            string (str, optional): file names without this string will be ignored.
                Use '*' for matching anything. Default is '*'.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``# ``. Attributes picked up
                from the header will be loaded too.

        Returns:
            None

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        # reset data
        self._data     = []
        self._restart_check_attr()
        self._calculated_shifts  = None
        self._calculated_calib  = None
        self._calculated_factors = None
        self._calculated_offsets = None

        # if dirpath is str
        if isinstance(dirpath, str) or isinstance(dirpath, Path):
            dirpath = Path(dirpath)
            if dirpath.is_file():
                if verbose: print('dirpath is a file')
                if verbose: print('Loading...')
                self.append(Spectrum(filepath=dirpath))
                if verbose: print('Done!')
                return
            elif dirpath.is_dir():
                dirpath = [dirpath, ]
            else:
                raise ValueError(f'cannot read dirpath.\ndirpath: {dirpath}')

        # if dirpath is iterable
        if isinstance(dirpath, Iterable):
            if verbose: print('Loading...')
            for j, filepath in enumerate(dirpath):
                if verbose: print(f'{j+1}/{len(dirpath)}: {filepath}')

                # if Path(filepath).is_dir():
                #     fl = filelist(dirpath=filepath, string=string)
                #     for i, f in enumerate(fl):
                #         if verbose: print(f'    {i+1}/{len(fl)}: {f}')
                #         self.append(Spectrum(filepath=f))

                elif Path(filepath).is_file():
                    self.append(Spectrum(filepath=filepath))
                else:
                    raise ValueError(f'cannot read dirpath.\dirpath: {dirpath}')
        
        # save dirpath
        self.dirpath = dirpath

        if verbose: print('Done!')
# %%
