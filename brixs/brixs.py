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
# """

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
from collections.abc import Iterable, MutableMapping
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

# backpack
from .backpack.filemanip import load_Comments, filelist
from .backpack.arraymanip import index, moving_average, extract, shifted, sort, get_attributes
from .backpack.arraymanip import is_integer, all_equal, factors, flatten, peak_fit, check_monotonicity, fix_monotinicity
from .backpack.figmanip import n_digits, n_decimal_places, figure, set_window_position
from .peaks import Peak, Peaks, _peak_function_creator
from .backpack.model_functions import voigt_fwhm, dirac_delta

# common definitions ===========================================================
# from .__init__ import settings
from .config import settings

cc = ['cross-correlation', 'cc']
roll = ['roll', 'rotate', 'r', 'rot']
hard = ['hard', 'x', 'h', 'Hard']
soft = ['soft', 'Soft', 'interp', 'y', 's']

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

    axis = c, col, column, columns, hor, horizontal, x, ... will return 1
    axis = r, row, rows, v, ver, vertical, y, ... will return 0

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
    elif type(axis)==str and (axis=='x' or axis.startswith('c') or axis.startswith('h')):
        return 1
    elif axis == 0:
        return 0
    elif type(axis)==str and (axis=='y' or axis.startswith('r') or axis.startswith('v')):
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
        if row_size is None:    row_size =    shape[1]
        if column_size is None: column_size = shape[0]
        if row_size <= 0 or column_size <= 0 or is_integer(row_size)==False or is_integer(column_size)==False:
            raise ValueError("Size of bins must be a positive integer.")
        else:
            if factor_check:
                assert shape[1] % column_size == 0, f"The {shape[1]} pixels in a row is not evenly divisible by {column_size}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[1])))}"
                assert shape[0] % row_size == 0,    f"The {shape[0]} pixels in a column is not evenly divisible by {row_size}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[0])))}"

            nbins = np.array((round(shape[0]/row_size), round(shape[1]/column_size)))
            # bins_size = np.array((shape[0]//_nbins[0], shape[1]//_nbins[1]))
            bins_size = np.array((shape[0]/_nbins[0], shape[1]/_nbins[1]))

    return nbins, bins_size

# BRIXS ========================================================================
class Image(metaclass=_Meta):
    """Image object.

    Args:
        data (2D array, optional): Image.
        filepath (string or path object, optional): filename or file handle.
            If the filename ends in .gz, the file is automatically saved in
            compressed gzip format. This is overwriten by data.

    Attributes:
        data (2D array): This is where we store the Image.
        x, y (1D array): x and y axis values.
        vmin, vmax (number): Minimum value in data.
        shape (tuple): Shape of data (vertical size, horizontal size).
        histogram (brixs.Spectrum): Data intensity histogram.

        nbins (tuple): Number of bins (number of rows, number of columns).
        bins_size (tuple): Bins size (size of rows, size of columns).
        reduced (brixs.Image): Binned image.

        calculated_shift (brixs.Spectrum): Calculated shifts.
        shifts_v, shifts_h (1D array): Shift values in the vertical and
            horizontal direction.

        spectrum_v, spectrum_h (brixs.Spectrum): Spectrum obtained by integrating
            pixels in the vertical and horizontal direction.
        columns, rows (brixs.Spectra): Spectra obtained from each pixel columns
            or row.

    Methods:
        save()
        load()
        plot()
        imshow()
        binning()
        calculate_histogram()
        calculate_spectrum()
        floor()
        calculate_shift()
        set_shifts()
        fix_curvature()


    """
    _read_only = ['shape', 'vmin', 'vmax', 'reduced', 'calculated_shift', 'p', 'f']

    def __init__(self, *args, **kwargs):
        # argument parsing
        data, filepath = self._sort_args(args, kwargs)

        # besic attr
        self._data = None
        self._vmin = None
        self._vmax = None
        self._x_centers = None
        self._y_centers = None
        self._x_edges = None
        self._y_edges = None
        self._shape = None

        # binning attr
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None

        # shifts
        self._calculated_shift = None
        self._shifts_v = None
        self._shifts_h = None
        self._p = None
        self._f = None

        # set data
        if data is not None:
            self.data = copy.deepcopy(data)
        elif filepath is not None:
            self.load(filepath)

    @property
    def data(self):
        return copy.deepcopy(self._data)
    @data.setter
    def data(self, value):
        # basic attr
        self._data  = np.array(value)
        self._vmin  = min([min(x) for x in self.data])
        self._vmax  = max([max(x) for x in self.data])
        self._shape = (self.data.shape[0], self.data.shape[1])
        self.x_centers = None
        self.y_centers = None
        # self._x_edges = None
        # self._y_edges = None

        # binning attr
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None

        # shift attr
        self._calculated_shift = None
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
        self._x_centers = np.array(value)
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
        self._y_centers = np.array(value)
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
        self._x_edges   = np.array(value)
        centers = moving_average(value, n=2)
        # self._x_centers = Spectrum(x=centers, y=np.zeros(len(centers)))
        self._x_centers = np.array(centers)
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
        self._y_edges   = np.array(value)
        centers = moving_average(value, n=2)
        # self._y_centers = Spectrum(x=centers, y=np.zeros(len(centers)))
        self._y_centers = np.array(centers)
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
        binning(self, bins_size=value)
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

    def __len__(self):
        if self._data is None:
            return 0
        else:
            return len(self._data)

    def save(self, filepath, only_data=False,  **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.

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
            kwargs['header'] += '==== Image attributes ===='  + '\n'
            for n in dict:
                if n not in ['_data', '_reduced', '_calculated_shift', '_histogram', '_spectrum_h', '_spectrum_v',]:
                    if isinstance(dict[n], Iterable):
                        kwargs['header'] += f'{n}: {list(dict[n])}'  + '\n'
                    else:
                        kwargs['header'] += f'{n}: {dict[n]}'  + '\n'
            np.savetxt(Path(filepath), self._data, **kwargs)

    def load(self, filepath, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first decompressed.

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
            if '==== Image attributes ====' in line:
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


    def fastplot(self, ax=None, colorbar=False, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.pcolorfast()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
            **kwargs: kwargs are passed to `matplotlib.pyplot.pcolorfast()`_.

        If not specified, the following parameters are passed to `matplotlib.pyplot.pcolorfast()`_:

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
            `matplotlib.image.AxesImage`_

        .. _matplotlib.pyplot.pcolorfast(): https://matplotlib.org/stable/api/_as_gen/matplotlib.axes.Axes.pcolorfast.html
        .. _matplotlib.image.AxesImage: https://matplotlib.org/3.5.0/api/image_api.html#matplotlib.image.AxesImage
        """
        # initialization
        if ax is None:
            ax = plt
            if settings.ALWAYS_PLOT_NEW_WINDOW:
                figure()
                if  settings.FIGURE_POSITION is not None:
                    try:
                        set_window_position(settings.FIGURE_POSITION)
                    except:
                        pass

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

        # # x and y
        # assert check_monotonicity(self.x) == 1, f'x axis (Image.x) must be increasingly monotonic.\nData cannot be plotted by Image.fastplot().\nPlease, use Image.plot() or set Image.x = None\nx: {self.x}'
        # assert check_monotonicity(self.y) == 1, f'y axis (Image.y) must be increasingly monotonic.\nData cannot be plotted by Image.fastplot().\nPlease, use Image.plot() or set Image.y = None\ny: {self.y}'
        # x = np.linspace(self.x[0], self.x[-1], len(self.x))
        # y = np.linspace(self.y[0], self.y[-1], len(self.y))

        # fix monotonic of labels x
        if check_monotonicity(self.x_centers) != 1:
            x, ordering = fix_monotinicity(self.x_centers, np.arange(len(self.x_centers)), mode='increasing')
            assert len(x)==len(self.x_centers), f'Cannot plot when Image.x have repeated elements.\nEither fix Image.x or set it to None.\nx: {self.x_centers}'
            ordering = [int(i) for i in ordering]
            data = copy.deepcopy(self.data)
            for i in ordering:
                if i != ordering[i]:
                    data[:, i] = self.data[:, ordering[i]]
        else:
            x = self.x_centers
            data = self.data

        # fix monotonic of labels y
        if check_monotonicity(self.y_centers) != 1:
            y, ordering = fix_monotinicity(self.y_centers, np.arange(len(self.y_centers)), mode='increasing')
            assert len(y)==len(self.y_centers), f'Cannot plot when Image.y have repeated elements.\nEither fix Image.x or set it to None.\ny: {self.y_centers}'
            ordering = [int(i) for i in ordering]
            data2 = copy.deepcopy(data)
            for i in ordering:
                if i != ordering[i]:
                    data2[i, :] = data[ordering[i], :]
        else:
            y = self.y_centers
            data2 = data

        # plot
        x = np.linspace(x[0], x[-1], len(x))
        dx  = np.mean(np.diff(x))
        x = np.linspace(x[0]-dx/2, x[-1]+dx/2, len(x))

        y = np.linspace(y[0], y[-1], len(y))
        dy  = np.mean(np.diff(y))
        y = np.linspace(y[0]-dy/2, y[-1]+dy/2, len(y))
        pos = ax.pcolorfast(x, y, data2, **kwargs)

        # colorbar
        if colorbar:
            plt.colorbar(pos, aspect=50)#, extend='both')

        return pos

    def meshplot(self, ax=None, colorbar=False, **kwargs):
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
            if settings.ALWAYS_PLOT_NEW_WINDOW:
                figure()
                if  settings.FIGURE_POSITION is not None:
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

        # fix monotonic of labels x
        if check_monotonicity(self.x_centers) != 1:
            x, ordering = fix_monotinicity(self.x_centers, np.arange(len(self.x_centers)), mode='increasing')
            assert len(x)==len(self.x_centers), f'Cannot plot when Image.x have repeated elements.\nEither fix Image.x or set it to None.\nx: {self.x_centers}'
            ordering = [int(i) for i in ordering]
            data = copy.deepcopy(self.data)
            for i in ordering:
                if i != ordering[i]:
                    data[:, i] = self.data[:, ordering[i]]
        else:
            # check if x_edges is defined
            if self.x_edges is not None:
                x = self.x_edges
            else:
                x = self.x_centers
            data = self.data

        # fix monotonic of labels y
        if check_monotonicity(self.y_centers) != 1:
            y, ordering = fix_monotinicity(self.y_centers, np.arange(len(self.y_centers)), mode='increasing')
            assert len(y)==len(self.y_centers), f'Cannot plot when Image.y have repeated elements.\nEither fix Image.x or set it to None.\ny: {self.y_centers}'
            ordering = [int(i) for i in ordering]
            data2 = copy.deepcopy(data)
            for i in ordering:
                if i != ordering[i]:
                    data2[i, :] = data[ordering[i], :]
        else:
            # check if y_edges is defined
            if self.y_edges is not None:
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

    def plot(self, ax=None, colorbar=False, verbose=True, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.imshow()`_.

        Warning:
            Pixels are always square. For irregular pixel row/columns, see Image.meshplot()

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
            verbose (bool, optional): if True, a warning will show up if data
                has iregular pixel sizes. Default is true.
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
            if settings.ALWAYS_PLOT_NEW_WINDOW:
                figure()
                if  settings.FIGURE_POSITION is not None:
                    try:
                        set_window_position(settings.FIGURE_POSITION)
                    except:
                        pass

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
            data2fit = np.array([[i, j] for i, j in zip(x2fit, y2fit) if j > max(y2fit)*0.0001])  # clean zeros
            kwargs['vmax'] = data2fit[-1, 0]

        # # x and y
        # assert check_monotonicity(self.x) == 1, f'x axis (Image.x) must be increasingly monotonic.\nData cannot be plotted by Image.imshow().\nPlease, use Image.plot() or set Image.x = None\nx: {self.x}'
        # assert check_monotonicity(self.y) == 1, f'y axis (Image.y) must be increasingly monotonic.\nData cannot be plotted by Image.imshow().\nPlease, use Image.plot() or set Image.y = None\ny: {self.y}'
        # x = np.linspace(self.x[0], self.x[-1], len(self.x))
        # y = np.linspace(self.y[0], self.y[-1], len(self.y))
        # if 'extent' not in kwargs:
        #     kwargs['extent'] = [x[0], x[-1], y[0], y[-1]]

        # check irregular spacing
        if verbose:
            sx = Spectrum(x=self.x_centers, y=self.x_centers)
            sy = Spectrum(x=self.y_centers, y=self.y_centers)
            try:
                sx.check_step_x()
                sy.check_step_x()
            except ValueError:
                print('Data seems to have different irregular pixel size. Maybe plot it using Image.meshplot().\nTo turn off this warning set verbose to False.')

        # fix monotonic of labels x
        if check_monotonicity(self.x_centers) != 1:
            x, ordering = fix_monotinicity(self.x_centers, np.arange(len(self.x_centers)), mode='increasing')
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
            y, ordering = fix_monotinicity(self.y_centers, np.arange(len(self.y_centers)), mode='increasing')
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

        # plot
        pos = ax.imshow(data2, **kwargs)

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
            None
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
            return Spectrum(x=self.y_centers, y=np.sum(self._data, axis=0))
        elif axis == 1:
            return Spectrum(x=self.x_centers, y=np.sum(self._data, axis=1))

    def floor(self, x=0, y=0, n=30, nx=None, ny=None):
        """Set background intensity to zero.

        Args:
            x, y (int, optional): x and y position to average background intensity.
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

        self._data -= np.mean(self._data[int(x_start):int(x_stop), int(y_start):int(y_stop)])
        self._vmin = min([min(x) for x in self.data])
        self._vmax = max([max(x) for x in self.data])

    def calculate_shift(self, axis=0, limit_size=1000):
        """Calculate intensity misalignments via cross-correlation.

        Args:
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images.
                Default is 1000. Set to False to bypass this limit.

        Returns:
            None
        """
        axis = _axis_interpreter(axis)

        mode = 'cross-correlation'
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

        # calculate
        ss.calculate_shift(mode=mode, peak=peak, bkg_check=bkg_check)
        self._calculated_shift = ss.calculated_shift
        self._calculated_shift.x = centers

    def set_shift(self, value=None, p=None, f=None, axis=0):
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

        # OLD... shift is absolute!!!!!!!!! Works ==============================
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

        if self.reduced is not None:
            self.binning(nbins=self.nbins)

    def fix_curvature(self, deg=2, axis=0):
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
        self.reduced.calculate_shift(axis=axis)

        p, f = self.reduced.calculated_shift.polyfit(deg=deg)
        self._p = p
        self._f = f

        self.set_shifts(p=p, axis=axis)
        self._calculated_shift = self.reduced.calculated_shift


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

        calculated_shift (brixs.Spectrum): Calculated shifts.
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
        calculate_shift()
        set_shifts()
        fix_curvature()

    """

    _read_only = ['reduced', 'calculated_shift', 'p', 'f']

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
        self._calculated_shift = None
        self._p      = None
        self._f      = None
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
            self._data = np.array([event for event in np.array(value) if not event[2]<=0])
        elif value.shape[1] == 2:
            data = np.array(value)
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
        self._calculated_shift = None
        self._p      = None
        self._f      = None
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

    def save(self, filepath, only_data=False,  **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.

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
            kwargs['header'] += '==== PhotonEvents attributes ===='  + '\n'
            for n in dict:
                if n not in ['_data', '_x', '_y', '_I', '_reduced', '_calculated_shift', '_f', '_shifts']:
                    if isinstance(dict[n], Iterable):
                        kwargs['header'] += f'{n}: {list(dict[n])}'  + '\n'
                    else:
                        kwargs['header'] += f'{n}: {dict[n]}'  + '\n'
            np.savetxt(Path(filepath), self._data, **kwargs)

    def load(self, filepath, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first decompressed.

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
            if '==== PhotonEvents attributes ====' in line:
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
            if settings.ALWAYS_PLOT_NEW_WINDOW:
                figure()
                if  settings.FIGURE_POSITION is not None:
                    try:
                        set_window_position(settings.FIGURE_POSITION)
                    except:
                        pass

        # kwargs
        if 's' not in kwargs:
            kwargs['s'] = 0.1

        # plot
        ax.scatter(self.data[:, 0], self.data[:, 1], **kwargs)

        # set limits
        ax.xlim(0, self.shape[1])
        ax.ylim(0, self.shape[0])

        return ax

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
        self.reduced._x = moving_average(_x_edges, n=2)
        self.reduced._y = moving_average(_y_edges, n=2)
        # self._x_centers = moving_average(self.x_edges, n=2)
        # self._y_centers = moving_average(self.y_edges, n=2)
        self._nbins     = _nbins
        self._bins_size = _bins_size

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

    def calculate_shift(self, axis=0):
        """For now, the only mode tested is cc"""
        axis = _axis_interpreter(axis)
        mode = 'cc'
        assert self.reduced is not None, 'Image was not binned yet.\nPlease, use Image.binning()'

        self.reduced.calculate_shift(axis=axis, mode=mode)

        if mode in cc:
            if axis == 0:
                self.reduced.calculated_shift._y = self.reduced.calculated_shift.y*self.bins_size[0]
            elif axis == 1:
                self.reduced.calculated_shift._y = self.reduced.calculated_shift.y*self.bins_size[1]
        self._calculated_shift = self.reduced.calculated_shift

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
        """mode = 'cc'"""
        axis = _axis_interpreter(axis)

        assert self.reduced is not None, 'Image was not binned yet.\nPlease, use Image.binning()'

        self.reduced.floor()
        self.calculate_shift(axis=axis)

        p, f = self.calculated_shift.polyfit(deg=deg)
        self._p = p
        self._f = f

        # set shifts
        self.set_shifts(p=p)

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

    _read_only = ['step', 'monotonicity',  'R2', 'fit', 'residue', 'guess']
    _non_removable = []

    def __init__(self, *args, **kwargs):
        # argument parsing
        data, x, y, filepath = self._args_checker(args, kwargs)

        # sorting data
        if data is not None:
            self.data = data
            # _restart_attr is builtin inside data
        elif filepath is not None:
            self.load(filepath)
            # _restart_attr is builtin inside load
        elif y is not None:
            self._restart_attr()
            if x is None:
                self._x = np.arange(0, len(y))
                self._step = 1
                self._monotonicity = 'increasing'
            else:
                if len(x) == len(y):
                    self._x = np.array(x)
                    self._step = None
                    self._monotonicity = None
                else:
                    raise ValueError('x and y data are not the same length.')
            self._y = np.array(y)
            self._data = np.vstack((self.x, self.y)).transpose()
        else:
            self._data = None
            self._x = None
            self._y = None
            self._restart_attr()


        # default values
        # self.calib = settings.DEFAULT_SPECTRUM_CALIB

    def __len__(self):
        if self.x is None:
            return 0
        else:
            return len(self.x)

    def __add__(self, s):
        ss = Spectra([self, s])
        try:
            ss.check_same_x()
        except ValueError:
            raise ValueError('Cannot add spectra. x axis is different.')
        return Spectrum(x=self.x, y=self.y + s.y)

    def __sub__(self, s):
        ss = Spectra([self, s])
        try:
            ss.check_same_x()
        except ValueError:
            raise ValueError('Cannot subtract spectra. x axis is different.')
        return Spectrum(x=self.x, y=self.y - s.y)

    def __mul__(self, s):
        ss = Spectra([self, s])
        try:
            ss.check_same_x()
        except ValueError:
            raise ValueError('Cannot subtract spectra. x axis is different.')
        return Spectrum(x=self.x, y=self.y * s.y)

    def __div__(self, s):
        ss = Spectra([self, s])
        try:
            ss.check_same_x()
        except ValueError:
            raise ValueError('Cannot subtract spectra. x axis is different.')
        return Spectrum(x=self.x, y=self.y / s.y)

    def _args_checker(self, args, kwargs):
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
                self._y = np.array(value)
                self._data = np.vstack((self.x, self.y)).transpose()
                self._step = 1
                self._monotonicity = 'increasing'
                self._restart_attr()
            else:
                self._data = np.array(value)
                self._x = self.data[:, 0]
                self._y = self.data[:, 1]
                self._step = None
                self._monotonicity = None
                self._restart_attr()
        except IndexError:  # one column data (list)
            self._x = np.arange(0, len(value))
            self._y = np.array(value)
            self._data = np.vstack((self.x, self.y)).transpose()
            self._step = 1
            self._monotonicity = 'increasing'
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
            raise ValueError('length of x you are trying to set is not compatible with length of y.')
        else:
            self._x = np.array(value)
            self._data[:, 0] = np.array(value)
            self._step = None
            self._monotonicity = None
            self._restart_attr()
            if self.calib != 1:
                print('Calibration not 1. Applying calibration to new data.')
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
            raise ValueError('length of y you are trying to set is not compatible with length of x.')
        else:
            self._y = np.array(value)
            self._data[:, 1] = np.array(value)
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
        # raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shift(value, mode='roll').")
        self.set_shift(value, mode='roll')
    @shift_roll.deleter
    def shift_roll(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_interp(self):
        return self._shift_interp
    @shift_interp.setter
    def shift_interp(self, value):
        # raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shift(value, mode='interp').")
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
        # return copy.deepcopy(self._peaks)
        return self._peaks
    @peaks.setter
    def peaks(self, value):
        if isinstance(value, Peaks):
            self._peaks = copy.deepcopy(value)
        elif isinstance(value, dict):
            self._peaks = Peaks(value)
        else:
            raise ValueError(f'Peaks must be a dictionary.')
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


    def _restart_attr(self):
        """Set relevant attributes to initial value."""
        self._step = None
        self._monotonicity = None

        self._calib = 1
        self._shift = 0
        self._shift_roll = 0
        self._shift_interp = 0
        self._offset = 0
        self._factor = 1

        self._fit = None
        self._residue = None
        self._guess = None
        self._R2 = None

        self._peaks = Peaks()

    def _check_ranges(self, ranges):
        """check if ranges is the right format.

        If any item of ranges is None, this item is replaced by the min or max
            value of the data.
        """
        text = 'Ranges should be a pair (x_init1, x_final1) or a list of pairs like this: ((x_init1, x_final1), (x_init2, x_final2), ...)\nUse None to indicate the minimum or maximum x value of the data.'
        # check format
        if ranges is None:
            ranges = ((min(self.x), max(self.x)),)
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
                         ranges[i][0] = min(self.x)
                    if ranges[i][1] is None:
                         ranges[i][1] = max(self.x)
            else:
                raise ValueError(f'Ranges pair {r} is not a valid pair.\n'+text)
        return tuple(ranges)


    def save(self, filepath, only_data=False,  **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.

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
            kwargs['header'] += '==== Spectrum attributes ===='  + '\n'
            for n in dict:
                if n not in ['_x', '_y', '_data']:
                    if isinstance(dict[n], Iterable):
                        kwargs['header'] += f'{n}: {list(dict[n])}'  + '\n'
                    else:
                        kwargs['header'] += f'{n}: {dict[n]}'  + '\n'
            np.savetxt(Path(filepath), self._data, **kwargs)

    def load(self, filepath, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first decompressed.

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

        # find where attributes listing starts
        for i, line in enumerate(header):
            if '==== Spectrum attributes ====' in line:
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
                if name in ['_peaks']:
                    if value == []:
                        self._peaks = Peaks()
                    else:
                        print('warning: peaks cannot be loaded yet')
                        self._peaks = Peaks()
                ### DEALING WITH OTHER ATTRIBUTES ###
                elif name not in []:  # except these attrs
                    try:
                        setattr(self, name, value)
                    except Exception as e:
                        print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')


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
        if self.monotonicity is None:
            self.check_monotonicity()

        # check step uniformity
        d = np.diff(self.x)
        if abs((max(d) - min(d))*100/np.mean(np.diff(self.x))) > max_error:
            self._step = None
            raise ValueError(f"Step in the x-coordinate seems not to be uniform.")

        self._step = np.mean(d)

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
            raise ValueError('x array is not monotonic. Use Spectrum.fix_monotinicity()')

    def fix_monotinicity(self, mode='increasing'):
        """Rearrange x, y such as x array is monotonically increasing or decreasing.

        Args:
            mode (str, optional): increasing or decreasing.

        Returns:
            None
        """
        if mode != 'increasing' and mode != 'decreasing':
            raise ValueError('mode should be "decreasing" or "increasing".')
        if self.monotonicity is None:
            try:
                self.check_monotonicity()
            except ValueError:
                unqa, ID, counts = np.unique(self.x, return_inverse=True, return_counts=True)
                self.data = np.column_stack(( unqa , np.bincount(ID,self.y)/counts ))
                self.check_monotonicity()
        if self.monotonicity != mode:
            temp = sort(self.x, self.x, self.y)
            if mode == 'increasing':
                self._x = np.array(temp[0])
                self._y = np.array(temp[1])
            if mode == 'decreasing':
                self._x = np.fliplr(temp[0])
                self._y = np.fliplr(temp[1])
            self._calib = -self.calib
            if self.step is not None:
                self._step = -self.step

            self.check_monotonicity()


    def set_calib(self, value):
        """Set calibration value.

        Args:
            value (number): calibration value (x-coordinates will be multiplied
                by this value).

        Returns:
            None
        """
        if value == 0:
            raise ValueError('cannot set calib = 0.')
        if self.calib != value:
            if self.calib != 1:
                self._x = self.x*self.calib**-1
                self.data[:, 0] = self._x
                if self.step is not None:
                    self._step *= self.calib**-1
            if value != 1:
                self._x = self.x*value
                self.data[:, 0] = self._x
                if self.step is not None:
                    self._step *= value
            self._calib = value

        # reset monotonicity ======================================
        # step is changed above
        self._monotonicity = None

        # pass value to child data =================================
        if self.peaks is not None:
            self.peaks.calib = value
        if self.guess is not None:
            self.guess.calib   = value
        if self.fit is not None:
            self.fit.calib     = value
        if self.residue is not None:
            self.residue.calib = value

    def set_shift(self, value, mode):
        """Set shift value.

        Args:
            value (float or int): shift value (value will be added to x-coordinates).
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

        if mode in roll and self.step is None:
            try:
                self.check_step_x()
            except ValueError:
                raise ValueError(f'Cannot shift data using mode = {mode}, because x-coordinates are not uniform.')

        if mode in roll:
            if is_integer(value) == False:
                raise ValueError('shift must be an integer for mode = `roll`.')

            if self.shift_roll != value:
                if self.shift_roll != 0:  # undo last shift
                    self._x, self._y = shifted(self.x, self.y, value=-self.shift_roll, mode='roll')

                    if self.peaks is not None:
                        self.peaks.shift += self.step*-self.shift_roll
                if value != 0:
                    self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                    if self.peaks is not None:
                        self.peaks.shift += self.step*value
                self._shift_roll = value
                self._data[:, 0] = self._x
                self._data[:, 1] = self._y
        elif mode in hard:
            if self.shift != value:
                if self.shift != 0:
                    self._x, self._y = shifted(self.x, self.y, value=-self._shift, mode='x')
                    if self.peaks is not None:
                        self.peaks.shift += -self._shift
                if value != 0:
                    self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                    if self.peaks is not None:
                        self.peaks.shift += value
                self._shift = value
                self._data[:, 0] = self._x
                self._data[:, 1] = self._y
        elif mode in soft:
            if self.shift_interp != value:
                if self.monotonicity is None:
                    self.check_monotonicity()
                if self.monotonicity != 'increasing':
                    raise ValueError('x array must be monotonicaly increasing.\nTip: use Spectrum.fix_monotinicity()')

                if self.shift_interp != 0:
                    self._x, self._y = shifted(self.x, self.y, value=-self.shift_interp, mode='interp')
                    if self.peaks is not None:
                        self.peaks.shift += -self.shift_interp
                if value != 0:
                    self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                    if self.peaks is not None:
                        self.peaks.shift += value
                self._shift_interp = value
                self._data[:, 0] = self._x
                self._data[:, 1] = self._y
        else:
            raise ValueError(f'Invalid mode. Valid options are `roll`, `x`, `interp`.')

        # pass value to child data =================================
        # peaks is changed above
        if self.guess is not None:
            self.guess.set_shift(value, mode=mode)
        if self.fit is not None:
            self.fit.set_shift(value, mode=mode)
        if self.residue is not None:
            self.residue.set_shift(value, mode=mode)

    def set_offset(self, value):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).

        Returns:
            None
        """
        if self.offset != value:
            if self.offset != 0:
                self._y = self.y - self.offset
                self.data[:, 1] = self._y
            if value != 0:
                self._y = self.y + value
                self.data[:, 1] = self._y
            self._offset = value

        # pass value to child data =================================
        if self.peaks is not None:
            self.peaks.offset = value
        if self.guess is not None:
            self.guess.offset = value
        if self.fit is not None:
            self.fit.offset = value
        if self.residue is not None:
            self.residue.offset = value

    def set_factor(self, value):
        """Set y multiplicative factor.

        Args:
            value (number): multiplicative factor (y-coordinates will be
                multiplied by this value).

        Returns:
            None
        """
        if value == 0:
            raise ValueError('cannot set factor = 0.')
        if self.factor != value:
            if self.factor != 1:
                self._y = self.y*self.factor**-1
                self.data[:, 1] = self._y
                # self._multiplicative_y_fix(self.factor**-1)
            if value != 1:
                self._y = self.y*value
                self.data[:, 1] = self._y
            self._factor = value

        # pass value to child data =================================
        if self.peaks is not None:
            self.peaks.factor = value
        if self.guess is not None:
            self.guess.factor = value
        if self.fit is not None:
            self.fit.factor = value
        if self.residue is not None:
            self.residue.factor = value


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
            If none arguments is given, this function will adust x to have regular spacing.

        Returns:
            None

        Warning:
            Spectrum.interp() will change data for good. There is no way to
                recover the previous state of the array.
        """
        self.check_monotonicity()
        if self.monotonicity != 'increasing':
            raise ValueError('monotonicity of data must be strictly increasing. Use Spectrum.fix_monotinicity().')

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
        self._data = np.vstack((self._x, self._y)).transpose()

        # reset
        self._step = None
        self._monotonicity = None
        self._fit     = None
        self._residue = None
        self._guess   = None
        self._R2      = None

    def extract(self, ranges):
        """Return the extracted data range from full data.

        Args:
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            :py:attr:`Spectrum`
        """
        ranges = self._check_ranges(ranges)
        x, y  = extract(self.x, self.y, ranges)
        s = Spectrum(x=x, y=y)
        s._offset       = self.offset
        s._factor       = self.factor
        s._calib        = self.calib
        s._shift        = self.shift
        s._shift_roll   = self.shift_roll
        s._shift_interp = self.shift_interp
        s._peaks        = self.peaks

        for attr in self.__dict__:
            if attr.startswith('_') == False:
                s.__setattr__(attr, self.__dict__[attr])
        return s

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
        self._data = np.vstack((self._x, self._y)).transpose()

    def floor(self, value=None, n=20, ranges=None):
        """Sets zero value for y-coordinates (shifts data verticaly).

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

            if value is None:
                i = index(self.x, min(self.x))
            else:
                i = index(self.x, value)

            if self.monotonicity == 'increasing' or self.monotonicity == 'decreasing':
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
                self.offset = self.offset+value2
            else:
                value2 = -self.y[i]
                self.offset = self.offset+value2
        else:
            ranges = self._check_ranges(ranges)
            _, y= extract(self.x, self.y, ranges)
            value2 = -np.mean(y)
            self.offset = self.offset+value2

    def flip(self):
        """Flip x axis.

        Returns:
            None
        """
        self.calib = -self.calib
        # self.x = -self.x

    def calculate_area(self):
        """Returns the calculated area under the curve. Wrapper for `numpy.trapz()`_.

        Returns:
            Number

        .. _numpy.trapz(): https://numpy.org/doc/stable/reference/generated/numpy.trapz.html
        """
        return np.trapz(y=self.y, x=self.x)

    def normalize(self, value, ranges):
        """Set a factor such as the average y between ranges is equal to value.

        It uses Spectra.set_factor().

        Args:
            value (number): value.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_factor`
        """
        s = self.extract(ranges=ranges)
        self.set_factor(self.factor*value/np.mean(s.y))

    def zero(self, mode='max', peak=0):
        """Uses Spectrum.set_shift() to shift a specified position to zero.

        It shifts data using shift mode: 'hard'.

        Args:
            mode (string, optional): method for identifing the point which the x-
                coordinate is zero. Default is 'max'. Options are: 'max', 'min',
                'fitted peaks', 'peak'.
            peak (int, optional): if mode='peak' or mode='fitted peaks', this variable selects which peak
                to align data.

        Returns:
            None
        """
        shift_mode = 'x'
        if mode == 'max':
            self.set_shift(-self.x[np.argmax(self.y)], mode=shift_mode)
        if mode == 'min':
            self.set_shift(-self.x[np.argmin(self.y)], mode=shift_mode)
        elif mode == 'fitted peaks':
            assert self.fit != None, 'fit is not defined for this spectrum.'
            assert self.fit.peaks != None, 'peaks are not defined for the fitted curve.'
            self.set_shift(-self.fit.peaks[peak]['c'], mode=shift_mode)
        elif mode == 'peak':
            assert self.peaks != None, 'peaks are not defined for this spectrum.'
            self.set_shift(-self.fit.peaks[peak]['c'], mode=shift_mode)
        else:
            raise ValueError('mode not valid.\nValid modes: max, min, fitted peaks, peak.')


    def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, **kwargs):
        """Plot spectrum. Wrapper for `matplotlib.pyplot.plot()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number, optional): defines a vertical offset. Default is 0.
            shift (number, optional): horizontal shift value. Default is 0.
            factor (number, optional): multiplicative factor on the y axis.
                Default is 1.
            calib (number, optional): multiplicative factor on the x axis.
                Default is 1.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            `Line2D`_

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        if ax is None:
            ax = plt
            if settings.ALWAYS_PLOT_NEW_WINDOW:
                figure()
                if  settings.FIGURE_POSITION is not None:
                    try:
                        set_window_position(settings.FIGURE_POSITION)
                    except:
                        pass

        return ax.plot((self.x*calib) + shift, self.y*factor + offset, **kwargs)


    def find_peaks(self, prominence=5, width=4, moving_average_window=8):
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
        # check monotonicity
        if self.monotonicity is None:
            self.check_monotonicity()
        # if self.monotonicity != 'increasing':
        #     raise ValueError('x axis must be monotonicaly increasing.\nUse Spectrum.fix_monotinicity().')


        # check points and moving_average_window
        if width < 1:
            raise ValueError('width must be 1 or higher.')
        if isinstance(width, int) == False:
            if width.is_integer() == False:
                raise ValueError('width must be an integer.')
        if moving_average_window < 1:
            raise ValueError('moving_average_window must be 1 or higher.')
        if isinstance(moving_average_window, int) == False:
            if moving_average_window.is_integer() == False:
                raise ValueError('moving_average_window must be an integer.')

        # step
        if self.step is None:
            self.check_step_x()

        # data smoothing
        if moving_average_window > 1:
            y2 = moving_average(self.y, moving_average_window)
            x2 = moving_average(self.x, moving_average_window)
        else:
            x2 = x[:]
            y2 = y[:]

        # parameters
        if prominence is None:
            prominence = (max(y2)-min(y2))*0.05
        else:
            prominence = (max(y2)-min(y2))*prominence/100

        # find peaks
        try:
            peaks, d = find_peaks(y2, prominence=prominence, width=width)
            # print(d)
            # self._peaks = {i: {'amp': d['width_heights'][i], 'c': x2[peaks[i]],  'fwhm': abs(d['widths'][i]*self.step)} for i in range(len(peaks))}
            # self._peaks = {i: {'amp': d['prominences'][i]+(y2[d['right_bases'][i]]+y2[d['left_bases'][i]])/2, 'c': x2[peaks[i]],  'fwhm': abs(d['widths'][i]*self.step)} for i in range(len(peaks))}
            shift = self.shift_roll*self.step + self.shift_interp + self._shift
            self._peaks = Peaks([{'amp': d['prominences'][i]+max([y2[d['right_bases'][i]], y2[d['left_bases'][i]]]), 'c': x2[peaks[i]],  'fwhm': abs(d['widths'][i]*self.step)} for i in range(len(peaks))],
                                      shift=shift, offset=self.offset, calib=self.calib, factor=self.factor)
        except IndexError:
            self._peaks = Peaks({}, shift=self.shift, offset=self.offset, calib=self.calib, factor=self.factor)

    def fit_peaks(self, asymmetry=False, fixed_m=0, offset=True, amp_bounds=(0, 3), c_bounds=(-2, 2), fwhm_bounds=(0, 2), verbose=False):
        """Fit peaks. Wrapper for `scipy.optimize.curve_fit()`_.

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
        self.check_monotonicity()

        # check if peaks is defined ============================================
        if len(self.peaks) == 0:
            raise ValueError('No peaks to fit. Run Spectrum.find_peaks().')

        # check asymmetry ======================================================
        if isinstance(asymmetry, dict):
            for i in range(len(self.peaks)):
                if i not in asymmetry:
                    asymmetry[i] = False
                else:
                    if type(asymmetry[i]) != bool:
                        raise ValueError(f'the asymmetry must be bool (True or False) for all peaks.\nasymmetry = {asymmetry}')
            for i in asymmetry:
                if i not in range(len(self.peaks)):
                    raise ValueError(f'Peak {i} not in peaks.\nDefined peaks: {self.peak}\nasymmetry = {asymmetry}).')
        else:
            if isinstance(asymmetry, Iterable):
                raise ValueError('asymmetry must be a number or dictionary.')
            asymmetry = {i: asymmetry for i in range(len(self.peaks))}

        # check fixed_m ========================================================
        if isinstance(fixed_m, dict):
            for i in range(len(self.peaks)):
                if i not in fixed_m:
                    fixed_m[i] = 0
                else:
                    if (isinstance(fixed_m[i], bool) and fixed_m[i] == True) or (fixed_m[i] < 0 or fixed_m[i] > 1):
                        raise ValueError(f'fixed_m must be False or a number from 0 to 1 for all peaks.\nfixed_m = {fixed_m}')
            for i in fixed_m:
                if i not in range(len(self.peaks)):
                    raise ValueError(f'Peak {i} not in peaks.\nDefined peaks: {self.peak}\nfixed_m = {fixed_m}).')
        else:
            if isinstance(fixed_m, Iterable):
                raise ValueError('fixed_m must be a number from 0 to 1, False, or dictionary.')
            fixed_m = {i: fixed_m for i in range(len(self.peaks))}

        # build model ==========================================================
        model_str = ''
        args_str = ''
        bounds_min = []
        bounds_max = []
        p0 = []
        peaks_reconstruction = {}

        if verbose: print('Building model...')
        for i in range(len(self.peaks)):
            if verbose: print(f'Peak {i}')
            f, f_str, a_str = _peak_function_creator(asymmetry=asymmetry[i], fixed_m=fixed_m[i], idx=i)
            p, bmin, bmax = self.peaks[i].build_guess(asymmetry=asymmetry[i], fixed_m=fixed_m[i], amp_bounds=amp_bounds, c_bounds=c_bounds, fwhm_bounds=fwhm_bounds)
            p0.extend(p)
            bounds_min.extend(bmin)
            bounds_max.extend(bmax)
            model_str += f_str + ' + '
            args_str += a_str + ', '
            peaks_reconstruction[i] = {'asymmetry':asymmetry[i], 'fixed_m':fixed_m[i], '_func':f, 'n_args':len(p)}

        if offset == True and type(offset) == bool:
            if verbose: print(f'Offfset')
            model_str = f'lambda x, {args_str}offset: {model_str}offset'
            p0 = np.append(p0, min(self.y))
            bounds_min = np.append(bounds_min, -np.inf)
            bounds_max = np.append(bounds_max, np.inf)
        elif offset == False and type(offset) == bool:
            model_str = f'lambda x, {args_str[:-2]}: {model_str[:-3]}'
        else:
            if verbose: print(f'Fixed offfset')
            model_str = f'lambda x, {args_str}: {model_str}{offset}'
        model = eval(model_str)

        # guess ================================================================
        self._guess = Spectrum(x=self.x, y=model(self.x, *p0))
        self.guess._offset       = self.offset
        self.guess._factor       = self.factor
        self.guess._calib        = self.calib
        self.guess._shift        = self.shift
        self.guess._shift_roll   = self.shift_roll
        self.guess._shift_interp = self.shift_interp
        self.guess.peaks         = copy.deepcopy(self.peaks)

        # curve_fit ============================================================
        if verbose: print(f'Fitting data...')
        popt, pcov = curve_fit(model, self.x, self.y, p0=p0, bounds=(bounds_min, bounds_max))
        if verbose: print(f'Done!')

        # fit ==================================================================
        x_temp = np.linspace(self.x[0], self.x[-1], len(self.x)*2)
        self._fit = Spectrum(x=x_temp, y=model(x_temp, *popt))
        self.fit._offset       = self.offset
        self.fit._factor       = self.factor
        self.fit._calib        = self.calib
        self.fit._shift        = self.shift
        self.fit._shift_roll   = self.shift_roll
        self.fit._shift_interp = self.shift_interp
        self.fit.peaks = Peaks()

        # residue ==============================================================
        if offset:
            self._residue = Spectrum(x=self.x, y=self.y-model(self.x, *popt)+popt[-1])
        else:
            self._residue = Spectrum(x=self.x, y=self.y-model(self.x, *popt))
        self.residue._offset       = self.offset
        self.residue._factor       = self.factor
        self.residue._calib        = self.calib
        self.residue._shift        = self.shift
        self.residue._shift_roll   = self.shift_roll
        self.residue._shift_interp = self.shift_interp

        # R2 ===================================================================
        self._R2 =  1- (sum((self.y-model(self.x, *popt))**2)/sum((self.y-np.mean(self.y))**2))

        # save fitted peaks parameters =========================================
        stop = 0
        for i in range(len(self.peaks)):
            peak = {}
            start = stop
            stop  = start + peaks_reconstruction[i]['n_args']
            peak['amp'] = popt[start]
            peak['c'] = popt[start+1]
            if peaks_reconstruction[i]['fixed_m'] == False and type(peaks_reconstruction[i]['fixed_m']) == bool:  # variable m
                if peaks_reconstruction[i]['asymmetry']:
                    peak['fwhm1'] = popt[start+2]
                    peak['fwhm2'] = popt[start+4]
                    peak['m1'] = popt[start+3]
                    peak['m2'] = popt[start+5]
                else:
                    peak['fwhm'] = popt[start+2]
                    peak['m'] = popt[start+3]
            else:
                if peaks_reconstruction[i]['asymmetry']:
                    peak['fwhm1'] = popt[start+2]
                    peak['fwhm2'] = popt[start+3]
                else:
                    peak['fwhm'] = popt[start+2]
                peak['m'] = peaks_reconstruction[i]['fixed_m']
            self.fit.peaks.append(peak)
            if offset:
                self.fit.peaks.offset = popt[-1]

    def polyfit(self, deg=2):
        """Fit data with a polynomial. Wrapper for `numpy.polyfit()`_.

        Args:
            deg (int, optional): degree of the fitting polynomial. Default is 2.

        Returns:
            Polynomial coefficients, highest power first.
            Model function f(x).

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        popt = np.polyfit(self.x, self.y, deg=deg)
        model = lambda x: np.polyval(popt, x)

        # fit ==================================================================
        x_temp = np.linspace(self.x[0], self.x[-1], len(self.x)*2)
        self._fit = Spectrum(x=x_temp, y=model(x_temp))
        self.fit._offset       = self.offset
        self.fit._factor       = self.factor
        self.fit._calib        = self.calib
        self.fit._shift        = self.shift
        self.fit._shift_roll   = self.shift_roll
        self.fit._shift_interp = self.shift_interp
        self.fit.peaks = Peaks()

        # guess ================================================================
        self._guess = None

        # residue ==============================================================
        self._residue = Spectrum(x=self.x, y=self.y-model(self.x))
        self.residue._offset       = self.offset
        self.residue._factor       = self.factor
        self.residue._calib        = self.calib
        self.residue._shift        = self.shift
        self.residue._shift_roll   = self.shift_roll
        self.residue._shift_interp = self.shift_interp

        return popt, model


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

        self._restart_attr()


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
                      'calculated_calib', 'calculated_factor',
                      'calculated_offset', 'calculated_shift']
    _non_removable = []


    def __init__(self, *args, **kwargs):
        # argument parsing
        data, dirpath, n = self._args_checker(args, kwargs)

        # sorting data
        if data is not None:
            self.data = data
            # _restart_attr is builtin inside data
        elif dirpath is not None:
            self.load(dirpath)
            # _restart_attr is builtin inside load
        elif n is not None:
            self._data     = [-1]*n
            self._restart_check_attr()
            self._calculated_shift  = None
            self._calculated_calib  = None
            self._calculated_factor = None
            self._calculated_offset = None
        else:
            self._data = []
            self._restart_check_attr()
            self._calculated_shift  = None
            self._calculated_calib  = None
            self._calculated_factor = None
            self._calculated_offset = None

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
         return self.data[item]

    def __setitem__(self, item, value):
        if isinstance(value, Spectrum) == False:
            raise ValueError(f'value must be of type brixs.spectrum, not {type(value)}')
        self.data[item] = value
        self._restart_check_attr()
        self._calculated_shift  = None
        self._calculated_calib  = None
        self._calculated_factor = None
        self._calculated_offset = None

    def __len__(self):
        return len(self.data)

    def __delitem__(self, item):
        self.remove(item)

    def _args_checker(self, args, kwargs):
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
            raise AttributeError(f'invalid attributes.\nValid atributes are `data`, `n`, and `filepaths`\nInput attributes: {kwargs.keys()}')

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
                self._data = copy.deepcopy(value)
            else:
                raise ValueError('data must be a list.')
        self._restart_check_attr()
        self._calculated_shift  = None
        self._calculated_calib  = None
        self._calculated_factor = None
        self._calculated_offset = None
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].shift
        return temp
    @shift.setter
    def shift(self, value):
        self.set_shift(value=value, mode='x')
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
        self.set_shift(value=value, mode='roll')
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
        self.set_shift(value=value, mode='interp')
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
        self.set_factor(value)
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
        self.set_offset(value)
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
    def R2(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].R2
        return temp
    @R2.setter
    def R2(self, value):
        raise AttributeError('Cannot be set.')
    @R2.deleter
    def R2(self):
        raise AttributeError('Cannot delete object.')

    @property
    def fit(self):
        ss = Spectra(n=len(self))
        for i in range(len(self)):
            ss[i] = self[i].fit
        return ss
    @fit.setter
    def fit(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @fit.deleter
    def fit(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def guess(self):
        ss = Spectra(n=len(self))
        for i in range(len(self)):
            ss[i] = self[i].guess
        return ss
    @guess.setter
    def guess(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @guess.deleter
    def guess(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def residue(self):
        ss = Spectra(n=len(self))
        for i in range(len(self)):
            ss[i] = self[i].residue
        # if None in temp:
        #     print(f'Some spectra do not have fit defined.\nfit={temp}')
        return ss
    @residue.setter
    def residue(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @residue.deleter
    def residue(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def peaks(self):
        return self.get_peaks(key='peaks')
    @peaks.setter
    def peaks(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @peaks.deleter
    def peaks(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def map(self):
        return self.calculate_map()
    @map.setter
    def map(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @map.deleter
    def map(self):
        raise AttributeError('Attribute cannot be deleted.')


    def append(self, s=None, x=None, y=None, data=None):
        """Append spectrum to the spectrum list.

        Args:
            s (Spectrum obj or list): Spectrum object to be added or
                list of Spectrum.
            data (list or array, optional): two column list (or array).
            x (list or array, optional): x values (1D list/array). Overwrites `data`.
            y (list or array, optional): y values (1D list/array). Overwrites `data`.

        Returns:
            None

        See Also:
            :py:func:`Spectra.remove`.
        """
        if s is None:
            s = Spectrum(x=x, y=y, data=data)
            self.append(s=s)
        elif s is not None:
            if isinstance(s, Iterable):
                for i, temp in enumerate(s):
                    if isinstance(temp, Spectrum) == False:
                        raise ValueError(f'All entries must be of type brixs.spectrum.\nEntry {i} is of type {type(temp)}')
                self._data += copy.deepcopy(s)
            else:
                if isinstance(s, Spectrum) == False:
                    raise ValueError('Spectrum must be of type brixs.Spectrum.')
                self._data += [copy.deepcopy(s)]
        else:
            raise ValueError('No data to append.')
        self._restart_check_attr()
        self._calculated_shift  = None
        self._calculated_calib  = None
        self._calculated_factor = None
        self._calculated_offset = None

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
        self._restart_check_attr()
        self._calculated_shift  = None
        self._calculated_calib  = None
        self._calculated_factor = None
        self._calculated_offset = None


    def save(self, dirpath='./', prefix='spectrum_', suffix='.dat', zfill=None, only_data=False, verbose=False, **kwargs):
        r"""Save Spectra. Wrapper for `numpy.savetxt()`_.

        Args:
            dirpath (string or pathlib.Path, optional): folderpath or folder handle.
                Default is the current directory.
            prefix (string, optional): prefix used for naming the files.
            suffix (string, optional): suffix used for naming the files. If the
                filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            zfill (int, optional): number of digits for file numbering. If `None`,
                zfill will be determined.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            verbose (bool, optional): turn verbose on and off. Default is `False`.

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
        dirpath = Path(dirpath)

        # check if dirpath is a directory
        assert dirpath.exists(), f'dirpath does not exists.\ndirpath: {dirpath}'
        assert dirpath.is_dir(), f'dirpath is not a directory.\ndirpath: {dirpath}'

        # set filenames
        if zfill is None:
            zfill = n_digits(len(self)-1)[0]

        # saving
        if verbose: print('saving files...')
        for i, s in enumerate(self.data):
            filename = f'{prefix}' + f'{i}'.zfill(zfill) + f'{suffix}'
            if verbose:  print(f':{i}/{len(self)-1}: {filename}')
            s.save(filepath=dirpath/filename, only_data=only_data, **kwargs)
        if verbose: print('Done!')

    def load(self, dirpath, string='*', verbose=False, **kwargs):
        """Load data from text files. Wrapper for `numpy.genfromtxt()`_.

        Args:
            dirpath (string, path object, or list): folderpath, folder handle,
                or a list of filepaths. The list can be comprised of file and
                folderpaths. All files inside a folder where string can be found
                in the filename are imported. If the filename extension is .gz or .bz2,
                the file is first decompressed.
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
        self._calculated_shift  = None
        self._calculated_calib  = None
        self._calculated_factor = None
        self._calculated_offset = None

        # if dirpath is str
        if type(dirpath) == str:
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

                if Path(filepath).is_dir():
                    fl = filelist(dirpath=filepath, string=string)
                    for i, f in enumerate(fl):
                        if verbose: print(f'        {j+1}/{len(fl)}: {f}')
                        self.append(Spectrum(filepath=f))

                elif Path(filepath).is_file():
                    self.append(Spectrum(filepath=filepath))
                else:
                    raise ValueError(f'cannot read filepath.\nfilepath: {dirpath}')
        if verbose: print('Done!')


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

    def _restart_check_attr(self):
        """restart check attributes."""
        # check x attrs
        self._length       = None
        self._step         = None
        self._x            = None
        self._monotonicity = None

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

    def fix_monotinicity(self, mode='increasing'):
        """Rearrange x, y such as x array is monotonically increasing or decreasing.

        Args:
            mode (str, optional): increasing or decreasing.

        Returns:
            None
        """
        for s in self.data:
            s.fix_monotinicity(mode=mode)

    def check_length(self):
        """Checks if all spectra has the same length.

        If all spectra have the same length, it sets :py:attr:`Spectra.length` = length.
        Otherwise, it raises an error.

        Returns:
            None

        Raises:
            ValueError: spectra does not have the same length.

        See Also:
            :py:func:`Spectra.check_step_x`, :py:func:`Spectra.check_same_x`.
        """
        for i, s in enumerate(self.data):
            try:
                if len(s.data) != len(self.data[i+1].data):
                    self._length = None
                    raise ValueError(f"Spectrum {i} and {i+1} have the different length.\nSpectrum {i}: {len(s)}\nSpectrum {i+1}: {len(self[i+1])}")
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

            Two checks are performed:

            1. Checks if the step between two data points is the same
            through out the x vector for each spectrum (vector uniformity).

            2. Checks if the step is the same between different spectra.

            Returns:
                None

            Raises:
                ValueError: If condition 1, or 2 are not satisfied.

            See Also:
                :py:func:`Spectra.check_length`, :py:func:`Spectra.check_same_x`.
        """
        if self.x is None:
            try:
                self.check_same_x(max_error=max_error)
            except ValueError:
                # 1) check step uniformity
                steps = np.zeros(len(self))
                for idx, s in enumerate(self.data):
                    if s.step is None:
                        try:
                            s.check_step_x()
                        except ValueError:
                            raise ValueError(f"Step in the x-coordinate of spectrum {idx} seems not to be uniform.")
                    steps[idx] = s.step

                # 2) check step between spectra
                avg_step = np.mean(steps)
                if sum([abs(steps[i]-steps[i+1]) > abs(avg_step*max_error/100) for i in range(len(self)-1)]) > 0:
                    self._step = None
                    raise ValueError(f"Spectra seems to have different step size. Calculated step sizes = {steps}")
                self._step = avg_step
                return
        else:
            # check step uniformity
            temp = Spectrum(x=self.x, y=self.x)
            try:
                temp.check_step_x(max_error=max_error)
            except ValueError:
                raise ValueError(f"Step in the x-coordinate of spectrum {idx} seems not to be uniform.")
            self._step = temp.step
            return

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
        # if self.step is None:  # calculate step if it hasn't yet
        #     self.check_step_x(max_error=max_error)

        # average step
        if self.step is None:
            step = 0
            for s in self:
                step += np.mean(np.diff(s.x))
            step = step/len(self)
        else:
            step = self.step

        # check x between spectra
        for idx, s in enumerate(self):
            try:
                if max(abs(s.x - self[idx+1].x))*100/abs(step) > max_error:
                    self._x = None
                    raise ValueError(f"x axis of spectrum {idx} and {idx+1} seem to be different.\nUse brixs.Spectra.interp() to interpolate the data and make the x axis for different spectra match.")
            except IndexError:
                pass
        self._x = self[0].x

        # update length
        self._length = len(self.x)


    def _check_fit(self):
        """Check if all data has fit object. Raises an error."""
        temp = {i: None for i in range(len(self))}
        for i in range(len(self)):
                temp[i] = self[i].fit
        if None in temp:
            raise RuntimeError(f'Fit is not defined for some spectra.\nList of fits: {temp}')

    def _check_number_of_peaks_per_spectrum(self, fitted_peaks=False):
        """Check if number of peaks for each spectrum is the same."""
        if fitted_peaks:
            peaks = self.fit.get_peaks('s')
        else:
            peaks = self.get_peaks('s')

        temp = [len(peaks[i]) for i in range(len(peaks))]
        if np.all(np.diff(temp) != 0) == True:
            temp = str({i: temp[i] for i in range(len(temp))}).replace(', ', '\n').replace('{', '').replace('}', '')
            print(f'number of peaks is different for different spectrum.\nCalculated values based on peaks might be wrong.\nNumber of peaks per spectrum:\n{temp}')


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

        None:
            If none arguments is given, this function will adust all spectra to
            have the same x.

        Returns:
            None

        Warning:
            Spectra.interp() will change data for good. There is no way to
                recover the previous state of the array.
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
        self._restart_check_attr()

    def extract(self, ranges):
        """Return the extracted data range from full data.

        Args:
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            :py:attr:`Spectra`
        """
        if isinstance(ranges, Iterable):
            temp = Spectra(n=len(self))
            for i, s in enumerate(self.data):
                temp[i] = s.extract(ranges=ranges)
        else:
            raise ValueError(f'Ranges is not a valid value.\nRanges should be a pair (x_init1, x_final1) or a list of pairs like this: ((x_init1, x_final1), (x_init2, x_final2), ...)\nUse None to indicate the minimum or maximum x value of the data.')


        for attr in self.__dict__:
            if attr.startswith('_') == False:
                temp.__setattr__(attr, self.__dict__[attr])
        return temp

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
        # reset attributes =====================================================
        # self._restart_check_attr()
        # self._calculated_shift  = None
        # self._calculated_calib  = None
        self._calculated_factor = None
        self._calculated_offset = None

    def flip(self):
        """Flip x axis.

        Returns:
            None
        """
        for s in self:
            s.flip()
        self._restart_check_attr()
        self._calculated_shift  = None
        self._calculated_calib  = None

    def concatenate(self):
        """Return spectrum of concatenate spectra."""
        x = np.concatenate([s.x for s in self.data])
        y = np.concatenate([s.y for s in self.data])
        return Spectrum(x=x, y=y)

    def align(self, mode=None, peak=0, bkg_check=True):
        """Uses Spectra.calculate_shift() and Spectra.set_shift() to align spectra.

        Args:
            mode (string, optional): method used to calculate the shifts. If
                None, it will try to use 'cross-correlation'. If it fails, it
                will use 'max'. See :py:func:`Spectra.calculate_shift()` for more info.
            bkg_check (bool, optional): if mode is 'cc', bkg_check=true prevents
                from cross correlating spectra with significant background (10 % of
                the maximum y) which can lead to wrong results.
                See :py:func:`Spectra.calculate_shift()` for more info.
            peak (int, optional): if mode='peak' or mode='fitted peaks', this variable selects which peak
                to align data.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_shift`
        """
        if mode is None:
            if self.x is None:
                try:
                    self.check_same_x()
                    mode = 'cc'
                except ValueError:
                    mode = 'max'
            else:
                mode = 'cc'
        self.calculate_shift(mode=mode, peak=peak, bkg_check=bkg_check)
        self.set_shift()

    def normalize(self, mode=None, peak=0, bkg_check=True):
        """Uses Spectra.calculate_factor() and Spectra.set_factor() to normalize spectra.

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
            :py:func:`Spectra.calculate_factor`
        """
        if mode is None:
            mode = 'delta'
        self.calculate_factor(mode=mode, peak=peak)
        self.set_factor()


    def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, **kwargs):
        """Plot spectra. Wrapper for `matplotlib.pyplot.plot()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number or list, optional): defines a vertical offset. Default is 0.
            shift (number or list, optional): horizontal shift value. Default is 0.
            factor (number or list, optional): multiplicative factor on the y axis.
                Default is 1.
            calib (number or list, optional): multiplicative factor on the x axis.
                Default is 1.
            vertical_increment or vi (number): vertical increment between curves
                in terms of percentage of the maximum delta y.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.


        Returns:
            `Line2D`_ list, offsets list

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        if ax is None:
            ax = plt
            if settings.ALWAYS_PLOT_NEW_WINDOW:
                figure()
                if  settings.FIGURE_POSITION is not None:
                    try:
                        set_window_position(settings.FIGURE_POSITION)
                    except:
                        pass

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

        # plot
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self.data[i].plot(ax=ax, offset=offset[i], shift=shift[i], factor=factor[i], calib=calib[i], **kwargs)

        return temp, offset


    def set_shift(self, value=None, mode=None, type='relative'):
        """Shift data recursively.

        Args:
            value (number or list, optional): value will be added to x-coordinates.
                If None, it will look for calculated values from Spectra.calculate_shift().
            mode (string, optional): If ``mode='x'`` or ``mode='hard'``, y is fully preserved
                while x is shifted. If ``mode='y'``, ``'interp'``, or ``'soft'``, x is preserved
                while y is interpolated with a shift. If ``mode='roll'`` (or rotate or r), x is also preserved
                and y elements are rolled along the array (``shift`` value must be an integer).
                The form of y-coordinates is fully preserved, but the edge of y-coordinates are lost.
            type (str, optional): set values 'relative' or in 'absolute' units.
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
            :py:func:`Spectra.calculate_shift()`
        """
        # if verbose: print(f'Setting shifts...')

        # check value and mode
        if value is None:
            if self.calculated_shift is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_shift().')

            value = self.calculated_shift.y

            # if mode is not defined, auto pick the right one
            if mode is None:
                if self.calculated_shift.mode in cc:
                    mode = 'roll'
                elif self.calculated_shift.mode == 'max':
                    mode = 'x'
                elif self.calculated_shift.mode == 'peak':
                    mode = 'x'
                elif self.calculated_shift.mode == 'fitted peaks':
                    mode = 'x'
                elif self.calculated_shift.mode == 'user-defined':
                    text = 'shift mode not define (please, select "roll", "soft", "hard").' +\
                            'Note that shifts were calculated elsewhere by the user, ' +\
                            'the script can has no way of knowing the most appropriate shift mode.'
                    raise ValueError(text)
                elif self.calculated_shift.mode in None:
                    raise ValueError('calculated shifts not defined. Use Spectra.calculate_shift().')

            # if shift were calculated by cc, one should use roll
            if self.calculated_shift.mode in cc and mode not in roll:
                # one could in principle do a roll even though the calculation mode
                # is cc by multiplying the shift value by the step, however,
                # the whole point of using cc is to be able to do a roll shift.
                # Therefore, when cc is used, roll shift is enforced.
                raise ValueError('shifts were calculated using cross-correlation. Only shift mode possible is "roll".')

            # if shift were calculated by peak or max, one can opt to shift via roll
            # as long as self.step is defined
            if self.calculated_shift.mode in ['peak', 'max'] and mode in roll:
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
                raise ValueError('Arrays must be monotonicaly increasing.\nUse Spectra.fix_monotinicity()')

        # if mode is roll, steps must be defined for all spectra
        if mode in roll:
            if self.step is None:
                try:
                    self.check_step_x()
                except ValueError:
                    temp = [self[i].step for i in  range(len(self))]
                    raise ValueError(f'Cannot apply roll because some spectra have no step defined.\nData must be homogeneous.\nSteps: {temp}\nUse Spectra.check_step_x()')

        # set values ===========================================================
        for i in range(len(self)):
            if type == 'relative':
                # check for previous shifts
                if mode in roll:
                    v = value[i] + self[i].shift_roll
                elif mode in hard:
                    v = value[i] + self[i].shift
                elif mode in soft:
                    v = value[i] + self[i].shift_interp
            elif type == 'absolute':
                v = value[i]
            else:
                raise ValueError("type has to be either 'relative' or 'absolute'.")
            self.data[i].set_shift(value=v, mode=mode)

        # if mode not interp and shifts are different, reset x
        if mode not in soft:
            if all_equal(value) == False:
                self._x = None
        # reset attributes =====================================================
        # self._restart_check_attr()
        self._calculated_shift  = None
        self._calculated_calib  = None
        # self._calculated_factor = None
        # self._calculated_offset = None

    def set_factor(self, value=None, type='relative'):
        """Apply multiplicative y factor recursively.

        Args:
            value (number or list, optional): value will be multiplied to y-coordinates.
                If None, it will look for calculated values from Spectra.calculate_factor().
            type (str, optional): set values 'relative' or in 'absolute' units.
                If 'relative', modifications will be done on top of previous
                modifications. Default is 'relative'.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_factor()`
        """
        if value is None:
            if self.calculated_factor is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_factor().')

            value = self.calculated_factor.y
        else:
            # check if value is a number
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'


        # set values ===========================================================
        for i in range(len(self)):
            if type == 'relative':
                v = value[i]*self[i].factor
            elif type == 'absolute':
                v = value[i]
            else:
                raise ValueError("type has to be either 'relative' or 'absolute'.")
            self[i].set_factor(value=v)

        # reset attributes =====================================================
        # self._restart_check_attr()
        # self._calculated_shift  = None
        # self._calculated_calib  = None
        self._calculated_factor = None
        # self._calculated_offset = None

    def set_calib(self, value, type='relative'):
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
            if type == 'relative':
                v = value[i] * self[i].calib
            elif type == 'absolute':
                v = value[i]
            else:
                raise ValueError("type has to be either 'relative' or 'absolute'.")
            self[i].set_calib(value=v)

        # reset attributes =====================================================
        self._restart_check_attr()
        self._calculated_shift  = None
        self._calculated_calib  = None
        # self._calculated_factor = None
        # self._calculated_offset = None

    def set_offset(self, value=None, type='relative'):
        """Apply additive y factor recursively.

        Args:
            value (number or list, optional): value will be added to y-coordinates.
                If None, it will look for calculated values from Spectra.calculate_offset().
            type (str, optional): set values 'relative' or in 'absolute' units.
                If 'relative', modifications will be done on top of previous
                modifications. Default is 'relative'.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_offset()`, :py:func:`Spectra.floor()`
        """
        if value is None:
            if self.calculated_offset is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_offset().')

            value = self.calculated_offset.y
        else:
            # check if value is a number
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'


        # set values ===========================================================
        for i in range(len(self)):
            if type == 'relative':
                v = value[i] + self[i].offset
            elif type == 'absolute':
                v = value[i]
            else:
                raise ValueError("type has to be either 'relative' or 'absolute'.")
            self[i].set_offset(value=v)

        # reset attributes =====================================================
        # self._restart_check_attr()
        # self._calculated_shift  = None
        # self._calculated_calib  = None
        # self._calculated_factor = None
        self._calculated_offset = None


    def calculate_shift(self, mode='cross-correlation', peak=0, bkg_check=True):
        """Calculate how much spectra must be shifted to align with the first spectrum.

        Args:
            mode (string, optional): method used to calculate the shifts.
                The current options are: 'cross-correlation' ('cc'), 'max',
                'fitted peaks', or 'peak'.
                For mode='peak' or mode='fitted peaks', data must be fitted via
                Spectrum.fit_peak(), i. e.,
                spectra must have a Spectrum.fit.peaks object defined.
            bkg_check (bool, optional): if mode is 'cc', bkg_check=true prevents
                from cross correlating spectra with significant background (10 % of
                the maximum y) which
                can lead to wrong results. The algorithm uses find_peaks to find
                approximately whats is data and what is background, Then the
                background is said to be big if average(bkg) > max(y)*0.1. If
                background region cannot be found clearly, the criteria
                becames average(y) > max(y)*0.1.
            peak (int, optional): if mode='peak' or mode='fitted peaks', this variable selects which peak
                to align data.

        Returns:
            None

        Note:
            For ``mode = 'cc'``, spectra must have the same x-coordinates (this
            is checked before execution).

        See Also:
            :py:func:`Spectra.set_shifts`
        """
        # spectra will be alined to the first spectrum
        ref = 0

        # try checking if x is the same (execution is much faster if it is) ====
        if self.x is None:
            try:
                self.check_same_x()
            except ValueError:
                pass

        # common variables =====================================================
        values = np.array([0.0]*len(self))

        # CALCULATION ==========================================================
        if mode in cc:
            if bkg_check:  # check if background is zero
                temp = [False]*len(self)
                for i in range(len(self)):
                    peaks, d = find_peaks(self[i].y, prominence=(max(self[i].y)-np.mean(self[i].y))*0.05, width=1)
                    bkg = list(self[i].y)
                    for j in range(len(peaks)):
                        del bkg[peaks[j]-int(d['widths'][j]*10):peaks[j]+int(d['widths'][j])*10]
                    if len(bkg) > 10:
                        temp[i] = np.mean(bkg) > max(self[i].y)*0.1
                    else:
                        temp[i] = np.mean(self[i].y)>max(self[i].y)*0.1
                if True in temp:
                    raise ValueError(f'some spectra have a background that seems to be too big for a reliable cross-corelation.\nSpectra with big bkg: {temp}\nYou can try using Spectra.floor()')

            # x must be the same for cc
            x, ys = self._gather_ys()

            # x must be uniform (same step between each data point)
            self.check_step_x()

            # calculate cross-correlation
            for i in range(len(self)):
                cross_correlation = np.correlate(ys[:, ref], ys[:, i], mode='full')
                values[i] = np.argmax(cross_correlation)
            ref_value = values[ref]
            values -= ref_value
        elif mode == 'max':

            j_ref = np.argmax(self[ref].y)
            ref_value = self[ref].x[j_ref]
            for i in range(len(self)):
                j = np.argmax(self[i].y)
                values[i] = -(self[i].x[j] - ref_value)
        elif mode == 'fitted peaks':
            # check if all data has fit
            self._check_fit()
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum(fitted_peaks=True)

            # check if peaks are defined
            assert len(self[0].fit.peaks) > 0, 'Spectra does not have defined fitted peaks.\nMaybe use Spectra.fit_peaks(), Spectra.find_peaks().'

            # calculate peak
            peaks     = self.fit.get_peaks('p')
            ref_value = peaks[peak]['c'][0]
            values    = [c for c in peaks[peak]['c']]
            if None in values:
                raise ValueError(f'cannot calculate shifts.\npeak {peak} is not defined for all spectra.\ncenter of peak {peak}: {values}')
            values    = ref_value-np.array(values)
        elif mode == 'peak':
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum()
            # check if peaks are defined
            assert len(self[0].peaks) > 0, 'Spectra does not have defined peaks.\nMaybe use Spectra.find_peaks().'

            # calculate peak
            peaks     = self.get_peaks('p')
            ref_value = peaks[peak]['c'][0]
            values    = [c for c in peaks[peak]['c']]
            if None in values:
                raise ValueError(f'cannot calculate shifts.\npeak {peak} is not defined for all spectra.\ncenter of peak {peak}: {values}')
            values    = ref_value-np.array(values)
        else:
            raise ValueError('mode not valid.\nValid modes: cross-correlation, max, fitted peaks, peak.')

        # save calculated values ===============================================
        self._calculated_shift           = Spectrum(y=values)
        self._calculated_shift.mode      = mode
        self._calculated_shift.ref_value = ref_value

    def calculate_factor(self, mode='peak', peak=0, bkg_check=True):
        """Calculate mult. factor for spectra to be same height as the first spectrum.

        Args:
            mode (string, optional): method used to calculate the multiplicative factors.
                The current options are: 'max', 'delta', 'area', 'fitted peaks',
                'fitted peaks area', 'peak', or 'peak area'.
                For peak related modes, data must be fitted via
                Spectrum.fit_peak(), i. e., spectra must have a
                Spectrum.fit.peaks object defined.
            bkg_check (bool, optional): if mode is 'area', bkg_check=true prevents
                from comparing spectra with significant background (10 % of
                the maximum y) which can lead to wrong results.
                The algorithm uses find_peaks to find
                approximately whats is data and what is background, Then the
                background is said to be big if average(bkg) > max(y)*0.1. If
                background region cannot be found clearly, the criteria
                becames average(y) > max(y)*0.1.
            peak (int, optional): if peak related mode, this
                variable selects which peak to consider.

        Returns:
            None

        See Also:
            :py:func:`Spectra.set_factor`
        """
        # spectra will be alined to the first spectrum
        ref = 0

        # common variables =====================================================
        values = np.array([0.0]*len(self))

        # CALCULATION ==========================================================
        if mode == 'max':
            ref_value = max(self.data[ref].y)
            for i in range(len(self)):
                values[i] = ref_value/max(self.data[i].y)
        elif mode == 'delta':
            ref_value = max(self.data[ref].y) - max(self.data[ref].y)
            for i in range(len(self)):
                values[i] = ref_value/(max(self.data[i].y) - min(self.data[i].y))
        elif mode == 'area':
            if bkg_check:  # check if background is zero
                temp = [False]*len(self)
                for i in range(len(self)):
                    peaks, d = find_peaks(self[i].y, prominence=(max(self[i].y)-np.mean(self[i].y))*0.05, width=1)
                    bkg = list(self[i].y)
                    for j in range(len(peaks)):
                        del bkg[peaks[j]-int(d['widths'][j]*10):peaks[j]+int(d['widths'][j])*10]
                    if len(bkg) > 10:
                        temp[i] = np.mean(bkg) > max(self[i].y)*0.1
                    else:
                        temp[i] = np.mean(self[i].y)>max(self[i].y)*0.1
                if True in temp:
                    raise ValueError(f'some spectra have a background that seems to be too big for a reliable cross-corelation.\nSpectra with big bkg: {temp}\nYou can try using Spectra.floor()')

            ref_value = self.data[ref].area
            for i in range(len(self)):
                values[i] = ref_value/self.data[i].area
        elif mode == 'fitted peaks':
            # check if all data has fit
            self._check_fit()
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum(fitted_peaks=True)
            # check if peaks are defined
            assert len(self[0].fit.peaks) > 0, 'Spectra does not have defined fitted peaks.\nMaybe use Spectra.fit_peaks(), Spectra.find_peaks().'

            # calculate peak
            peaks     = self.fit.get_peaks('p')
            ref_value = peaks[peak]['amp'][0]
            values    = [c for c in peaks[peak]['amp']]
            if None in values:
                raise ValueError(f'cannot calculate multiplicative factors.\npeak {peak} is not defined for all spectra.\ncenter of peak {peak}: {values}')
            values    = ref_value/np.array(values)
        elif mode == 'fitted peaks area':
            # check if all data has fit
            self._check_fit()
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum(fitted_peaks=True)
            # check if peaks are defined
            assert len(self[0].fit.peaks) > 0, 'Spectra does not have defined fitted peaks.\nMaybe use Spectra.fit_peaks(), Spectra.find_peaks().'

            # calculate peak
            peaks     = self.fit.get_peaks('p')
            ref_value = peaks[peak]['area'][0]
            values    = [c for c in peaks[peak]['area']]
            if None in values:
                raise ValueError(f'cannot calculate multiplicative factors.\npeak {peak} is not defined for all spectra.\ncenter of peak {peak}: {values}')
            values    = ref_value/np.array(values)
        elif mode == 'peak':
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum()
            # check if peaks are defined
            assert len(self[0].peaks) > 0, 'Spectra does not have defined peaks.\nMaybe use Spectra.find_peaks().'

            # calculate peak
            peaks     = self.get_peaks('p')
            ref_value = peaks[peak]['amp'][0]
            values    = [c for c in peaks[peak]['amp']]
            if None in values:
                raise ValueError(f'cannot calculate multiplicative factors.\npeak {peak} is not defined for all spectra.\ncenter of peak {peak}: {values}')
            values    = ref_value/np.array(values)
        elif mode == 'peak area':
            # check if number of peaks for each spectrum is the same
            self._check_number_of_peaks_per_spectrum()
            # check if peaks are defined
            assert len(self[0].peaks) > 0, 'Spectra does not have defined peaks.\nMaybe use Spectra.find_peaks().'

            # calculate peak
            peaks     = self.get_peaks('p')
            ref_value = peaks[peak]['area'][0]
            values    = [c for c in peaks[peak]['area']]
            if None in values:
                raise ValueError(f'cannot calculate multiplicative factors.\npeak {peak} is not defined for all spectra.\ncenter of peak {peak}: {values}')
            values    = ref_value/np.array(values)
        else:
            raise ValueError('mode not valid.\nValid modes: max, delta, area, fitted peaks, fitted peaks area, peak, peak area.')

        # save calculated values ===============================================
        self._calculated_factor           = Spectrum(y=values)
        self._calculated_factor.mode      = mode
        self._calculated_factor.ref_value = ref_value

    def calculate_offset(self, ranges=None):
        """Calculate additive factor for spectra to be verticaly aligned with the first spectrum.

        Args:
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                The additive factor will be calculated for spectra such the average value of
                data inside ranges is the same as for the first spectrum.


        Returns:
            None

        See Also:
            :py:func:`Spectra.set_offset`, :py:func:`Spectra.floor()`,
        """
        # spectra will be alined to the first spectrum
        ref = 0

        # common variables =====================================================
        values = np.array([0.0]*len(self))

        # CALCULATION ==========================================================
        s_ref = self[ref].extract(ranges)
        ref_value = np.mean(s_ref.y)
        for i in range(len(self)):
            s = self[i].extract(ranges)
            values[i] = ref_value-np.mean(s.y)

        # save calculated values ===============================================
        self._calculated_offset           = Spectrum(y=values)
        # self._calculated_offset.mode      = mode
        self._calculated_offset.ref_value = ref_value

    def calculate_calib(self, start=None, stop=None, values=None, mode='cross-correlation', peak=0, bkg_check=True, deg=1):
        """Calculate calibration factor via :py:func:`Spectra.calculate_shift()`.

        Args:
            start (number): value used for measuring the first spectrum.
            stop (number): value used for measuring the last spectrum.
            values (list): value list. It overwrites start and stop. Must be the
                same length as number of spectra.
            mode (string, optional): method used to calculate the shifts.
                Default is 'cross-correlation'. This
                parameter is passed to :py:func:`Spectra.calculate_shift()`.
                see :py:func:`Spectra.calculate_shift()` for more.
            bkg_check (bool, optional): see :py:func:`Spectra.calculate_shift()` for more.
            peak (int, optional): see :py:func:`Spectra.calculate_shift()` for more.
            deg (int, ooptional): degree of the fitting polynomial. Default is 1.

        Returns:
            calibration factor (number)
        """
        # check number of values matches the numbre of spectra
        if values is None:
            values = np.linspace(start, stop, len(self))
        if len(self) != len(values):
            raise ValueError(f'number of values ({len(values)}) do not match the number of spectra ({len(self)})')

        # CALCULATION ==========================================================
        self.calculate_shift(mode=mode, peak=peak, bkg_check=bkg_check)
        centers = -(self.calculated_shift.y + self.calculated_shift.ref_value)

        if mode in cc:
            centers = centers*self.step

        # save calculated values ===============================================
        self._calculated_calib      = Spectrum(x=values, y=centers)
        self._calculated_calib.mode = mode
        # self._calculated_calib.ref_value = ref_value
        popt, model = self.calculated_calib.polyfit(deg=deg)
        self._calculated_calib.popt = popt

        return 1/popt[0]


    def calculate_sum(self):
        """Returns Spectrum object with the sum of all spectra.

        Returns:
            :py:class:`Spectra` object.

        Note:
            All spectra have to have the same x-coordinates. This is verified
            before summing up the spectra.
        """
        # check x
        if self.x is None:
            self.check_same_x()

        # calculate sum
        y = np.zeros(len(self.x))
        for i in range(len(self)):
            y += self.data[i].y

        s = Spectrum(x=self.x, y=y)

        # copy user defined attributes
        for attr in self.__dict__:
            if attr.startswith('_') == False:
                s.__setattr__(attr, self.__dict__[attr])

        return s

    def calculate_map(self, axis=0):
        """Return image.

        Args:
            axis (int or string, optional): Axis along which x axis will be laid
                down. By default, x axis will run across the vertical (0) direction.

        Returns:
            :py:class:`Image`.
        """
        axis = _axis_interpreter(axis)

        y, ys = self._gather_ys()
        x = None

        if axis == 1:
            ys = ys.transpose()
            x = y
            y = None

        im = Image(data=ys)
        im.x_centers = x
        im.y_centers = y

        # copy user defined attributes
        for attr in self.__dict__:
            if attr.startswith('_') == False:
                im.__setattr__(attr, self.__dict__[attr])

        return im


    def get_peaks(self, key='peak'):
        """Return a list with peak data from all spectra.

        Args:
            key (string, optional): if 'peaks', list is organized by peak
                indices. If 'spectrum', list is organized by spectrum
                indices.

        Returns:
            dictionary
        """
        if key.startswith('p'):
            # peaks_attr = ['amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2']
            # n_peaks = max([len(s.peaks) for s in self])
            # if n_peaks == 0:
            #     return {}
            # return {peak: {k: [s.peaks[peak][k] if peak in s.peaks else None for s in self] for k in peaks_attr} for peak in range(n_peaks)}

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

    def find_peaks(self, prominence=None, width=4, moving_average_window=8):
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
            s.find_peaks(prominence=prominence, width=width, moving_average_window=moving_average_window)

    def fit_peaks(self, asymmetry=False, fixed_m=0, offset=True, amp_bounds=(0, 3), c_bounds=(-2, 2), fwhm_bounds=(0, 2), verbose=False):
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
        if verbose: print('Fitting...\n')
        for i in range(len(self)):
            if verbose: print(f'spectrum {i} ===============================')
            self[i].fit_peaks(asymmetry=asymmetry, fixed_m=fixed_m, offset=offset, amp_bounds=amp_bounds, c_bounds=c_bounds, fwhm_bounds=fwhm_bounds, verbose=verbose)
            if verbose: print('='*20)
            if verbose: print('\n')


    def polyfit(self, deg=2):
        """Fit data recursively with a polynomial. Wrapper for `numpy.polyfit()`_.

        Args:
            deg (int, optional): degree of the fitting polynomial. Default is 2.

        Returns:
            list with polynomial coefficients, highest power first.
            list with Model function f(x).

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        popt = [0]*len(self)
        model = [0]*len(self)
        for i in range(len(self)):
            popt[i], model[i] = self[i].polyfit(deg=deg)

        return popt, model




# %%
