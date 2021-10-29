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
import decimal
import warnings
from scipy.optimize import curve_fit
from scipy.signal import find_peaks

def n_decimal_places(number):
    """Return the number of decimal places of a number.

    Args:
        number (float or int): number.

    Returns:
        number of decimal places in number.
    """

    return abs(decimal.Decimal(str(number)).as_tuple().exponent)

def n_digits(number):
    """Return the number of digits of a number.

    Args:
        number (float or int): number.

    Returns:
        a tuple with number of digits and number of decimal places.
    """
    if n_decimal_places(number) != 0:
        return (len(str(int(np.around(number))) ), n_decimal_places(number))
    else:
        return (len(str(int(np.around(number))) ), 0)

def index(x, value):
    """Returns the index of the element in array which is closest to value.

    Args:
        x (list or array): 1D array.
        value (float or int): value.

    Returns:
        index (int)
    """
    return np.argmin(np.abs(np.array(x)-value))

def movingaverage(x, window_size):
    """Returns the moving average of an array.

    Args:
        x (list or array): 1D array.
        window_size (int): number of points to average.

    Returns:
        array of lenght given by (len(x)-window_size+1).

    Example:
        >>> x = [0,1,2,3,4,5,6,7,8,9]
        >>> print(movingaverage(x, 1))
        [0. 1. 2. 3. 4. 5. 6. 7. 8. 9.]
        >>> print(movingaverage(x, 2))
        [0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]
        >>> print(movingaverage(x, 3))
        [1. 2. 3. 4. 5. 6. 7. 8.]
        >>> print(movingaverage(x, 4))
        [1.5 2.5 3.5 4.5 5.5 6.5 7.5]

    """
    if window_size < 1:
        raise ValueError('window_size must be a positive integer (> 1).')

    x = np.array(x)
    window = np.ones(int(window_size))/float(window_size)

    return np.convolve(x, window, 'valid')

def extract(x, y, ranges):
    """Returns specifc data ranges from x and y.

    Args:
        x (list or array): 1D reference vector.
        y (list or array): 1D y-coordinates or list of several data sets sharing
            the same x-coordinates.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from ref.

    Returns:
        x and y arrays.


    Examples:

        >>> x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
        >>> y = np.array(x)**2
        >>> ranges = ((0, 3), (7.5, 9))
        >>> x_sliced, y_sliced = am.extract(x, y, ranges)
        >>> print(x_sliced)
        [0 1 2 3 8 9]
        >>> print(y_sliced)
        [0 1 4 9 64 81]

    """
    x = np.array(x)
    y = np.array(y)

    try:
        choose_range = []
        for x_init, x_final in ranges:
            choose_range.append(np.logical_and(x>=x_init, x<=x_final))
        choose_range = [True if x == 1 else False for x in sum(choose_range)]
    except TypeError:
        choose_range = []
        x_init, x_final = ranges
        choose_range.append(np.logical_and(x>=x_init, x<=x_final))
        choose_range = [True if x == 1 else False for x in sum(choose_range)]

    if len(y.shape) > 1:
        s = []
        for i in range(y.shape[0]):
            s.append(np.hstack(y[i][choose_range]))
        return np.hstack(x[choose_range]), s
    else:
        return np.hstack(x[choose_range]), np.hstack(y[choose_range])

def shift(x, y, value, mode='hard'):
    """Shift (x, y) data.

    Args:
        x (list or array): 1D array.
        y (list or array): 1D array.
        value (float or int): shift value.
        mode (string, optional):
            #. ``mode='x'`` or ``mode='hard'``
                y is fully preserved while x is shifted.
            #. ``mode='y'``, ``'interp'``, or ``'soft'``
                x is preserved while y is interpolated with a shift
            #. ``mode='roll'``,
                x and y are preserved and y elements are just rolled along the
                array (in this case ``shift`` value must be an integer).

    Returns:
        Shifted x and y.

    Warning:
        It is always better to use ``mode='hard'`` or ``'roll'`` since the form of y is fully
        preserved (no interpolation). After applying a shift using the ``mode='interp'``,
        one can apply a
        'inverse' shift to retrieve the original data. The diference between the retrieved
        y data and the original data will give an ideia of the information loss
        caused by the interpolation.
    """
    x = copy.deepcopy(np.array(x))
    y = copy.deepcopy(np.array(y))

    if mode == 'y' or mode == 'interp' or mode=='soft':
        y = np.interp(x, x + value, y)

    elif mode == 'x' or mode == 'hard':
        x = np.array(x) + value

    elif mode == 'roll' or mode == 'rotate':
        try:
            if value.is_integer():
                y = np.roll(y, int(value))
            else:
                raise ValueError("value must be an interger for mode='roll'.")
        except AttributeError:
            y = np.roll(y, int(value))
        if value > 0:
            y[:int(value)] = 0
        elif value < 0:
            y[int(value):] = 0
    else:
        raise ValueError('mode not recognized')

    return x, y

def fwhmVoigt(x, amp, c, w, m):
    r"""Pseudo-voigt curve.

    .. math:: y(x) = A \left[ m  \frac{w^2}{w^2 + (x-c)^2}   + (1-m) e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM
    :param m: Factor from 1 to 0 of the lorentzian amount
    :return: :math:`y(x)`
    """
    lorentz = fwhmLorentz(x, 1, c, w)
    gauss = fwhmGauss(x, 1, c, w)

    return amp*(m*lorentz + (1-m)*gauss)

def Lorentz(x, gamma, c):
    r"""Cauchy–Lorentz distribution.

    .. math:: y(x) = \frac{1}{\pi \gamma} \frac{\gamma^2}{\gamma^2 + (x-c)^2}

    where,

    .. math:: \text{Amplitude }= \frac{1}{\pi \gamma}

    and,

    .. math:: \text{fwhm }= 2 \gamma

    and,

    .. math:: \text{Area }= 1

    :param x: x array
    :param gamma: Scale factor
    :param c: Center
    :return: :math:`y(x)`
    """
    return (1/(np.pi*gamma))*((gamma**2)/(gamma**2 + (x-c)**2))

def fwhmLorentz(x, amp, c, w):
    r"""Cauchy–Lorentz distribution.

    .. math:: y(x) = \text{amp } \frac{w^2}{w^2 + (x-c)^2}

    where,

    .. math:: \text{Area }= \text{amp } \pi w

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """
    return amp*(np.pi*w) * Lorentz(x, gamma=w, c=c)
    # return A*((w**2)/(w**2 + 4* (x-c)**2))

def Gauss(x, amp, c, sigma):
    r"""Gaussian distribution.

    .. math:: y(x) = \text{amp } e^{-\frac{(x-c)^2}{2 \sigma^2}}

    where,

    .. math:: \text{Area }= \sqrt{2 \pi} \text{ amp } |\sigma|

    and,

    .. math:: \text{fwhm }= 2 \sqrt{2 \ln(2)} \sigma


    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param sigma: standard deviation
    :return: :math:`y(x)`
    """
    return amp*np.exp(-(x-c)**2/(2*sigma**2))

def fwhmGauss(x, amp, c, w):
    r"""Gaussian distribution.

    .. math:: y(x) = \text{amp } e^{-\frac{4 \ln(2) (x-c)^2}{w^2}}

    where,

    .. math:: \text{Area }= \frac{\sqrt{\pi} \text{ amp } w}{2 \sqrt{\ln(2)}}

    :param x: x array
    :param A: Amplitude
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """
    return Gauss(x, amp, c, w/(2*np.sqrt(2*np.log(2))))
    # return A*np.exp((-4*np.log(2)*((x-c)**2))/(w**2))

def peak_fit(x, y, guess_c=None, guess_A=None, guess_w=None, guess_offset=0, fixed_m=False, asymmetry=False):
    r"""Simple peak fit function. Data is fitted with a pseudo-voigt curve.

    .. math:: y(x) = A \left[ m \frac{w^2}{w^2 + (x-c)^2}   + (1-m) e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

    Args:
        x (list or array): 1D array x-coordinates.
        y (list or array): 1D array y-coordinates.
        guess_c (float or int, optional): guess Center. If None, it will be guessed by the
            position of ``guess_A``.
        guess_A (float or int, optional): guess Amplitude. If None, it will be guessed by the
            the maximum y-coordinate.
        guess_w (float or int, optional): guess FWHM. If None, it will be guessed as 10% of
            ``guess_c``
        guess_offset (float or int, optional): guess Offset. If None, it will be guessed as zero [0].
        fixed_m (False or number): `Factor from 1 to 0 of the lorentzian amount.
            If False, ``m`` will be a fitting parameter. If
            ``fixed_m=<number>``, ``<number>`` will be used for ``m``.
        asymmetry (Boolean, , optional). If True, peak asymmetry is taken into account by fiting first
            half of the peak with a different ``w`` and ``m`` than the second half. The optimal ``w`` parameter
            returned will be the sum of the ``w`` of the first and second halfs.

    Returns:
        1) 2 column (x, y) array with the fitted peak.
        2) 2 column (x, y) array with "Smoothed" fitted peak. This is just the
            fitted peak array with a linear interpolation with 100 times more data points.
        3) An array with the optimized parameters for Amplitude, Center, FWHM and offset.
        4) One standard deviation errors on the parameters

    See Also:
        :py:func:`backpack.model_functions.fwhmVoigt`

    Example:
        >>> import matplotlib.pyplot as plt
        >>> from backpack.model_functions import fwhmGauss
        >>> x = np.linspace(0, 100, 1000)
        >>> amp = 1
        >>> w = 10
        >>> c = 25
        >>> y = fwhmGauss(x, amp, c, w) + np.random.normal(-0.1, 0.1, 1000)
        >>> _, smooth, popt, err = am.peak_fit(x, y)
        >>> print(f'A = {popt[0]} +/- {err[0]}')
        A = 1.0187633097698565 +/- 0.015832212669397393
        >>> print(f'c = {popt[1]} +/- {err[1]}')
        c = 24.940397238212615 +/- 0.0674408860045828
        >>> print(f'w = {popt[2]} +/- {err[2]}')
        w = 9.952479117885305 +/- 0.24269432941302957
        >>> print(f'offset = {popt[3]} +/- {err[3]}')
        offset = -0.10529818707281721 +/- 0.032115907889677234
        >>> print(f'm = {popt[4]} +/- {err[4]}')
        m = 2.4612897990598044e-14 +/- 0.005287100795322343
        >>> plt.scatter(x, y)
        >>> plt.plot(smooth[:, 0], smooth[:, 1], color='r', lw=3)
        >>> plt.show()

        .. image:: _static/peak_fit.png
            :width: 600
            :align: center
    """
    start = x[0]
    stop = x[-1]

    if guess_A is None:
        guess_A = max(y)

    if guess_c is None:
        guess_c = x[index(y, guess_A)]

    if guess_w is None:
        guess_w = 0.1*guess_c

    if not fixed_m or fixed_m != 0:  # variable m
        if asymmetry:
            p0 = [guess_A, guess_c, guess_w, 0.5, guess_w, 0.5, guess_offset]
            def function2fit(x, A, c, w1, m1, w2, m2, offset):
                f = np.heaviside(x-c, 0)*fwhmVoigt(x, A, c, w1, m1) + offset +\
                    np.heaviside(c-x, 0)*fwhmVoigt(x, A, c, w2, m2)
                return f
            bounds=[[0,      start,   0,      0, 0,      0, -np.inf],
                    [np.inf, stop,  np.inf, 1, np.inf, 1, np.inf]]
        else:
            p0 = [guess_A, guess_c, guess_w, 0.5, guess_offset]
            def function2fit(x, A, c, w, m, offset):
                return fwhmVoigt(x, A, c, w, m) + offset
            bounds=[[0,      start,   0,      0, -np.inf],
                    [np.inf, stop,  np.inf, 1, np.inf]]

    else:
        if fixed_m > 1:
            fixed_m = 1
        elif fixed_m < 0:
            fixed_m = 0
        if asymmetry:
            p0 = [guess_A, guess_c, guess_w, guess_w, guess_offset]
            def function2fit(x, A, c, w1, w2, offset):
                f = np.heaviside(x-c, 0)*fwhmVoigt(x, A, c, w1, fixed_m) + offset +\
                    np.heaviside(c-x, 0)*fwhmVoigt(x, A, c, w2, fixed_m)
                return f
            bounds=[[0,      start,   0,      0,  -np.inf],
                    [np.inf, stop,  np.inf, 1,   np.inf]]
        else:
            p0 = [guess_A, guess_c, guess_w, guess_offset]
            def function2fit(x, A, c, w, offset):
                return fwhmVoigt(x, A, c, w, fixed_m) + offset
            bounds=[[0,      start,   0,      -np.inf],
                    [np.inf, stop,  np.inf, np.inf]]

    # Fit data
    popt, pcov = curve_fit(function2fit, x, y, p0,  # sigma = sigma,
                           bounds=bounds)
    err = np.sqrt(np.diag(pcov))  # One standard deviation errors on the parameters

    # smooth data
    arr100 = np.zeros([100*len(x), 2])
    arr100[:, 0] = np.linspace(x[0], x[-1], 100*len(x))
    arr100[:, 1] = function2fit(arr100[:, 0],  *popt)

    if fixed_m:
        if asymmetry:
            popt_2 = (popt[0], popt[1], popt[2]/2+popt[4]/2, popt[-1])
        else:
            popt_2 = (popt[0], popt[1], popt[2], popt[-1])
    else:
        if asymmetry:
            popt_2 = (popt[0], popt[1], popt[2]/2+popt[4]/2, popt[-1], popt[-2])
        else:
            popt_2 = (popt[0], popt[1], popt[2], popt[-1], popt[-2])
    return function2fit(x, *popt), arr100, popt_2, err


def save_data(obj, filepath='./untitled.txt', add_labels=True, data_format='% .10e', header='', footer='', delimiter=', ', comment_flag='# ', newline='\n', checkOverwrite=False):
    r"""Save an array or a dictionary in a txt file.

    Args:
        obj (dict, list, or numpy.array): data to be saved to a file. If obj is
        a dictonary, use ``*`` in front of a key to do not save it to the file.
        filepath (str or pathlib.Path, optional): path to save file.
        add_labels (bool, optional): When obj is a dictonary, ``add_labels=True``
            makes the dict keys to be added to the header as label for each data column.
        data_format (string, or list, optional): If obj is a list, fmt can also
            be a list where each fmt element is associated with a column. If
            obj is a dict, fmt can also be a dict with same keys of obj. Then,
            each fmt value is associated with the corresponding column.

            See `np.savetxt <https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html?highlight=savetxt#numpy.savetxt>`_ documentation::

                fmt = (%[flag]width[.precision]specifier)

            * flag can be: '-' for left justify, '+', whch forces to precede

            * result with + or -, or '0' to Left pad the number with zeros
              instead of space (see width).

            * width is the minimum number of characters to be printed.

            * precision is tipically the number of significant digits

            * specifier is the type of notation. Tipically, either 'e' for
              scientific notation of 'f' for decimal floating point.

            * a common fmt strings is: '%.3f' for 3 decimal places.

        header (str, oprional): string that will be written at the beggining of
            the file (comment flag is added automatically).
        footer (str, oprional): string that will be written at the end of the
            file (comment flag is added automatically).
        delimiter (str, optional): The string used to separate values.
        comment_flag (str, optional): string that flag comments.
        newline (str, optional): string to indicate new lines.
        checkOverwrite (bool, optional): if True, it will check if file exists
            and ask if user want to overwrite file.

    See Also:
        :py:func:`load_data`
    """
    filepath = Path(filepath)

    if checkOverwrite:
        if filepath.exists() == True:
            if filepath.is_file() == True:
                if query_yes_no('File already exists!! Do you wish to ovewrite it?', 'yes') == True:
                    pass
                else:
                    warnings.warn('File not saved.')
                    return
            else:
                warnings.warn('filepath is pointing to a folder. Saving file as Untitled.txt')
                filepath = filepath/'Untitled.txt'

    if type(obj) == dict:
        # remove keys that start with star (*)
        obj2 = {key: obj[key] for key in obj if str(key).startswith('*') is False}

        # col labels
        if add_labels:
            if not header == '' and not header.endswith('\n'):
                header += '\n'
            for key in obj2:
                header += str(key) + f'{delimiter}'
            header = header[:-(len(delimiter))]

        # dict to array
        obj = []
        for key in obj2:
            obj.append(obj2[key])
        obj = np.array(obj).transpose()

    np.savetxt(filepath, obj, fmt=data_format, delimiter=delimiter, newline=newline, header=header, footer=footer, comments=comment_flag)



class photon_events():
    """Creates a ``photon_event`` class type object to deal with photon events lists.

    Args:
        events (list or array): two (x, y) or three (x, y, intensity) list or array with
            photon events.
        x_max (float, optional): maximum x value. If ``None``, it will be infered
            by the data.
        y_max (float, optional): maximum y value. If ``None``, it will be infered
            by the data.

    """


    def __init__(self, events, x_max=None, y_max=None):
        # basic attr
        self.events     = None
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

        self.load(events, x_max=x_max, y_max=y_max)


    def load(self, events, x_max=None, y_max=None):
        """Load photon events data.

        args:
            events (list or array): two (x, y) or three (x, y, intensity) list or array with
                photon events.
            x_max (float, optional): maximum x value. If ``None``, it will be infered
                by the data.
            y_max (float, optional): maximum y value. If ``None``, it will be infered
                by the data.

        returns:
            None

        See Also:
            :py:func:`photon_events.save`.
        """
        # check if data is the right Format
        if events.shape[1] == 2 or events.shape[1] == 3:
            self.events = events
        else: raise ValueError("data must have 2 or 3 columns.")

        # infer x_max and y_max if necessary
        if x_max is None:
            self.x_max = max(self.events[:, 0])
        else:
            self.x_max = copy.deepcopy(x_max)

        if y_max is None:
            self.y_max = max(self.events[:, 1])
        else:
            self.y_max = copy.deepcopy(y_max)


    def save2file(self, filepath, delimiter=','):
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
        if self.events.shape[1] == 3:
            header += f'x y I'
        else:
            header += f'x y'
        save_data(self.events, filepath=Path(filepath), delimiter=delimiter, header=header)


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

        self.hist, self.x_edges, self.y_edges = np.histogram2d(self.events[:, 0],
                                                                 self.events[:, 1],
                                                                 bins=(x_bins, y_bins),
                                                                 # weights=self.data[:, 2],
                                                                 range=((0, self.x_max), (0, self.y_max))
                                                                )

        self.x_centers = movingaverage(self.x_edges, window_size=2)
        self.y_centers = movingaverage(self.y_edges, window_size=2)


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
        y_centers, ref_column = extract(self.y_centers, self.hist[ref], ranges=ranges)
        for i in range(self.hist.shape[0]):

            y_centers, column     = extract(self.y_centers, self.hist[i], ranges=ranges)

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

        Example:
            f = lambda x, y: (x, y**2)

        args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        returns:
            None
        """
        self.events[:, 0], self.events[:, 1] = f(self.events[:, 0], self.events[:, 1])
        self.x_max, self.y_max  = f(self.x_max, self.y_max)
        self.set_binning(bins=self.bins)


    def calculate_spectrum(self, y_bins=None, y_bins_size=None, x_type='bins'):
        """Sum the photon events in the x direction.

        args:
            y_bins (int, optional): number of y bins. If None, the current binning is used.
            bins_size (int or tuple, optional): size of the y bins. This overwrites
                the argument ``y_bins``. If None, the current binning is used.
            x_type (string, optional): Type of the x axis. If 'bins', the x axis
                is numbered from 0 to the number of y_bins. If 'lenght', the
                detector x values are used.

        returns:
            spectrum type
        """


        if y_bins_size is None and y_bins is None:
            if x_type == 'bins':
                x = np.arange(0, y_bins)
            elif x_type == 'lenght':
                x = self.y_centers
            else:
                raise ValueError('x_type can only be `bins` or `lenght`.')

            self.spectrum = spectrum(data=np.vstack((x, sum(self.hist))).transpose())

        else:
            temp = photon_events(events=self.events)
            if y_bins_size is None:
                temp.set_binning(bins=(1, y_bins))
            else:
                temp.set_binning(bins_size=(self.x_max+1, y_bins_size))

            if x_type == 'bins':
                x = np.arange(0, y_bins)
            elif x_type == 'lenght':
                x = temp.y_centers
            else:
                raise ValueError('x_type can only be `bins` or `lenght`.')

            self.spectrum = spectrum(data=np.vstack((x, sum(temp.hist))).transpose())

        return self.spectrum


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
            fig = plt.figure()
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

        ax.plot(self.events[:, 0], self.events[:, 1],
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
            fig = plt.figure()
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
            fig = plt.figure()
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
            fig = plt.figure()
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


class spectrum():
    """Creates a ``spectrum`` class type object to deal with (x, y) data types.

    Args:
        data (list or array, optional): three column list (or array) with photon
            events. Column order should be x, y, and intensity.
        x (list or array, optional): x values (1D list/array). Overwrites `data`.
        y (list or array, optional): y values (1D list/array). Overwrites `data`.

    """

    def __init__(self, data=None, x=None, y=None):
        self.data = None

        self.load(data=data, x=x, y=y)

    def load(self, data=None, x=None, y=None):
        """Load spectrum..

        Data must have two columns, x (energy or distance) and y (intensity).

        args:
            data (list or array, optional): two column list (or array).
            x (list or array, optional): x values (1D list/array). Overwrites `data`.
            y (list or array, optional): y values (1D list/array). Overwrites `data`.

        returns:
            None

        See Also:
            :py:func:`spectrum.save`.
        """
        if x is None and y is None:
            if data is None:
                warnings.warn('No data to load.', stacklevel=2)
                return
            else:
                self.data = copy.deepcopy(data)
        else:
            if x is None:
                x = np.arange(0, len(y))
            self.data = np.vstack((x, y)).transpose()

    def set_x(self, x):
        if len(x) != len(data[:, 0]):
            raise ValueError('Lenght of x is not compatible with data.')
        self.data[:, 0] = x

    def get_x(self):
        return self.data[:, 0]

    def set_y(self, y):
        if len(y) != len(data[:, 1]):
            raise ValueError('Lenght of y is not compatible with data.')
        self.data[:, 1] = y

    def get_y(self):
        return self.data[:, 1]

    def save2file(self, filepath, delimiter=',', header=None):
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
        save_data(self.data, filepath=Path(filepath), delimiter=delimiter, header=header)

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
        x, y = shift(self.data[:, 0], self.data[:, 1], value=value, mode=mode)
        self.data = np.column_stack((x, y))

    def calib(self, dispersion=1, zero=0):
        """Calibrate data. Dispersion if adjusted before adjusting the zero value.

        args:
            dispersion (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].
            zero (number or string, optional): if number, the spectrum is shifted
                so the value given in `zero` is now set to 0 (zero). If `zero=elastic_peak`,
                the position of the elastic peak is guess by `spectrum.guess_elastic_peak()`
                and the center of the elastic peak is set to 0 (zero).

        returns:
            None
        """
        if type(zero) == int or type(zero) == float:
            offset = -zero
            f = lambda x, y: (x*dispersion + offset, y)
        elif zero == 'elastic_peak' and dispersion == 1:
            _, popt, _ = self.guess_elastic_peak()
            offset = -popt[1]
            f = lambda x, y: (x*dispersion + offset, y)
        elif zero == 'elastic_peak':
            f = lambda x, y: (x*dispersion, y)
            self.apply_correction(f=f)
            _, popt, _ = self.guess_elastic_peak()
            offset = -popt[1]
            f = lambda x, y: (x + offset, y)
        else:
            raise ValueError('zero must be a number.')

        self.apply_correction(f=f)

    def guess_elastic_peak(self, min_points=8, tail_factor=(1, 8), verbose=False):
        height_min = np.mean(self.data[:, 1])*1.5
        peaks, d = find_peaks(self.data[:, 1],
                              height=None,
                              threshold=None,
                              distance=None,
                              prominence=(height_min, None),
                              width=min_points,
                              wlen=None,
                              rel_height=0.5,
                              plateau_size=None)
        guess_c = self.data[:, 0][peaks[-1]]
        guess_w = d['widths'][-1]*(self.data[1, 0]-self.data[0, 0])
        guess_A = d['prominences'][-1]

        try:
            len(tail_factor) == 2
        except TypeError:
            tail_factor = (tail_factor, tail_factor)
        x_init      = index(self.data[:, 0], self.data[:, 0][peaks[-1]]-guess_w*(1+tail_factor[0]))
        x_final     = index(self.data[:, 0], self.data[:, 0][peaks[-1]]+guess_w*(1+tail_factor[1]))

        try:
            _, smooth, popt, err = peak_fit(self.data[x_init:x_final, 0], self.data[x_init:x_final, 1], guess_c=guess_c, guess_A=guess_A, guess_w=guess_w, guess_offset=0, fixed_m=False, asymmetry=True)
        except RuntimeError:
            if verbose:
                warnings.warn('Cannot find elastic peak. Reducing tail parameters.', stacklevel=2)
            tail_factor_new = (tail_factor[0]-1, tail_factor[0]-1)
            if tail_factor_new[0] < 0 or tail_factor_new[1] < 0:
                tail_factor_new = tail_factor
                min_points = min_points-1
                if verbose:
                    warnings.warn('Cannot reduce tail parameters. Reducing min_points', stacklevel=2)
                if min_points <= 0:
                    raise RuntimeError('cannot find elastic peak')
            smooth, popt, err = self.guess_elastic_peak(min_points=min_points, tail_factor=tail_factor_new)
        return smooth, popt, err

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
        self.data[:, 0], self.data[:, 1] = f(self.data[:, 0], self.data[:, 1])

    def floor(x_value=None):
        if x is None:
            i = 0
        else:
            i = am.index(self.data[:, 0], x_value)
        f = lambda x, y: (x, y - y[i])
        self.apply_correction(f)

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

    def peak_fit(ranges=None, guess_c=None, guess_A=None, guess_w=None, guess_offset=0, fixed_m=False, asymmetry=False):
            r"""Simple peak fit function. Data is fitted with a pseudo-voigt curve.

            .. math:: y(x) = A \left[ m \frac{w^2}{w^2 + (x-c)^2}   + (1-m) e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

            Args:
                ranges (list): a pair of values or a list of pairs. Each pair represents
                    the start and stop of a data range from ref. If `None`, full data is used.
                guess_c (float or int, optional): guess Center. If None, it will be guessed by the
                    position of ``guess_A``.
                guess_A (float or int, optional): guess Amplitude. If None, it will be guessed by the
                    the maximum y-coordinate.
                guess_w (float or int, optional): guess FWHM. If None, it will be guessed as 10% of
                    ``guess_c``
                guess_offset (float or int, optional): guess Offset. If None, it will be guessed as zero [0].
                fixed_m (False or number): `Factor from 1 to 0 of the lorentzian amount.
                    If False, ``m`` will be a fitting parameter. If
                    ``fixed_m=<number>``, ``<number>`` will be used for ``m``.
                asymmetry (Boolean, , optional). If True, peak asymmetry is taken into account by fiting first
                    half of the peak with a different ``w`` and ``m`` than the second half. The optimal ``w`` parameter
                    returned will be the sum of the ``w`` of the first and second halfs.

            Returns:
                1) 2 column (x, y) array with the fitted peak.
                2) 2 column (x, y) array with "Smoothed" fitted peak. This is just the
                    fitted peak array with a linear interpolation with 100 times more data points.
                3) An array with the optimized parameters for Amplitude, Center, FWHM and offset.
                4) One standard deviation errors on the parameters
            """
            if ranges is None:
                ranges = [[min(self.data[ref].data[:, 0]), max(self.data[ref].data[:, 0])]]
            x, y = extract(self.data[:, 0], self.data[:, 1], ranges=ranges)
            return peak_fit(x=x, y=y, guess_c=guess_c, guess_A=guess_A, guess_w=guess_w, guess_offset=guess_offset, fixed_m=fixed_m, asymmetry=asymmetry)

    def plot(self, ax=None, normalized=False, vertical_increment=0, shift=0, factor=1, **kwargs):
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

        # if 'marker' not in kwargs:
        #     kwargs['marker'] = 'o'
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
        data (list or array, optional): list of :py:class:`spectrum` objects.
    """

    def __init__(self, data=None):
        # basic attr
        self.data = None

        # shift attr
        self.shift_mode = None
        self.shifts = None
        self.shift_ranges = None

        # sum attr
        self.sum = None

        self.load(data=data)

    def load(self, data=None):
        """Load spectra.

        args:
            data (list or array, optional): list of :py:class:`spectrum` objects.

        returns:
            None

        See Also:
            :py:func:`spectra.save`.
        """
        if data is None:
            # warnings.warn('No filepath or data to load.', stacklevel=1)
            self.data = []
            return
        else:
            for s in data:
                print(type(s))
                if isinstance(s, spectrum) == False:
                    raise ValueError('all entries must be of type brixs.spectrum.')
            # if sum(isinstance(s, spectrum) for s in data) != len(data):
            #     raise ValueError('all entries must be of type brixs.spectrum.')
            self.data = copy.deepcopy(data)

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
            if isinstance(s, spectrum):
                self.data.append(s)
                return
            else:
                raise ValueError('spectrum must be of type brixs.spectrum.')
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
        n_digits = n_digits(self.get_spectra_count()-1)[0]
        for idx, s in enumerate(self.data):
            filename = f'{prefix}' + f'{idx}'.zfill(n_digits) + f'{suffix}'
            s.save(filepath=folderpath/filename, delimiter=delimiter, header=header)

    def calculate_shifts(self, ref=0, mode='cross-correlation', output='int', ranges=None, verbose=False, **kwargs):
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
        self.check_step_x()
        self.shifts = np.zeros(len(self.data))

        if ranges is None:
            ranges = [[min(self.data[ref].data[:, 0]), max(self.data[ref].data[:, 0])]]

        if mode == 'cross-correlation':
            _, y_ref = extract(self.data[ref].data[:, 0], self.data[ref].data[:, 1], ranges=ranges)
            for i, s in enumerate(self.data):
                _, y = extract(s.data[:, 0], s.data[:, 1], ranges=ranges)
                cross_correlation = np.correlate(y_ref, y, mode='Same')
                self.shifts[i] = np.argmax(cross_correlation)

                if verbose:
                    print(f'spectrum {i} shift calculated!')

            if output == 'float':
                step = self.data[ref].data[1, 0], self.data[ref].data[0, 1]
                for j, value in enumerate(self.shifts):
                    self.shifts[j] = value*step

        elif mode == 'max':
            _, y_ref = extract(self.data[ref].data[:, 0], self.data[ref].data[:, 1], ranges=ranges)
            j_ref = np.argmax(y_ref)
            for i, s in enumerate(self.data):
                _, y = extract(s.data[:, 0], s.data[:, 1], ranges=ranges)
                self.shifts[i] = j_ref - np.argmax(y)

                if verbose:
                    print(f'spectrum {i} shift calculated!')

            if output == 'float':
                step = self.data[ref].data[1, 0], self.data[ref].data[0, 1]
                for j, value in enumerate(self.shifts):
                    self.shifts[j] = value*step


        elif mode == 'elastic_peak':
            x_ref, y_ref = extract(self.data[ref].data[:, 0], self.data[ref].data[:, 1], ranges=ranges)
            s_ref = spectrum(x=x_ref, y=y_ref)
            _, popt_ref, _ = s_ref.guess_elastic_peak(**kwargs)
            for i, s in enumerate(self.data):
                x, y = extract(s.data[:, 0], s.data[:, 1], ranges=ranges)
                s_temp = spectrum(x=x, y=y)
                _, popt, _ = s_temp.guess_elastic_peak(**kwargs)
                self.shifts[i] = popt_ref[1]-popt[1]

                if verbose:
                    print(f'spectrum {i} shift calculated!')

            if output == 'int':
                step = x_ref[1] - x_ref[0]
                for j, value in enumerate(self.shifts):
                    self.shifts[j] = int(round(value/step))

        self.shifts -= self.shifts[ref]
        self.shift_mode = mode
        self.shift_ranges = ranges

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
        for i in range(self.get_spectra_count()):
            self.data[i].apply_shift(value=self.shifts[i], mode=mode)

    def calculate_sum(self):
        """Sum all spectra."""

        self.check_same_x()

        temp = copy.deepcopy(self.data[0])
        for i in range(1, self.get_spectra_count()):
            temp.data[:, 1] += self.data[i].data[:, 1]
        self.sum = spectrum(data=temp.data)
        return self.sum

    def check_step_x(self, max_error=0.001):
        """Check the step between data point in the x-coordinates

            args:
                max_error (number, optional): percentage value of the max error.

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
            d = np.diff(s.data[:, 0])
            if (max(d) - min(d))*100/np.mean(np.diff(s.data[:, 0])) > max_error:
                raise ValueError(f"Step in the x-coordinate of spectrum {idx} seems not to be uniform.")

        # check step between spectra
        for idx, s in enumerate(self.data):
            try:
                avg_step = (s.data[1, 0] - s.data[0, 0])
                avg_step_2 = (self.data[idx+1].data[1, 0] - self.data[idx+1].data[0, 0])
                if  (avg_step - avg_step_2)*100/avg_step > max_error:
                    raise ValueError(f"Spectrum {idx} ({avg_step}) and {idx+1} ({avg_step_2}) seems to have different step size.")
            except IndexError:
                pass

    def check_same_x(self, max_error=0.001):
        """Compare spectra to see if they have same x-coordinates.

            args:
                max_error (number, optional): percentage value of the max error.

            raises:
                ValueError: If any x-coodinate of any two spectrum is different.
        """
        self.check_step_x(max_error=max_error)

        # check x between spectra
        for idx, s in enumerate(self.data):
            try:
                step = s.data[1, 0] - s.data[0, 0]
                if max(abs(s.data[:, 0] - self.data[idx+1].data[:, 0]))*100/step > max_error:
                    raise ValueError(f"Spectrum {idx} and {idx+1} seems to be different.")
            except IndexError:
                pass

    def calib(self, dispersion=1, zero=0):
        """Calibrate data (from length to energy).

        args:
            dispersion (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].
            zero (number or string, optional): if number, the spectrum is shifted
                so the value given in `zero` is now set to 0 (zero). If `zero=elastic_peak`,
                the position of the elastic peak is guess by `spectrum.guess_elastic_peak()`
                and the center of the elastic peak is set to 0 (zero).

        returns:
            None
        """
        for s in self.data:
            s.calib(dispersion=dispersion, zero=zero)

    #
    # def interp(self, x=None, start=None, stop=None, num=1000, step=None):
    #     """Interpolate data.
    #
    #     args:
    #         x (list or array, optional): The x-coordinates at which to
    #             evaluate the interpolated values. This overwrites all other arguments.
    #         start (number, optional): The starting value of the sequence. If None,
    #             the minium x value will be used.
    #         stop (number, optional): The end value of the sequence. If None,
    #             the maximum x value will be used.
    #         num (int, optional): Number of samples to generate.
    #         step (number, optional): Spacing between values. This overwrites ``num``.
    #
    #     returns:
    #         None
    #     """
    #     if x is None:
    #         if start is None:
    #             start = max([min(s.data[:, 0]) for s in self.data])
    #         if stop is None:
    #             stop = min([max(s.data[:, 0]) for s in self.data])
    #
    #     for s in self.data:
    #         s.interp(x=x, start=start, stop=stop, num=num, step=step)
    #
    # def crop(self, start=None, stop=None):
    #     """Crop spectra ends.
    #
    #     args:
    #         start (number, optional): The starting value. If None,
    #             the minium x value will be used.
    #         stop (number, optional): The end value. If None,
    #             the maximum x value will be used.
    #
    #     returns:
    #         None
    #     """
    #     if start is None:
    #         start = max([min(s.data[:, 0]) for s in self.spectrum])
    #     if stop is None:
    #         stop = min([max(s.data[:, 0]) for s in self.spectrum])
    #
    #     for spectrum in self.spectrum:
    #         step = spectrum.data[1, 0] - spectrum.data[0, 0]
    #         data_range = (start + step/2, stop - step/2)
    #         a, b = am.extract(spectrum.data[:, 0], spectrum.data[:, 1], ranges=(data_range,) )
    #
    #         temp = np.zeros((len(a), 2))
    #         temp[:, 0] = a
    #         temp[:, 1] = b
    #         spectrum.data = copy.deepcopy(temp)

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

        # if 'marker' not in kwargs:
        #     kwargs['marker'] = 'o'
        if 'ms' not in kwargs:
            kwargs['ms'] = 5

        if idx == 'all':
            for i in range(len(self.data)):
                self.data[i].plot(ax=ax, normalized=normalized, vertical_increment=-vertical_increment*i, shift=shift, factor=factor, label=i, **kwargs)
        else:
            try:
                if len(idx) == 1:
                    self.data[idx[0]].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, label=idx[0], **kwargs)
                else:
                    for i in range(len(idx)):
                        self.data[i].plot(ax=ax, normalized=normalized, vertical_increment=-vertical_increment*i, shift=shift, factor=factor, label=i, **kwargs)
            except TypeError:
                self.data[idx].plot(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, label=idx, **kwargs)

        plt.legend()


        if show_ranges:
            if self.shifts is None:
                warnings.warn('Shift range not defined. Use calculate_shifts().', stacklevel=2)
            else:
                for r in self.shift_ranges:
                    plt.axvline(r[0], color='green', linewidth=1.2, zorder=10)
                    plt.axvline(r[1], color='red', linewidth=1.2, zorder=10)

        return ax
