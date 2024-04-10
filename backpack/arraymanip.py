#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> array and iterables"""

# %% ------------------------- Standard Imports --------------------------- %% #
from scipy.optimize import curve_fit
from collections.abc import Iterable
import numpy as np

# %% ======================== calculation/modification ==================== %% #
def sort(ref, *args):
    """Returns sorted arrays based on a reference array.

    Args:
        ref (list or array): 1D array.
        *args: array to sort.

    Returns:
        sorted arrays."""

    s = []
    for x in args:
        s.append( [x1 for (y,x1) in sorted(zip(ref,x), key=lambda pair: pair[0])])
    if len(args) == 1:
        s = s[0]
    return s

def moving_average(x, n):
    """Returns the moving average of an array.

    Args:
        x (list or array): 1D array.
        n (int): number of points to average.

    Returns:
        array of length given by (len(x)-n+1).

    Example:
        >>> x = [0,1,2,3,4,5,6,7,8,9]
        >>> print(am.moving_average(x, 1))
        [0. 1. 2. 3. 4. 5. 6. 7. 8. 9.]
        >>> print(am.moving_average(x, 2))
        [0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]
        >>> print(am.moving_average(x, 3))
        [1. 2. 3. 4. 5. 6. 7. 8.]
        >>> print(am.moving_average(x, 4))
        [1.5 2.5 3.5 4.5 5.5 6.5 7.5]

    """
    if n < 1:
        raise ValueError('n must be a positive integer (> 1).')
    if isinstance(n, int) == False:
        if n.is_integer() == False:
            raise ValueError('n must be a positive integer (> 1).')

    x = np.array(x)
    window = np.ones(int(n))/float(n)

    return np.convolve(x, window, 'valid')

def derivative(x, y, order=1):
    """Returns the derivative of y-coordinates as a function of x-coordinates.

        Args:
            x (list or array): 1D array x-coordinates.
            y (list or array): 1D array y-coordinates.
            order (number, optional): derivative order.

        Returns:
            x and y arrays.
    """

    # check order   
    if order<0:
        raise ValueError('order must be a positive integer.')
    
    # check monotonicity
    if check_monotonicity(x) < 0:
        raise ValueError('x array must be monotonically increasing. Please use brixs.fix_monotonocity().')


    x = np.array(x)
    y = np.array(y)

    x_diff = np.diff(x)
    y_diff = np.diff(y)/x_diff
    for i in range(order-1):
        y_diff = np.diff(y_diff)/x_diff[:len(x_diff)-(i+1)]

    for i in range(order):
        x = moving_average(x, n=2)

    return x, y_diff

def shifted(x, y, value, mode='hard'):
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
            #. ``mode='roll'` or ``r``,
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
    x = np.array(x)
    y = np.array(y)

    if value == 0:
        return x, y

    if mode == 'y' or mode == 'interp' or mode=='soft':
        if check_monotonicity(x) != 1:
            raise ValueError('x array must be increasingly monotonic.')
        y = np.interp(x, x + value, y)

    elif mode == 'x' or mode == 'hard':
        x = np.array(x) + value

    elif mode == 'roll' or mode == 'r':
        if is_integer(value):
            y = np.roll(y, int(value))
        else:
            raise ValueError("value must be an interger for mode='roll'.")
        # try:
        #     if value.is_integer():
        #         y = np.roll(y, int(value))
        #     else:
        #         raise ValueError("value must be an interger for mode='roll'.")
        # except AttributeError:
        #     y = np.roll(y, int(value))
        # if value > 0:
        #     y[:int(value)] = y[int(value)]
        # elif value < 0:
        #     # print(y[int(value):])
        #     # print( y[int(value-1)])
        #     y[int(value):] = y[int(value-1)]
    else:
        raise ValueError("mode not recognized (valid: 'y', 'x', 'roll').")

    return x, y

def transpose(array):
    """Returns transposed lists/arrays.

    Differently from numpy.transpose(), this function also transposes 1D arrays/lists.
    """
    try:
        row_count, col_count = np.shape(array)
        return np.transpose(array)
    except ValueError:
        return [[x] for x in array]

def mask(x, mask):
    """Returns a reduced array based on mask.

    Usage:
        mask('ABCDEF', [1, 0, 1, 0, 1, 1]) --> A C E F

    Args:
        x (Iterable): array
        mask (list): list with bools (True, False, ...)
        
    Returns:
        reduced array
    """
    return [d for d, s in zip(x, mask) if s]

# %% ============================= array check ============================ %% #
def index(x, value, closest=True):
    """Returns the first index of the element in array.

    Args:
        x (list or array): 1D array.
        value (float or int): value.
        closest (book, optional): if True, returns the index of the element in 
            array which is closest to value.

    Returns:
        index (int)
    """
    if closest:
        return int(np.argmin(np.abs(np.array(x)-value)))
    else:
        return np.where(x == value)[0]

def all_equal(array):
    """Returns True if all elements of an array are equal."""
    iterator = iter(array)
    try:
        first = next(iterator)
    except StopIteration:
        return True
    return all(first == x for x in iterator)

def has_duplicates(array):
    """Returns True if given array contains any duplicates."""
    if len(array) == len(set(array)):
        return False
    else:
        return True
    
def remove_duplicates(array):
    """Returns the sorted unique elements of an array. Wrapper for np.unique()"""
    temp = np.unique(array)
    if isinstance(array, list):
        temp = list(temp)
    return temp

# %% ============================ monotonicity ============================ %% #
def check_monotonicity(array):
    """return 1 (-1) if increas. (decre.) monotonic or 0 if not monotonic."""
    if np.all(np.diff(array) > 0) == True:
        return 1
    elif np.all(np.diff(array) < 0) == True:
        return -1
    else:
        return 0

def fix_monotonicity(x, y, mode='increasing'):
    """return x, y where the x array is monotonically increasing or decreasing."""
    if mode != 'increasing' and mode != 'decreasing':
        raise ValueError('mode should be "decreasing" or "increasing".')

    if mode == 'increasing':
        if check_monotonicity(x) == 1:
            return x, y
        else:
            unqa, ID, counts = np.unique(x, return_inverse=True, return_counts=True)
            return unqa, np.bincount(ID, y)/counts
    if mode == 'decreasing':
        if check_monotonicity(x) == -1:
            return x, y
        else:
            unqa, ID, counts = np.unique(x, return_inverse=True, return_counts=True)
            return np.fliplr(unqa), np.fliplr(np.bincount(ID, y)/counts)



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

# %% ============================== extract =============================== %% #
def choose(x, ranges):
    """Return a mask of x values inside range pairs.

    Args:
        x (list or array): 1d array.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. [[xi, xf], [xi_2, xf_2], ...]

    Returns:
        1d list.
    """
    assert isinstance(x, Iterable), 'input must be a iterable'
    x = np.array(x)

    try:  # ((start, end), )
        choose_range = [None]*len(ranges)
        for i, (x_init, x_final) in enumerate(ranges):
            choose_range[i] = np.logical_and(x>=x_init, x<=x_final)
        choose_range = [True if x == 1 else False for x in np.sum(choose_range, axis=0)]
    except TypeError:  # (start, end)
        x_init, x_final = ranges
        choose_range = np.logical_and(x>=x_init, x<=x_final)
    return choose_range

def extract(x, y, ranges, invert=False):
    """Returns specific data ranges from x and y.

    Args:
        x (list or array): 1D reference vector.
        y (list or array): 1D y-coordinates or list of several data sets.
        ranges (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x.
        invert (bool, optional): if inverted is True, data outside of the data 
            will be returned. Default is False.

    Returns:
        x and y arrays. If `y` is 1d, the returned `y` is 1d. If `y` is
        a multicolumn array then the returned `y` is also multicolumn


    Examples:

        if `y` is 1d, the returned `y` is 1d:

        >>> x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        >>> y = np.array(x)**2
        >>> ranges = ((0, 3), (7.5, 9))
        >>> x_sliced, y_sliced = am.extract(x, y, ranges)
        >>> print(x_sliced)
        [0 1 2 3 8 9]
        >>> print(y_sliced)
        [0 1 4 9 64 81]

        if `y` is multicolumn, the returned `y` is also multicolumn:

        >>> x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        >>> y = np.zeros((10, 2))
        >>> y[:, 0] = x**2
        >>> y[:, 0] = x**3
        >>> ranges = ((0, 3), (7.5, 9))
        >>> x_sliced, y_sliced = am.extract(x, y, ranges)
        >>> print(x_sliced)
        [0. 1. 2. 3. 8. 9.]
        >>> print(y_sliced)
        [[  0.   0.]
         [  1.   0.]
         [  8.   0.]
         [ 27.   0.]
         [512.   0.]
         [729.   0.]]


    """
    x = np.array(x)
    y = np.array(y)

    # if data is all inside ranges, then nothing is done
    if len(ranges) == 1:
        if ranges[0][0] <= min(x) and ranges[0][1] >= max(x):
            if invert:
                return np.array([]), np.array([])
            else:
                return x, y

    choose_range = choose(x, ranges)
    # print(choose_range[0:10])
    if invert:
        choose_range = np.invert(choose_range)
    # temp = np.compress(choose_range, np.c_[y.transpose(), x], axis=0)
    # print(choose_range[0:10])
    temp = np.compress(choose_range, np.c_[y, x], axis=0)
    # print(temp)
    if temp.any():
        if len(temp[0]) > 2:
            return temp[:, -1], temp[:, :-1]#.transpose()
        else:
            # print('here')
            return temp[:, -1], temp[:, 0]
    else:
        raise RuntimeError('No data points within the selected range.')

# %% ========================== data fitting ============================== %% #
def _gaussian(x, amp, c, sigma):
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

def _gaussian_fwhm(x, amp, c, w):
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
    return _gaussian(x, amp, c, w/(2*np.sqrt(2*np.log(2))))
    # return A*np.exp((-4*np.log(2)*((x-c)**2))/(w**2))


def _lorentzian(x, gamma, c):
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

def _lorentzian_fwhm(x, amp, c, w):
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
    return amp*(np.pi*w) * _lorentzian(x, gamma=w, c=c)
    # return A*((w**2)/(w**2 + 4* (x-c)**2))

def _voigt_fwhm(x, amp, c, w, m):
    r"""Pseudo-voigt curve.

    .. math:: y(x) = A \left[ m  \frac{w^2}{w^2 + (x-c)^2}   + (1-m) e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM
    :param m: Factor from 1 to 0 of the lorentzian amount
    :return: :math:`y(x)`
    """



    lorentz = _lorentzian_fwhm(x, 1, c, w)
    gauss   = _gaussian_fwhm(x, 1, c, w)

    return amp*(m*lorentz + (1-m)*gauss)

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
        fixed_m (False or number): Factor from 1 to 0 of the lorentzian amount.
            If False, ``m`` will be a fitting parameter. If
            ``fixed_m=<number>``, ``<number>`` will be used for ``m``.
        asymmetry (Bool, optional). If True, peak asymmetry is taken into account by fitting first
            half of the peak with a different ``w`` and ``m`` than the second half. The optimal ``w`` parameter
            returned will be the sum of the ``w`` of the first and second half.

    Returns:
        1) 2 column (x, y) array with "Smoothed" fitted peak (array length 100 bigger than input x, y).
        2) An array with the optimized parameters.
            if asymmetry=True, fixed_m=False: amp, c, fwhm1, m1, fwhm2, m2, offset
            if asymmetry=True, fixed_m=True: amp, c, fwhm1, fwhm2, offset
            if asymmetry=False, fixed_m=False: amp, c, fwhm, m, offset
            if asymmetry=False, fixed_m=True: amp, c, fwhm, offset

        3) One standard deviation errors on the parameters
        4) Peak function

    See Also:
        :py:func:`backpack.model_functions.voigt_fwhm`

    Example:
        >>> import matplotlib.pyplot as plt
        >>> from backpack.model_functions import fwhmGauss
        >>> x = np.linspace(0, 100, 1000)
        >>> amp = 1
        >>> w = 10
        >>> c = 25
        >>> y = fwhmGauss(x, amp, c, w) + np.random.normal(-0.1, 0.1, 1000)
        >>> smooth, popt, err, f = am.peak_fit(x, y)
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
    """
        # .. image:: _static/peak_fit.png
        #     :width: 600
        #     :align: center

    start = x[0]
    stop = x[-1]

    if guess_A is None:
        guess_A = max(y)

    if guess_c is None:
        guess_c = x[index(y, guess_A)]

    if guess_w is None:
        guess_w = 0.1*guess_c

    if fixed_m == False and type(fixed_m)==bool:  # variable m
        if asymmetry:
            p0 = [guess_A, guess_c, guess_w, 0.5, guess_w, 0.5, guess_offset]
            def function2fit(x, A, c, w1, m1, w2, m2, offset):
                f = np.heaviside(c-x, 0)*_voigt_fwhm(x, A, c, w1, m1) + offset +\
                    np.heaviside(x-c, 0)*_voigt_fwhm(x, A, c, w2, m2)
                return f
            bounds=[[-np.inf, start,   0,    0,   0,    0, -np.inf],
                    [ np.inf, stop,  np.inf, 1, np.inf, 1,  np.inf]]
        else:
            p0 = [guess_A, guess_c, guess_w, 0.5, guess_offset]
            def function2fit(x, A, c, w, m, offset):
                return _voigt_fwhm(x, A, c, w, m) + offset
            bounds=[[-np.inf, start,   0,    0, -np.inf],
                    [ np.inf, stop,  np.inf, 1,  np.inf]]

    else:
        if fixed_m > 1:
            fixed_m = 1
        elif fixed_m < 0:
            fixed_m = 0
        if asymmetry:
            p0 = [guess_A, guess_c, guess_w/2, guess_w/2, guess_offset]
            def function2fit(x, A, c, w1, w2, offset):
                f = np.heaviside(c-x, 0)*_voigt_fwhm(x, A, c, w1, fixed_m) + offset +\
                    np.heaviside(x-c, 0)*_voigt_fwhm(x, A, c, w2, fixed_m)
                return f
            bounds=[[-np.inf, start,   0,    0, -np.inf],
                    [ np.inf, stop,  np.inf, np.inf,  np.inf]]
        else:
            p0 = [guess_A, guess_c, guess_w, guess_offset]
            def function2fit(x, A, c, w, offset):
                return _voigt_fwhm(x, A, c, w, fixed_m) + offset
            bounds=[[-np.inf, start,   0,   -np.inf],
                    [ np.inf, stop,  np.inf, np.inf]]

    # Fit data
    popt, pcov = curve_fit(function2fit, x, y, p0,  # sigma = sigma,
                           bounds=bounds)
    err = np.sqrt(np.diag(pcov))  # One standard deviation errors on the parameters

    # smooth data
    arr100 = np.zeros([100*len(x), 2])
    arr100[:, 0] = np.linspace(x[0], x[-1], 100*len(x))
    arr100[:, 1] = function2fit(arr100[:, 0],  *popt)

    # if fixed_m == False and type(fixed_m)==bool:
    #     if asymmetry:
    #         popt_2 = (popt[0], popt[1], popt[2], popt[3], popt[4], popt[5])
    #     else:
    #         popt_2 = (popt[0], popt[1], popt[2], popt[3], popt[4])
    # else:
    #
    #     if asymmetry:
    #         popt_2 = (popt[0], popt[1], popt[2], popt[4], popt[-1])
    #     else:
    #         popt_2 = (popt[0], popt[1], popt[2], popt[-1])

    return arr100, popt, err, lambda x: function2fit(x, *popt)

# %% ========================== Experimental ============================== %% #
def flatten(x):
    """[EXPERIMENTAL] Return a copy of the array collapsed into one dimension."""
    assert isinstance(x, Iterable), 'input must be Iterable (list, tuple, array)'

    # empty array
    if len(x) == 0:
        return x

    if len(np.array(x).shape) == 1:
        return x
    return np.concatenate(x).ravel()