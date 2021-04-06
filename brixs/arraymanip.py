#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Everyday use functions for array manipulation."""

import numpy as np
import copy
from scipy.optimize import curve_fit
from .model_functions import fwhmVoigt


def index(x, value):
    """Returns the index of the element in array which is closest to value.

    Args:
        x (list or array): 1D array.
        value (float or int): value.

    Returns:
        index (int)
    """
    return np.argmin(np.abs(np.array(x)-value))


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
    return s


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


def movingaverage(x, window_size):
    """Returns the moving average of an array.

    Args:
        x (list or array): 1D array.
        window_size (int): number of points to average.

    Returns:
        array.

    Example:
        >>> x = [0,1,2,3,4,5,6,7,8,9]
        >>> print(am.movingaverage(x, 1))
        [0. 1. 2. 3. 4. 5. 6. 7. 8. 9.]
        >>> print(am.movingaverage(x, 2))
        [0.5 1.5 2.5 3.5 4.5 5.5 6.5 7.5 8.5]
        >>> print(am.movingaverage(x, 3))
        [1. 2. 3. 4. 5. 6. 7. 8.]
        >>> print(am.movingaverage(x, 4))
        [1.5 2.5 3.5 4.5 5.5 6.5 7.5]

    """
    if window_size < 1:
        raise ValueError('window_size must be a positive integer (> 1).')

    x = np.array(x)
    window = np.ones(int(window_size))/float(window_size)

    return np.convolve(x, window, 'valid')


def derivative(x, y, order=1):
    """Returns the derivative of y-coordinates as a function of x-coodinates.

        Args:
            x (list or array): 1D array x-coordinates.
            y (list or array): 1D array y-coordinates.
            order (number, optional): derivative order.

        Returns:
            x and y arrays.
    """

    if order<0:
        raise ValueError('order must be a positive integer.')

    x = np.array(x)
    y = np.array(y)

    x_diff = np.diff(x)
    y_diff = np.diff(y)/x_diff
    for i in range(order-1):
        y_diff = np.diff(y_diff)/x_diff[:len(x_diff)-(i+1)]

    for i in range(order):
        x = movingaverage(x, window_size=2)

    return x, y_diff


def shift(x, y, shift, mode='hard'):
    """Shift (x, y) data.

    Args:
        x (list or array): 1D array.
        y (list or array): 1D array.
        shift (float or int): shift value.
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
        y = np.interp(x, x + shift, y)

    elif mode == 'x' or mode == 'hard':
        x = np.array(x) + shift

    elif mode == 'roll' or mode == 'rotate':
        if shift.is_integer():
            y = np.roll(y, int(shift))
        else:
            raise ValueError("shift must be an interger for mode='roll'.")
        if shift > 0:
            y[:int(shift)] = 0
        elif shift < 0:
            y[int(shift):] = 0
    else:
        raise ValueError('mode not recognized')

    return x, y


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
