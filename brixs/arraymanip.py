#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Everyday use functions for array manipulation."""

import numpy as np
import copy
from scipy.optimize import curve_fit
from .model_functions import fwhmVoigt


def index(array, value):
    """Returns the index of the element in array which is closest to value.

    Args:
        array (list or array): 1d array.
        value (float or int): value.

    Returns:
        index
    """
    return np.argmin(np.abs(np.array(array)-value))


def extract(x, y, ranges):
    """Returns x and y elements that fall whithin x intervals.

    Args:
        x (list or array): 1d array
        y (list or array): 1d array
        ranges (list): a pair of x values or a list of pairs. Each pair represents
            the start and stop of a data range.

    Notes:
        y is reduced based on x array.

    Examples:
        Use::

            ranges=[[2, 7], [9, 14]]

        to extract data from 2 to 7 and 9 to 14.

    Warning:
        If data ranges intersept with each other, the returned data with have repeated elements.

    Returns:
        Reduced x and y arrays
    """
    x_clean = []
    y_clean = []

    for xinit, xfinal in ranges:
        choose_range = np.logical_and(x>xinit, x<xfinal)
        x_clean.append(x[choose_range])
        y_clean.append(y[choose_range])

    return np.hstack(x_clean), np.hstack(y_clean)


def peak_fit(x, y, guess_c, guess_A, guess_w, guess_offset=0, fixed_m=False, start=None, stop=None, asymmetry=True):
    """Fit a peak with a pseudo-voigt curve.


    Args:
        x (list or array): 1d array
        y (list or array): 1d array
        guess_c (float or int): guess Center
        guess_A (float or int): guess Amplitude
        guess_w (float or int): guess FWHM
        guess_offset (float or int, optional): guess Offset [0]
        fixed_m (False or number): if false, ``m`` will be a fitting parameter. If
            ``fixed_m=number``, ``number`` will be used for ``m``.
        start (float or int): start x value to fit the peak. If ``None``, full
            data is used.
        stop (float or int): final x value to fit the peak.  If ``None``, full
            data is used.
        asymmetry: Bool value. If ``asymmetry=True``, peak asymmetry is taken into account by fiting first
        half of the peak with a different FHWM and m than the second half (m is the
        factor from 1 to 0 of the lorentzian amount).

    Returns:
        1) 2 column (x,y) array with the fitted peak.
        2) 2 column (x,y) array with "Smoothed" fitted peak. This is just the
            fitted peak array with a linear interpolation with 100 times more data points.
        3) An array with the optimized parameters for Amplitude, Center, FWHM and offset.
    """
    if start is None: start=x[0]
    if stop is None: stop=x[-1]

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
    x2fit, y2fit = extract(x, y, [[start, stop],])
    popt, pcov = curve_fit(function2fit, x2fit, y2fit, p0,  # sigma = sigma,
                           bounds=bounds)

    # smooth data
    arr100 = np.zeros([100*len(x2fit), 2])
    arr100[:, 0] = np.linspace(x2fit[0], x2fit[-1], 100*len(x2fit))
    arr100[:, 1] = function2fit(arr100[:, 0],  *popt)

    if asymmetry:
        popt_2 = (popt[0], popt[1], popt[2]/2+popt[4]/2, popt[-1])
    else:
        popt_2 = (popt[0], popt[1], popt[2], popt[-1])

    return function2fit(x, *popt), arr100, popt_2


def shift(x, y, shift, mode='hard'):
    """Shift (x, y) data.

    Args:
        x (list or array): 1D array.
        y (list or array): 1D array.
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
        Shifted x and y.
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

def movingaverage(array, window_size, remove_boundary_effects=True):
    """Returns the moving average of an array.

    The returned array has the same length of the original array.

    Example:
        >>> print(manip.movingaverage([0,1,2,3,4,5,6,7,8,9], 1))
        [0. 1. 2. 3. 4. 5. 6. 7. 8. 9.]
        >>> print(manip.movingaverage([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 2))
        [0. , 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5]
        >>> print(manip.movingaverage([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 3))
        [0.33333333 1. 2. 3. 4. 5. 6.  7. 8. 5.66666667]
        >>> print(manip.movingaverage([0, 1, 2, 3, 4, 5, 6, 7, 8, 9], 4))
        [0.25 0.75 1.5  2.5  3.5  4.5  5.5  6.5  7.5  6.  ]

    Warning:
        Note by the example that the resulting array contains boundary effects.

    Args:
        array (list or np.array): array.
        window_size (int): number of points to average.

    Returns:
        array.
    """
    if window_size < 1:
        raise ValueError('window_size must be a positive integer (> 1).')

    array = np.array(array)
    window = np.ones(int(window_size))/float(window_size)
    final = np.convolve(array, window, 'same')

    if remove_boundary_effects:
        final = final[int(window_size/2):]
    return final


def derivative(x, y, order=1, window_size=1):
    """ if window_size > 1, boundary effects might be """

    if order<0:
        raise ValueError('order must be a positive integer.')

    x = np.array(x)
    y = np.array(y)

    x_diff = np.diff(x)
    y_diff = np.diff(y)/x_diff
    for i in range(order-1):
        y_diff = np.diff(y_diff)/x_diff[:len(x_diff)-(i+1)]


    i = int(order/2)
    f = - int(order/2)
    if (order % 2) != 0:
        f -= 1

    if window_size>1:
        x = movingaverage(x, window_size=window_size)
        y_diff = movingaverage(y_diff, window_size=window_size)

    return x[i: f], y_diff




# def increasing_monotonicity(dataX, dataY):
#     """Returns an array sorted and monotonic.
#
#     The sorting is based on dataX and the monotonicity is done by averaging
#     dataY for same dataX values.
#
#     If you need decreasing monotonicity just run this function and invert the
#     returned arrays.
#
#     :param dataX: list of numbers
#     :param dataY: list of numbers
#     :return: two numpy arrays (x_monotonic, y_monotonic)
#
#     Example:    dataX= [1, 2, 4, 2, 1]
#                 dataY = [5, 6, 9, 8, 9]
#
#                 x_return = [1, 2, 4]
#                 y_return = [7, 7, 9]
#     """
#     # sort increasingly
#     data2sort = np.array(np.transpose([dataX, dataY]))
#     data_sorted = data2sort[data2sort[:, 0].argsort()]
#
#     done = False
#     data_sorted_clean = np.copy(data_sorted)
#     i = 0
#
#     while not done:
#         val = data_sorted_clean[i, 0]
#
# #        print(i, val)
#         if i == len(data_sorted_clean)-1:
# #            print('aqui')
#             done = True
#             return data_sorted_clean[:, 0], data_sorted_clean[:, 1]
#
#         #Find how many duplicates there is
#         number_of_duplicates = 0
#         k = np.copy(i)
#         while val == data_sorted_clean[k+1, 0]:
#             k = k+1
#             number_of_duplicates = number_of_duplicates+1
#             if k==(len(data_sorted_clean)-1):
#                 done = True
#                 break
# #        print(i, val, k, number_of_duplicates)
#
#         #Mean
#         if number_of_duplicates>=1:
#             data_sorted_clean[i, 1] = np.mean(data_sorted_clean[i:(i+number_of_duplicates+1), 1], dtype=np.float64)
# #            print(data_sorted_clean)
#
#             for j in range(number_of_duplicates):
# #                print('e')
#                 data_sorted_clean = np.delete(data_sorted_clean, i+1, axis=0)
# #                print(data_sorted_clean)
#         i = i + 1
#
#         if done:
#             return data_sorted_clean[:, 0], data_sorted_clean[:, 1]
