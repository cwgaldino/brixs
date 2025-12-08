#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""basic peak fitting functions"""

# %% ------------------------- Standard Imports --------------------------- %% #
from scipy.optimize import curve_fit
import numpy as np
# %%

# %% -------------------------- brixs Imports ----------------------------- %% #
import brixs as br
# %%

# %% ------------------------ supporting functions ------------------------ %% #
def _index(x, value, closest=True):
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

# %% ====================== basic data fitting ============================ %% #
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

    .. math:: y(x) = \text{amp } \frac{w^2}{w^2 + 4 * (x-c)^2}

    where,

    .. math:: \text{Area }= \text{amp } \pi w \frac{1}{2}

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """
    return amp*((w**2)/(w**2 + 4* (x-c)**2))

def _voigt_fwhm(x, amp, c, w, m):
    r"""Pseudo-voigt curve.

    .. math:: y(x) = \text{amp } \left[ m  \frac{w^2}{w^2 + 4*(x-c)^2}   + (1-m) e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

    where,

    .. math:: \text{Area } = \frac{\text{amp } w}{2}     \left[ m \pi +   (1-m)      \frac{\sqrt{\pi}}{\sqrt{\ln(2)} }  \right]  

    
    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM
    :param m: Factor from 1 to 0 of the lorentzian amount
    :return: :math:`y(x)`
    """
    lorentz = _lorentzian_fwhm(x, 1, c, w)
    gauss = _gaussian_fwhm(x, 1, c, w)

    return amp*(m*lorentz + (1-m)*gauss)

def fit_peak(x, y, guess_c=None, guess_A=None, guess_w=None, guess_offset=0, fixed_m=False, asymmetry=False):
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
        dictionary {fit, popt, sigma, model}

        fit: 2 column (x, y) array with "Smoothed" fitted peak (array length 100 bigger than input x, y).
        popt: An array with the optimized parameters.
            if asymmetry=True, fixed_m=False: amp, c, fwhm1, m1, fwhm2, m2, offset
            if asymmetry=True, fixed_m=True: amp, c, fwhm1, fwhm2, offset
            if asymmetry=False, fixed_m=False: amp, c, fwhm, m, offset
            if asymmetry=False, fixed_m=True: amp, c, fwhm, offset
        sigma: One standard deviation errors on the parameters
        model: Peak function -> model(x)
    """
        # .. image:: _static/peak_fit.png
        #     :width: 600
        #     :align: center

    start = x[0]
    stop = x[-1]

    if guess_A is None:
        guess_A = max(y)

    if guess_c is None:
        guess_c = x[_index(y, guess_A)]

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

    return {'fit':arr100, 'popt':popt, 'sigma':err, 'model':lambda x: function2fit(x, *popt)}

# %% ====================== Spectrum peak fitting ========================= %% #
def _fit_peak(self, guess_c=None, guess_A=None, guess_w=None, guess_offset=0, fixed_m=False, asymmetry=False, moving_average_window=1, limits=None):     
    r"""Simple peak fit function. Data is fitted with a pseudo-voigt curve.

    .. math:: y(x) = A \left[ m \frac{w^2}{w^2 + (x-c)^2}   + (1-m) e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

    Usage:
        >>> import brixs as br
        >>> import brixs.addons.fitting
        >>>
        >>> x = np.linspace(0, 100, 1000)
        >>> amp = 100
        >>> w = 10
        >>> c = 25
        >>> y = br.gaussian_fwhm(x, amp, c, w)
        >>>
        >>> s = br.Spectrum(x, y)
        >>>
        >>> smooth, popt, err, f = s.fit_peak()
        >>>
        >>> print(f'A = {popt[0]} +/- {err[0]}')

    Args:
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
        moving_average_window (int, optional): window size for smoothing the
            data for finding the peak. Default is 1.
        limits (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. Use None to indicate
            the minimum or maximum x value of the data.

    Returns:
        Dictionary {fit, popt, sigma, model}

        1) fit: Spectrum with fitted peak (array length 100 bigger than original spectrum).
        2) popt: An array with the optimized parameters.
            if asymmetry=True, fixed_m=False: amp, c, fwhm1, m1, fwhm2, m2, offset
            if asymmetry=True, fixed_m=True: amp, c, fwhm1, fwhm2, offset
            if asymmetry=False, fixed_m=False: amp, c, fwhm, m, offset
            if asymmetry=False, fixed_m=True: amp, c, fwhm, offset
        3) sigma: One standard deviation errors on the parameters
        4) model: Peak function
    """
    if limits is None:
        x0 = self.x
        y0 = self.y
    else:
        limits = self._check_limits(limits)
        s0 = self._copy(limits)
        x0 = s0.x
        y0 = s0.y

    # smoothing
    assert moving_average_window > 0, f'moving_average_window must be positive different than 0, not {moving_average_window}'
    if moving_average_window == 1:
        x = x0
        y = y0
    else:
        x = br.moving_average(x0, moving_average_window)
        y = br.moving_average(y0, moving_average_window)

    # guess amp and c
    amp = max(y)
    c = x[np.argmax(y)]

    # guess fwhm
    try:
        w1 = x[np.argmax(y)] - x[:np.argmax(y)][::-1][_index(y[:np.argmax(y)][::-1], max(y)/2)]
    except ValueError:
        w1 = x[np.argmax(y):][_index(y[np.argmax(y):], max(y)/2)] - x[np.argmax(y)]
    try:
        w2 = x[np.argmax(y):][_index(y[np.argmax(y):], max(y)/2)] - x[np.argmax(y)]
    except ValueError:
        w2 = x[np.argmax(y)] - x[:np.argmax(y)][::-1][_index(y[:np.argmax(y)][::-1], max(y)/2)]
    w = w1 + w2

    # x, y = derivative(moving_average(self.x, 10), moving_average(self.y, 10))
    # w = np.abs(x[np.argmin(y)] - x[np.argmax(y)])
    if w == 0:
        w = 0.1*(max(self.x)-min(self.x))

    result = fit_peak(x, y, guess_c=c, guess_A=amp, guess_w=w, guess_offset=0, fixed_m=fixed_m, asymmetry=asymmetry)
    s = br.Spectrum(x=result['fit'][:, 0], y=result['fit'][:, 1])
    s.copy_attrs_from(self)
    return {'fit':s, 'popt':result['popt'], 'sigma':result['sigma'], 'model':result['model']}
br.Spectrum.fit_peak = _fit_peak

# %% ======================= Spectra peak fitting ========================= %% #
def _fit_peak_spectra(self, guess_c=None, guess_A=None, guess_w=None, guess_offset=0, fixed_m=False, asymmetry=False, moving_average_window=1, limits=None):     
    r"""Simple peak fit function. Data is fitted with a pseudo-voigt curve.

    .. math:: y(x) = A \left[ m \frac{w^2}{w^2 + (x-c)^2}   + (1-m) e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

    Usage:
        >>> import brixs as br
        >>> import brixs.addons.fitting
        >>>
        >>> x = np.linspace(0, 100, 1000)
        >>> amp = 1
        >>> w = 10
        >>> c = 25
        >>> y = br.fwhmGauss(x, amp, c, w)
        >>>
        >>> s = br.Spectrum(x, y)
        >>>
        >>> smooth, popt, err, f = s.fit_peak()
        >>>
        >>> print(f'A = {popt[0]} +/- {err[0]}')

    Args:
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
        moving_average_window (int, optional): window size for smoothing the
            data for finding the peak. Default is 1.
        limits (list): a pair of values or a list of pairs. Each pair represents
            the start and stop of a data range from x. Use None to indicate
            the minimum or maximum x value of the data.

    Returns:
        Dictionary {fit, popt, sigma, model}

        1) fit: Spectra with every fitted peak (arrays have length 100 bigger than input).
        2) popt: List of arrays with the optimized parameters.
            if asymmetry=True, fixed_m=False: amp, c, fwhm1, m1, fwhm2, m2, offset
            if asymmetry=True, fixed_m=True: amp, c, fwhm1, fwhm2, offset
            if asymmetry=False, fixed_m=False: amp, c, fwhm, m, offset
            if asymmetry=False, fixed_m=True: amp, c, fwhm, offset
        3) sigma: list of one-standard deviation errors on the parameters
        4) model: Peak functions
    """
    ss   = br.Spectra()
    ss.copy_attrs_from(self)
    popt  = []
    sigma = []
    model = []
    for s in self:
        result = s.fit_peak(guess_c=guess_c, guess_A=guess_A, guess_w=guess_w, guess_offset=guess_offset, fixed_m=fixed_m, moving_average_window=moving_average_window, asymmetry=asymmetry, limits=limits)
        ss.append(result['fit'])
        popt.append(result['popt'])
        sigma.append(result['sigma'])
        model.append(result['model'])
    return {'fit':ss, 'popt':popt, 'sigma':sigma, 'model':model}
br.Spectra.fit_peak = _fit_peak_spectra

# %% ====================== basic model functions ========================= %% #
br.gaussian          = _gaussian
br.gaussian_fwhm     = _gaussian_fwhm
br.lorentzian        = _lorentzian
br.lorentzian_fwhm   = _lorentzian_fwhm
br.voigt_fwhm        = _voigt_fwhm


