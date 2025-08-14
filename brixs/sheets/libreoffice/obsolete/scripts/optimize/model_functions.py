#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Everyday use mathematical functions and distributions."""

import numpy as np
from scipy.special import erf


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


def areaGauss(x, A, c, sigma):
    r"""Gaussian distribution.

    .. math:: y(x) = \frac{\text{Area}}{\sqrt{2\pi} w} e^{-\frac{(x-c)^2}{2 w^2}}

    :param x: x array
    :param A: Area
    :param c: Center
    :param sigma: standard deviation
    :return: :math:`y(x)`
    """
    return Gauss(x, A/(np.sqrt(2*np.pi)*sigma), c, sigma)
    # return A/(np.sqrt(2*np.pi)*abs(w))  *np.exp(-(x-c)**2/(2*w**2))


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


def fwhmAreaGauss(x, A, c, w):
    r"""Gaussian distribution.

    .. math:: y(x) = \frac{2 \sqrt{\ln(2)} A}{w \sqrt{\pi}}  e^{-\frac{4 \ln(2) (x-c)^2}{w^2}}

    :param x: x array
    :param A: Area
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """
    w = w/(2*np.sqrt(2*np.log(2)))
    return Gauss(x, A/(np.sqrt(2*np.pi)*w), c, w)
    # return (A/(w*np.sqrt(np.pi/4*np.log(2))))*np.exp((-4*np.log(2)*((x-c)**2))/(w**2))


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


def fwhmAreaLorentz(x, A, c, w):
    r"""Cauchy–Lorentz distribution.

    .. math:: y(x) = A \frac{1}{\pi w} \frac{w^2}{w^2 + (x-c)^2}

    :param x: x array
    :param A: Area
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """
    return A * Lorentz(x, gamma=w, c=c)
    # return ((2*A)/(np.pi))*((w)/(w**2 + 4*(x-c)**2))


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


def fwhmAreaVoigt(x, A, c, w, m):
    r"""Pseudo-voigt curve.

    .. math:: y(x) = A \left[ m \frac{1}{\pi w} \frac{w^2}{w^2 + 4 (x-c)^2}   + (1-m) \frac{2 \sqrt{\ln(2)}}{w \sqrt{\pi}}  e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

    where,

    .. math:: \text{amp }= \frac{A}{w} \left[ m\frac{1}{\pi} + (1-m)\frac{2\sqrt{\ln(2)}}{\sqrt{\pi}} \right]

    :param x: x array
    :param A: is the Area
    :param c: Center
    :param w: FWHM
    :param m: Factor from 1 to 0 of the lorentzian amount
    :return: :math:`y(x)`
    """
    lorentz = fwhmAreaLorentz(x, 1, c, w)
    gauss = fwhmAreaGauss(x, 1, c, w)

    return A*(m*lorentz + (1-m)*gauss)


def fwhmArctan(x, amp, c, w):
    r"""Arctangent function.

    .. math:: y(x) =   \frac{A}{\pi} \left[ \arctan(\frac{1}{w}(x-c)) + \frac{\pi}{2} \right]

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM (it will take fwhm units to go from amp/4 to (3amp)/4)
    :return: :math:`y(x)`
    """

    return amp * (np.arctan((w**-1)*(x - c)) + (np.pi/2))/np.pi


def square(x, amp, c, w):
    r"""Square step function.

    .. math::

        y(x) =   \begin{cases}
                    0, & \text{ for    $x < c-w/2$}\\
                    \text{amp}, & \text{ for  }  c-w/2 < x < c+w/2\\
                    0, & \text{ for  }  x > c+w/2\\
                    \end{cases}

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """
    return - np.heaviside(x-c-w/2, amp)*amp + np.heaviside(x-c+w/2, amp)*amp


def fwhmErr(x, amp, c, w):
    r"""Error function. Integal of gaussian function calculated by ``scipy.special.erf()``.

    .. math:: y(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM (it will take roughly fwhm units to go from amp/4 to (3amp)/4)
    :return: :math:`y(x)`
    """
    return amp/2 * (erf((w**-1)/2*(x - c))+1)
