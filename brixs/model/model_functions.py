#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Mathematical functions and distributions."""

# %% ------------------------- Standard Imports --------------------------- %% #
from scipy.special import erf
import numpy as np

# %% ========================== model functions =========================== %% #
def gaussian(x, amp, c, sigma):
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


def gaussian_area(x, A, c, sigma):
    r"""Gaussian distribution.

    .. math:: y(x) = \frac{\text{Area}}{\sqrt{2\pi} w} e^{-\frac{(x-c)^2}{2 w^2}}

    where,

    .. math:: \text{amp } = \frac{\text{Area }}{\sqrt{2 \pi} |\sigma| }

    :param x: x array
    :param A: Area
    :param c: Center
    :param sigma: standard deviation
    :return: :math:`y(x)`
    """
    return gaussian(x, A/(np.sqrt(2*np.pi)*sigma), c, sigma)


def gaussian_fwhm(x, amp, c, w):
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
    return gaussian(x, amp, c, w/(2*np.sqrt(2*np.log(2))))


def gaussian_area_fwhm(x, A, c, w):
    r"""Gaussian distribution.

    .. math:: y(x) = \frac{2 \sqrt{\ln(2)} A}{w \sqrt{\pi}}  e^{-\frac{4 \ln(2) (x-c)^2}{w^2}}

    where,

    .. math:: \text{ amp }  = \frac{ \text{Area } 2 \sqrt{\ln(2)}} {\sqrt{\pi}w }

    :param x: x array
    :param A: Area
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """
    w = w/(2*np.sqrt(2*np.log(2)))
    return gaussian(x, A/(np.sqrt(2*np.pi)*w), c, w)


def lorentzian(x, gamma, c):
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


def lorentzian_fwhm(x, amp, c, w):
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


def lorentzian_area_fwhm(x, A, c, w):
    r"""Cauchy–Lorentz distribution.

    .. math:: y(x) = A \frac{1}{\pi w} \frac{w^2}{w^2 + (x-c)^2}

    where,

    .. math:: \text{amp }  = \frac{2 \text{Area }}{\pi w }

    :param x: x array
    :param A: Area
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """
    amp = 2*A/(np.pi*w)
    return lorentzian_fwhm(x, amp, c, w)


def voigt_fwhm(x, amp, c, w, m):
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
    lorentz = lorentzian_fwhm(x, 1, c, w)
    gauss = gaussian_fwhm(x, 1, c, w)

    return amp*(m*lorentz + (1-m)*gauss)


def voigt_area_fwhm(x, A, c, w, m):
    r"""Pseudo-voigt curve.

    .. math:: y(x) = A \left[ m \frac{1}{\pi w} \frac{w^2}{w^2 + 4 (x-c)^2}   + (1-m) \frac{2 \sqrt{\ln(2)}}{w \sqrt{\pi}}  e^{-\frac{4 \ln(2) (x-c)^2}{w^2}} \right]

    where,

    .. math:: \text{amp }= \frac{2}{w}  \frac{\text{Area }}{m \pi + (1-m) \frac{\sqrt{\pi}}{\sqrt{\ln(2)} } } 

    :param x: x array
    :param A: is the Area
    :param c: Center
    :param w: FWHM
    :param m: Factor from 1 to 0 of the lorentzian amount
    :return: :math:`y(x)`
    """
    lorentz = lorentzian_fwhm(x, 1, c, w)
    gauss   = gaussian_fwhm(x, 1, c, w)

    amp = 2/w * A/(m*np.pi + ((1-m)*np.sqrt(np.pi))/(np.sqrt(np.log(2))))

    return amp*(m*lorentz + (1-m)*gauss)

def asymmetric_voigt_fwhm(x, amp, c, w1, w2, m1, m2):
    """Asymmetric pseudo-voigt curve 
    
    The first half to the bell curve (lowest x) has a different FWHM as the second
    half the bell curve (highest x).

    The total FHWM of given by

    .. math:: w = \frac{w1 + w2}{2}

    :param x: x array
    :param amp: peak amplitude
    :param c: Center
    :param w1: FWHM of the first half of the bell curve (lowest x)
    :param w2: FWHM of the second half of the bell curve (highest x)
    :param m1: Factor from 1 to 0 of the lorentzian amount of the first half of the bell curve (lowest x)
    :param m2: Factor from 1 to 0 of the lorentzian amount of the second half of the bell curve (highest x)
    :return: :math:`y(x)`
    """
    return np.heaviside(c-x, 0)*voigt_fwhm(x, amp, c, w1, m1) + np.heaviside(x-c, 0)*voigt_fwhm(x, amp, c,  w2,  m2) + dirac_delta(x, amp, c)

def arctan_fwhm(x, amp, c, w):
    r"""Arctangent function.

    .. math:: y(x) =   \frac{A}{\pi} \left[ \arctan(\frac{1}{w}(x-c)) + \frac{\pi}{2} \right]

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """

    return amp * (np.arctan((x - c)/w*2) + (np.pi/2))/np.pi


def square_pulse(x, amp, c, w):
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


def heaviside(x, amp, c, flip=False):
    r"""Heaviside step function.

    .. math::

        y(x) =   \begin{cases}
                    0, & \text{ for    $x < c$}\\
                    \text{amp}, & \text{ for  }  $x < c$\\
                    \end{cases}

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param flip: flip the step. Default is False (1 if x > c)
    :return: :math:`y(x)`
    """
    if flip:
        return np.heaviside(c-x, amp)
    else:
        return np.heaviside(x-c, amp)


def dirac_delta(x, amp, c):
    r"""Crude implementation of Dirac delta function.

    .. math::

        y(x) =   \begin{cases}
                    \text{amp}, & \text{ for  }  $x = c$\\
                    0, & \text{ otherwise}\\
                    \end{cases}

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :return: :math:`y(x)`
    """
    x = list(x)
    final = np.zeros(len(x))
    if c in x:
        final[x.index(c)] = amp
        return final
    else:
        return final
        

def erf_fwhm(x, amp, c, w):
    r"""Error function. Integral of gaussian function calculated by ``scipy.special.erf()``.

    :param x: x array
    :param amp: Amplitude
    :param c: Center
    :param w: FWHM
    :return: :math:`y(x)`
    """
    #     .. math:: y(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} dt
    # return amp/2 * (erf((w**-1)/2*(x - c))+1)
    return amp/2 * (erf((2*np.sqrt(2*np.log(2))+1)/w/2*(x - c)) + 1)