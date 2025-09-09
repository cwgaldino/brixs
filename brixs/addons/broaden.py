#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Spectrum broadening functions"""

# %% ------------------------- Standard Imports -------------------------- %% #
import numpy as np
# %%

# %% ------------------------- Special Imports ---------------=----------- %% #
try:
    from scipy.ndimage import gaussian_filter1d
except ModuleNotFoundError:
    pass
# %%

# %% -------------------------- brixs Imports ---------------------------- %% #
import brixs as br
# %%

# %% ============================ broadening ============================= %% #
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

def _broaden(self, w, m=0, mode='same'):
    """Return broadened spectrum (gaussian or lorentzian convolution).

    Warning: 
        Convolution may cause boundary effects at the edges of the spectrum.
 
    Args:
        w (number): fwhm value.
        m (number, optional): Lorentzian amount from 0 to 1. Default is 0.
        mode (str, optional): if `same`, length of spectrum is preserved, but 
            boundary effects are still visible. If `valid`, the edges of the
            spectrum are cropped to remove boundary effects.

    Return:
        Convolved spectrum        
    """
    assert mode in ['same', 'valid'], 'mode can only be `same` or `valid`.'
    assert w >= 0, 'w must be a positive number'
    assert m >=0 and m<=1, 'm must be a number between 0 and 1'

    # if broadening is zero, does nothing
    if w == 0:
        return self.copy()

    # check step uniformity
    self.check_step()

    # get x and y
    step = self.step
    x = self.x
    y = self.y

    # error message
    error = 'Broadening seems to be too large for this data range ' +\
        f'({min(x)}, {max(x)}) and boundary effects will dominate the broadened spectrum.' +\
        'Please, reduce `w`.' +\
        'If possible, one can also increase the data range of the spectrum' +\
        'One can `fake` a bigger x range by interpolating the spectrum with a larger `x` array [ s.interp(x=x) ]'
    
    # gaussian
    x1 = np.arange(-4*w, 4*w, step)
    y1 = _voigt_fwhm(x1, 1/np.exp(1), 0, w, m)
    if len(x1) > len(x):
        raise ValueError(error)
    
    # check if w is too large
    # if len(x1) > 0.3*len(x):
    #     raise ValueError(error)
    # Raises:
    #     ValueError: if Broadening is so large for the data range (min(x), max(x)) 
    #         so that boundary effects dominate the broadened spectrum.

    # get new x and y
    if mode == 'same':
        x2 = x
        y2 = np.convolve(y, y1, mode='same')
        y2[:-1] = y2[1:]  # padding necessary because of np.convolve
    else:
        j = (len(y1) + 1)/2  # (max(M, N) - min(M, N) + 1)
        if br.is_integer(j):
            x2 = x[int(j):-int(j)]
        else:
            x2 = x[int(j)-1:-int(j)]
        y2 = np.convolve(y, y1, mode='valid')

    # final spectrum
    s = self.copy()
    a1 = s.calculate_area()
    s._x = x2
    s._y = y2
    a2 = s.calculate_area()
    s._y = y2*a1/a2

    return s
br.Spectrum.broaden = _broaden

def _broaden_x(self, w, m=0):
    """Return broadened image in the x direction (gaussian or lorentzian convolution).
    
    Note: 
        Convolution may cause boundary effects at the edges of the spectrum.

    Raises:
        ValueError: if Broadening is so large for the data range (min(x), max(x)) 
            so that boundary effects dominate the broadened spectrum.
        
    Args:
        w (number): fwhm value.
        m (number, optional): Lorentzian amount from 0 to 1. Default is 0.

    Return:
        Convolved image        
    """
    # get rows
    _ss1 = self.get_rows(max_number_of_rows=self.shape[0])   

    # broaden
    _ss2 = br.Spectra()
    for _s in _ss1:
        _ss2.append(_s.broaden(w=w, m=m, mode='same'))
    _im2 = _ss2.stack_spectra_as_rows()

    # final image
    im = self.copy()
    im._data = _im2.data
    # im.x_centers = self.x_centers
    # im.y_centers = self.y_centers

    return im
br.Image.broaden_x = _broaden_x

def _broaden_y(self, w, m=0):
    """Return broadened image in the y direction (gaussian or lorentzian convolution).
    
    Note: 
        Convolution may cause boundary effects at the edges of the spectrum.

    Raises:
        ValueError: if Broadening is so large for the data range (min(y), max(y)) 
            so that boundary effects dominate the broadened spectrum.
        
    Args:
        w (number): fwhm value.
        m (number, optional): Lorentzian amount from 0 to 1. Default is 0.

    Return:
        Convolved image        
    """
    # get rows
    _ss1 = self.get_columns(max_number_of_columns=self.shape[1])   

    # broaden
    _ss2 = br.Spectra()
    for _s in _ss1:
        _ss2.append(_s.broaden(w=w, m=m, mode='same'))
    _im2 = _ss2.stack_spectra_as_columns()

    # final image
    im = self.copy()
    im._data = _im2.data
    # im.x_centers = self.x_centers
    # im.y_centers = self.y_centers

    return im
br.Image.broaden_y = _broaden_y

def _broaden_spectra(self, w, m=0, mode='same'):
    """Return broadened spectra (gaussian or lorentzian convolution).

    Note: 
        Convolution may cause boundary effects at the edges of the spectrum.

    Raises:
        ValueError: if Broadening is so large for the data range (min(x), max(x)) 
            so that boundary effects dominate the broadened spectrum.
 
    Args:
        w (number): fwhm value.
        m (number, optional): Lorentzian amount from 0 to 1. Default is 0.
        mode (str, optional): if `same`, length of spectrum is preserved, but 
            boundary effects are still visible. If `valid`, the edges of the
            spectrum are cropped to remove boundary effects.

    Return:
        Convolved spectrum        
    """
    # check same X
    ss = self.copy()
    for i, _s in enumerate(self):
        ss[i] = _s.broaden(w=w, m=m, mode=mode)

    return ss
br.Spectra.broaden = _broaden_spectra
# %%

# %% ===================== broadening via scipy filter =================== %% #
def _broaden2(self, w):
    """Gaussian broadening via scipy gaussian_filter1d

    Args:
         w (number): fwhm value.

    Return:
        broadened spectrum
    """
    assert w >= 0, 'w must be a positive number'

    # check step uniformity
    self.check_step()

    # get x, y, and step
    step = self.step
    x    = self.x
    y    = self.y

    # final spectrum
    s = self.copy()
    try:
        s._y = gaussian_filter1d(y, w/step/(2*np.sqrt(2*np.log(2))))
    except NameError:
        raise ModuleNotFoundError('broaden2() function requires the scipy package installed.\nPlease, install the scipy package or use one of the other broaden methods such as s.broaden() and s.broadenfft()')

    return s
br.Spectrum.broaden2 = _broaden2
# %%