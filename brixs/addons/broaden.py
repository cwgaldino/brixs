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

    # gaussian
    x1 = np.arange(-4*w, 4*w, step)
    y1 = _voigt_fwhm(x1, 1/np.exp(1), 0, w, m)

    # check if w is too large
    error = 'Broadening seems to be too large for this data range ' +\
            f'({min(x)}, {max(x)}) and boundary effects will dominate the broadened spectrum.' +\
            'Please, reduce `w`.' +\
            'If possible, one can also increase the data range of the spectrum' +\
            'One can `fake` a bigger x range by interpolating the spectrum with a larger `x` array [ s.interp(x=x) ]'
    if len(x1) > 0.3*len(x):
        raise ValueError(error)

    # get new x and y
    if mode == 'same':
        x2 = x
        y2 = np.convolve(y, y1, mode='same')
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
        _ss1.append(_s.broaden(w=w, m=m, mode='same'))
    _im2 = _ss2.stack_spectra_as_rows()

    # final image
    im = self.copy()
    im._data = _im2.data

    return _im2
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
    _ss1 = self.get_columns(max_number_of_rows=self.shape[1])   

    # broaden
    _ss2 = br.Spectra()
    for _s in _ss1:
        _ss1.append(_s.broaden(w=w, m=m, mode='same'))
    _im2 = _ss2.stack_spectra_as_columns()

    # final image
    im = self.copy()
    im._data = _im2.data

    return _im2
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

# %% =========== broadening via fft (copied from crispy 0.7.3) =========== %% #
MIN_KERNEL_SUM = 1e-8

def _gaussian_kernel1d(sigma=None, truncate=6):
    size = int(2 * truncate * sigma)
    if size % 2 == 0:
        size = size + 1
    x = np.arange(size)
    # print('The size of the kernel is: {}'.format(size))
    mu = np.median(x)
    # The prefactor 1 / (sigma * np.sqrt(2 * np.pi))
    # drops in the normalization.
    kernel = np.exp(-0.5 * ((x - mu)**2 / sigma**2))
    if kernel.sum() < MIN_KERNEL_SUM:
        raise Exception(
            'The kernel can\'t be normalized, because its sum is close to '
            'zero. The sum of the kernel is < {0}'.format(MIN_KERNEL_SUM))
    kernel /= kernel.sum()
    return kernel

def _gaussian_kernel2d(sigma=None, truncate=(6, 6)):
    if sigma.size != 2 or len(truncate) != 2:
        raise Exception('Sigma and the truncation parameter don\'t have the '
                        'required dimenstion.')
    kernel_x = _gaussian_kernel1d(sigma[0], truncate[0])
    kernel_y = _gaussian_kernel1d(sigma[1], truncate[1])
    kernel = np.outer(kernel_y, kernel_x)
    return kernel

def _convolve_fft(array, kernel):
    """
    Convolve an array with a kernel using FFT.
    Implemntation based on the convolve_fft function from astropy.

    https://github.com/astropy/astropy/blob/master/astropy/convolution/convolve.py
    """

    array = np.asarray(array, dtype=complex)
    kernel = np.asarray(kernel, dtype=complex)

    if array.ndim != kernel.ndim:
        raise ValueError("Image and kernel must have same number of "
                         "dimensions")

    array_shape = array.shape
    kernel_shape = kernel.shape
    new_shape = np.array(array_shape) + np.array(kernel_shape)

    array_slices = []
    kernel_slices = []
    for (new_dimsize, array_dimsize, kernel_dimsize) in zip(
            new_shape, array_shape, kernel_shape):
        center = new_dimsize - (new_dimsize + 1) // 2
        array_slices += [slice(center - array_dimsize // 2,
                         center + (array_dimsize + 1) // 2)]
        kernel_slices += [slice(center - kernel_dimsize // 2,
                          center + (kernel_dimsize + 1) // 2)]

    array_slices = tuple(array_slices)
    kernel_slices = tuple(kernel_slices)

    if not np.all(new_shape == array_shape):
        big_array = np.zeros(new_shape, dtype=complex)
        big_array[array_slices] = array
    else:
        big_array = array

    if not np.all(new_shape == kernel_shape):
        big_kernel = np.zeros(new_shape, dtype=complex)
        big_kernel[kernel_slices] = kernel
    else:
        big_kernel = kernel

    array_fft = np.fft.fftn(big_array)
    kernel_fft = np.fft.fftn(np.fft.ifftshift(big_kernel))

    rifft = np.fft.ifftn(array_fft * kernel_fft)

    return rifft[array_slices].real

def _broadenfft(self, w):
    """Gaussian broadening via FFT

    Args:
         w (number): fwhm value.

    Return:
        broadened spectrum
    """
    assert w >= 0, 'w must be a positive number'

    # check step uniformity
    self.check_step()  # not 100% if uniform step is required. Check!

    # get x and y
    step = self.step
    x = self.x
    y = self.y


    sigma = w / (2 * np.sqrt(2 * np.log(2)))
    kernel = _gaussian_kernel1d(sigma)
    # elif w.size == 2:
    #     kernel = gaussian_kernel2d(sigma)

    # final spectrum
    s = self.copy()
    s._y = _convolve_fft(y, kernel)

    return s
br.Spectrum.broadenfft = _broadenfft
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