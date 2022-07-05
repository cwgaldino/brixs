#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Peak objects"""

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
from collections.abc import Iterable, MutableMapping

# backpack
from .backpack.arraymanip import  sort
from .backpack.model_functions import voigt_fwhm, dirac_delta

# BRIXS
from . import brixs as br



# %% suport functions ==========================================================
def _peak_function_creator(asymmetry, fixed_m, idx=0):
    """Returns instructions for building peak functions.

    Args:
        assymetry (bool): if True, the returned peak function will require two
            fwhm values, one for each half of the peak.
        fixed_m (False or number): m is the amount of lorentzian
            contribution for a peak. If False, m will be a required as an input
            argument on the returned peak function.
        idx (int, optional): number to be inprinted in the string

    Returns:
        function f(x), function f(x) as string, argument list as string
    """
    if fixed_m == False and type(fixed_m) == bool:  # variable m
        if asymmetry:
            def function2fit(x, amp, c, w1, m1, w2, m2):
                f = np.heaviside(x-c, 0)*voigt_fwhm(x, amp, c, w1, m1) +\
                    np.heaviside(c-x, 0)*voigt_fwhm(x, amp, c, w2, m2) +\
                    dirac_delta(x, amp, c)
                return f
            f_str = f'np.heaviside(x-c_{idx}, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w1_{idx}, m1_{idx}) + np.heaviside(c_{idx}-x, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w2_{idx}, m2_{idx})'
            args_str = f'amp_{idx}, c_{idx}, w1_{idx}, m1_{idx}, w2_{idx}, m2_{idx}'
        else:
            def function2fit(x, amp, c, w, m):
                return voigt_fwhm(x, amp, c, w, m)
            f_str = f'voigt_fwhm(x, amp_{idx}, c_{idx}, w_{idx}, m_{idx})'
            args_str = f'amp_{idx}, c_{idx}, w_{idx}, m_{idx}'
    else:
        if fixed_m > 1:
            m = 1
        elif fixed_m < 0:
            m = 0
        else:
            m = fixed_m
        if asymmetry:
            def function2fit(x, A, c, w1, w2):
                f = np.heaviside(x-c, 0)*voigt_fwhm(x, A, c, w1, fixed_m) +\
                    np.heaviside(c-x, 0)*voigt_fwhm(x, A, c, w2, fixed_m) +\
                    dirac_delta(x, amp, c)
                return f
            f_str = f'np.heaviside(x-c_{idx}, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w1_{idx}, {m}) + np.heaviside(c_{idx}-x, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w2_{idx}, {m})'
            args_str = f'amp_{idx}, c_{idx}, w1_{idx}, w2_{idx}'
        else:
            def function2fit(x, A, c, w):
                return voigt_fwhm(x, A, c, w, fixed_m)
            f_str = f'voigt_fwhm(x, amp_{idx}, c_{idx}, w_{idx}, {m})'
            args_str = f'amp_{idx}, c_{idx}, w_{idx}'
    return function2fit, f_str, args_str

def peak_function(amp, c, fwhm, m=0, fwhm1=None, fwhm2=None, m1=None, m2=None,):
    """Returns a function f(x) of a peak.

    Args:
        amp (number): peak Amplitude (maximum value).
        c (number): peak position.
        fwhm (number): width of the peak (full width half maximum). If peak is
            assymetric (fwhm1 and fwhm2 are defined), fwhm will be set to fwhm1+fwhm2.
        fwhm1 (number, optional): fwhm of the first half of the peak (for
            assimetric peaks).
        fwhm2(number, optional): fwhm of the second half of the peak (for
            assimetric peaks).
        m (number, optional): factor from 0 to 1 of the lorentzian amount.
            Default is 0 (fully gaussian peak). If peak is assymetric (m1 and m2
            are defined), m will be set the average m value.
        m1 (number, optional): m of the first half of the peak (for
            assimetric peaks). Default is 0.
        m2 (number, optional): m of the second half of the peak (for
            assimetric peaks). Default is 0.

    Returns:
        function f(x)
    """
    if fwhm1 is not None or fwhm2 is not None:
        if fwhm2 is None or fwhm2 is None:
            raise ValueError('fwhm1 and fwhm2 must be defined.')
        else:
            asymmetry = True
    else:
        asymmetry = False

    function2fit, _, _ = _peak_function_creator(asymmetry, fixed_m=False, idx=0)

    if asymmetry:
        return lambda x: function2fit(x, amp, c, fwhm1, m1, fwhm2, m2)
    else:
        return lambda x: function2fit(x, amp, c, fwhm, m)

# %% Peaks =====================================================================
class Peak(MutableMapping):
    """A dictionary for saving peaks parameters.

    Only keyword arguments can be used.

    Args:
        amp (number): peak Amplitude (maximum value).
        c (number): peak position.
        fwhm (number): width of the peak (full width half maximum). If peak is
            assymetric (fwhm1 and fwhm2 are defined), fwhm will be set to fwhm1+fwhm2.
        fwhm1 (number, optional): fwhm of the first half of the peak (for
            assimetric peaks).
        fwhm2(number, optional): fwhm of the second half of the peak (for
            assimetric peaks).
        m (number, optional): factor from 0 to 1 of the lorentzian amount.
            Default is 0 (fully gaussian peak). If peak is assymetric (m1 and m2
            are defined), m will be set the average m value.
        m1 (number, optional): m of the first half of the peak (for
            assimetric peaks). Default is 0.
        m2 (number, optional): m of the second half of the peak (for
            assimetric peaks). Default is 0.
        shift (number): initial shift value (does not have initial effect over the data).
        calib (number): initial calibration value (does not have initial effect over the data).
        offset( number): initial offset value (does not have initial effect over the data).
        factor (number): initial multiplicative factor (does not have initial effect over the data).


    Attributes:
        area (number): peak area.
        assymetry (bool): True if peak is assymetric. This is defined when the
            object is created. If assymetry=False, and fhwm1 (or fwhm2) is set,
            then assymetry is set to True. If assymetry=True, it cannot reverse
            to False, and another object must be created from scratch.
        store (dictionary):
        shift( number): shift value (value will be added to x-coordinates).
        calib (number): calibration value (x-coordinates will be multiplied by this value).
        offset( number): offset value (value will be added to y-coordinates).
        factor (number): multiplicative factor (y-coordinates will be multiplied by this value).
    """

    def __init__(self, *args, **kwargs):
        self.store = {'amp':None,
                      'area':None,
                      'c':None,
                      'fwhm':None,
                      'fwhm1':None,
                      'fwhm2':None,
                      'm':None,
                      'm1':None,
                      'm2':None,
                      }
        self._asymmetry = False
        if 'shift' in kwargs:
            self._shift = kwargs.pop('shift')
        else:
            self._shift = 0
        if 'calib' in kwargs:
            self._calib = kwargs.pop('calib')
        else:
            self._calib = 1
        if 'offset' in kwargs:
            self._offset = kwargs.pop('offset')
        else:
            self._offset = 0
        if 'factor' in kwargs:
            self._factor = kwargs.pop('factor')
        else:
            self._factor = 1

        # check keys
        if any([item not in self.store for item in kwargs.keys()]):
            for key in kwargs:
                if key not in self.store:
                    raise AttributeError(f'{key} is not a valid attribute.')

        # minimum requirements to initialize peak
        if 'amp' not in kwargs:
            raise ValueError("Cannot initialize peak without 'amp'.")
        else:
            self._check_amp(kwargs['amp'])
        if 'area' in kwargs:
            raise ValueError("Area is a read-only parameter.")
        if 'c' not in kwargs:
            raise ValueError("Cannot initialize peak without 'c'.")
        if 'fwhm' not in kwargs:
            if 'fwhm1' not in kwargs and 'fwhm2' not in kwargs:
                raise ValueError("Cannot initialize peak without 'fwhm' or ('fwhm1' and 'fhwm2').")
            elif 'fwhm1' in kwargs and 'fwhm2' in kwargs:
                self._asymmetry = True
                self._check_fwhm(kwargs['fwhm1'])
                self._check_fwhm(kwargs['fwhm2'])
                kwargs['fwhm'] = kwargs['fwhm1']+kwargs['fwhm2']
            else:
                raise ValueError(f'When "fwhm1" or "fwhm2" is defined, both must be defined.')
        else:
            if 'fwhm1' in kwargs or 'fwhm2' in kwargs:
                raise ValueError("Cannot initialize peak with 'fwhm', ('fwhm1' and 'fhwm2').\nUse either 'fwhm' or (fwhm1 and fwhm2).")
            self._check_fwhm(kwargs['fwhm'])
        if self.asymmetry:
            if 'm1' not in kwargs:
                kwargs['m1'] = 0
            if 'm2' not in kwargs:
                kwargs['m2'] = 0
            kwargs['m'] = (kwargs['m1']+kwargs['m2'])/2
        else:
            if 'm' not in kwargs:
                kwargs['m'] = 0

        # fix area
        if self.asymmetry:
            c1 = kwargs['m1']/np.pi + (1-kwargs['m1'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            c2 = kwargs['m2']/np.pi + (1-kwargs['m2'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            kwargs['area'] = (kwargs['amp']*kwargs['fwhm1']*c1 + kwargs['amp']*kwargs['fwhm2']*c2)/2
        else:
            c = kwargs['m']/np.pi + (1-kwargs['m'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            kwargs['area'] = kwargs['amp']*kwargs['fwhm']*c
        self.store.update(dict(**kwargs))


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
        self.set_shift(value)
    @shift.deleter
    def shift(self):
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
    def asymmetry(self):
        return self._asymmetry
    @asymmetry.setter
    def asymmetry(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @asymmetry.deleter
    def asymmetry(self):
            raise AttributeError('Cannot delete object.')

    @property
    def spectrum(self):
        return self.calculate_spectrum()
    @spectrum.setter
    def spectrum(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @spectrum.deleter
    def spectrum(self):
            raise AttributeError('Cannot delete object.')

    def __str__(self):
        # return str({name:self.store[name] for name in self.store if self.store[name] is not None})[1:-1].replace(', ', '\n')
        # return str({name:self.store[name] for name in self.store if self.store[name] is not None})[1:-1].replace(', ', '\n')
        return str(self.store)[1:-1].replace(', ', '\n')
        # return str(self.store).replace('}, ', '\n ')

    def __repr__(self):
        return str({name:self.store[name] for name in self.store if self.store[name] is not None})
        # return str(self.store)[1:-1].replace(', ', '\n')
        # return str({i:val for i, val in enumerate(self.store)})[1:-1].replace('}, ', '}\n')
        # return str(self.store).replace('}, ', '\n ')

    def __getitem__(self, name):
        return self.store[self._check_name(name)]

    def __setitem__(self, name, value):
        if name not in self.store:
            raise KeyError(f'{name} is not defined and cannot be created.')
        if name == 'amp':
            self._check_amp(value)
            self.store['amp'] = value
        if name == 'c':
            self.store['c'] = value
        if name == 'area':
            raise ValueError('Area is read-only parameter. Change amp instead.')
        if name == 'fwhm':
            if self.asymmetry:
                raise ValueError('Cannot change fwhm value when asymmetry = True.\nChange fwhm1 and fhwm2 instead.')
            else:
                self.store['fwhm'] = value
        if name == 'fwhm1':
            if self.asymmetry:
                if value is None:
                    print('Setting asymmetry of the peak to False. Done!')
                    self.store['fwhm1'] = None
                    self.store['fwhm2'] = None
                    self.store['m1'] = self.store['m']
                    self.store['m2'] = self.store['m']
                else:
                    self._check_fwhm(value)
                    self.store['fwhm1'] = value
                    self.store['fwhm'] = self.store['fwhm1']+self.store['fwhm2']
            else:
                if value is None:
                    pass
                else:
                    self._check_fwhm(value)
                    print('Setting asymmetry of the peak to True. Done!')
                    self._asymmetry = True
                    self.store['fwhm1'] = value
                    if self.store['fwhm']-self.store['fwhm1'] >= 0:
                        self.store['fwhm2'] = self.store['fwhm']-self.store['fwhm1']
                    else:
                        self.store['fwhm2'] = 0
                        self.store['fwhm'] = self.store['fwhm1']+self.store['fwhm2']
                    self.store['m1'] = self.store['m']
                    self.store['m2'] = self.store['m']
        if name == 'fwhm2':
            if self.asymmetry:
                if value is None:
                    print('Setting asymmetry of the peak to False. Done!')
                    self.store['fwhm1'] = None
                    self.store['fwhm2'] = None
                    self.store['m1'] = None
                    self.store['m2'] = None
                else:
                    self._check_fwhm(value)
                    self.store['fwhm2'] = value
                    self.store['fwhm'] = self.store['fwhm1']+self.store['fwhm2']
            else:
                if value is None:
                    pass
                else:
                    self._check_fwhm(value)
                    print('Setting asymmetry of the peak to True. Done!')
                    self._asymmetry = True
                    self.store['fwhm2'] = value
                    if self.store['fwhm']-self.store['fwhm2'] >= 0:
                        self.store['fwhm1'] = self.store['fwhm']-self.store['fwhm2']
                    else:
                        self.store['fwhm1'] = 0
                        self.store['fwhm'] = self.store['fwhm1']+self.store['fwhm2']
                    self.store['m1'] = 0
                    self.store['m2'] = 0
        if name == 'm':
            if self.asymmetry:
                raise ValueError('cannot change m value when asymmetry = True.\nChange m1 and m2 instead.')
            else:
                self.store['m'] = value
        if name == 'm2':
            if self.asymmetry:
                self.store['m2'] = value
                self.store['m'] = (self.store['m1']+self.store['m2'])/2
            else:
                raise ValueError('cannot set m2 because asymmetry is False.\nChange asymmetry to True first.\nTo change asymmetry just define fwhm1 or fwhm2.')
        if name == 'm1':
            if self.asymmetry:
                self.store['m1'] = value
                self.store['m'] = (self.store['m1']+self.store['m2'])/2
            else:
                raise ValueError('cannot set m1 because asymmetry is False.\nChange asymmetry to True first.\nTo change asymmetry just define fwhm1 or fwhm2.')

        # fix area
        self.calculate_area()

    def __delitem__(self, key):
        raise AttributeError('itens cannot be deleted')

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        # return len(self.store)
        return len({name:self.store[name] for name in self.store if self.store[name] is not None})


    def _check_name(self, name):
        if name not in self.store:
            raise KeyError(f'{name} is not a valid peak parameter.')
        return name

    def _check_amp(self, value):
        pass
        # if value <= 0:
        #     raise ValueError('amp must be a number higher than 0.')

    def _check_area(self, value):
        pass
        # if value <= 0:
        #     raise ValueError('area must be a number higher than 0.')

    def _check_fwhm(self, value):

        if value <= 0:
            raise ValueError(f'fwhm must be a number higher than 0.\nInvalid fwhm:{value}')

    def _check_m(self, value):
        if value < 0 or value > 1:
            raise ValueError(f'm must be a number between 0 and 1.\nInvalid m:{value}')


    def set_calib(self, value):
        """Set calibration value.

        Args:
            value (number): calibration value (x-coordinates will be multiplied
                by this value).

        Returns:
            None
        """
        if self.calib != value:
            if self.calib != 1:
                self.store['c'] = self.store['c']*self.calib**-1
                self.store['fwhm'] = abs(self.store['fwhm']*self.calib**-1)
                if self.asymmetry:
                    self.store['fwhm1'] = abs(self.store['fwhm1']*self.calib**-1)
                    self.store['fwhm2'] = abs(self.store['fwhm2']*self.calib**-1)
            if value != 1:
                self.store['c'] = self.store['c']*value
                self.store['fwhm'] = abs(self.store['fwhm']*value)
                if self.asymmetry:
                    self.store['fwhm1'] = abs(self.store['fwhm1']*value)
                    self.store['fwhm2'] = abs(self.store['fwhm2']*value)
            self._calib = value

            # fix area
            self.calculate_area()

    def set_shift(self, value):
        """Set shift value.

        Args:
            value (float or int): shift value (value will be added to x-coordinates).

        Returns:
            None
        """
        if self.shift != value:
            if self.shift != 0:
                self.store['c'] = self.store['c']-self.shift
            if value != 0:
                self.store['c'] = self.store['c']+value
            self._shift = value

    def set_offset(self, value):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).

        Returns:
            None
        """
        if self.offset != value:
            if self.offset != 0:
                self.store['amp'] = self.store['amp']-self.offset
            if value != 0:
                self.store['amp'] = self.store['amp']+value
            self._offset = value

    def set_factor(self, value):
        """Set y multiplicative factor.

        Args:
            value (number): multiplicative factor (y-coordinates will be
                multiplied by this value).

        Returns:
            None
        """
        if self.factor != value:
            if self.factor != 1:
                self.store['amp'] = self.store['amp']*self.factor**-1
            if value != 1:
                self.store['amp'] = self.store['amp']*value
            self._factor = value

            # fix area
            self.calculate_area()


    def calculate_area(self):
        """Calculate peak area.

        Returns:
            None
        """
        if self.asymmetry:
            c1 = self.store['m1']/np.pi + (1-self.store['m1'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            c2 = self.store['m2']/np.pi + (1-self.store['m2'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            self.store['area'] = (self.store['amp']*self.store['fwhm1']*c1 + self.store['amp']*self.store['fwhm2']*c2)/2
        else:
            c = self.store['m']/np.pi + (1-self.store['m'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            self.store['area'] = self.store['amp']*self.store['fwhm']*c

    def calculate_spectrum(self, x=None):
        """Return peak curve.

        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be constructed.

        Returns:
            :py:class:`Spectrum`.
        """
        if x is None:
            x = np.arange(self['c']-10*self['fwhm'], self['c']+10*self['fwhm'], self['fwhm']/20)
        # non_none = {name:self.store[name] for name in self.store if self.store[name] is not None}
        # f = peak_function(**non_none)
        f = peak_function(amp=self['amp'], c=self['c'],  fwhm=self['fwhm'], m=self['m'], fwhm1=self['fwhm1'], fwhm2=self['fwhm2'], m1=self['m1'], m2=self['m2'])
        s = br.Spectrum(x=x, y=f(x))
        s._shift  = self.shift
        s._factor = self.factor
        s._calib  = self.calib
        s._offset = self.offset
        return s

    def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, **kwargs):
        """Place a marker at the maximum of every peak position. Wrapper for `matplotlib.pyplot.errorbar()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number, optional): defines a vertical offset. Default is 0.
            shift (number, optional): horizontal shift value. Default is 0.
            factor (number, optional): multiplicative factor on the y axis.
                Default is 1.
            calib (number, optional): multiplicative factor on the x axis.
                Default is 1.
            **kwargs: kwargs are passed to `matplotlib.pyplot.errorbar()`_ that plots the data.

        Returns:
            `ErrorbarContainer`_

        .. matplotlib.pyplot.errorbar(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
        .. ErrorbarContainer: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
        """
        if ax is None:
            ax = plt

        # data
        c = self['c']
        amp = self['amp']
        xerr = self['fwhm']/2

        if 'lw' not in kwargs and 'linewidth' not in kwargs:
            kwargs['lw'] = 0
        if 'elinewidth' not in kwargs :
            kwargs['elinewidth'] = 2
        if 'marker' not in kwargs :
            kwargs['marker'] = 'o'
        if 'markersize' not in kwargs and 'ms' not in kwargs:
            kwargs['markersize'] = 5

        return ax.errorbar((c*calib)+shift, amp*factor+offset, xerr=xerr*calib, **kwargs)

    def build_guess(self, asymmetry, fixed_m, amp_bounds=(0, 3), c_bounds=(-2, 2), fwhm_bounds=(0, 2)):
        """Returns initial guess parameter for fitting.

        Args:
            assymetry (bool): if True, the returned peak function will require two
                fwhm values, one for each half of the peak.
            fixed_m (False or number): m is the amount of lorentzian
                contribution for a peak. If False, m will be a required as an input
                argument on the returned peak function.
            amp_bounds, c_bounds, fwhm_bounds (tuple, optional): minimum and
                maximum multiplication factor for boundary values. For amp, the
                bounds are set between amp*amp_bounds[0] and amp*amp_bounds[1].
                For c, bounds are set c+fwhm*c_bound[0] and c+fwhm*c_bound[1].
                Finaly, for fwhm, the bounds are set between fwhm1*fwhm_bounds[0]
                and fwhm1*fwhm_bounds[0]. Note that if fwhm1*fwhm_bounds[0] is
                less than zero, the minimum fwhm boundary is set to zero as it
                cannot be negative.

        Returns:
            three lists: p0, bounds_min, bounds_max
        """
        # check values
        # assert amp_bounds[1]  > amp_bounds[0],  f'Minimum amp boundary cannot be bigger than the maximum.\namp_bounds = {amp_bounds}'
        # assert c_bounds[1]    > c_bounds[0],    f'Minimum c boundary cannot be bigger than the maximum.\nc_bounds = {c_bounds}'
        # assert fwhm_bounds[1] > fwhm_bounds[0], f'Minimum fwhm boundary cannot be bigger than the maximum.\nfwhm_bounds = {fwhm_bounds}'
        assert amp_bounds[1]  != amp_bounds[0],  f'Minimum amp boundary cannot be equal to the maximum.\namp_bounds = {amp_bounds}'
        assert c_bounds[1]    != c_bounds[0],    f'Minimum c boundary cannot be equal to the maximum.\nc_bounds = {c_bounds}'
        assert fwhm_bounds[1] != fwhm_bounds[0], f'Minimum fwhm boundary cannot be equal to the maximum.\nfwhm_bounds = {fwhm_bounds}'


        p0         = [       self['amp'],                    self['c']]
        bounds_min = [self['amp']*amp_bounds[0], self['c']+self['fwhm']*c_bounds[0]]
        bounds_max = [self['amp']*amp_bounds[1], self['c']+self['fwhm']*c_bounds[1]]
        assert self['amp'] > self['amp']*amp_bounds[0]                   and self['amp'] < self['amp']*amp_bounds[1],        f'guess amp is out of bounds\namp = {p0[-2]}\nbounds = {(bounds_min[-2], bounds_max[-2])}'
        assert self['c'] > self['c']+self['fwhm']*c_bounds[0] and self['c'] < self['c']+self['fwhm']*c_bounds[1], f'guess c is out of bounds\nc = {p0[-1]}\nbounds = {(bounds_min[-1], bounds_max[-1])}'

        if fixed_m == False and type(fixed_m) == bool:  # variable m
            if asymmetry:
                if self['fwhm1'] is not None:
                    p0.extend(        [      self['fwhm1'],       0.5,       self['fwhm2'],         0.5])
                    bounds_min.extend([self['fwhm1']*fwhm_bounds[0], 0, self['fwhm2']*fwhm_bounds[0], 0])
                    bounds_max.extend([self['fwhm1']*fwhm_bounds[1], 1, self['fwhm2']*fwhm_bounds[1], 1])
                    assert self['fwhm1'] > self['fwhm1']*fwhm_bounds[0] and self['fwhm1'] < self['fwhm1']*fwhm_bounds[1], f'guess fwhm1 is out of bounds\namp = {p0[-4]}\nbounds = {(bounds_min[-4], bounds_max[-4])}'
                    assert self['fwhm2'] > self['fwhm2']*fwhm_bounds[0] and self['fwhm2'] < self['fwhm2']*fwhm_bounds[1], f'guess fwhm2 is out of bounds\namp = {p0[-2]}\nbounds = {(bounds_min[-2], bounds_max[-2])}'
                else:
                    p0.extend(        [      self['fwhm'],       0.5,       self['fwhm'],         0.5])
                    bounds_min.extend([self['fwhm']*fwhm_bounds[0], 0, self['fwhm']*fwhm_bounds[0], 0])
                    bounds_max.extend([self['fwhm']*fwhm_bounds[1], 1, self['fwhm']*fwhm_bounds[1], 1])
                    assert self['fwhm'] > self['fwhm']*fwhm_bounds[0] and self['fwhm'] < self['fwhm']*fwhm_bounds[1], f'guess fwhm is out of bounds\namp = {p0[-2]}\nbounds = {(bounds_min[-2], bounds_max[-2])}'
                if bounds_min[-2] < 0:
                    bounds_min[-2] = 0
                if bounds_min[-4] < 0:
                    bounds_min[-4] = 0
            else:
                p0.extend(        [       self['fwhm'],       0.5])
                bounds_min.extend([self['fwhm']*fwhm_bounds[0], 0])
                bounds_max.extend([self['fwhm']*fwhm_bounds[1], 1])
                assert self['fwhm'] > self['fwhm']*fwhm_bounds[0] and self['fwhm'] < self['fwhm']*fwhm_bounds[1], f'guess fwhm is out of bounds\namp = {p0[-2]}\nbounds = {(bounds_min[-2], bounds_max[-2])}'
                if bounds_min[-2] < 0:
                    bounds_min[-2] = 0
        else:
            if asymmetry:
                if self['fwhm1'] is not None:
                    p0.extend(        [        self['fwhm1'],              self['fwhm2']         ])
                    bounds_min.extend([self['fwhm1']*fwhm_bounds[0], self['fwhm2']*fwhm_bounds[0]])
                    bounds_max.extend([self['fwhm1']*fwhm_bounds[1], self['fwhm2']*fwhm_bounds[1]])
                    assert self['fwhm1'] > self['fwhm1']*fwhm_bounds[0] and self['fwhm1'] < self['fwhm1']*fwhm_bounds[1], f'guess fwhm1 is out of bounds\namp = {p0[-2]}\nbounds = {(bounds_min[-2], bounds_max[-2])}'
                    assert self['fwhm2'] > self['fwhm2']*fwhm_bounds[0] and self['fwhm2'] < self['fwhm2']*fwhm_bounds[1], f'guess fwhm2 is out of bounds\namp = {p0[-1]}\nbounds = {(bounds_min[-1], bounds_max[-1])}'
                else:
                    p0.extend(        [        self['fwhm'],              self['fwhm']         ])
                    bounds_min.extend([self['fwhm']*fwhm_bounds[0], self['fwhm']*fwhm_bounds[0]])
                    bounds_max.extend([self['fwhm']*fwhm_bounds[1], self['fwhm']*fwhm_bounds[1]])
                assert self['fwhm'] > self['fwhm']*fwhm_bounds[0] and self['fwhm'] < self['fwhm']*fwhm_bounds[1], f'guess fwhm is out of bounds\namp = {p0[-1]}\nbounds = {(bounds_min[-1], bounds_max[-1])}'
                if bounds_min[-1] < 0:
                    bounds_min[-1] = 0
                if bounds_min[-2] < 0:
                    bounds_min[-2] = 0
            else:
                p0.extend(        [       self['fwhm']        ])
                bounds_min.extend([self['fwhm']*fwhm_bounds[0]])
                bounds_max.extend([self['fwhm']*fwhm_bounds[1]])
                assert self['fwhm'] > self['fwhm']*fwhm_bounds[0] and self['fwhm'] < self['fwhm']*fwhm_bounds[1], f'guess fwhm is out of bounds\namp = {p0[-1]}\nbounds = {(bounds_min[-1], bounds_max[-1])}'
                if bounds_min[-1] < 0:
                    bounds_min[-1] = 0

        return p0, bounds_min, bounds_max

class Peaks(MutableMapping):
    """A special dictionary for saving peaks.

    Data can be passed as positional argument. The other arguments can only be
        passed as keyword arguments.

    Example:
        .. code-block:: python

            p1 = br.Peak(amp=1, c=0, fwhm=10)
            p2 = br.Peak(amp=1, c=1, fwhm1=10, fwhm2=5)

            ps = br.Peaks(p1, p2)
            ps = br.Peaks([p1, p2])
            ps = br.Peaks([p1, p2], shift=10)

    Args:
        data (Peak or list): can be passed as positional arguments
        shift( number): initial shift value (does not have initial effect over the data).
        calib (number): initial calibration value (does not have initial effect over the data).
        offset( number): initial offset value (does not have initial effect over the data).
        factor (number): initial multiplicative factor (does not have initial effect over the data).
        *args: Peak objects
    """

    def __init__(self, *args, **kwargs):
        self.store = []
        if 'shift' in kwargs:
            self._shift = kwargs.pop('shift')
        else:
            self._shift = 0
        if 'calib' in kwargs:
            self._calib = kwargs.pop('calib')
        else:
            self._calib = 1
        if 'offset' in kwargs:
            self._offset = kwargs.pop('offset')
        else:
            self._offset = 0
        if 'factor' in kwargs:
            self._factor = kwargs.pop('factor')
        else:
            self._factor = 1
        if 'data' in kwargs:
            if isinstance(kwargs['data'], dict) or isinstance(kwargs['data'], Peak):#(id(Peak) == id(kwargs['data'].__class__)):#
                self.append(kwargs['data'])
            elif isinstance(kwargs['data'], Iterable):
                for p in kwargs['data']:
                    self.append(p)
            else:
                raise ValueError('data must be a list of peaks (dictionaries).')
        else:
            if len(args) == 1:
                if isinstance(args[0], dict) or isinstance(args[0], Peak): #(id(Peak) == id(args[0].__class__)):#
                    self.append(args[0])
                elif isinstance(args[0], Iterable):
                    for p in args[0]:
                        self.append(p)
                else:
                    raise ValueError('data must be a list of peaks (dictionaries).')
            elif len(args) > 1:
                for p in args:
                    self.append(p)


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
        self.set_shift(value)
    @shift.deleter
    def shift(self):
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


    def __str__(self):
        return str({i:val for i, val in enumerate(self.store)})[1:-1].replace('}, ', '}\n')

    def __repr__(self):
        return str({i:val for i, val in enumerate(self.store)})[1:-1].replace('}, ', '}\n')

    def __getitem__(self, key):
        return self.store[self._check_key(key)]

    def __setitem__(self, key, value):

        if isinstance(value, Peak):
            self.store[self._check_key(key)] = value
        elif isinstance(value, dict):
            self.store[self._check_key(key)] = Peak(**value)
        else:
            raise ValueError('valueaaaaa must be a dict or a peak object')
        self._fix_order()

    def __delitem__(self, key):
        del self.store[self._check_key(key)]
        self._fix_order()

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)


    def _check_key(self, key):
        """Check if key exists. Allows for minus (-) assignment."""
        if key > len(self)-1 or key < -len(self):
            raise KeyError('key out of range of defined peaks.\n')
        return key

    def _fix_order(self):
        """Returns another PeakDict where peak number is numbered according to its position c."""
        if self == []:
            return
        else:
            l = len(self)
            c = [0]*l
            for i, peak in enumerate(self):
                c[i] = peak['c']
            if sum(1 for test in np.diff(c) if test < 0) > 0:
                self.store = sort(c, self.store)


    def append(self, value):
        """Append peak. Peak order is reassigned.

        Args:
            value (Peak or dict): peak to be appended.

        Returns:
            None
        """
        if isinstance(value, Peak):
        # if id(Peak) == id(value.__class__):
            self.store.append(value)
        elif isinstance(value, dict):
            self.store.append(Peak(**value))
        else:
            raise ValueError('value must be a dict or a peak object')
        self._fix_order()

    def remove(self, key):
        """Remove peak. Peak order is reassigned.

        Args:
            key (int or list): peaks to be removed.

        Returns:
            None
        """
        if isinstance(key, Iterable):
            key = [self._check_key(k) for k in key]
            key = sort(key)[::-1]
            if has_duplicates(key):
                raise ValueError('list has duplicated peaks')
            for k in key:
                del self.store[k]
        else:
            del self.store[key]


    def set_calib(self, value):
        """Set calibration value.

        Args:
            value (number): calibration value (x-coordinates will be multiplied
                by this value).

        Returns:
            None
        """
        for peak in self.store:
            peak.calib = value
        self._calib = value
        self._fix_order()

    def set_shift(self, value):
        """Set shift value.

        Args:
            value (float or int): shift value (value will be added to x-coordinates).

        Returns:
            None
        """
        for peak in self.store:
            peak.shift = value
        self._shift = value

    def set_offset(self, value):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).

        Returns:
            None
        """
        for peak in self.store:
            peak.offset = value
        self._offset = value

    def set_factor(self, value):
        """Set y multiplicative factor.

        Args:
            value (number): multiplicative factor (y-coordinates will be
                multiplied by this value).

        Returns:
            None
        """
        for peak in self.store:
            peak.factor = value
        self._factor = value


    def split(self, key, n=1):
        """Splits a peak in n symmetric peaks.

        Amplitude and fwhm are divided between the peaks. Position of the peaks
            are equaly spread between c-fwhm and c+fwhm.

        Args:
            key (int or list): key of the peak.
            n (int, optional): number of extra peaks.

        Returns:
            None
        """
        fwhm = self[key]['fwhm']/(n+1)
        amp  = self[key]['amp']/(n+1)
        cs = np.linspace(self[key]['c']-self[key]['fwhm']/2, self[key]['c']+self[key]['fwhm']/2, n+1)
        m = self[key]['m']

        self.remove(key)
        for c in cs:
            self.append(Peak(amp=amp, fwhm=fwhm, c=c, m=m))


    def add_near(self, key):
        """Appends an extra peak 1/4 of the fwhm away from a peak.

        Args:
            key (int or list): key of the peak. If list, use repeated elements
                to append multiple peaks.

        Returns:
            None
        """
        if isinstance(key, Iterable):
            temp = []
            for i in range(len(self.store)):
                count = len([k2 for k2 in key if k2 == i])   # counts same key
                if count > 0:
                    for n in range(count):
                        temp = copy.deepcopy(self.store[self._check_key(key)])
                        temp['c']+=temp['fwhm']/4
                        self.append(temp)
            for value in temp:
                self.append(value)
        else:
            temp = copy.deepcopy(self.store[self._check_key(key)])
            self.store[self._check_key(key)]['c'] -= self.store[self._check_key(key)]['fwhm']/4
            temp['c'] += temp['fwhm']/4
            self.append(temp)


    def calculate_spectrum(self, x=None):
        """Return the spectrum with all peaks.

        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be constructed.

        Returns:
            :py:class:`Spectrum`, :py:class:`Spectra`.
        """
        if len(self) > 0:
            if x is None:
                # find minimum value
                vmin = self[0]['c']-10*self[0]['fwhm']
                for i in range(len(self)):
                    v = self[i]['c']-10*self[i]['fwhm']
                    if v < vmin:
                        vmin = v
                # find maximum value
                vmax = self[-1]['c']+10*self[-1]['fwhm']
                for i in range(len(self)):
                    v = self[i]['c']+10*self[i]['fwhm']
                    if v > vmax:
                        vmax = v
                # find smallest fwhm
                fwhm_min = self[0]['fwhm']
                for i in range(len(self)):
                    fwhm = self[i]['fwhm']
                    if fwhm < fwhm_min:
                        fwhm_min = fwhm

                x = np.arange(vmin, vmax, fwhm_min/20)

            ss = br.Spectra(n=len(self))
            for i in range(len(self)):
                ss[i] = self[i].calculate_spectrum(x=x)

            s = br.Spectrum(x=x, y=np.zeros(len(x)))
            for s1 in ss:
                s += s1

            return s, ss
        else:
            raise ValueError('No peaks defined.')


    def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, **kwargs):
        """Place a marker at the maximum of every peak position. Wrapper for `matplotlib.pyplot.errorbar()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number, optional): defines a vertical offset. Default is 0.
            shift (number, optional): horizontal shift value. Default is 0.
            factor (number, optional): multiplicative factor on the y axis.
                Default is 1.
            calib (number, optional): multiplicative factor on the x axis.
                Default is 1.
            **kwargs: kwargs are passed to `matplotlib.pyplot.errorbar()`_ that plots the data.

        Returns:
            dict with `ErrorbarContainer`_

        .. matplotlib.pyplot.errorbar(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
        .. ErrorbarContainer: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
        """
        if ax is None:
            ax = plt

        r = {}
        for i in range(len(self)):
            r[i] = self[i].plot(ax=ax, offset=offset, shift=shift, factor=factor, **kwargs)

        return r
