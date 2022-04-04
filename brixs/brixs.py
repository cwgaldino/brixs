#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Module for analysis of RIXS spectra.

.. autosummary::

    Image
    PhotonEvents
    Spectrum
    Spectra

IDEAS:
- when you crop the data in Spectra, you modify the spectrum. There's no way to
go back (like we do with the shift and offset). Maybe figure a away to go back.


Todo:
    * crop() method
    * flip() or transpose() method
    * noise_filter() method (remove cosmic rays)
    * save image as tiff
"""

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
from collections.abc import Iterable, MutableMapping
from scipy.signal import find_peaks
from scipy.optimize import curve_fit

# backpack
from .backpack.filemanip import save_data, save_obj, load_obj, load_data, load_Comments, filelist
from .backpack.arraymanip import index, moving_average, extract, shifted, sort, get_attributes
from .backpack.arraymanip import is_integer, all_equal, factors, flatten, peak_fit
from .backpack.figmanip import n_digits, n_decimal_places
from .backpack.model_functions import gaussian_fwhm, voigt_fwhm

# common definitions
cc = ['cross-correlation', 'cc']
roll = ['roll', 'rotate', 'r', 'rot']
hard = ['hard', 'x', 'h', 'Hard']
soft = ['soft', 'Soft', 'interp', 'y', 's']


class _Meta(type):
    """Metaclass to facilitate creation of read-only and non-removable attributes."""
    def __new__(self, class_name, bases, attrs):

        def lazy_read_only(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    raise AttributeError('Attribute is "read only". Cannot set attribute.')
                def deleter(self):
                    raise AttributeError('Attribute is "read only". Cannot delete object.')
            return getter, setter, deleter, 'read only attribute'

        def lazy_non_removable(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    return setattr(self, variable, value)
                def deleter(self):
                    raise AttributeError('Attribute cannot be deleted.')
            return getter, setter, deleter, 'non-removable attribute'

        new_attrs = {}
        for name, value in attrs.items():
            if name == '_read_only':
                for attr in value:
                    _property = property(*lazy_read_only(attr))
                    new_attrs[attr] = _property
            elif name == '_non_removable':
                for attr in value:
                    _property = property(*lazy_non_removable(attr))
                    # print(_property.__doc__)
                    new_attrs[attr] = _property
            else:
                new_attrs[name] = value

        return type(class_name, bases, new_attrs)


class _Peak(MutableMapping):
    """A special dictionary for saving peaks parameters.

    Keys must be integer. Value of the keys are always numbered according to the
        center `c` of each peak.


    minimum peak attributes:
    amp
    c
    fwhm (or fhwm1 and fwhm2)

    Elements must be dictionaries with the following items:
    amp
    c
    fwhm
    fwhm1
    fwhm2
    m
    m1
    m2

    shift
    calib
    factor
    offset

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
        # if 'amp' not in kwargs and 'area' not in kwargs:
        #     raise ValueError("Cannot initialize peak without at least 'amp' or 'area'.")
        # elif 'amp' in kwargs and 'area' in kwargs:
        #     raise ValueError("Cannot initialize peak with both 'amp' and 'area' defined, pick one.")
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
        # self.update(dict(*args, **kwargs))  # use the free update to set keys


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
            raise ValueError('area is read-only parameter. Change amp instead.')
        if name == 'fwhm':
            if self.asymmetry:
                raise ValueError('cannot change fwhm value when asymmetry = True.\nChange fwhm1 and fhwm2 instead.')
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
        if value <= 0:
            raise ValueError('amp must be a number higher than 0.')

    def _check_area(self, value):
        if value <= 0:
            raise ValueError('area must be a number higher than 0.')

    def _check_fwhm(self, value):
        # print(value)
        if value <= 0:
            raise ValueError('fwhm must be a number higher than 0.')

    def _check_m(self, value):
        if value < 0 or value > 1:
            raise ValueError('m must be a number between 0 and 1.')


    def set_calib(self, value):
        """Calibrate data.

        Args:
            value (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].

        Returns:
            None
        """
        if self.calib != value:
            if self.calib != 1:
                self.store['c'] = self.store['c']*self.calib**-1
                self.store['fwhm'] = self.store['fwhm']*self.calib**-1
                if self.asymmetry:
                    self.store['fwhm1'] = self.store['fwhm1']*self.calib**-1
                    self.store['fwhm2'] = self.store['fwhm2']*self.calib**-1
            if value != 1:
                self.store['c'] = self.store['c']*value
                self.store['fwhm'] = self.store['fwhm']*value
                if self.asymmetry:
                    self.store['fwhm1'] = self.store['fwhm1']*value
                    self.store['fwhm2'] = self.store['fwhm2']*value
            self._calib = value

            # fix area
            self.calculate_area()

    def set_shift(self, value):
        """Shift data.
        """
        if self.shift != value:
            if self.shift != 0:
                self.store['c'] = self.store['c']-self.shift
            if value != 0:
                self.store['c'] = self.store['c']+value
            self._shift = value

    def set_offset(self, value):
        if self.offset != value:
            # if self.offset != 0:
            #     self.store['amp'] = self.store['amp']-self.offset
            # if value != 0:
            #     self.store['amp'] = self.store['amp']+value
            self._offset = value

    def set_factor(self, value):

        if self.factor != value:
            if self.factor != 1:
                self.store['amp'] = self.store['amp']*self.factor**-1
            if value != 1:
                self.store['amp'] = self.store['amp']*value
            self._factor = value

            # fix area
            self.calculate_area()

    def calculate_area(self):
        if self.asymmetry:
            c1 = self.store['m1']/np.pi + (1-self.store['m1'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            c2 = self.store['m2']/np.pi + (1-self.store['m2'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            self.store['area'] = (self.store['amp']*self.store['fwhm1']*c1 + self.store['amp']*self.store['fwhm2']*c2)/2
        else:
            c = self.store['m']/np.pi + (1-self.store['m'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            self.store['area'] = self.store['amp']*self.store['fwhm']*c

    def get_curve(self, x):
        non_none = {name:self.store[name] for name in self.store if self.store[name] is not None}
        peak = _peak(**non_none)
        s = Spectrum(x=x, y=peak(x))
        s._shift = self.shift
        s._factor = self.factor
        s._calib = self.calib
        s.offset = self.offset
        return s


class _PeaksDict(MutableMapping):
    """A special dictionary for saving peaks parameters.

    Keys must be integer. Value of the keys are always numbered according to the
        center `c` of each peak.


    minimum peak attributes:
    amp
    c
    fwhm (or fhwm1 and fwhm2)

    Elements must be dictionaries with the following items:
    amp
    c
    fwhm
    fwhm1
    fwhm2
    m
    m1
    m2
    offset

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
            if isinstance(kwargs['data'], dict) or isinstance(kwargs['data'], _Peak):
                self.append(kwargs['data'])
            elif isinstance(kwargs['data'], Iterable):
                for peak in kwargs['data']:
                    self.append(peak)
            else:
                # self.append(kwargs['data'])
                raise ValueError('data must be a list of peaks (dictionaries).')
        else:
            if len(args) == 1:
                if isinstance(args[0], dict) or isinstance(args[0], _Peak):
                    self.append(args[0])
                elif isinstance(args[0], Iterable):
                    for peak in args[0]:
                        self.append(peak)
                else:
                    # self.append(kwargs['data'])
                    raise ValueError('data must be a list of peaks (dictionaries).')
            elif len(args) > 1:
                for peak in args:
                    self.append(peak)
        # self.update(dict(*args, **kwargs))  # use the free update to set keys


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
        # return str(self.store).replace('}, ', '\n ')

    def __repr__(self):
        return str({i:val for i, val in enumerate(self.store)})[1:-1].replace('}, ', '}\n')
        # return str(self.store).replace('}, ', '\n ')

    def __getitem__(self, key):
        return self.store[self._check_key(key)]

    def __setitem__(self, key, value):
        # if isinstance(key, int) == False:
        #     raise ValueError(f'Dict keys must be integers (Error key: {key}).')
        if isinstance(value, _Peak):
            self.store[self._check_key(key)] = value
        elif isinstance(value, dict):
            self.store[self._check_key(key)] = _Peak(**value)
        else:
            # self.store[self._check_key(key)] = value
            raise ValueError('value must be a dict or a peak object')
        self._fix_order()
        #     if all(item in list(value.keys()) for item in ['amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2']):
        #         raise ValueError(f"Invalid keys for peak {key}). The valid keys are 'amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2'.")
        #     if 'amp' not in value:
        #         raise ValueError(f'Cannot find "amp" for peak: {key}.')
        #     else:
        #         if value['amp'] <= 0:
        #             raise ValueError('amp must be a number higher than 0.')
        #     if 'fwhm' not in value:
        #         if 'fwhm1' not in value and 'fwhm2' not in value:
        #             raise ValueError(f'Cannot find "fwhm" for peak: {key}.')
        #         elif 'fwhm1' in value and 'fwhm2' in value:
        #             if value['fwhm1'] <= 0:
        #                 raise ValueError('fwhm1 must be a number higher than 0.')
        #             if value['fwhm2'] <= 0:
        #                 raise ValueError('fwhm2 must be a number higher than 0.')
        #             value['fwhm'] = value['fwhm1']+value['fwhm2']
        #         else:
        #             raise ValueError(f'Cannot find "fwhm1" or "fwhm2" for peak: {key}.')
        #     if value['fwhm'] <= 0:
        #         raise ValueError('fwhm must be a number higher than 0.')
        #     if 'm' not in value:
        #         if 'm1' in value and 'm2' in value:
        #             if value['m1'] > 1 or value['m1']<0:
        #                 raise ValueError('m1 must be a number between 0 and 1.')
        #             if value['m2'] > 1 or value['m2']<0:
        #                 raise ValueError('m2 must be a number between 0 and 1.')
        #             if (value['m1']+value['m2'])/2 > 1 or (value['m1']+value['m2'])/2 < 0:
        #                 raise ValueError('m = m1+m2 must be a number between 0 and 1.')
        #             value['m'] = (value['m1']+value['m2'])/2
        #         elif 'm1' in value or 'm2' in value:
        #             raise ValueError(f'Cannot find "m1" or "m2" for peak: {key}.')
        #     else:
        #         if value['m'] > 1 or value['m']<0:
        #             raise ValueError('m must be a number between 0 and 1.')
        #     if 'c' not in value:
        #         raise ValueError(f'Cannot find "c" for peak: {key}.')
        # else:
        #     raise ValueError('input must be a dictionary with keys `amp`, `fwhm`, and `c`. Optional keys are `fwhm1`, `fwhm2`, `m`, `m1`, and `m2`. Input is not a dictionary.')
        #
        # self.store[self._check_key(key)] = value
        # self._fix_order()

    def __delitem__(self, key):
        del self.store[self._check_key(key)]
        self._fix_order()

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)


    def _check_key(self, key):
        if key > len(self)-1 or key < -len(self):
            raise KeyError('key out of range of defined peaks.\n')
        # if key < 0:
        #     key = self.__len__() + key
        # if key < 0:
        #     raise KeyError('key out of range of peaks.')
        return key

    def _fix_order(self):
        """Returns another PeakDict where peak number is numbered according to its position."""
        if self == []:
            return
        else:
            l = len(self)
            c = [0]*l
            for i, peak in enumerate(self):
                c[i] = peak['c']
            if sum(1 for test in np.diff(c) if test < 0) > 0:
                self.store = sort(c, self.store)
        # if self == {}:
        #     return
        # else:
        #     l = len(self)
        #     c = [0]*l
        #     for i, key in enumerate(self):
        #         c[i] = self[key]['c']
        #     if sum(1 for test in np.diff(c) if test < 0) > 0 or any([x for x in np.diff(list(self.keys())) if x!=1]) or 0 not in self.keys():
        #         i_ordered = sort(c, [i for i in list(self.keys())])
        #         self.store = _PeaksDict({i:self[i_ordered[i]] for i in range(l)})

    def append(self, value):
        if isinstance(value, _Peak):
            self.store.append(value)
        elif isinstance(value, dict):
            # print(value)
            peak = _Peak(**value)
            self.store.append(peak)
        else:
            # self.store[len(self)-1] = _Peak(**value)
            raise ValueError('value must be a dict or a peak object')
        self._fix_order()
        # self.store.append({})
        # self.__setitem__(len(self)-1, value)
        # self._fix_order()

    def remove(self, key):
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
        """Calibrate data.

        Args:
            value (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].

        Returns:
            None
        """
        for peak in self.store:
            peak.calib = value
        self._calib = value

    def set_shift(self, value):
        """Shift data.
        """
        for peak in self.store:
            peak.shift = value
        self._shift = value

    def set_offset(self, value):
        for peak in self.store:
            peak.offset = value
        self._offset = value

    def set_factor(self, value):
        for peak in self.store:
            peak.factor = value
        self._factor = value



    def add_near(self, key):
        """add another peak 1/4 of the fwhm away from the the peak.
        Can also accept a list.

        List can have multiple keys from the same peak for multiple sholders.
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
                        # temp.append({'amp':self.store[k]['amp'], 'c':self.store[k]['c']+self.store[k]['fwhm']*(n+1)/4, 'fwhm':self.store[k]['fwhm']})
                    # self.store[k]['amp'] = self.store[k]['amp']/(count+1)
            for value in temp:
                self.append(value)
        else:
            temp = copy.deepcopy(self.store[self._check_key(key)])
            self.store[self._check_key(key)]['c'] -= self.store[self._check_key(key)]['fwhm']/4
            temp['c']+=temp['fwhm']/4
            self.append(temp)
            # self.store[key]['amp'] = self.store[key]['amp']/2

    def get_curve(self, key, x):
        return self.store[key].get_curve(x)
        # peak = _peak(**self.store[key])
        # s = Spectrum(x=x, y=peak(x))
        # s.offset = self.offset
        # return s


def _peak_function_creator(asymmetry, fixed_m, idx=0):
    if fixed_m == False and type(fixed_m) == bool:  # variable m
        if asymmetry:
            def function2fit(x, amp, c, w1, m1, w2, m2):
                f = np.heaviside(x-c, 0)*voigt_fwhm(x, amp, c, w1, m1) +\
                    np.heaviside(c-x, 0)*voigt_fwhm(x, amp, c, w2, m2)
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
                    np.heaviside(c-x, 0)*voigt_fwhm(x, A, c, w2, fixed_m)
                return f
            f_str = f'np.heaviside(x-c_{idx}, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w1_{idx}, {m}) + np.heaviside(c_{idx}-x, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w2_{idx}, {m})'
            args_str = f'amp_{idx}, c_{idx}, w1_{idx}, w2_{idx}'
        else:
            def function2fit(x, A, c, w):
                return voigt_fwhm(x, A, c, w, fixed_m)
            f_str = f'voigt_fwhm(x, amp_{idx}, c_{idx}, w_{idx}, {m})'
            args_str = f'amp_{idx}, c_{idx}, w_{idx}'
    return function2fit, f_str, args_str


def _peak(x, amp, c, fwhm, m=0, fwhm1=None, fwhm2=None, m1=None, m2=None,):
    """returns a function f(x) of a peak."""
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

def _axis_interpreter(axis):
    """Allows for flexibility when signaling axis direction.

    It follows numpy's convention. ``0`` is the vertical axis, while ``1`` is the horizontal one.

    axis = c, col, column, columns, hor, horizontal, x, ... will return 1
    axis = r, row, rows, v, ver, vertical, y, ... will return 0

    Returns:
        0 or 1
    """
    if (axis < 0 or axis > 1) and type(axis)!=str:
        raise ValueError("axis must be 0 ('y') or 1 ('x')")
    elif axis == 1:
        return 1
    elif type(axis)==str and (axis=='x' or axis.startswith('c') or axis.startswith('h')):
        return 1
    elif axis == 0:
        return 0
    elif type(axis)==str and (axis=='y' or axis.startswith('r') or axis.startswith('v')):
        return 0
    else:
        raise ValueError("axis must be 0 ('y') or 1 ('x')")

def _xaxis_interpreter(xaxis):
    if xaxis.startswith('c'):
        return 'centers'
    elif xaxis.startswith('p') or xaxis.startswith('b'):
        return 'bins'

def _bins_interpreter(*args, **kwargs):
    """Allows for flexible binning parameters.

    Keyword arguments (kwargs) cannot be mixed with positional arguments. If
    positional arguments are used, they are assumend to be the number of bins.
    Size of bins cannot be assigned via positional arguments.

    Args:
        shape (tuple): shape of data. MUST BE PASSED AS KEYWORD argument.
        factor_check (bool, optional): if True, it will assert that nbins or
            bins_size are a factor of shape (divisible). MUST BE PASSED AS
            KEYWORD argument.

    kwargs:
        nbins (int or tuple, optional): number of bins. If one value is given,
            this is used for both vertical (number of rows) and horizontal
            (number of columns) directions. If two values are given,
            they are used separetely for the number of rows and number of
            columns, respectively. If the number of pixels in the image
            cannot be divided by the selected number of bins, it will raise an error.
            ``bins_size`` is calculated. nbins overwrites all other arguments.
        nrows (int, optional): number of rows.
        ncolumns (int, optional): number of columns.
        bins_size (int or tuple, optional): size of the bins. If one value is given,
            this is used for both vertical (number of rows) and horizontal
            (number of columns) directions. If two values are given,
            they are used separetely for rows and columns size, respectively.
            If number of pixels cannot be divided by ``bins_size``, the value of
            ``bins_size`` will be recalculated to the closest possible value.
            ``nbins`` is calculated.
        rows_size (int or float, optional): size of the rows.
        columns_size (int or float, optional): size of the columns.

    Returns:
        bins, bins_size
    """
    shape         = kwargs['shape']
    factor_check  = kwargs['factor_check']
    del kwargs['factor_check']
    del kwargs['shape']

    # initial check
    if kwargs != {} and args != ():
        raise AttributeError('Cannot mix keyword arguments with positional arguments.')

    # initialization
    nbins        = None
    nrows        = None
    ncolumns     = None
    bins_size    = None
    rows_size    = None
    columns_size = None

    # keyword arguments
    if 'nbins' in kwargs:
        nbins = kwargs['nbins']
        if isinstance(nbins, Iterable):
            if len(nbins) == 1:
                nrows, ncolumns = nbins[0], nbins[0]
            else:
                nrows, ncolumns = nbins[0], nbins[1]
        else:
            nrows, ncolumns = nbins, nbins
    else:
        if 'nrows' in kwargs:
            nrows = kwargs['nrows']
        if 'ncolumns' in kwargs:
            ncolumns = kwargs['ncolumns']

    if 'bins_size' in kwargs:
        bins_size = kwargs['bins_size']
        if isinstance(bins_size, Iterable):
            if len(bins_size) == 1:
                rows_size, columns_size = bins_size[0], bins_size[0]
            else:
                rows_size, columns_size = bins_size[0], bins_size[1]
        else:
            rows_size, columns_size = bins_size, bins_size
    else:
        if 'rows_size' in kwargs:
            rows_size = kwargs['rows_size']
        if 'columns_size' in kwargs:
            columns_size = kwargs['columns_size']

    # positional arguments
    if len(args) == 1:
        if isinstance(args[0], Iterable):
            if len(args[0]) == 1:
                nrows, ncolumns = args[0][0], args[0][0]
            elif len(args[0]) > 1:
                nrows, ncolumns = args[0][0], args[0][1]
            else:
                raise AttributeError(f'Cannot understand the arguments passed to the function.\nArguments passed: {args[0]}.\nExamples:\n10\n(10)\n10, 5\n(10, 5)\n')
        else:
            nrows, ncolumns = args[0], args[0]
        nbins = np.array((nrows, ncolumns))
    elif len(args) == 2:
        if isinstance(args[0], Iterable) or isinstance(args[1], Iterable):
            raise AttributeError(f'Cannot understand the arguments passed to the function.\nArguments passed: {args[0]}.\nExamples:\n10\n(10)\n10, 5\n(10, 5)\n')
        else:
            nrows, ncolumns = args[0], args[1]
        nbins = np.array((nrows, ncolumns))
    elif len(args) > 2:
        raise AttributeError(f'Cannot understand the arguments passed to the function.\nArguments passed: {args[0]}.\nExamples:\n10\n(10)\n10, 5\n(10, 5)\n')


    # after this part, either nbins, nrows, and ncolumns are defined or bins_size, row_size, and column_size are defined

    # final
    if bins_size is None and nbins is None:
        raise AttributeError('Number of bins or bin size must be passed as an argument.')
    elif nbins is not None:
        if nrows is None:    nrows    = shape[0]
        if ncolumns is None: ncolumns = shape[1]
        if nrows <= 0 or ncolumns <= 0 or is_integer(nrows)==False or is_integer(ncolumns)==False:
            raise ValueError("Number of bins must be a positive integer.")
        else:
            if factor_check:
                assert shape[1] % ncolumns == 0, f"The {shape[1]} pixels in a row is not evenly divisible by {ncolumns}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[1])))}"
                assert shape[0] % nrows    == 0, f"The {shape[0]} pixels in a column is not evenly divisible by {nrows}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[0])))}"

            nbins = np.array((nrows, ncolumns))
            # bins_size = np.array((shape[0]//nbins[0], shape[1]//nbins[1]))
            bins_size = np.array((shape[0]/nbins[0], shape[1]/nbins[1]))

    else:
        if row_size is None:    row_size =    shape[1]
        if column_size is None: column_size = shape[0]
        if row_size <= 0 or column_size <= 0 or is_integer(row_size)==False or is_integer(column_size)==False:
            raise ValueError("Size of bins must be a positive integer.")
        else:
            if factor_check:
                assert shape[1] % column_size == 0, f"The {shape[1]} pixels in a row is not evenly divisible by {column_size}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[1])))}"
                assert shape[0] % row_size == 0,    f"The {shape[0]} pixels in a column is not evenly divisible by {row_size}\nPlease, pick one of the following numbers: {np.sort(list(factors(shape[0])))}"

            nbins = np.array((round(shape[0]/row_size), round(shape[1]/column_size)))
            # bins_size = np.array((shape[0]//_nbins[0], shape[1]//_nbins[1]))
            bins_size = np.array((shape[0]/_nbins[0], shape[1]/_nbins[1]))

    return nbins, bins_size


class Image(metaclass=_Meta):
    """Image object.

    Args:
        data (2D array, optional): Image.
        filepath (string or path object, optional): filename or file handle.
            If the filename ends in .gz, the file is automatically saved in
            compressed gzip format. This is overwriten by data.

    Attributes:
        data (2D array): This is where we store the Image.
        x, y (1D array): x and y axis values.
        vmin, vmax (number): Minimum value in data.
        shape (tuple): Shape of data (vertical size, horizontal size).
        histogram (brixs.Spectrum): Data intensity histogram.

        nbins (tuple): Number of bins (number of rows, number of columns).
        bins_size (tuple): Bins size (size of rows, size of columns).
        reduced (brixs.Image): Binned image.

        calculated_shifts (brixs.Spectrum): Calculated shifts.
        shifts_v, shifts_h (1D array): Shift values in the vertical and
            horizontal direction.

        spectrum_v, spectrum_h (brixs.Spectrum): Spectrum obtained by integrating
            pixels in the vertical and horizontal direction.
        columns, rows (brixs.Spectra): Spectra obtained from each pixel columns
            or row.

    Methods:
        save()
        load()
        plot()
        imshow()
        binning()
        calculate_histogram()
        calculate_spectrum()
        floor()
        calculate_shifts()
        set_shifts()
        fix_curvature()


    """
    _read_only = ['shape', 'vmin', 'vmax',
                  'reduced',
                  'calculated_shifts']

    def __init__(self, *args, **kwargs):
        # argument parsing
        data, filepath = self._sort_args(args, kwargs)

        # besic attr
        self._data = None
        self._vmin = None
        self._vmax = None
        self._x    = None
        self._y    = None
        self._shape = None

        # binning attr
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None
        # self._x_edges    = None
        # self._y_edges    = None
        # self._x_centers  = None
        # self._y_centers  = None

        # shifts
        self._calculated_shifts = None
        self._shifts_v = None
        self._shifts_h = None

        # set data
        if data is not None:
            self.data = copy.deepcopy(data)
        elif filepath is not None:
            self.load(filepath)

    @property
    def data(self):
        return copy.deepcopy(self._data)
    @data.setter
    def data(self, value):
        # basic attr
        self._data  = np.array(value)
        self._vmin  = min([min(x) for x in self.data])
        self._vmax  = max([max(x) for x in self.data])
        self._shape = (self.data.shape[0], self.data.shape[1])
        self._x     = np.arange(0, self.data.shape[1])
        self._y     = np.arange(0, self.data.shape[0])

        # binning attr
        # if self.nbins[0] > 0 and self.nbins[1] > 0:
        #     try:
        #         self.binning(nbins=self.nbins)
        #     except:
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None
        # self._x_edges    = None
        # self._y_edges    = None
        # self._x_centers  = None
        # self._y_centers  = None

        # shift attr
        self._calculated_shifts = None
        self._shifts_v = np.zeros(self.data.shape[1])
        self._shifts_h = np.zeros(self.data.shape[0])
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, value):
        if value is None:
            self._x = np.arange(0, self.data.shape[1])
        elif isinstance(value, Iterable):
            if len(value) == self.shape[1]:
                self._x = np.array(value)
            else:
                raise ValueError(f"Length of x must be the same as the number of pixels in the horizontal direction.\nNumber of pixels = {self.data.shape[1]}\nLength of the array = {len(value)}")
        else:
            raise ValueError(f"x must be None or an iterable (list, tuple, or 1D array)")
    @x.deleter
    def x(self):
        self._x  = np.arange(0, self.data.shape[1])

    @property
    def y(self):
        return self._y
    @y.setter
    def y(self, value):
        if value is None:
            self._y = np.arange(0, self.data.shape[0])
        elif isinstance(value, Iterable):
            if len(value) == self.shape[0]:
                self._y = np.array(value)
            else:
                raise ValueError(f"Length of y must be the same as the number of pixels in the vertical direction.\nNumber of pixels = {self.data.shape[0]}\nLength of the array = {len(value)}")
        else:
            raise ValueError(f"y must be None or an iterable (list, tuple, or 1D array)")
    @y.deleter
    def y(self):
        self._y  = np.arange(0, self.data.shape[0])

    @property
    def shifts_v(self):
        return copy.deepcopy(self._shifts_v)
    @shifts_v.setter
    def shifts_v(self, value):
        self.set_shifts(value=value, axis=0)
    @shifts_v.deleter
    def shifts_v(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shifts_h(self):
        return copy.deepcopy(self._shifts_h)
    @shifts_h.setter
    def shifts_h(self, value):
        self.set_shifts(value=value, axis=1)
    @shifts_h.deleter
    def shifts_h(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nbins(self):
        return self._nbins
    @nbins.setter
    def nbins(self, value):
        if type(value) != str:
            binning(self, nbins=value)
        elif value == 'guess':
            raise NotImplementedError('not implemented yet')
        else:
            raise ValueError("Not a valid option of `nbins`.\nValid options are: 'guess', a number, a tuple, or a list.")
    @nbins.deleter
    def nbins(self):
        raise AttributeError('Cannot delete object.')

    @property
    def bins_size(self):
        return self._bins_size
    @bins_size.setter
    def bins_size(self, value):
        binning(self, bins_size=value)
    @bins_size.deleter
    def bins_size(self):
        raise AttributeError('Cannot delete object.')

    @property
    def histogram(self):
        return self.calculate_histogram()
    @histogram.setter
    def histogram(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @histogram.deleter
    def histogram(self):
        raise AttributeError('Cannot delete object.')

    @property
    def columns(self):
        ss = Spectra(n=self.shape[1])
        for i in range(self.shape[1]):
            ss[i] = Spectrum(x=self.y, y=self.data[:, i])
        return ss
    @columns.setter
    def columns(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @columns.deleter
    def columns(self):
        raise AttributeError('Cannot delete object.')

    @property
    def rows(self):
        ss = Spectra(n=self.shape[0])
        for i in range(self.shape[0]):
            ss[i] = Spectrum(x=self.x, y=self.data[i, :])
        return ss
    @rows.setter
    def rows(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @rows.deleter
    def rows(self):
        raise AttributeError('Cannot delete object.')

    @property
    def spectrum_v(self):
        return self.calculate_spectrum(axis=0)
    @spectrum_v.setter
    def spectrum_v(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @spectrum_v.deleter
    def spectrum_v(self):
        raise AttributeError('Cannot delete object.')

    @property
    def spectrum_h(self):
        return self.calculate_spectrum(axis=1)
    @spectrum_h.setter
    def spectrum_h(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @spectrum_h.deleter
    def spectrum_h(self):
        raise AttributeError('Cannot delete object.')

    def _sort_args(self, args, kwargs):
        """checks initial arguments.

        Keyword arguments (kwargs) cannot be mixed with positional arguments.

        For positional arguments, if one data set is passed, it assumes it is
            `data`. If this one argument is of type string or Pathlib.Path, it
            assumes it is a filepath.

        Raises:
            AttributeError: if kwargs and args cannot be read.

        Returns:
            data, filepath
        """
        # initial check
        if kwargs != {} and args != ():
            raise AttributeError('Cannot mix keyword arguments with positional arguments. Keyword arguents are `data`, and `filepath`.')

        # initialization
        data = None
        filepath = None

        # keyword arguments
        if 'data' in kwargs:
            data = kwargs['data']
        if 'filepath' in kwargs:
            filepath = kwargs['filepath']

        # positional arguments
        if len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepath = args[0]
            elif isinstance(args[0], Iterable):
                data = args[0]
        elif len(args) > 2:
            raise AttributeError('brixs.Image() cannot figure out the data out of the arguments passed. Maybe use keyword arguments (data, filepath).')

        return data, filepath

    def __len__(self):
        if self._data is None:
            return 0
        else:
            return len(self._data)

    def save(self, filepath, only_data=False,  **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.

        If not specified, the following parameters are passed to `numpy.savetxt()`_:

        Args:
            fmt (string, or list, optional): A single format (like ``%10.5f``), or a
                sequence of formats. See numpy's documentation for more information.
                If not specified, best fmt is calculated based on the
                number of decimal places of the data. Specifing fmt makes the code
                runs a little faster (not much tough, but it might make a difference
                if saving a lot of files).
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (", ").
            newline (str, optional): String or character separating lines.
                Default is ``\n``.
            header (bool, optional): String that will be written at the beginning of the file.
                Note that, attributes are always saved at the beginning of the file.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """

        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([n_decimal_places(x) for x in flatten(self._data)])
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # save
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
            np.savetxt(Path(filepath), self._data, **kwargs)
        else:
            if 'header' not in kwargs:
                kwargs['header'] = ''
            else:
                kwargs['header'] += '\n'
            dict = get_attributes(self)
            kwargs['header'] += '==== Image attributes ===='  + '\n'
            for n in dict:
                if n not in ['_data', '_reduced', '_calculated_shifts', '_histogram', '_spectrum_h', '_spectrum_v',]:
                    if isinstance(dict[n], Iterable):
                        kwargs['header'] += f'{n}: {list(dict[n])}'  + '\n'
                    else:
                        kwargs['header'] += f'{n}: {dict[n]}'  + '\n'
            np.savetxt(Path(filepath), self._data, **kwargs)

    def load(self, filepath, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first decompressed.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``# ``. Attributes picked up
                from the header will be loaded too.

        Returns:
            None

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # read data
        data = np.genfromtxt(Path(filepath), **kwargs)
        self.data = data

        # read header
        header = load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
        attr_start = 0
        binning_flag = False  # check if data needs binning

        # find where attributes listing starts
        for i, line in enumerate(header):
            if '==== Image attributes ====' in line:
                attr_start = i
                break
            attr_start = -1

        # read attributes
        if attr_start != -1:
            for i, line in enumerate(header[attr_start+1:-1]):

                # extract name and value
                name = line[1:-1].split(':')[0].strip()
                value = eval(line[1:-1].split(':')[1].strip())


                ### DEALING WITH ATTRIBUTES THAT NEED TO RUN SOMETHING ###
                if name in ['_nbins']:
                    if value[0] > 0 and value[1] > 0:
                        binning_flag = True
                elif name in ['_bins_size']:
                    if value[0] > 0 and value[1] > 0:
                        binning_flag = True
                ### DEALING WITH OTHER ATTRIBUTES ###
                elif name not in ['_vmin', '_vmax', '_shape']:  # except these attrs
                    try:
                        setattr(self, name, value)
                    except Exception as e:
                        print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')

        ### RUN ###
        if binning_flag:
            try:
                self.binning(nbins=self.nbins, bins_size=self.bins_size)
            except:
                pass

    def plot(self, ax=None, colorbar=False, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.imshow()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
            **kwargs: kwargs are passed to `matplotlib.pyplot.imshow()`_.

        If not specified, the following parameters are passed to `matplotlib.pyplot.imshow()`_:

        Args:
            cmap: The Colormap instance. Default is 'jet'.
            aspect: The aspect ratio of the Axes. Default is 'auto'.
            interpolation: The interpolation method used. Default is 'none'.
            origin: Place the [0, 0] index. Default is 'lower' This is necessary
                to make image plots comparable to photon event plots.
            vmin: Minimum intensity that the colormap covers. The intensity histogram is
                calculated and vmin is set on the position of the maximum.
            vmax: Maximmum intensity that the colormap covers.  The intensity histogram is
                calculated and vmax is set to the value where the :del:` integral of the
                intensity histogram is 0.998 of the total integral after vmin.`
                intensity drops below 0.01 % of the maximum.
                **THIS MIGHT CHANGE IN THE FUTURE**.

        Returns:
            `matplotlib.image.AxesImage`_

        .. _matplotlib.pyplot.imshow(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.imshow.html
        .. _matplotlib.image.AxesImage: https://matplotlib.org/3.5.0/api/image_api.html#matplotlib.image.AxesImage
        """
        # initialization
        if ax is None:
            ax = plt

        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'
        if 'vmin' not in kwargs:
            kwargs['vmin'] = self.histogram.x[np.argmax(self.histogram.y)]
        if 'vmax' not in kwargs:
            # vamx is set when histogram drops below 0.01% of the max
            x2fit = self.histogram.x[np.argmax(self.histogram.y)+1:]
            y2fit = self.histogram.y[np.argmax(self.histogram.y)+1:]
            data2fit = np.array([[i, j] for i, j in zip(x2fit, y2fit) if j > max(y2fit)*0.0001])  # clean zeros
            kwargs['vmax'] = data2fit[-1, 0]

        # plot
        X, Y = np.meshgrid(self.x, self.y)
        pos = plt.pcolormesh(X, Y, self.data, **kwargs)

        # colorbar
        if colorbar:
            plt.colorbar(pos, aspect=50, extend='both')

        return pos

    def binning(self, *args, **kwargs):
        """Compute the 2D histogram of the data (binning of the data).

        Args:
            nbins (int or tuple, optional): number of bins. If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the number of rows and number of
                columns, respectively. If the number of pixels in the image
                cannot be divided by the selected number of bins, it will raise an error.
                ``bins_size`` is calculated.
            bins_size (int or tuple, optional): size of the bins. This overwrites
                the argument ``bins``. If one value is given,
                it is used for both x and y directions. If two values are given,
                they are used separetely for rows and columns size, respectively.
                If number of pixels cannot be divided by ``bins_size``, the value of
                ``bins_size`` will be recalculated to the closest possible value.
                ``nbins`` is calculated.

        Example:
            .. code-block:: python

                nbins = 10       # (10 rows, 10 columns)
                nbins = (10)     # (10 rows, 10 columns)
                nbins = (10, 5)  # (10 rows, 5 columns)

        Returns:
            None
        """
        kwargs['shape']        = self.shape
        kwargs['factor_check'] = True
        _nbins, _bins_size = _bins_interpreter(*args, **kwargs)
        reduced = Image(np.add.reduceat(np.add.reduceat(self._data, np.arange(0, self.shape[0], _bins_size[0]), axis=0), np.arange(0, self.shape[1], _bins_size[1]), axis=1))
        _y_edges   = np.arange(0, self.shape[1], _bins_size[1])
        _x_edges   = np.arange(0, self.shape[0], _bins_size[0])
        _y_centers = moving_average(self.y_edges, n=2)
        _x_centers = moving_average(self.x_edges, n=2)

        # saving
        self._nbins     = _nbins
        self._bins_size = _bins_size
        self._reduced   = reduced
        self.reduced._x = _x_centers
        self.reduced._y = _y_centers
        # self._y_edges   = np.arange(0, self.shape[1], _bins_size[1])
        # self._x_edges   = np.arange(0, self.shape[0], _bins_size[0])
        # self._y_centers = moving_average(self.y_edges, n=2)
        # self._x_centers = moving_average(self.x_edges, n=2)

    def calculate_histogram(self, **kwargs):
        """Compute the histogram of data. Wrapper for `numpy.histogram()`_.

        If not specified, the following parameters are passed to `numpy.histogram()`_:

        Args:
            bins (int or tuple, optional): number of bins. If not specified, it will
                be set to 1000. If 1000 is too high (maximum value of the histogram
                is less than 5 % of the total integrated intensity) bins will be
                reduced 10 by 10 until the criteria is satisfied.

        Returns
            brixs.Spectrum

        .. _numpy.histogram(): https://numpy.org/doc/stable/reference/generated/numpy.histogram.html
        """
        if 'bins' not in kwargs:
            kwargs['bins'] = 1000
            hist, bin_edges = np.histogram(flatten(self._data), **kwargs)
            while max(hist) < self.shape[0]*self.shape[1]*0.05:
                kwargs['bins'] -= 10
                hist, bin_edges = np.histogram(flatten(self._data), **kwargs)
                if kwargs['bins'] < 50:
                    break
        else:
            hist, bin_edges = np.histogram(flatten(self._data), **kwargs)
        x = moving_average(bin_edges, 2)
        return Spectrum(x=x, y=hist)

    def calculate_spectrum(self, axis=1):
        """Integrate data in one direction (sum columns or rows).

        Args:
            axis (int or string, optional): Axis along which elements are integrated.
                By default, data is integrated in the horizontal direction.

        Returns:
            :py:class:`Spectrum`.
        """
        axis = _axis_interpreter(axis)

        if axis == 0:
            return Spectrum(x=self.y, y=np.sum(self._data, axis=0))
        elif axis == 1:
            return Spectrum(x=self.x, y=np.sum(self._data, axis=1))

    def floor(self, x=0, y=0, n=30, nx=None, ny=None):
        """Set background intensity to zero.

        Args:
            x, y (int, optional): x and y position to average background intensity.
            n, nx, ny (int, optional): size of the pixel window around x, y.

        Returns:
            None
        """
        assert x >= 0 and is_integer(x) and x<self.shape[1], f'x must be a positive integer smaller than {self.shape[1]}.'
        assert y >= 0 and is_integer(y) and y<self.shape[0], f'y must be a positive integer smaller than {self.shape[0]}.'

        # sorting n
        if nx is None: nx = n
        if ny is None: ny = n

        # check if x falls inside the image
        if x-nx/2 >= 0 and x+nx/2 < self.shape[1]:
            x_start = x-nx/2
            x_stop  = x+nx/2
        elif x-nx/2 < 0:
            x_start = 0
            x_stop = x+nx/2-(x-nx/2)
        elif x+n/2 >= self.shape[1]:
            x_start = x-nx/2-(x+nx/2-self.shape[1])
            x_stop  = self.shape[1]-1
        else:
            raise ValueError('Averaging range falls outside of the image. Please, change x or n.')

        # check if y falls inside the image
        if y-ny/2 >= 0 and y+ny/2 < self.shape[0]:
            y_start = y-ny/2
            y_stop  = y+ny/2
        elif y-nx/2 < 0:
            y_start = 0
            y_stop = y+ny/2-(y-ny/2)
        elif y+n/2 >= self.shape[0]:
            y_start = y-ny/2-(y+ny/2-self.shape[0])
            y_stop  = self.shape[0]-1
        else:
            raise ValueError('Averaging range falls outside of the image. Please, change y or n.')

        self._data -= np.mean(self._data[int(x_start):int(x_stop), int(y_start):int(y_stop)])
        self._vmin = min([min(x) for x in self.data])
        self._vmax = max([max(x) for x in self.data])

    def calculate_shifts(self, ref=0, axis=0, mode='cross-correlation', ranges=None, verbose=False, idx=0, bypass=False, **kwargs):

        axis = _axis_interpreter(axis)
        if axis == 0:
            ss = self.columns
            ss.calculate_shifts(ref=ref, mode=mode, ranges=ranges, verbose=verbose, idx=idx, bypass=bypass, **kwargs)
            self._calculated_shifts = ss.calculated_shifts
            self._calculated_shifts.x = self.x
        elif axis == 1:
            ss = self.rows
            ss.calculate_shifts(ref=ref, mode=mode, ranges=ranges, verbose=verbose, idx=idx, bypass=bypass, **kwargs)
            self._calculated_shifts = ss.calculated_shifts
            self._calculated_shifts.x = self.x

    def set_shifts(self, value, axis=0):
        """Roll array of pixels along a given axis.

        Elements that roll beyond the edge are NOT re-introduced at the first.
        They are lost. This might change in the future.

        Args:
            value (int or tuple): The number of pixels by which the data are shifted.
                If a tuple, then it must be of the same size as the number of
                columns (``for axis=0``) or rows (``for axis=1``). If elements are
                not int, it will rounded to an integer value.
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical direction.

        Returns:
            None
        """
        axis = _axis_interpreter(axis)

        if axis == 0:
            if isinstance(value, Iterable):
                assert len(value) == self.shape[1], f'Number of values must be the same as the number of columns ({self.shape[1]})'
                value = [round(k) for k in value]
            else:
                value = [int(round(value))]*self.shape[1]

            # undo last shift
            for i, v in enumerate(self._shifts_v):
                if v != 0:
                    temp = Spectrum(self._data[:, i])
                    temp.shift_roll = -v
                    self._data[:, i] = temp.y
                    self._shifts_v[i] = 0

            # apply shift
            for i, v in enumerate(value):
                if v != 0:
                    temp = Spectrum(self._data[:, i])
                    temp.shift_roll = v
                    self._data[:, i] = temp.y
                    self._shifts_v[i] = v
        elif axis == 1:
            if isinstance(value, Iterable):
                assert len(value) == self.shape[0], f'Number of values must be the same as the number of rows ({self.shape[0]})'
                value = [round(v) for v in value]
            else:
                value = [int(round(value))]*self.shape[0]

            # undo last shift
            for i, h in enumerate(self._shifts_h):
                if h != 0:
                    temp = Spectrum(self._data[i, :])
                    temp.shift_roll = -h
                    self._data[i, :] = temp.y
                    self._shifts_h[i] = 0

            # apply shift
            for i, h in enumerate(value):
                if h != 0:
                    temp = Spectrum(self._data[i, :])
                    temp.shift_roll = h
                    self._data[i, :] = temp.y
                    self._shifts_h[i] = h

        if self.reduced is not None:
            self.binning(nbins=self.nbins)

    def fix_curvature(self, deg=2, axis=0):
        """Fix curvature.

        Args:
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical direction.

        Returns:
            None
        """
        axis = _axis_interpreter(axis)

        assert self.reduced is not None, 'Image was not binned yet.\nPlease, use Image.binning()'

        # calculate shifts
        temp = Image(data=self.reduced.data)
        temp._x = self.x
        temp._y = self.y
        temp.floor()
        temp.calculate_shifts(axis=axis)

        # set shift
        self.set_shifts(temp.calculated_shifts.y, axis=axis)
        self._calculated_shifts = temp.calculated_shifts


class PhotonEvents(metaclass=_Meta):
    """Photon events object.

    Args:
        data (list or array): two (x, y) or three (x, y, intensity) list or array with
            photon events.
        filepath (string or path object, optional): filename or file handle.
            If the filename ends in .gz, the file is automatically saved in
            compressed gzip format. This is overwriten by data.
        shape (tuple, optional): Shape of data (maximum vertical size, maximum horizontal size).
        x, y, I (list or array): data columns.

    Attributes:
        data (2D array): This is where we store the Image.
        shape (tuple): Shape of data (vertical size, horizontal size).
        x, y, I (1D array): x, y, and intensity arrays.

        nbins (tuple): Number of bins (number of rows, number of columns).
        bins_size (tuple): Bins size (size of rows, size of columns).
        reduced (brixs.Image): Binned image.

        calculated_shifts (brixs.Spectrum): Calculated shifts.
        p (list): Polynomial coefficients of the fitted shift values.
        f (function): Function f(x) of the fitted shift values.
        shifts (brixs.Spectrum): Curve of the fitted shift values

    Methods:
        save()
        load()
        plot()
        binning()
        calculate_histogram()
        calculate_spectrum()
        floor()
        calculate_shifts()
        set_shifts()
        fix_curvature()

    """

    _read_only = ['reduced', 'calculated_shifts', 'p', 'f', 'shifts']

    def __init__(self, *args, **kwargs):
        # argument parsing
        data, filepath, shape = self._sort_args(args, kwargs)

        # basic attr
        self._data  = None
        self._shape = [None, None]

        # binning attr
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None

        # shifts
        self._calculated_shifts = None
        self._p      = None
        self._f      = None
        self._shifts = None

        # set data
        if data is not None:
            self.data = copy.deepcopy(data)
        elif filepath is not None:
            self.load(filepath)
        if shape != (None, None):
            self.shape = shape

    @property
    def data(self):
        return copy.deepcopy(self._data)
    @data.setter
    def data(self, value):
        # basic attr
        if value.shape[1] == 3:
            self._data = np.array([event for event in np.array(value) if not event[2]<=0])
        elif value.shape[1] == 2:
            data = np.array(value)
            self._data = np.c_[data, np.ones(data.shape[0])]
        else:
            raise ValueError("Data must have 2 or 3 columns.")

        # set shape
        if self.shape[0] is None:
            self._shape[0] = max(self.data[:, 1])
        elif max(self.data[:, 1]) > self._shape[0]:
            self._shape[0] = max(self.data[:, 1])
        if self.shape[1] is None:
            self._shape[0] = max(self.data[:, 0])
        elif max(self.data[:, 0]) > self._shape[1]:
            self._shape[1] = max(self.data[:, 0])

        # binning attr
        self._nbins      = np.array((-1, -1))
        self._bins_size  = np.array((-1, -1))
        self._reduced    = None

        # shift attr
        self._calculated_shifts = None
        self._p      = None
        self._f      = None
        self._shifts = None
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shape(self):
        return copy.deepcopy(self._shape)
    @shape.setter
    def shape(self, value):
        if isinstance(value, Iterable):
            if len(value) == 2:
                if value[0] is None:
                    value[0] = max(self.data[:, 1])
                if value[1] is None:
                    value[1] = max(self.data[:, 0])
                if value[0] < 0 or value[1] < 0:
                    raise ValueError('Shape cannot be negative.')
                self._shape = value
            else:
                raise ValueError(f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}')
        else:
            raise ValueError(f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}')
    @shape.deleter
    def shape(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x(self):
        return copy.deepcopy(self._data[:, 0])
    @x.setter
    def x(self, value):
        raise AttributeError('Cannot set object. This might change in the future.')
    @x.deleter
    def x(self):
        raise AttributeError('Cannot delete object.')

    @property
    def y(self):
        return copy.deepcopy(self._data[:, 1])
    @y.setter
    def y(self, value):
        raise AttributeError('Cannot set object. This might change in the future.')
    @y.deleter
    def y(self):
        raise AttributeError('Cannot delete object.')

    @property
    def I(self):
        return copy.deepcopy(self._data[:, 2])
    @I.setter
    def I(self, value):
        raise AttributeError('Cannot set object. This might change in the future.')
    @I.deleter
    def I(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nbins(self):
        return self._nbins
    @nbins.setter
    def nbins(self, value):
        if type(value) != str:
            self.binning(nbins=value)
        # elif value == 'guess':
        #     guess_bins = self.guess_bins()
        #     self.bins = guess_bins
        else:
            raise ValueError("Not a valid option of `nbins`.\nValid options are: a number, a tuple, or a list.")
    @nbins.deleter
    def nbins(self):
        raise AttributeError('Cannot delete object.')

    @property
    def bins_size(self):
        return self._bins_size
    @bins_size.setter
    def bins_size(self, value):
        self.binning(bins_size=value)
    @bins_size.deleter
    def bins_size(self):
        raise AttributeError('Cannot delete object.')

    def _sort_args(self, args, kwargs):
        """checks initial arguments.

        Keyword arguments (kwargs) cannot be mixed with positional arguments.

        For positional arguments, if one data set is passed, it assumes it is
            `data`. If this one argument is of type string or Pathlib.Path, it
            assumes it is a filepath.

        Raises:
            AttributeError: if kwargs and args cannot be read.

        Returns:
            data, filepath, shape
        """
        # initial check
        if kwargs != {} and args != ():
            raise AttributeError('Cannot mix keyword arguments with positional arguments. Keyword arguents are `x`, `y`, `I`, `shape`, `data`, and `filepath`.')

        # initialization
        data     = None
        filepath = None
        x        = None
        y        = None
        I        = None
        shape    = (None, None)
        error = 'brixs.PhotonEvents() cannot figure out the data out of the arguments passed.\nMaybe use keyword arguments.\nValid arguments: x, y, I, shape, data, filepath.'

        # keyword arguments
        if 'data' in kwargs:
            data = kwargs['data']
        if 'x' in kwargs:
            x = kwargs['x']
        if 'y' in kwargs:
            y = kwargs['y']
        if 'I' in kwargs:
            I = kwargs['I']
        if 'shape' in kwargs:
            shape = kwargs['shape']
        if 'data' in kwargs:
            data = kwargs['data']
        if 'filepath' in kwargs:
            filepath = kwargs['filepath']

        if x is not None and y is not None and data is None:
            assert len(x) == len(y), f'x and y must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}.'
            data = np.c_[np.array(x), np.array(y)]
            if I is not None:
                assert len(x) == len(I), f'x, y and I must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}\nLength of I = {len(I)}.'
                data = np.c_[data, np.array(I)]


        # positional arguments
        if len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepath = Path(args[0])
            else:
                data = args[0]
        elif len(args) == 2:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepath = Path(args[0])
                shape = list(args[1])
                assert len(shape) == 2, f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}'
            else:
                if len(args[0]) == len(args[1]):
                    x = args[0]
                    y = args[1]
                    assert len(x) == len(y), f'x and y must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}.'
                    data = np.c_[np.array(x), np.array(y)]
                else:
                    data = args[0]
                    shape = list(args[1])
                    assert len(shape) == 2, f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}'
        elif len(args) == 3:
            if len(args[0]) == len(args[1]):
                x = args[0]
                y = args[1]
                assert len(x) == len(y), f'x and y must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}.'
                data = np.c_[np.array(x), np.array(y)]
                if len(args[0]) == len(args[2]):
                    I = args[2]
                    data = np.c_[data, np.array(I)]
                else:
                    shape = args[2]
                    assert len(shape) == 2, f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}'
        elif len(args) == 4:
            x = args[0]
            y = args[1]
            assert len(x) == len(y), f'x and y must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}.'
            data = np.c_[np.array(x), np.array(y)]
            I = args[2]
            assert len(x) == len(I), f'x, y and I must have the same length.\nLength of x = {len(x)}\nLength of y = {len(y)}\nLength of I = {len(I)}.'
            data = np.c_[data, np.array(I)]
            shape = args[3]
            assert len(shape) == 2, f'Shape must be a list or tuple (vertical size, horizontal size).\nInvalid shape: {value}'
        elif len(args) > 4:
            raise AttributeError(error)

        return data, filepath, shape

    def __len__(self):
        if self._data is None:
            return 0
        else:
            return len(self._data[:, 0])

    def save(self, filepath, only_data=False,  **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.

        If not specified, the following parameters are passed to `numpy.savetxt()`_:

        Args:
            fmt (string, or list, optional): A single format (like ``%10.5f``), or a
                sequence of formats. See numpy's documentation for more information.
                If not specified, best fmt is calculated based on the
                number of decimal places of the data. Specifing fmt makes the code
                runs a little faster (not much tough, but it might make a difference
                if saving a lot of files).
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (", ").
            newline (str, optional): String or character separating lines.
                Default is ``\n``.
            header (bool, optional): String that will be written at the beginning of the file.
                Note that, attributes are always saved at the beginning of the file.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """

        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([n_decimal_places(x) for x in flatten(self._data)])
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # save
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
            np.savetxt(Path(filepath), self._data, **kwargs)
        else:
            if 'header' not in kwargs:
                kwargs['header'] = ''
            else:
                kwargs['header'] += '\n'
            dict = get_attributes(self)
            kwargs['header'] += '==== PhotonEvents attributes ===='  + '\n'
            for n in dict:
                if n not in ['_data', '_x', '_y', '_I', '_reduced', '_calculated_shifts', '_f', '_shifts']:
                    if isinstance(dict[n], Iterable):
                        kwargs['header'] += f'{n}: {list(dict[n])}'  + '\n'
                    else:
                        kwargs['header'] += f'{n}: {dict[n]}'  + '\n'
            np.savetxt(Path(filepath), self._data, **kwargs)

    def load(self, filepath, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first decompressed.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``# ``. Attributes picked up
                from the header will be loaded too.

        Returns:
            None

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # read data
        data = np.genfromtxt(Path(filepath), **kwargs)
        self.data = data

        # read header
        header = load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
        attr_start = 0
        binning_flag = False  # check if data needs binning

        # find where attributes listing starts
        for i, line in enumerate(header):
            if '==== PhotonEvents attributes ====' in line:
                attr_start = i
                break
            attr_start = -1

        # read attributes
        if attr_start != -1:
            for i, line in enumerate(header[attr_start+1:-1]):

                # extract name and value
                name = line[1:-1].split(':')[0].strip()
                value = eval(line[1:-1].split(':')[1].strip())


                ### DEALING WITH ATTRIBUTES THAT NEED TO RUN SOMETHING ###
                if name in ['_nbins']:
                    if value[0] > 0 and value[1] > 0:
                        binning_flag = True
                elif name in ['_bins_size']:
                    if value[0] > 0 and value[1] > 0:
                        binning_flag = True
                ### DEALING WITH OTHER ATTRIBUTES ###
                elif name not in []:  # except these attrs
                    try:
                        setattr(self, name, value)
                    except Exception as e:
                        print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')

        ### RUN ###
        if binning_flag:
            try:
                self.binning(nbins=self.nbins, bins_size=self.bins_size)
            except:
                pass

    def plot(self, ax=None, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.scatter()`_.

        The limits of the plot is also set to the size of the detector (given in
            the attribute ``shape``.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.

        If not specified, the following parameters are passed to `matplotlib.pyplot.scatter()`_:

        Args:
            s: The marker size in points**2. Default is 0.1.

        Returns:
            `matplotlib.image.AxesImage`_

        .. _matplotlib.pyplot.scatter(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.scatter.html
        .. _matplotlib.image.AxesImage: https://matplotlib.org/3.5.0/api/image_api.html#matplotlib.image.AxesImage
        """
        if ax is None:
            ax = plt

        # kwargs
        if 's' not in kwargs:
            kwargs['s'] = 0.1

        # plot
        ax.scatter(self.data[:, 0], self.data[:, 1], **kwargs)

        # set limits
        ax.xlim(0, self.shape[1])
        ax.ylim(0, self.shape[0])

        return ax

    def binning(self, *args, **kwargs):
        """Compute the 2D histogram of the data (binning of the data).

        Args:
            nbins (int or tuple, optional): number of bins. If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the number of rows and number of
                columns, respectively. ``bins_size`` is calculated.
            bins_size (int or tuple, optional): size of the bins. This overwrites
                the argument ``bins``. If one value is given,
                it is used for both x and y directions. If two values are given,
                they are used separetely for rows and columns size, respectively.
                ``nbins`` is calculated.

        Example:
            .. code-block:: python

                nbins = 10       # (10 rows, 10 columns)
                nbins = (10)     # (10 rows, 10 columns)
                nbins = (10, 5)  # (10 rows, 5 columns)


        Returns:
            None
        """
        kwargs['shape']        = self.shape
        kwargs['factor_check'] = False
        _nbins, _bins_size = _bins_interpreter(*args, **kwargs)

        temp, _x_edges, _y_edges = np.histogram2d(self.data[:, 0],
                                                            self.data[:, 1],
                                                            bins=_nbins[::-1],
                                                             weights=self.data[:, 2],
                                                             range=((0, self.shape[1]), (0, self.shape[0]))
                                                            )
        self._reduced   = Image(temp.transpose())
        self.reduced._x = moving_average(_x_edges, n=2)
        self.reduced._y = moving_average(_y_edges, n=2)
        # self._x_centers = moving_average(self.x_edges, n=2)
        # self._y_centers = moving_average(self.y_edges, n=2)
        self._nbins     = _nbins
        self._bins_size = _bins_size

    def calculate_spectrum(self, nbins=None, bins_size=None, axis=1, xaxis='bins'):
        """Integrate data in one direction (sum columns or rows).

        Args:
            nbins (int, optional): number of bins. nbins overwrites bins_size. The
                nbins of the opposite axis is automatically set to 1.
            bins_size (int or tuple, optional): size of the bins.
            axis (int or string, optional): Axis along which elements are integrated.
                By default, data is integrated in the horizontal direction.

        Returns:
            :py:class:`Spectrum`.
        """
        axis = _axis_interpreter(axis)

        temp = PhotonEvents(data=self.data, shape=self.shape)
        temp.binning(nbins=nbins, bins_size=bins_size)
        s = temp.reduced.calculate_spectrum(axis=axis)
        return s

    def calculate_shifts(self, ref=0, axis=0, mode='cross-correlation', ranges=None, verbose=False, idx=0, bypass=False, **kwargs):
        axis = _axis_interpreter(axis)
        assert self.reduced is not None, 'Image was not binned yet.\nPlease, use Image.binning()'

        self.reduced.calculate_shifts(ref=0, axis=0, mode='cross-correlation', ranges=None, verbose=False, idx=0, bypass=False, **kwargs)
        self._calculated_shifts = self.reduced.calculated_shifts
        if axis == 0:
            self.reduced.calculated_shifts._y = self.reduced.calculated_shifts.y*self.bins_size[0]
        elif axis == 1:
            self.reduced.calculated_shifts._y = self.reduced.calculated_shifts.y*self.bins_size[1]

    def transform(self, f):
        """Changes the values of x, y based on a function.

        Args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        Example:
            f = lambda x, y: (x, y**2)

        Returns:
            None
        """
        self._data[:, 0], self._data[:, 1] = f(self.data[:, 0], self.data[:, 1])

    def set_shifts(self, value=None, p=None, f=None, axis=0):
        """
        value
        p
        f
        axis=0 or 1
        """
        axis = _axis_interpreter(axis)

        if value is not None:
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)
            assert len(value) == len(self), f'value must have the same length as the data\nLength of the data = {len(self)}.'

            if axis == 0:
                self._data[:, 1] = self.data[:, 1] + value
                self._shifts = Spectrum(x=self.data[:, 0], y=value)
            elif axis == 1:
                self._data[:, 0] = self.data[:, 0] + value
                self._shifts = Spectrum(x=self.data[:, 1], y=value)
        elif p is not None:
            if axis == 0:
                self._data[:, 1] = self.data[:, 1] + np.polyval(p, self.data[:, 0])
            elif axis == 1:
                self._data[:, 0] = self.data[:, 0] + np.polyval(p, self.data[:, 1])
            self._p = p
        elif f is not None:
            if axis == 0:
                self._data[:, 1] = self.data[:, 1] + f(self.data[:, 0])
            elif axis == 1:
                self._data[:, 0] = self.data[:, 0] + f(self.data[:, 1])
            self._f = f

    def fix_curvature(self, deg=2, axis=0):
        axis = _axis_interpreter(axis)

        assert self.reduced is not None, 'Image was not binned yet.\nPlease, use Image.binning()'

        self.reduced.floor()
        self.reduced.calculate_shifts()
        if axis == 0:
            self.reduced.calculated_shifts._y = self.reduced.calculated_shifts.y*self.bins_size[0]
        elif axis == 1:
            self.reduced.calculated_shifts._y = self.reduced.calculated_shifts.y*self.bins_size[1]

        p, f, sfit = self.reduced.calculated_shifts.polyfit(deg=deg)
        self._p = p
        self._f = f
        self._shifts = sfit

        # set shifts
        self.set_shifts(p=p)

    #
    # def guess_bins(self, bins_initial=(9, 100), bins_step=(1, 50), max_x_iter=3, max_iter=1000, mode='cross-correlation', ranges=None):
    #     """
    #     Args:
    #         max_x_iter (int, optional): every iteration the algorithm calculates
    #             the offsets for different x_bins. Sometimes, subsequant x_bins
    #             will give the same result (i.e., same number of equal offsets).
    #             The algorithm will keep trying max_x_iter different x_bins until
    #             finaly change to a different y_bins.
    #         max_iter (int, optional): maximum number of different binnigs to try.
    #     """
    #     # offset parameters
    #     ref = int(bins_initial[0]/2)
    #
    #     # inital iteration
    #     # temp = PhotonEvents(self.data, x_max=self.x_max, y_max=self.y_max)
    #     # temp.bins = bins_initial
    #     # temp.calculate_offsets(mode=mode, ranges=ranges, ref=ref)
    #     diff = [0, 0]
    #     diff_sum = sum(diff)
    #
    #     # bins_initial[0] -= 1
    #     x_counter = 0
    #     total_counter = 0
    #     bins = (bins_initial[0], bins_initial[1])
    #     while 0 in diff and total_counter < max_iter:
    #         temp = PhotonEvents(self.data, x_max=self.x_max, y_max=self.y_max)
    #         temp.bins = bins
    #         temp.calculate_offsets(mode=mode, ranges=ranges, ref=ref)
    #         diff = np.diff(temp.offsets)
    #
    #         if x_counter < max_x_iter-1:
    #             if diff_sum == sum(diff):
    #                 bins = (bins[0]+bins_step[0], bins[1])
    #                 x_counter += 1
    #         else:
    #             bins = (bins_initial[0], bins[1]+bins_step[1])
    #             x_counter = 0
    #
    #         diff_sum = sum(diff)
    #         total_counter += 1
    #     if total_counter >= max_iter:
    #         raise RuntimeError('cannot find solution.')
    #     else:
    #         return temp.bins
    #
    # def best_bins(self, y_bins=None, x_bins=None, deg=2, spectrum_bins=6000, **kwargs):
    #
    #     if y_bins is None:
    #         y_bins = [100, 150, 300, 400, 600, 1000, 1500, 2500, 5000, 6000, 8000, 12000, 15000]
    #
    #     if x_bins is None:
    #         x_bins = [3, 5, 7, 9, 12, 15, 20, 30, 60]
    #
    #     # bins_old = self.bins
    #
    #     fwhm = {y_bin: {x_bin:None for x_bin in x_bins} for y_bin in y_bins}
    #     best_bins = (0, 0)
    #     best = 10000
    #     print(f'Starting iteration...')
    #     for i, y_bin in enumerate(y_bins):
    #         print(f'({i}/{len(y_bins)-1}) Trying y_bin = {y_bin}.')
    #         for x_bin in x_bins:
    #             temp = PhotonEvents(self.data, x_max=self.x_max, y_max=self.y_max)
    #             temp.bins = (x_bin, y_bin)
    #             # print(self.bins)
    #             temp.calculate_offsets(ref=int(temp.bins[0]/2))
    #             temp.fit_offsets(deg=deg)
    #             temp.offsets_correction()
    #             s = temp.calculate_spectrum(bins=spectrum_bins)
    #             try:
    #                 s.guess_elastic_peak(**kwargs)
    #                 fwhm[y_bin][x_bin] = s.elastic_fwhm
    #                 if s.elastic_fwhm < best:
    #                     best_bins = (x_bin, y_bin)
    #                     best = s.elastic_fwhm
    #             except RuntimeError:
    #                 pass
    #
    #
    #     print(f'Done. Best bins is {best_bins}.')
    #
    #     return best_bins, best, fwhm
    #
    # def best_bins_for_spectrum(self, y_bins=None, **kwargs):
    #
    #     if y_bins is None:
    #         y_bins = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 15000, 20000]
    #
    #     fwhm = {y_bin: None for y_bin in y_bins}
    #     best = 10000
    #     print(f'Starting iteration...')
    #     for i, y_bin in enumerate(y_bins):
    #         print(f'({i}/{len(y_bins)-1}) Trying y_bin = {y_bin}')
    #         s = self.calculate_spectrum(bins=y_bin, mode='length')
    #         try:
    #             s.guess_elastic_peak(**kwargs)
    #             fwhm[y_bin] = s.elastic_fwhm
    #             if s.elastic_fwhm < best:
    #                 best_bins = (1, y_bin)
    #                 best = s.elastic_fwhm
    #         except RuntimeError:
    #             fwhm[y_bin] = -1
    #
    #
    #     return best_bins, best, fwhm


class Spectrum(metaclass=_Meta):
    """Creates a ``spectrum`` class type object to deal with (x, y) data types.

    The hierarchy for Keyword arguments is: 1) data, 2) y (and x), and finaly
        3) filepath. For example, if `data` and `filepath` are passed as
        arguments, `filepath` is ignored.

    Keyword arguments (kwargs) cannot be mixed with positional arguments.

    For positional arguments, if one data set is passed, it assumes it is
        `data`. If this one argument is of type string or Pathlib.Path, it
        assumes it is a filepath. If two data sets are passed, it will
        assume one is the x coordinates and the next one is the y coordinates.

    Args:
        data (list or array, optional): two column list (or array).
        x (list or array, optional): x values (1D list/array). Overwrites `data`.
        y (list or array, optional): y values (1D list/array). Overwrites `data`.

    Attributes:
        data (array): three column array (x, y, and intensity)
        x (array): vector with x-coordinate values
        y (array): vector with y-coordinate values
        step (number, read only): step size for x-coordinates.
        shift( number): shift value (value will be added to x-coordinates).
        calib (number): calibration value (x-coordinates will be multiplied by this value).
        offset( number): offset value (value will be added to y-coordinates).
        factor (number): ultiplicative factor (y-coordinates will be multiplied by this value).
        peaks (dict): each entry represents a peak.
        fit_data (dict): fit data for each peak.
        fit_popt (list): optimized parameters of the full fit (only for inspection).
        fit_func (function): function of the fitted peaks.
        fit_ranges (list): data range used to fit the peaks.

    Note, peaks can be copyed around wihtout need of deepcopy.
    """

    _read_only = ['step', 'monotonicity',  'R2'] # 'rolled_out',
    _non_removable = ['fit', 'residue', #'fit_data', #'fit_popt', #'fit_func', 'shift_mode',
                      'fit_ranges', 'guess']

    def __init__(self, *args, **kwargs):
        # argument parsing
        data, x, y, filepath = self._args_checker(args, kwargs)

        # sorting data
        if data is not None:
            self.data = data
        elif filepath is not None:
            self.load(filepath)
        elif y is not None:
            self._restart_attr()
            if x is None:
                self._x = np.arange(0, len(y))
                self._step = 1
                self._monotonicity = 'increasing'
            else:
                if len(x) == len(y):
                    self._x = np.array(x)
                    self._step = None
                    self._monotonicity = None
                else:
                    raise ValueError('x and y data are not the same length.')
            self._y = np.array(y)
            self._data = np.vstack((self.x, self.y)).transpose()
        else:
            self._data = None
            self._x = None
            self._y = None
            self._restart_attr()
    #
    # def __str__(self):
    #     return str(self.data)
    #
    # def __repr__(self):
    #     return str(self.data)

    def __len__(self):
        if self.x is None:
            return 0
        else:
            return len(self.x)


    def _args_checker(self, args, kwargs):
        """checks initial arguments.

        Keyword arguments (kwargs) cannot be mixed with positional arguments.

        For positional arguments, if one data set is passed, it assumes it is
            `data`. If this one argument is of type string or Pathlib.Path, it
            assumes it is a filepath. If two data sets are passed, it will
            assume one is the x coordinates and the next one is the y coordinates.

        Raises:
            AttributeError: if kwargs and args cannot be read.

        Returns:
            data, x, y, filepath
        """
        # initial check
        if kwargs != {} and args != ():
            raise AttributeError('cannot mix key word arguments with positional arguments. Key word arguents are `x`, `y`, `data`, and `filepath`.')

        # initialization
        data = None
        x    = None
        y    = None
        filepath = None

        # keyword arguments
        if 'data' in kwargs:
            data = kwargs['data']
        if 'y' in kwargs:
            if 'x' in kwargs:
                x = kwargs['x']
            else:
                x = None
            y = kwargs['y']
        elif 'x' in kwargs:
            raise AttributeError('cannot seem to find y data.')
        if 'filepath' in kwargs:
            filepath = kwargs['filepath']

        # positional arguments
        if len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepath = args[0]
            elif isinstance(args[0], Iterable):
                data = args[0]
        elif len(args) == 2:
            x = args[0]
            y = args[1]
        elif len(args) > 2:
            raise AttributeError('brixs.Spectrum() cannot figure out the data out of the arguments passed. Maybe use key word arguments (x, y, data, filepath).')

        return data, x, y, filepath

    def load(self, filepath, comments='#', delimiter=None):
        data = load_data(filepath, delimiter=delimiter, comment_flag=comments, force_array=True)
        self.data = data
        header = load_Comments(filepath, comment_flag=comments, stop_flag=comments)
        if 'file generated by brixs V0.1' in header[0]:
            if 'V0.1' in header[0]:
                if 'step' in header[1]:
                    self._step = eval(header[1].split(':')[-1])
                if 'calib' in header[2]:
                    self._calib = eval(header[2].split(':')[-1])
                if 'offset' in header[3]:
                    self._offset = eval(header[3].split(':')[-1])
                if 'factor' in header[4]:
                    self._factor = eval(header[4].split(':')[-1])
                if 'shift' in header[5]:
                    self._shift = eval(header[5].split(':')[-1])
                if 'shift_roll' in header[6]:
                    self._shift_roll = eval(header[6].split(':')[-1])
                if 'shift_interp' in header[7]:
                    self._shift_interp = eval(header[7].split(':')[-1])
                # if 'rolled_out' in header[8]:
                #     self._rolled_out = eval(header[8].split(':')[-1])


    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        try:
            if value.shape[1] != 2:
                raise ValueError('Data must have two columns (x, y).')
            elif value is None:
                raise ValueError('No data to load.')
            elif value.shape[1] == 1:   # one column data (list)
                self._x = np.arange(0, len(value))
                self._y = np.array(value)
                self._data = np.vstack((self.x, self.y)).transpose()
                self._step = 1
                self._monotonicity = 'increasing'
                self._restart_attr()
            else:
                self._data = np.array(value)
                self._x = self.data[:, 0]
                self._y = self.data[:, 1]
                self._step = None
                self._monotonicity = None
                self._restart_attr()
        except IndexError:  # one column data (list)
            self._x = np.arange(0, len(value))
            self._y = np.array(value)
            self._data = np.vstack((self.x, self.y)).transpose()
            self._step = 1
            self._monotonicity = 'increasing'
            self._restart_attr()
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, value):
        if len(value) != len(self.x):
            raise ValueError('length of x you are trying to set is not compatible with length of y.')
        else:
            self._x = np.array(value)
            self._data[:, 0] = np.array(value)
            self._step = None
            self._monotonicity = None
            self._restart_attr()
            if self.calib != 1:
                print('Calibration not 1. Applying calibration to new data.')
                f = lambda x, y: (x*self.calib, y)
                self.apply_correction(f=f)
    @x.deleter
    def x(self):
        raise AttributeError('Cannot delete object.')

    @property
    def y(self):
        return self._y
    @y.setter
    def y(self, value):
        if len(value) != len(self.y):
            raise ValueError('length of y you are trying to set is not compatible with length of x.')
        else:
            self._y = np.array(value)
            self._data[:, 1] = np.array(value)
            self._restart_attr()
    @y.deleter
    def y(self):
        raise AttributeError('Cannot delete object.')

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
        self.set_shift(value, mode='x')
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shifts(self):
        return {'hard': self.shift, 'interp': self.shift_interp, 'roll':self.shift_roll}
    @shifts.setter
    def shifts(self, value):
        raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shift(value, mode).")
    @shifts.deleter
    def shifts(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_roll(self):
        return self._shift_roll
    @shift_roll.setter
    def shift_roll(self, value):
        # raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shift(value, mode='roll').")
        self.set_shift(value, mode='roll')
    @shift_roll.deleter
    def shift_roll(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_interp(self):
        return self._shift_interp
    @shift_interp.setter
    def shift_interp(self, value):
        # raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shift(value, mode='interp').")
        self.set_shift(value, mode='interp')
    @shift_interp.deleter
    def shift_interp(self):
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
    def peaks(self):
        # return copy.deepcopy(self._peaks)
        return self._peaks
    @peaks.setter
    def peaks(self, value):
        if isinstance(value, _PeaksDict):
            self._peaks = copy.deepcopy(value)
        elif isinstance(value, dict):
            self._peaks = _PeaksDict(value)
        else:
            raise ValueError(f'Peaks must be a dictionary.')
        # value = copy.deepcopy(value)
        # if isinstance(value, dict):
        #     if value.keys() != []:
        #         for key in value.keys():
        #             if isinstance(key, int) == False:
        #                 raise ValueError(f'Dict keys must be integers (Error: {key}).')
        #             else:
        #                 if 'amp' not in value[key]:
        #                     raise ValueError(f'Cannot find "amp" for peak: {key}.')
        #                 if 'fwhm' not in value[key]:
        #                     raise ValueError(f'Cannot find "fwhm" for peak: {key}.')
        #                 if 'c' not in value[key]:
        #                     raise ValueError(f'Cannot find "c" for peak: {key}.')
        #     self._peaks = value
        # else:
        #     raise ValueError(f'Peaks must be a dictionary.')
    @peaks.deleter
    def peaks(self):
        raise AttributeError('Cannot delete object.')


    def set_x(self, x):
        """Set x-coordinates."""
        self.x = x

    def get_x(self):
        """Returns x-coordinates."""
        return self.x

    def set_y(self, y):
        """Set y-coordinates."""
        self.y = y

    def get_y(self):
        """Returns y-coordinates."""
        return self.y

    def get_deltay(self):
        return max(self.y) - min(self.y)

    def get_deltax(self):
        return max(self.x) - min(self.x)

    def _restart_attr(self):
        """Set relevant attributes to initial value."""
        self._step = None
        self._monotonicity = None

        self._calib = 1
        self._shift = 0
        self._shift_roll = 0
        self._shift_interp = 0
        self._offset = 0
        self._factor = 1

        self.fit = None
        self.residue = None
        self.guess0 = None
        self.fit_ranges   = None
        self._R2 = None

        self._peaks = _PeaksDict()

    def save(self, filepath, delimiter=', ', fmt='%.18e', header=True, comments=''):
        r"""Saves spectrum to a file.

        Args:
            filepath (string or pathlib.Path, optional): filepath to file.
            delimiter (str, optional): The string used to separate values.
                If whitespaces are used, consecutive whitespaces act as delimiter.
                Use ``\\t`` for tab. The default is comma (', ').
            header (bool, optional):if true, brixs.Spectrum parameters are
                saved at the beginning of each file.
            fmt (string, or list, optional): format for saving data.
                If string, the value is used for x- and y-coordinates. If tuple
                of strings, the first string is used for x-coordinates and the
                second for y-coordinates.

                    fmt = (%[flag]width[.precision]specifier)

                * `flag` can be:
                    1. '-' for left justify
                    2. '+', puts + or - in front of the numbers
                    3. '0' to Left pad the number with zeros instead of space (see width).

                * `width` is the minimum number of characters to be printed.

                * `precision` is the number of significant digits.

                * `specifier` is the type of notation. Tipically, either 'e' for
                scientific notation of 'f' for decimal floating point.

                * a common `fmt` strings is: '%.3f' for 3 decimal places.

                *  for more information see `np.savetxt <https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html?highlight=savetxt#numpy.savetxt>`_ documentation::

        Returns:
            None

        """
        if header:
            header = 'file generated by brixs V0.1' + '\n'
            header += f'step: {self.step}' + '\n'
            header += f'calib: {self.calib}' + '\n'
            header += f'offset: {self.offset}' + '\n'
            header += f'factor: {self.factor}' + '\n'
            header += f'shift: {self.shift}' + '\n'
            header += f'shift_roll: {self.shift_roll}' + '\n'
            header += f'shift_interp: {self.shift_interp}' + '\n'
            if comments != '':
                header += str(comments) + '\n'
            header += f'x{delimiter}y'
        save_data(self.data, filepath=Path(filepath), delimiter=delimiter, header=header, fmt=fmt)

    def check_step_x(self, max_error=0.1):
        """Checks vector uniformity of the x-coordinates.

            If the step between two data points is the same through out the
            x vector, it sets :py:attr:`step` with the value of the step size.

            Args:
                max_error (number, optional): percentage (of the x step) value of the max error.

            Returns:
                None

            Raises:
                ValueError: If x-ccordinates are not uniform.
        """
        if self.monotonicity is None:
            self.check_monotonicity()

        # check step uniformity
        d = np.diff(self.x)
        if (max(d) - min(d))*100/np.mean(np.diff(self.x)) > max_error:
            self._step = None
            raise ValueError(f"Step in the x-coordinate seems not to be uniform.")

        self._step = np.mean(d)

    def find_peaks(self, prominence=5, width=4, moving_average_window=8, ranges=None):
        """Find peaks using scipy.signal.find_peaks() function.

        Sets :py:attr:`peaks` attribute (sets it to None if no peak is found).

        Args:
            prominence (number, optional): minimum prominence of peaks. This
                parameter will be passed directly to the ``scipy.signal.find_peaks()``
                function. If ``None``, prominence is set as 5 % of the delta intensity
                value.
            width (number, optional): minimum number of data points defining a peak.
            moving_average_window (int, optional): window size for smoothing the
                data for finding peaks. Default is 4.
            ranges (list, optional): data range to find peaks. Pair of
                x-coordinate values or a list of pairs. Each pair represents the
                start and stop of a data range.

        Returns:
            None

        Raises:
            ValueError: if points or moving_average_window have an improper value.
        """

        if self.monotonicity is None:
            self.check_monotonicity()
        # if self.monotonicity is 'decreasing':
        #     raise ValueError('array must be monotonicaly increasing.\nTip: use Spectrum.fix_monotinicity().')


        # check points and moving_average_window
        if width < 1:
            raise ValueError('width must be 1 or higher.')
        if isinstance(width, int) == False:
            if width.is_integer() == False:
                raise ValueError('width must be an integer.')
        if moving_average_window < 1:
            raise ValueError('moving_average_window must be 1 or higher.')
        if isinstance(moving_average_window, int) == False:
            if moving_average_window.is_integer() == False:
                raise ValueError('moving_average_window must be an integer.')

        # step
        if self.step is None:
            self.check_step_x()

        # extracting data
        if ranges is None:
            ranges = [[min(self.x), max(self.x)]]
            x = self.x
            y = self.y
        else:
            ranges = self._check_ranges(ranges)
            x, y = extract(self.x, self.y, ranges=ranges)

        # data smoothing
        if moving_average_window > 1:
            y2 = moving_average(y, moving_average_window)
            x2 = moving_average(x, moving_average_window)
        else:
            x2 = x[:]
            y2 = y[:]

        # parameters
        if prominence is None:
            prominence = (max(y2)-min(y2))*0.05
        else:
            prominence = (max(y2)-min(y2))*prominence/100

        # find peaks
        try:
            peaks, d = find_peaks(y2, prominence=prominence, width=width)
            # print(d)
            # self._peaks = {i: {'amp': d['width_heights'][i], 'c': x2[peaks[i]],  'fwhm': abs(d['widths'][i]*self.step)} for i in range(len(peaks))}
            # self._peaks = {i: {'amp': d['prominences'][i]+(y2[d['right_bases'][i]]+y2[d['left_bases'][i]])/2, 'c': x2[peaks[i]],  'fwhm': abs(d['widths'][i]*self.step)} for i in range(len(peaks))}
            shift = self.shift_roll*self.step + self.shift_interp + self._shift
            self._peaks = _PeaksDict([{'amp': d['prominences'][i]+max([y2[d['right_bases'][i]], y2[d['left_bases'][i]]]), 'c': x2[peaks[i]],  'fwhm': abs(d['widths'][i]*self.step)} for i in range(len(peaks))],
                                      shift=shift, offset=self.offset, calib=self.calib, factor=self.factor)
        except IndexError:
            peaks = _PeaksDict({}, shift=self.shift, offset=self.offset, calib=self.calib, factor=self.factor)

    def fit_peaks(self, idx='all', asymmetry=False, fixed_m=0, ranges=None, offset=True, amp_max_error=3, fwhm_max_error=2, c_error=2):
        """

        Args:
            idx (number, list, or 'all'): idx of the peak to fit.
            ranges (list, optional) if None, ranges is automatically set based
                on the first and last peak to be fitted (given on `idx`). If the
                last (first) peak to be fitted (peaks in idx) is indeed the
                last (first) peak of the data, then the fit limit is set up to (from)
                the last (first) data point. If the spectrum has peaks not included in the
                fitting, the fit goes up to the mid point point
                between the last (or first) peak and next peak not included in
                the fitting. If the inclued and exclued peaks are too close
                (their fwhm overlap), them the fitted range goes up to half the
                fwhm of the last (or first) peak).Finally, if ranges is not None,
                only data inside ranges is considered in the fitting. Also, it
                raises an error if any peak in idx is outside of ranges.
            assymetry (bool or dict, optional): if True, fits each half of the
                with a different width.
            fixed_m (False, number, or dict, optional): m is the amount of lorentzian
                contribution for a peak. If False, m will be fitted for each peak.
                If a number (from 0 to 1), this will be used as the value of m
                (fixed m).
            offset (bool, optional): if True, a offset value will be fitted.
            amp_max_error (number, optional): must be >= 1. The maximum possible
                fitted amp value will be `amp_max_error*guessed_amp`
            fwhm_max_error (number, optional): must be >= 1. The maximum possible
                fitted fhwm value will be `fwhm_max_error*guessed_fwhm`.
            c_error (number, optional): must be > 0. The minimum possible fitted
                c value will be `guessed_x-(guessed_fwhm/2)*c_error`.The maximum
                possible fitted c value will be `guessed_x+(guessed_fwhm/2)*c_error`.


        Returns:
            None

        """
        self.check_monotonicity()

        # check if peaks is defined ============================================
        if self.peaks is None:
            raise ValueError('Spectrum.peaks is not defined. Run Spectrum.find_peaks().')
        if self.peaks == []:
            raise ValueError('Spectrum.peaks = {}. No peaks to fit.')

        # check errors =========================================================
        if amp_max_error < 1:
            raise ValueError('amp_max_error must be >= 1. The maximum possible fitted amp value will be `amp_max_error*guessed_amp`.')
        if fwhm_max_error < 1:
            raise ValueError('fwhm_max_error must be >= 1. The maximum possible fitted fhwm value will be `fwhm_max_error*guessed_fwhm`.')
        if c_error <= 0:
            raise ValueError('c_error must be > 0. \nThe minimum possible fitted c value will be `guessed_x-(guessed_fwhm/2)*c_error`.\nThe maximum possible fitted c value will be `guessed_x+(guessed_fwhm/2)*c_error`.')

        # fix idx ==============================================================
        if idx == 'all':
            idx = [k for k in range(len(self.peaks))]
        if isinstance(idx, Iterable) == False:
            idx = [idx, ]
        idx = np.unique([i if i>=0 else len(self.peaks)+i for i in idx])

        # fix ranges ===========================================================
        if ranges is None:    # try fitting assigned peaks (range is automatically set)
            ranges = [min(self.x),  max(self.x)]
            p_max = self.peaks[max(idx)]['c']+self.peaks[max(idx)]['fwhm']
            p_min = self.peaks[min(idx)]['c']-self.peaks[min(idx)]['fwhm']

            if max(idx)+1 < len(self.peaks):  # if last peak to fit isn't the last peak
                next_min = self.peaks[max(idx)+1]['c']-self.peaks[max(idx)+1]['fwhm']
                if p_max > next_min:
                    ranges[1] = self.peaks[max(idx)]['c'] + self.peaks[max(idx)]['fwhm']/2
                # elif p_max > next_min:
                else:
                    ranges[1] = p_max + (next_min-p_max)/2
                    # x_temp, y_temp = extract(self.x, self.y, (p_max, next_min))
                    # ranges[1] = x_temp[np.argmin(y_temp)]
            else:
                next_min = max(self.x)
                ranges[1] = max(self.x)

            if min(idx)-1 >= 0:
                previous_max = self.peaks[min(idx)-1]['c']-self.peaks[min(idx)-1]['fwhm']
                if p_min < previous_max:
                    ranges[0] = p_min - self.peaks[min(idx)]['fwhm']/2
                else:
                    x_temp, y_temp = extract(self.x, self.y, (previous_max, p_min))
                    ranges[0] = x_temp[np.argmin(y_temp)]
            else:
                previous_max = min(self.x)
                ranges[0] = min(self.x)

            ranges = (ranges, )
            x2fit, y2fit = extract(self.x, self.y, ranges=ranges)
        else:
            ranges = self._check_ranges(ranges)
            x2fit, y2fit = extract(self.x, self.y, ranges=ranges)

            # check if peaks are inside ranges
            for i in idx:
                flag = True
                for r in ranges:
                    if self.peaks[i]['c'] > r[0] and self.peaks[i]['c'] < r[1]:
                        flag = False
                if flag:
                    raise ValueError(f'peak {i} outside of range.')

        # check asymmetry ======================================================
        if isinstance(asymmetry, dict):
            for i in idx:
                if i not in asymmetry:
                    asymmetry[i] = False
                else:
                    if type(asymmetry[i]) != bool:
                        raise ValueError('asymmetry must be bool (True or False)')
            for i in asymmetry:
                if i not in idx:
                    raise ValueError(f'Peak {i} not in idx (idx = {idx}) (asymmetry = {asymmetry}).')
        else:
            if isinstance(asymmetry, Iterable):
                raise ValueError('asymmetry must be a number or dictionary.')
            asymmetry = {i: asymmetry for i in idx}

        # check fixed_m ========================================================
        if isinstance(fixed_m, dict):
            for i in idx:
                if i not in fixed_m:
                    fixed_m[i] = 0
                else:
                    if (isinstance(fixed_m[i], bool) and fixed_m[i] == True) or (fixed_m[i] < 0 or fixed_m[i] > 1):
                        raise ValueError('fixed_m must be a False or a number from 0 to 1.')
            for i in fixed_m:
                if i not in idx:
                    raise ValueError(f'Peak {i} not in idx (idx = {idx}) (fixed_m = {fixed_m}).')
        else:
            if isinstance(fixed_m, Iterable):
                raise ValueError('fixed_m must be a number from 0 to 1, False, or dictionary.')
            fixed_m = {i: fixed_m for i in idx}

        # fit function =====================================================
        f_temp = {}

        model_str = ''
        args_str = ''
        bounds_min = []
        bounds_max = []
        p0 = []

        for i in idx:
            f, f_str, a_str = _peak_function_creator(asymmetry=asymmetry[i], fixed_m=fixed_m[i], idx=i)
            if fixed_m[i] == False and type(fixed_m[i]) == bool:  # variable m
                if asymmetry[i]:
                    n_args = 6
                    p0 = np.append(p0,                 [self.peaks[i]['amp']-min(y2fit),                  self.peaks[i]['c'],                                self.peaks[i]['fwhm'],  0.5,          self.peaks[i]['fwhm'],     0.5])
                    bounds_min = np.append(bounds_min, [                  0,                self.peaks[i]['c']-self.peaks[i]['fwhm']/2*c_error,                        0,             0,                     0,                 0])
                    bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*amp_max_error, self.peaks[i]['c']+self.peaks[i]['fwhm']/2*c_error, self.peaks[i]['fwhm']*fwhm_max_error, 1, self.peaks[i]['fwhm']*fwhm_max_error, 1])
                else:
                    n_args = 4
                    p0 = np.append(p0,                 [self.peaks[i]['amp']-min(y2fit),                  self.peaks[i]['c'],                               self.peaks[i]['fwhm'],    0.5])
                    bounds_min = np.append(bounds_min, [                0,                  self.peaks[i]['c']-self.peaks[i]['fwhm']/2*c_error,                       0,               0])
                    bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*amp_max_error, self.peaks[i]['c']+self.peaks[i]['fwhm']/2*c_error,  self.peaks[i]['fwhm']*fwhm_max_error, 1])
            else:
                if asymmetry[i]:
                    n_args = 4
                    p0 = np.append(p0,                 [self.peaks[i]['amp']-min(y2fit),                   self.peaks[i]['c'],                           self.peaks[i]['fwhm']/2,                  self.peaks[i]['fwhm']/2])
                    bounds_min = np.append(bounds_min, [                0,                  self.peaks[i]['c']-self.peaks[i]['fwhm']/2*c_error,                     0,                                     0,             ])
                    bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*amp_max_error, self.peaks[i]['c']+self.peaks[i]['fwhm']/2*c_error, self.peaks[i]['fwhm']/2*fwhm_max_error, self.peaks[i]['fwhm']/2*fwhm_max_error])
                else:
                    n_args = 3
                    p0 = np.append(p0,                  [self.peaks[i]['amp']-min(y2fit),                  self.peaks[i]['c'],                           self.peaks[i]['fwhm']      ])
                    bounds_min = np.append(bounds_min, [                     0,             self.peaks[i]['c']-self.peaks[i]['fwhm']/2*c_error,                     0,              ])
                    bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*amp_max_error, self.peaks[i]['c']+self.peaks[i]['fwhm']/2*c_error, self.peaks[i]['fwhm']*fwhm_max_error])
            f_temp[i] = {'asymmetry':asymmetry[i], 'fixed_m':fixed_m[i], '_func':f, 'n_args':n_args}
            model_str += f_str + ' + '
            args_str += a_str + ', '

        if offset:
            model_str = f'lambda x, {args_str}offset: {model_str}offset'
            p0 = np.append(p0, min(y2fit))
            bounds_min = np.append(bounds_min, min(y2fit))
            bounds_max = np.append(bounds_max, max(y2fit))
        else:
            model_str = f'lambda x, {args_str[:-2]}: {model_str[:-3]}'
        model = eval(model_str)

        # fit peaks ============================================================
        # print(p0)
        # print(bounds_min)
        # print(bounds_max)
        # plt.plot(x2fit, y2fit)
        # plt.plot(x2fit, model(x2fit, *p0))
        self.guess = Spectrum(x=x2fit, y=model(x2fit, *p0))
        self.guess._offset = self.offset
        self.guess._factor = self.factor
        self.guess._calib = self.calib
        self.guess._shift = self.shift
        self.guess._shift_roll = self.shift_roll
        self.guess._shift_interp = self.shift_interp
        self.guess.peaks = copy.deepcopy(self.peaks)
        # print(x2fit)
        # print(y2fit)

        popt, pcov = curve_fit(model, x2fit, y2fit, p0=p0, bounds = (bounds_min, bounds_max))

        x_temp = np.linspace(self.x[0], self.x[-1], len(self.x)*2)
        self.fit = Spectrum(x=x_temp, y=model(x_temp, *popt))
        self.fit._offset = self.offset
        self.fit._factor = self.factor
        self.fit._calib = self.calib
        self.fit._shift = self.shift
        self.fit._shift_roll = self.shift_roll
        self.fit._shift_interp = self.shift_interp
        self.fit.peaks = _PeaksDict()

        if offset:
            self.residue = Spectrum(x=x2fit, y=y2fit-model(x2fit, *popt)+popt[-1])
        else:
            self.residue = Spectrum(x=x2fit, y=y2fit-model(x2fit, *popt))
        self.residue._offset = self.offset
        self.residue._factor = self.factor
        self.residue._calib = self.calib
        self.residue._shift = self.shift
        self.residue._shift_roll = self.shift_roll
        self.residue._shift_interp = self.shift_interp

        self._R2 =  1- (sum((y2fit-model(x2fit, *popt))**2)/sum((self.y-np.mean(self.y))**2))

        # separate peak parameters =============================================
        # x_temp = np.linspace(self.x[0], self.x[-1], len(x2fit)*3)
        stop = 0
        for i in idx:
            peak_temp = {}
            start = stop
            stop  = start + f_temp[i]['n_args']
            # f_temp[i]['popt'] = popt[start:stop]
            # f_temp[i]['pcov'] = pcov[start:stop]
            peak_temp['amp'] = popt[start]
            peak_temp['c'] = popt[start+1]
            if f_temp[i]['fixed_m'] == False and type(f_temp[i]['fixed_m']) == bool:  # variable m
                if f_temp[i]['asymmetry']:
                    peak_temp['fwhm1'] = popt[start+2]
                    peak_temp['fwhm2'] = popt[start+4]
                    peak_temp['m1'] = popt[start+3]
                    peak_temp['m2'] = popt[start+5]
                else:
                    peak_temp['fwhm'] = popt[start+2]
                    peak_temp['m'] = popt[start+3]
            else:
                if f_temp[i]['asymmetry']:
                    peak_temp['fwhm1'] = popt[start+2]
                    peak_temp['fwhm2'] = popt[start+3]
                else:
                    peak_temp['fwhm'] = popt[start+2]
                peak_temp['m'] = f_temp[i]['fixed_m']
            self.fit.peaks.append(peak_temp)
            self.fit.peaks.offset = popt[-1]
            # if fixed_m[i] != 0 and type(fixed_m[i]) != bool:
            #     print(fixed_m[i])
            #     print(f_temp[i]['fixed_m'])
            #     self.fit.peaks[i]['m'] = fixed_m[i]
            # if offset:
            #     def peak(*args):
            #         return lambda x: (f_temp[i]['_func'](x/self.calib - self.shift, *args) + popt[-1])*self.factor + self.offset
            #     f_temp[i]['func'] = peak(*f_temp[i]['popt'])
            # else:
            #     def peak(*args):
            #         return lambda x: (f_temp[i]['_func'](x/self.calib - self.shift, *args))*self.factor + self.offset
            #     f_temp[i]['func'] = peak(*f_temp[i]['popt'])
            # f_temp[i]['curve'] = Spectrum(x=x_temp, y=f_temp[i]['func'](x_temp))

        # works, but not usefull anymore
        # self.fit_data = f_temp
        # self.fit_popt = popt
        # self.fit_func = lambda x: (model(x/self.calib - self.shift, *popt))*self.factor + self.offset
        # self.fit_func = lambda x: (model(x, *popt))
        self.fit_ranges = ranges
        # self.fit = Spectrum(x=x_temp, y=self.fit_func(x_temp))
        # self.fit.calib = self.calib
        # self.fit.shift = self.shift

    def fit_peak(self, asymmetry=False, fixed_m=0, ranges=None, offset=True):
        self.check_monotonicity()

        # fix ranges ===========================================================
        if ranges is None:
            x2fit = self.x
            y2fit = self.y
        else:
            ranges = self._check_ranges(ranges)
            x2fit, y2fit = extract(self.x, self.y, ranges=ranges)

        # input parameters =============================================
        amp   = max(y2fit)-np.mean(y2fit)
        c     = x2fit[np.argmax(y2fit)]
        i2 = index([1 if x < amp/2+np.mean(y2fit) else 0 for x in y2fit[np.argmax(y2fit):]], 1)
        if i2 <= 0:
            print('trouble estimating width')
            width = (x2fit[np.argmax(y2fit)+4]-c)*2
        else:
            width = (x2fit[np.argmax(y2fit)+i2]-c)*2

        if fixed_m == False and type(fixed_m) == bool:
            if asymmetry:
                p0         = [ amp,     c,   width,  0.5, width, 0.5]
                bounds_min = [  0,   -np.inf,   0,    0,    0,    0]
                bounds_max = [np.inf, np.inf, np.inf, 1,  np.inf, 1]
            else:
                p0         = [ amp,     c,    width, 0.5]
                bounds_min = [  0,   -np.inf,   0,    0]
                bounds_max = [np.inf, np.inf, np.inf, 1]
        else:
            if asymmetry:
                p0         = [ amp,     c,     width, width]
                bounds_min = [  0,   -np.inf,    0,      0]
                bounds_max = [np.inf, np.inf, np.inf, np.inf]
            else:
                p0         = [ amp,     c,    width]
                bounds_min = [  0,   -np.inf,   0]
                bounds_max = [np.inf, np.inf, np.inf]

        # model ===================================
        # print(p0)
        f, f_str, a_str = _peak_function_creator(asymmetry=asymmetry, fixed_m=fixed_m, idx=0)
        model_str = f_str + ' + '
        args_str = a_str + ', '

        if offset:
            model_str = f'lambda x, {args_str}offset: {model_str}offset'
            p0 = np.append(p0, min(y2fit))
            bounds_min = np.append(bounds_min, -np.inf)
            bounds_max = np.append(bounds_max, np.inf)
        else:
            model_str = f'lambda x, {args_str[:-2]}: {model_str[:-3]}'
        model = eval(model_str)

        # fit =================================
        self.guess = Spectrum(x=x2fit, y=model(x2fit, *p0))
        self.guess._offset = self.offset
        self.guess._factor = self.factor
        self.guess._calib = self.calib
        self.guess._shift = self.shift
        self.guess._shift_roll = self.shift_roll
        self.guess._shift_interp = self.shift_interp

        try:
            popt, pcov = curve_fit(model, x2fit, y2fit, p0=p0, bounds = (bounds_min, bounds_max))
        except RuntimeError:
            raise RuntimeError('Optimal parameters not found: The maximum number of function evaluations is exceeded.\n IMPLEMENT A CHANGE IN p0 WHEN THIS HAPPENS.')
        # x_temp = np.linspace(self.x[0], self.x[-1], len(self.x)*2)

        self.fit = Spectrum(x=x2fit, y=model(x2fit, *popt))
        self.fit._offset = self.offset
        self.fit._factor = self.factor
        self.fit._calib = self.calib
        self.fit._shift = self.shift
        self.fit._shift_roll = self.shift_roll
        self.fit._shift_interp = self.shift_interp
        if fixed_m == False and type(fixed_m) == bool:
            if asymmetry:
                temp = {'amp':popt[0], 'c':popt[1], 'fwhm1':popt[2], 'm1':popt[3], 'fwhm2':popt[4], 'm2':popt[5]}
            else:
                temp = {'amp':popt[0], 'c':popt[1], 'fwhm':popt[2], 'm':popt[3]}

        else:
            if asymmetry:
                temp = {'amp':popt[0], 'c':popt[1], 'fwhm1':popt[2], 'fwhm2':popt[4]}
            else:
                temp = {'amp':popt[0], 'c':popt[1], 'fwhm':popt[2]}
        self.fit.peaks = _PeaksDict([temp], shift=self.shift, offset=self.offset, factor=self.factor, calib=self.calib)


        self.residue = Spectrum(x=x2fit, y=model(x2fit, *popt))
        self.residue._offset = self.offset
        self.residue._factor = self.factor
        self.residue._calib = self.calib
        self.residue._shift = self.shift
        self.residue._shift_roll = self.shift_roll
        self.residue._shift_interp = self.shift_interp

    def polyfit(self, deg=2):
        p = np.polyfit(self.x, self.y, deg=deg)
        f = lambda x: np.polyval(p, x)
        x = np.linspace(min(self.x), max(self.x), len(self.x)*20)
        return p, f, Spectrum(x=x, y=f(x))

    def check_monotonicity(self):
        if np.all(np.diff(self.x) > 0) == True:
            self._monotonicity = 'increasing'
        elif np.all(np.diff(self.x) < 0) == True:
            self._monotonicity = 'decreasing'
        else:
            self._monotonicity = None
            raise ValueError('x array is not monotonicaly. Use Spectrum.fix_monotinicity()')

    def fix_monotinicity(self, mode='increasing'):
        if mode != 'increasing' and mode != 'decreasing':
            raise ValueError('mode should be "decreasing" or "increasing".')
        if self.monotonicity is None:
            try:
                self.check_monotonicity()
            except ValueError:
                unqa, ID, counts = np.unique(self.x, return_inverse=True, return_counts=True)
                # self.x = unqa
                # self.y = np.bincount(ID, self.y)/counts
                # f = np.column_stack(( unqa , np.bincount(ID,self.y)/counts ))
                # print(f.shape)
                self.data = np.column_stack(( unqa , np.bincount(ID,self.y)/counts ))
                # print(self.x)
                self.check_monotonicity()
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

    def _multiplicative_x_fix(self, value):
        # peaks =================================
        if self.peaks is not None:
            self.peaks.calib = value

        # fit ==================================
        if self.fit is not None:
            self.guess.calib = value
            self.fit.calib = value
            self.residue.calib = value

        # # fit_ranges =============================
        # if self.fit_ranges is not None:
        #     self.fit_ranges = [(r[0]*value, r[1]*value) for r in self.fit_ranges]

    # def _additive_x_fix(self, value):
    #     # fit ==================================
    #     if self.fit is not None:
    #         self.guess.set_shift(value, mode=self.shift_mode)
    #         self.fit.set_shift(value, mode=self.shift_mode)
    #         self.residue.set_shift(value, mode=self.shift_mode)
    #
    #     # peaks =================================
    #     if self.peaks is not None:
    #         if self.shift_mode in roll:
    #             value = self.step*value
    #         self.peaks.shift = value

        # # fit_ranges =============================
        # if self.fit_ranges is not None:
        #     self.fit_ranges = [(r[0]+value, r[1]+value) for r in self.fit_ranges]

    def _additive_y_fix(self, value):
        # peaks =================================
        if self.peaks is not None:
            self.peaks.offset = value

        # fit ==================================
        if self.fit is not None:
            self.guess.offset = value
            self.fit.offset = value
            self.residue.offset = value

    def _multiplicative_y_fix(self, value):
        # peaks =================================
        if self.peaks is not None:
            self.peaks.factor = value
        # fit ==================================
        if self.fit is not None:
            self.guess.factor = value
            self.fit.factor = value
            self.residue.factor = value

    def _check_ranges(self, ranges):
        text = 'Ranges should be a pair (x_init1, x_final1) or a list of pairs like this: ((x_init1, x_final1), (x_init2, x_final2), ...)\nUse None to idicate the minimum or maximum x value of the data.'
        # check format
        if ranges is None:
            ranges = ((min(self.x), max(self.x)),)
        elif isinstance(ranges, Iterable) == True:
            if isinstance(ranges[0], Iterable) == False:
                    ranges = (ranges, )
        else:
            raise ValueError(text)
        # check pairs
        ranges = list(ranges)
        for i in range(len(ranges)):
            if len(ranges[i]) == 2:
                if None in ranges[i]:
                    ranges[i] = list(ranges[i])
                    if ranges[i][0] is None:
                         ranges[i][0] = min(self.x)
                    if ranges[i][1] is None:
                         ranges[i][1] = max(self.x)
                    # r = tuple(r)
            else:
                raise ValueError(f'Ranges pair {r} is not a valid pair.\n'+text)
        # print(ranges)
        return tuple(ranges)



    def set_calib(self, value):
        """Calibrate data.

        Args:
            value (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].

        Returns:
            None
        """
        if value == 0:
            raise ValueError('cannot set calib = 0.')
        if self.calib != value:
            if self.calib != 1:
                self._x = self.x*self.calib**-1
                self.data[:, 0] = self._x
                if self.step is not None:
                    self._step *= self.calib**-1
            if value != 1:
                self._x = self.x*value
                self.data[:, 0] = self._x
                if self.step is not None:
                    self._step *= value
            self._calib = value
        self._multiplicative_x_fix(value)

    def set_shift(self, value, mode):
        """Shift data.

        Args:
            value (float or int): shift value.
            mode (string, optional): If ``mode='x'`` or ``mode='hard'``, y is fully preserved
                while x is shifted. If ``mode='y'``, ``'interp'``, or ``'soft'``, x is preserved
                while y is interpolated with a shift. If ``mode='roll'`` (or rotate or r), x is also preserved
                and y elements are rolled along the array (``shift`` value must be an integer).
                The form of y-coordinates is fully preserved, but the edge of y-coordinates are lost.

        Returns:
            None

        Warning:
            It is always better to use ``mode='hard'`` or ``'roll'`` since the form of y is fully
            preserved (no interpolation). After applying a shift using the ``mode='interp'``,
            one can apply a
            'inverse' shift to retrieve the original data. The diference between the retrieved
            y data and the original data will give an ideia of the information loss
            caused by the interpolation. Also, multiple shifts using interp may cause
            built up of data loss. It's not drastic if data is suficiently well
            resolved, but it should be avoided.
        """

        if mode in roll and self.step is None:
            try:
                self.check_step_x()
            except ValueError:
                raise ValueError(f'Cannot shift data using mode = {mode}, because x-coordinates are not uniform.')

        if mode in roll:
            if is_integer(value) == False:
                raise ValueError('shift must be an integer for mode = `roll`.')

            if self.shift_roll != value:
                if self.shift_roll != 0:  # undo last shift
                    self._x, self._y = shifted(self.x, self.y, value=-self.shift_roll, mode='roll')
                    # if -self.shift < 0:
                    #     self._y[-self.shift_roll:] = self.rolled_out
                    # elif -self.shift > 0:
                    #     self._y[:-self.shift_roll] = self.rolled_out
                    # self._rolled_out = None
                    if self.peaks is not None:
                        self.peaks.shift += self.step*-self.shift_roll
                if value != 0:
                    # if value > 0:
                    #     self._rolled_out = self.y[-int(value):]
                    # elif value < 0:
                    #     self._rolled_out = self.y[:-int(value)]
                    self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                    if self.peaks is not None:
                        self.peaks.shift += self.step*value
                self._shift_roll = value
                self._data[:, 0] = self._x
                self._data[:, 1] = self._y

            if self.fit is not None:
                self.guess.set_shift(value, mode=mode)
                self.fit.set_shift(value, mode=mode)
                self.residue.set_shift(value, mode=mode)


        elif mode in hard:
            if self.shift != value:
                if self.shift != 0:
                    self._x, self._y = shifted(self.x, self.y, value=-self._shift, mode='x')
                    if self.peaks is not None:
                        self.peaks.shift += -self._shift
                if value != 0:
                    self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                    if self.peaks is not None:
                        self.peaks.shift += value
                self._shift = value
                self._data[:, 0] = self._x
                self._data[:, 1] = self._y
            if self.fit is not None:
                self.guess.set_shift(value, mode=mode)
                self.fit.set_shift(value, mode=mode)
                self.residue.set_shift(value, mode=mode)
        elif mode in soft:
            if self.shift_interp != value:
                if self.monotonicity is None:
                    self.check_monotonicity()
                if self.monotonicity != 'increasing':
                    raise ValueError('x array must be monotonicaly increasing.\nTip: use Spectrum.fix_monotinicity()')

                if self.shift_interp != 0:
                    self._x, self._y = shifted(self.x, self.y, value=-self.shift_interp, mode='interp')
                    if self.peaks is not None:
                        self.peaks.shift += -self.shift_interp
                if value != 0:
                    self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
                    if self.peaks is not None:
                        self.peaks.shift += value
                self._shift_interp = value
                self._data[:, 0] = self._x
                self._data[:, 1] = self._y
            if self.fit is not None:
                self.guess.set_shift(value, mode=mode)
                self.fit.set_shift(value, mode=mode)
                self.residue.set_shift(value, mode=mode)
        else:
            raise ValueError(f'Invalid mode. Valid options are `roll`, `x`, `interp`.')

        # if self.shift != value:
        #     if self.shift != 0:
        #         self._x, self._y = shifted(self.x, self.y, value=-self.shift, mode=self.shift_mode)
        #         if self.shift_mode in roll:
        #             # print(-self.shift)
        #             if -self.shift < 0:
        #                 self._y[int(-self.shift):] = self.rolled_out
        #             elif -self.shift > 0:
        #                 self._y[:int(-self.shift)] = self.rolled_out
        #             self._rolled_out = None
        #         # self._additive_x_fix(self.shift, mode=self.shift_mode)
        #     if value != 0:
        #         if mode in roll:
        #             if value > 0:
        #                 self._rolled_out = self.y[-int(value):]
        #             elif value < 0:
        #                 self._rolled_out = self.y[:-int(value)]
        #         self._x, self._y = shifted(self.x, self.y, value=value, mode=mode)
        #         # self._shift_mode = mode
        #     else:
        #         # self._shift_mode = None
        #     self._shift = value
        #     self._data[:, 0] = self._x
        #     self._data[:, 1] = self._y
        # self._additive_x_fix(value)

    def set_offset(self, value):
        if self.offset != value:
            if self.offset != 0:
                self._y = self.y - self.offset
                self.data[:, 1] = self._y
                # self._additive_y_fix(-self.offset)
            if value != 0:
                self._y = self.y + value
                self.data[:, 1] = self._y
            self._offset = value
        self._additive_y_fix(value)

    def set_factor(self, value):
        if value == 0:
            raise ValueError('cannot set factor = 0.')
        if self.factor != value:
            if self.factor != 1:
                self._y = self.y*self.factor**-1
                self.data[:, 1] = self._y
                # self._multiplicative_y_fix(self.factor**-1)
            if value != 1:
                self._y = self.y*value
                self.data[:, 1] = self._y
            self._factor = value
        self._multiplicative_y_fix(value)

    def apply_correction(self, f, verbose=True):
        """Changes the values of x, y based on a function.

        resets everything.

        Example:
            f = lambda x, y: (x, y**2)

        Args:
            f (function): function ``x, y = f(x, y)`` that takes as input the
                position of a photon event and returns its corrected values.

        Returns:
            None

        Warning:
            Spectra.calib is set to 1.
            Spectra.step is recalculated if possible.
            Fit information is lost.
        """
        # core ================================
        self._x, self._y = f(self.x, self.y)
        self._data[:, 0] = self._x
        self._data[:, 1] = self._y

        self._restart_attr()
        if self.step is not None:
            try:
                self.check_step_x()
            except ValueError:
                pass



    def interp(self, start=None, stop=None, num=None, step=None, x=None):
        """Interpolate data.

        Args:
            x (list or array, optional): The x-coordinates at which to
                evaluate the interpolated values. This overwrites all other arguments.
            start (number, optional): The starting value of the sequence. If None,
                the minium x value will be used.
            stop (number, optional): The end value of the sequence. If None,
                the maximum x value will be used.
            num (int, optional): Number of samples to generate.
            step (number, optional): Spacing between values. This overwrites ``num``.

        Returns:
            None

        Warning:
            Spectrum.interp() will change data for good. There is no way to
                recover the previous state of the array.
        """
        self.check_monotonicity()
        if self.monotonicity == 'decreasing':
            raise ValueError('monotonicity of data must be strictly increasing. Use Spectrum.fix_monotinicity().')

        if x is None:
            if start is None:
                start = min(self.x)
            if stop is None:
                stop = max(self.x)

            if step is not None:   # step overwrites num
                x = np.arange(start, stop, step=step)
            if num is not None:
                x = np.linspace(start, stop, num=num)
            else:
                temp, _ = extract(self.x, self.y, (start, stop))
                x = np.linspace(start, stop, num=len(temp))

        self._y = np.interp(x, self.x, self.y)
        self._x = x
        self._data = np.vstack((self._x, self._y)).transpose()

        # self._restart_attr()
        if self.step is not None:
            try:
                self.check_step_x()
            except ValueError:
                pass

    def compress(self, ranges):
        """Extract data range from full data."""
        ranges = self._check_ranges(ranges)
        self._x, self._y  = extract(self.x, self.y, ranges)
        self._data = np.vstack((self._x, self._y)).transpose()

        self._restart_attr()
        if self.step is not None:
            try:
                self.check_step_x()
            except ValueError:
                pass

    def crop(self, start=None, stop=None):
        """crop data."""
        if start is None and stop is None:
            return
        if start is None:
            start = min(self.x)
        if stop is None:
            stop = max(self.x)
        self._x, self._y  = extract(self.x, self.y, (start, stop))
        self._data = np.vstack((self._x, self._y)).transpose()

    def floor(self, value=None, n=20, ranges=None):
        """Sets zero value for y-coordinates (shifts data verticaly).

        Args:
            value (number, optional): x-coordinate to set y-coordinate value as
                zero. If `None`, the minium value of x is used.
            n (int, optional): number of data points to average and setting to
                zero. This only applies if data is monotonicaly.
            ranges (list, optional): Pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                The mean value of the y-coordinates inside the defined range
                will be set to zero. This overwrites value and n.

        Returns:
            None
        """
        if ranges is None:
            if self.monotonicity is None:
                try:
                    self.check_monotonicity()
                except ValueError:
                    pass

            if value is None:
                i = index(self.x, min(self.x))
            else:
                i = index(self.x, value)

            if self.monotonicity == 'increasing' or self.monotonicity == 'decreasing':
                start = i-int(n/2)
                end   = i+int(n/2)
                if start < 0:
                    end   = i+int(n/2)-start
                    start = 0
                if end > len(self.y)-1:
                    start = i-int(n/2) - (end - len(self.y)-1)
                    end   = len(self.y)-1
                if start == end:
                    end = start + 1
                value2 = -np.mean(self.y[start:end])
                self.offset = self.offset+value2
            else:
                value2 = -self.y[i]
                self.offset = self.offset+value2
        else:
            ranges = self._check_ranges(ranges)
            _, y= extract(self.x, self.y, ranges)
            value2 = -np.mean(y)
            self.offset = self.offset+value2
        self._additive_y_fix(value2)

    def flip(self):
        self.x = -self.x
        # self._y = self.y[::-1]
        # self._data = np.vstack((self._x, self._y)).transpose()


    def plot(self, ax=None, offset=0, shift=0, factor=1, ranges=None, **kwargs):
        """Plot spectrum.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            vertical_increment (float or int, optional): defines the vertical offset. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            ranges (list, optional): ranges.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            None
        """

        if ax is None:
            ax = plt

        if ranges is None:
            ax.plot(self.x + shift, self.y*factor + offset, **kwargs)
        else:
            ranges = self._check_ranges(ranges)
            x, y = extract(self.x, self.y, ranges=ranges)
            ax.plot(x + shift, y*factor + offset, **kwargs)

    def plot_detected_peaks(self, ax=None, offset=0, shift=0, factor=1, **kwargs):
        """Plut a marker where peaks were detected.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            **kwargs: kwargs are passed to ``plt.errorbar()`` or ``plt.scatter()``.

        Returns:
            None

        Raises:
            ValueError: if :py:attr:`Spectrum.peaks` is not defined.
        """
        if self.peaks is None:
            raise ValueError('No peaks defined.')

        if ax is None:
            ax = plt

        # data
        c = np.array([self.peaks[i]['c'] for i in range(len(self.peaks))])
        amp = np.array([self.peaks[i]['amp'] for i in range(len(self.peaks))])
        xerr = np.array([self.peaks[i]['fwhm']/2 for i in range(len(self.peaks))])

        if 'lw' not in kwargs and 'linewidth' not in kwargs:
            kwargs['lw'] = 0
        if 'elinewidth' not in kwargs :
            kwargs['elinewidth'] = 2
        if 'marker' not in kwargs :
            kwargs['marker'] = 'o'
        if 'markersize' not in kwargs and 'ms' not in kwargs:
            kwargs['markersize'] = 5

        # for i in self.peaks:
        ax.errorbar(c+shift, amp*factor+offset, xerr=xerr, **kwargs)

    # OBSOLETE =================================================================
    def plot_fitted_ranges(self, ax=None, **kwargs):
        if self.fit_ranges is None:
            raise ValueError('data was not fitted. Use fit_peak().')


        if ax is None:
            ax = plt


        if 'lw' not in kwargs and 'linewidth' not in kwargs:
            kwargs['lw'] = 2
        if 'zorder' not in kwargs:
            kwargs['zorder'] = 10

        kwargs1 = copy.deepcopy(kwargs)
        kwargs2 = copy.deepcopy(kwargs)
        if 'color' not in kwargs1 and 'c' not in kwargs1:
            kwargs1['color'] = 'green'
        if 'color' not in kwargs2 and 'c' not in kwargs2:
            kwargs2['color'] = 'red'

        for r in self.fit_ranges:
            ax.vlines(r[0], min(self.y), max(self.y), **kwargs1)
            ax.vlines(r[1], min(self.y), max(self.y), **kwargs2)

    def fit_peaks_old(self, idx='all', ranges=None, asymmetry=False, fixed_m=0, multiplicity=1, offset=True):
            """

            Args:
                idx (number, list, or 'all'): idx of the peak to fit.
                ranges (list, optional) if None, ranges is automatically set, based
                    on the first and last peak to be fitted (given on `idx`).
                    If ranges is not None, only data inside ranges is considered in
                    the fitting. Also, it raises an error if any peak in idx is
                    outside of ranges.
                assymetry (bool or list of bools, optional): if True, fits each half of the
                    with a different width.
                fixed_m (False, number, or list, optional): m is the amount of lorentzian
                    contribution for a peak. If False, m will be fitted for each peak.
                    If a number (from 0 to 1), this will be used as the value of m
                    (fixed m).
                multiplicity (int or list, optional): use this if peaks have a shoulder
                    or some spliting. The peak will be fitted with `multiplicity` peaks.
                offset (bool, optional): if True, a offset value will be fitted.

            Returns:
                None

            """
            self.check_monotonicity()

            # check if peaks is defined ============================================
            if self.peaks is None:
                raise ValueError('Spectrum.peaks is not defined. Run Spectrum.find_peaks().')
            if self.peaks == []:
                raise ValueError('Spectrum.peaks = []. No peaks to fit.')

            # fix idx ==============================================================
            if idx == 'all':
                idx = [k for k in self.peaks]
            if isinstance(idx, Iterable) == False:
                idx = [idx, ]
            idx = [i if i>=0 else len(self.peaks)+i for i in idx]

            # check multiplicity ======================================================
            if isinstance(multiplicity, dict):
                for i in idx:
                    if i not in multiplicity:
                        multiplicity[i] = 1
                    else:
                        if multiplicity[i] < 1 and type(multiplicity) != int:
                            raise ValueError('multiplicity must be a positive integer.')
                for i in multiplicity:
                    if i not in idx:
                        raise ValueError(f'Peak {i} not in idx (idx = {idx}) (multiplicity = {multiplicity}).')
            else:
                if isinstance(multiplicity, Iterable):
                    raise ValueError('multiplicity must be a number or dictionary.')
                multiplicity = {i: multiplicity for i in idx}

            # check asymmetry ======================================================
            if isinstance(asymmetry, dict):
                for i in idx:
                    if i not in asymmetry:
                        asymmetry[i] = False
                    else:
                        if type(asymmetry[i]) != bool:
                            raise ValueError('asymmetry must be bool (True or False)')
                for i in asymmetry:
                    if i not in idx:
                        raise ValueError(f'Peak {i} not in idx (idx = {idx}) (asymmetry = {asymmetry}).')
            else:
                if isinstance(asymmetry, Iterable):
                    raise ValueError('asymmetry must be a number or dictionary.')
                asymmetry = {i: asymmetry for i in idx}

            # check fixed_m ========================================================
            if isinstance(fixed_m, dict):
                for i in idx:
                    if i not in fixed_m:
                        fixed_m[i] = 0
                    else:
                        if type(fixed_m[i]) == False or (fixed_m[i] > 0 and fixed_m[i] < 1):
                            raise ValueError('fixed_m must be a False or a number from 0 to 1.')
                for i in fixed_m:
                    if i not in idx:
                        raise ValueError(f'Peak {i} not in idx (idx = {idx}) (fixed_m = {fixed_m}).')
            else:
                if isinstance(fixed_m, Iterable):
                    raise ValueError('fixed_m must be a number, False, or dictionary.')
                fixed_m = {i: fixed_m for i in idx}

            # fix ranges ===========================================================
            if ranges is None:    # try fitting assigned peaks (range is automatically set)
                ranges = [min(self.x),  max(self.x)]
                p_max = self.peaks[max(idx)]['c']+self.peaks[max(idx)]['fwhm']
                p_min = self.peaks[min(idx)]['c']-self.peaks[min(idx)]['fwhm']

                if max(idx)+1 < len(self.peaks):
                    next_min = self.peaks[max(idx)+1]['c']-self.peaks[max(idx)+1]['fwhm']
                    if p_max > next_min:
                        ranges[1] = p_max + self.peaks[max(idx)]['fwhm']/2
                    else:
                        x_temp, y_temp = extract(self.x, self.y, (p_max, next_min))
                        ranges[1] = x_temp[np.argmin(y_temp)]
                else:
                    next_min = max(self.x)
                    ranges[1] = max(self.x)

                if min(idx)-1 >= 0:
                    previous_max = self.peaks[min(idx)-1]['c']-self.peaks[min(idx)-1]['fwhm']
                    if p_min < previous_max:
                        ranges[0] = p_min - self.peaks[min(idx)]['fwhm']/2
                    else:
                        x_temp, y_temp = extract(self.x, self.y, (previous_max, p_min))
                        ranges[0] = x_temp[np.argmin(y_temp)]
                else:
                    previous_max = min(self.x)
                    ranges[0] = min(self.x)

                ranges = (ranges, )
                x2fit, y2fit = extract(self.x, self.y, ranges=ranges)
            else:
                if isinstance(ranges[0], Iterable) == False:
                    ranges = (ranges, )
                x2fit, y2fit = extract(self.x, self.y, ranges=ranges)

                # check if peaks are inside ranges
                for i in idx:
                    flag = True
                    for r in ranges:
                        if self.peaks[i]['c'] > r[0] and self.peaks[i]['c'] < r[1]:
                            flag = False
                    if flag:
                        raise ValueError(f'peak {i} outside of range.')

            # fit function =====================================================
            f_temp = {}

            model_str = ''
            args_str = ''
            bounds_min = []
            bounds_max = []
            p0 = []

            total_i = -1
            for i in idx:
                for extra in range(multiplicity[i]):
                    total_i += 1
                    f, f_str, a_str = _peak_function_creator(asymmetry=asymmetry[i], fixed_m=fixed_m[i], idx=total_i)
                    if fixed_m[i] == False and type(fixed_m[i]) == bool:  # variable m
                        if asymmetry[i]:
                            n_args = 6
                            p0 = np.append(p0, [self.peaks[i]['amp']-min(y2fit), self.peaks[i]['c']+self.peaks[i]['fwhm']/4*extra, self.peaks[i]['fwhm'], 0.5, self.peaks[i]['fwhm'], 0.5])
                            bounds_min = np.append(bounds_min, [            0,          self.peaks[i]['c']-self.peaks[i]['fwhm']*2,            0,          0,       0,               0])
                            bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*3, self.peaks[i]['c']+self.peaks[i]['fwhm']*2, self.peaks[i]['fwhm']*4, 1, self.peaks[i]['fwhm']*4, 1])
                        else:
                            n_args = 4
                            p0 = np.append(p0, [self.peaks[i]['amp']-min(y2fit), self.peaks[i]['c']+self.peaks[i]['fwhm']/4*extra, self.peaks[i]['fwhm'], 0.5])
                            bounds_min = np.append(bounds_min, [            0,          self.peaks[i]['c']-self.peaks[i]['fwhm']*2,            0,          0])
                            bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*3, self.peaks[i]['c']+self.peaks[i]['fwhm']*2, self.peaks[i]['fwhm']*4, 1])
                    else:
                        if asymmetry[i]:
                            n_args = 4
                            p0 = np.append(p0, [self.peaks[i]['amp']-min(y2fit), self.peaks[i]['c']+self.peaks[i]['fwhm']/4*extra, self.peaks[i]['fwhm'], self.peaks[i]['fwhm']])
                            bounds_min = np.append(bounds_min, [            0,          self.peaks[i]['c']-self.peaks[i]['fwhm']*2,            0,                0,             ])
                            bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*3, self.peaks[i]['c']+self.peaks[i]['fwhm']*2, self.peaks[i]['fwhm']*4, self.peaks[i]['fwhm']*4])
                        else:
                            n_args = 3
                            p0 = np.append(p0, [self.peaks[i]['amp']-min(y2fit), self.peaks[i]['c']+self.peaks[i]['fwhm']/4*extra, self.peaks[i]['fwhm']])
                            bounds_min = np.append(bounds_min, [            0,          self.peaks[i]['c']-self.peaks[i]['fwhm']*2,            0,        ])
                            bounds_max = np.append(bounds_max, [self.peaks[i]['amp']*3, self.peaks[i]['c']+self.peaks[i]['fwhm']*2, self.peaks[i]['fwhm']*4])
                    f_temp[total_i] = {'asymmetry':asymmetry[i], 'fixed_m':fixed_m[i], '_func':f, 'n_args':n_args}
                    model_str += f_str + ' + '
                    args_str += a_str + ', '

            if offset:
                model_str = f'lambda x, {args_str}offset: {model_str}offset'
                p0 = np.append(p0, min(y2fit))
                bounds_min = np.append(bounds_min, min(y2fit))
                bounds_max = np.append(bounds_max, np.mean(y2fit))
            else:
                model_str = f'lambda x, {args_str[:-2]}: {model_str[:-3]}'
            model = eval(model_str)

            # fit peaks ============================================================
            # print(p0)
            # print(bounds_min)
            # print(bounds_max)
            # plt.plot(x2fit, y2fit)
            # plt.plot(x2fit, model(x2fit, *p0))
            popt, pcov = curve_fit(model, x2fit, y2fit, p0=p0, bounds = (bounds_min, bounds_max))
            # plt.plot(x2fit, model(x2fit, *popt))

            # final ============================================================
            x_temp = np.linspace(self.x[0], self.x[-1], len(x2fit)*3)
            stop = 0
            for i in range(total_i+1):
                start = stop
                stop  = start + f_temp[i]['n_args']
                f_temp[i]['popt'] = popt[start:stop]
                f_temp[i]['pcov'] = pcov[start:stop]
                f_temp[i]['amp'] = popt[start]
                f_temp[i]['c'] = popt[start+1]
                if f_temp[i]['fixed_m'] == False and type(f_temp[i]['fixed_m']) == bool:  # variable m
                    if f_temp[i]['asymmetry']:
                        f_temp[i]['fwhm'] = popt[start+2]/2 + popt[start+4]/2
                        f_temp[i]['fwhm1'] = popt[start+2]
                        f_temp[i]['fwhm2'] = popt[start+4]
                        f_temp[i]['m1'] = popt[start+3]
                        f_temp[i]['m2'] = popt[start+5]
                    else:
                        f_temp[i]['fwhm'] = popt[start+2]
                        f_temp[i]['m'] = popt[start+3]
                else:
                    if f_temp[i]['asymmetry']:
                        f_temp[i]['fwhm'] = popt[start+2]/2 + popt[start+3]/2
                    else:
                        f_temp[i]['fwhm'] = popt[start+2]
                    f_temp[i]['m'] = f_temp[i]['fixed_m']

                if offset:
                    def peak(*args):
                        return lambda x: (f_temp[i]['_func'](x/self.calib - self.shift, *args) + popt[-1])*self.factor + self.offset
                    f_temp[i]['func'] = peak(*f_temp[i]['popt'])
                else:
                    def peak(*args):
                        return lambda x: (f_temp[i]['_func'](x/self.calib - self.shift, *args))*self.factor + self.offset
                    f_temp[i]['func'] = peak(*f_temp[i]['popt'])
                f_temp[i]['curve'] = Spectrum(x=x_temp, y=f_temp[i]['func'](x_temp))


            self.fit_data = f_temp
            self.fit_popt = popt
            # self.fit_func = lambda x: (model(x/self.calib - self.shift, *popt))*self.factor + self.offset
            self.fit_func = lambda x: (model(x, *popt))
            self.fit_ranges = ranges
            self.fit = Spectrum(x=x_temp, y=self.fit_func(x_temp))
            # self.fit.calib = self.calib
            # self.fit.shift = self.shift

    def _plot_old(self, ax=None, normalized=False, vertical_increment=0, shift=0, factor=1, show_fit_ranges=False, show_fit=False, show_fit_partial=False, mark_peaks=False, kwargs_fit={}, **kwargs):
        """Plot spectrum.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.
            vertical_increment (float or int, optional): defines the vertical offset. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            show_fit_ranges (bool, optional): show ranges used to fit peaks.
            show_fit (bool, optional): show fitted curve.
            show_fit_partial (bool, optional): same as `show_fit`, but it shows all
                peaks used to fit.
            mark_peaks (book, optional): put a marker where peaks were detected.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes
        """

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if normalized:
            factor = factor/max(self.y)


        ax.plot((self.x + shift), self.y*factor + vertical_increment, **kwargs)
        if show_fit:
            self.plot_fit(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, **kwargs_fit)
        if show_fit_partial:
            self.plot_fit_partial(ax=ax, normalized=normalized, vertical_increment=vertical_increment, shift=shift, factor=factor, color='black', lw=0.5, ls='--')

        if mark_peaks:
            self.mark_peaks(ax=ax)

        if show_fit_ranges:
            for r in self.fit_ranges:
                ax.vlines(r[0], min(self.y), max(self.y), color='green', linewidth=2, zorder=10)
                ax.vlines(r[1], min(self.y), max(self.y), color='red', linewidth=2, zorder=10)

        return ax

    def _plot_fit_old(self, ax=None, normalized=False, vertical_increment=0, shift=0, factor=1, **kwargs):
        """Plot spectrum.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.
            vertical_increment (float or int, optional): defines the vertical offset. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            show_ranges (bool, optional): show ranges in which offsets were calculated.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes

        Raises:
            ValueError: if :py:attr:`Spectrum.fit_func` is not defined.
        """
        if self.fit_func is None:
            raise ValueError('No fit data to plot.')

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if normalized:
            factor = factor/max(self.y)

        # x,  y = extract(self.x, self.y, self.elastic.ranges)
        x = np.linspace(min(self.x), max(self.x), len(self.x)*10)
        ax.plot(x + shift, self.fit_func(x)*factor + vertical_increment, **kwargs)

        return ax

    def _plot_fit_partial_old(self, ax=None, normalized=False, vertical_increment=0, shift=0, factor=1, **kwargs):
        """Plot spectrum.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            normalized (bool, optional): if True, spectrum is normalized by its
                maximum value.
            vertical_increment (float or int, optional): defines the vertical offset. Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            show_ranges (bool, optional): show ranges in which offsets were calculated.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes

        Raises:
            ValueError: if :py:attr:`Spectrum.fit_data` is not defined.
        """
        if self.fit_data is None:
            raise ValueError('No fit data to plot.')

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        if normalized:
            factor = factor/max(self.y)

        # x,  y = extract(self.x, self.y, self.elastic.ranges)
        x = np.linspace(min(self.x), max(self.x), len(self.x)*10)
        for i in self.fit_data:
            ax.plot(x + shift, self.fit_data[i]['func'](x)*factor + vertical_increment, **kwargs)

        return ax

    def _mark_peaks_old(self, ax=None, show_width=True, **kwargs):
            """Plut a marker where peaks were detected.

            Args:
                ax (matplotlib.axes, optional): axes for plotting on.
                show_width (bool, optional): if True, it will use ``plt.errorbar()`` to
                    plot a horizontal line with the marker the size of the estimated
                    width. If False, it uses ``plt.scatter()`` to draw the markers.
                **kwargs: kwargs are passed to ``plt.errorbar()`` or ``plt.scatter()``.

            Returns:
                matplotlib.axes

            Raises:
                ValueError: if :py:attr:`Spectrum.peaks` is not defined.
            """
            if self.peaks is None:
                raise ValueError('No peaks defined.')

            if ax is None:
                fig = plt.figure()
                ax = fig.add_subplot(111)

            # data
            c = [self.peaks[i]['c'] for i in self.peaks]
            amp = [self.peaks[i]['amp'] for i in self.peaks]
            xerr = [self.peaks[i]['fwhm']/2 for i in self.peaks]



            if show_width:
                if 'lw' not in kwargs and 'linewidth' not in kwargs:
                    kwargs['lw'] = 0
                if 'elinewidth' not in kwargs :
                    kwargs['elinewidth'] = 2
                if 'marker' not in kwargs :
                    kwargs['marker'] = 'o'
                if 'markersize' not in kwargs and 'ms' not in kwargs:
                    kwargs['markersize'] = 5

                for i in self.peaks:
                    ax.errorbar(c, amp, xerr=xerr, **kwargs)
            else:
                if 'lw' not in kwargs and 'linewidth' not in kwargs:
                    kwargs['lw'] = 0
                if 'marker' not in kwargs :
                    kwargs['marker'] = 'o'
                if 's' not in kwargs:
                    kwargs['s'] = 100

                for i in self.peaks:
                    ax.scatter(c, amp, **kwargs)

            return ax


class Spectra(metaclass=_Meta):
    """Returns a ``spectra`` class type object for dealing with many spectrum at a time.

    Args:
        n (int, optional): when the number of spectra `n` is known, one can preallocate
            the data array.
        data (list or array, optional): list of :py:class:`spectrum` objects.
            Overwrites ``n``.

    Attributes:
        data (list of :py:attr:`Spectrum`): list with :py:attr:`Spectrum` objects.
        step (number, read only): step size for x-coordinates.
        length (number, read only): length of x-coordinates vector.
        x (number, read only): x-coordinates (only applicable if x-coordinates
            are the same for all spectra.
        shifts (list of int): Calculated shifts as int. Must always have the same length as data.
        shift_ref (number, read only): index of the spectrum used as reference
            for shift calculation.
        shift_ref_value (number, read only): value of the shift for the reference
            spectrum, e.g.:
                * if shift_calc_mode is 'peak'
                    shift_ref_value returns the position of the elastic peak of the reference spectrum.
                * if shift_calc_mode is 'cross-correlation'
                    shift_ref_value returns the max value of the cross-correlation of the reference spectrum with itself.
                * if shift_calc_mode is 'max'
                    shift_ref_value returns the position of the maximum value of the reference spectrum.
        shift_calc_mode (string, read only): mode used to calculate shifts.
        shift_ranges (list, read only): data range used to calculate shifts.
        disp_values (list, read only): value of the 'relative spectrum position'
            for each spectrum. See :py:attr:`shift_ref_value`
        disp_energies (number, read only): energy value associated with each spectrum.
        disp_popt (number, read only): optimal parameters of the fit of
            disp_values vs disp_energies.
        disp_func (number, read only): function of the fitted curve.
        disp (number, read only): final calculated dispersion.
        sum (number, read only): spectrum resulting from the sumation of all spectra.
    """

    _read_only     = ['step', 'length', 'x', 'monotonicity',
                      'shift_calculated',# 'shifts_ref', 'shifts_ref_value', 'shifts_calc_mode', 'shifts_ranges',
                      'factor_calculated',# 'factors_ref', 'factors_ref_value', 'shift_calc_mode', 'shift_ranges',
                      'calib_calculated',
                      'offset_calculated',
                      'calculated_shifts',
                      'sum']
    _non_removable = []


    def __init__(self, *args, **kwargs):

        data, filepaths, n = self._args_checker(args, kwargs)

        if data is not None:
            self.data = data
        elif filepaths is not None:
            self.load(filepaths)
        elif n is not None:
            self._data     = [-1]*n
            self._restart_sum_attr()
            self._restart_check_attr()
            self._restart_shift_attr()
            self._restart_calib_attr()
            self._restart_factor_attr()
        else:
            self.data = data

    def _args_checker(self, args, kwargs):
        """checks initial arguments.

         Keyword arguments (kwargs) cannot be mixed with positional arguments.

        The hierarchy for Keyword arguments is: 1) data, 2) y (and x), and finaly
            3) filepath. For example, if `data` and `filepath` are passed as
            arguments, `filepath` is ignored.

        For positional arguments, if one data set is passed, it assumes it is
            `data`. If this one argument is of type string or Pathlib.Path, it
            assumes it is a filepath. If two data sets are passed, it will
            assume one is the x coordinates and the next one is the y coordinates.

        Raises:
            AttributeError: if kwargs and args cannot be read.

        Returns:
            data, x, y, filepath
        """
        # print(kwargs)
        # print(args)
        if kwargs != {} and args != ():
            raise AttributeError('cannot mix key word arguments with positional arguments. Key word arguents are `x`, `y`, `data`, and `filepath`.')
        if any([item not in ['data', 'n', 'filepaths'] for item in kwargs.keys()]):
            raise AttributeError(f'invalid attributes.\nValid atributes are `data`, `n`, and `filepaths`\nInput attributes: {kwargs.keys()}')

        data = None
        n    = None
        filepaths = None
        if 'data' in kwargs:
            data = kwargs['data']
        elif 'filepaths' in kwargs:
            filepaths = kwargs['filepaths']
        elif 'n' in kwargs:
            n = kwargs['n']
        elif len(args) == 1:
            if isinstance(args[0], int):
                n = args[0]
            elif isinstance(args[0], str) or isinstance(args[0], Path):
                filepaths = [args[0], ]
            elif isinstance(args[0], Iterable):
                if isinstance(args[0][0], str) or isinstance(args[0][0], Path):
                    filepaths = args[0]
                else:
                    data = args[0]
        elif len(args) > 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepaths = args
            elif isinstance(args[0], Iterable):
                data = args
        return data, filepaths, n

    def load(self, filepaths, comments='#', delimiter=None, string='*'):
        self._data     = []
        for filepath in filepaths:
            if Path(filepath).is_dir():
                fl = filelist(dirpath=filepath, string=string)
                for i, f in enumerate(fl):
                    print(f'Loading {i}/{len(fl)-1}: {f}')
                    self.append(Spectrum(filepath=f))
                print('done')
            elif Path(filepath).is_file():
                self.append(Spectrum(filepath=filepath))
            else:
                raise ValueError('filepath is not valid.')


    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        if value is None:
            self._data = []
        else:
            # print(type(value))
            if isinstance(value, Iterable):
                # for s in value:
                #     if isinstance(s, Spectrum) == False:
                #         raise ValueError('all entries must be of type brixs.spectrum.')
                self._data = copy.deepcopy(value)
            else:
                raise ValueError('data must be a list.')
        self._restart_sum_attr()
        self._restart_check_attr()
        self._restart_shift_attr()
        self._restart_calib_attr()
        self._restart_factor_attr()
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].shift
        return temp
    @shift.setter
    def shift(self, value):
        self.set_shift(value=value, mode='x')
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_roll(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].shift_roll
        return temp
    @shift_roll.setter
    def shift_roll(self, value):
        self.set_shift(value=value, mode='roll')
    @shift_roll.deleter
    def shift_roll(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_interp(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].shift_interp
        return temp
    @shift_interp.setter
    def shift_interp(self, value):
        self.set_shift(value=value, mode='interp')
    @shift_interp.deleter
    def shift_interp(self):
        raise AttributeError('Cannot delete object.')

    @property
    def calib(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].calib
        return temp
    @calib.setter
    def calib(self, value):
        self.set_calib(value)
    @calib.deleter
    def calib(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].factor
        return temp
    @factor.setter
    def factor(self, value):
        self.set_factor(value)
    @factor.deleter
    def factor(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].offset
        return temp
    @offset.setter
    def offset(self, value):
        self.set_offset(value)
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    @property
    def R2(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].R2
        if None in temp:
            print(f'Some spectra do not have R2 defined.\nR2={temp}')
        return temp
    @R2.setter
    def R2(self, value):
        raise AttributeError('Cannot be set.')
    @R2.deleter
    def R2(self):
        raise AttributeError('Cannot delete object.')

    @property
    def peaks(self):
        n = len(self)
        n_peaks = [None]*n
        for i in range(n):
            if self[i].peaks is not None:
                n_peaks[i] = len(self[i].peaks)
        if None in n_peaks:
            t = {i:n_peaks[i] for i in range(len(n_peaks))}
            raise RuntimeError(f'Peaks are not defined for some spectra. Use Spectra.find_peaks().\nnumber of peaks for each spectrum: {t}')
        if any(x!=n_peaks[0] for x in n_peaks):
            t = {i:n_peaks[i] for i in range(len(n_peaks))}
            raise RuntimeError(f'Peak data can only be retrieved if all spectra have the same number of peaks.\nnumber of peaks for each spectrum: {t}\n')

        temp = [0]*len(self[0].peaks)
        for i in range(len(self[0].peaks)):
                temp[i] = {key:np.array([s.peaks[i][key] for s in self]) for key in self[0].peaks[i].store}
        return temp
    @peaks.setter
    def peaks(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @peaks.deleter
    def peaks(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def fit(self):
        ss = Spectra(n=len(self))
        for i in range(len(self)):
            ss[i] = self[i].fit
        # if None in temp:
        #     print(f'Some spectra do not have fit defined.\nfit={temp}')
        return ss
    @fit.setter
    def fit(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @fit.deleter
    def fit(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def guess(self):
        ss = Spectra(n=len(self))
        for i in range(len(self)):
            ss[i] = self[i].guess
        # if None in temp:
        #     print(f'Some spectra do not have fit defined.\nfit={temp}')
        return ss
    @guess.setter
    def guess(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @guess.deleter
    def guess(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def residue(self):
        ss = Spectra(n=len(self))
        for i in range(len(self)):
            ss[i] = self[i].residue
        # if None in temp:
        #     print(f'Some spectra do not have fit defined.\nfit={temp}')
        return ss
    @residue.setter
    def residue(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @residue.deleter
    def residue(self):
        raise AttributeError('Attribute cannot be deleted.')

    def __str__(self):
        return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

    def __repr__(self):
        return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

    def __next__(self):
        result = self.data[self._index]
        self._index +=1
        return result
        # End of Iteration
        raise StopIteration

    def __getitem__(self, item):
         return self.data[item]

    def __setitem__(self, item, value):
        # if isinstance(value, Spectrum) == False:
        #     raise ValueError('value must be brixs.spectrum.')
        self.data[item] = value
        self._restart_sum_attr()
        self._restart_check_attr()
        self._restart_shift_attr()
        self._restart_calib_attr()

    def __len__(self):
        return self.get_spectra_count()

    def __delitem__(self, item):
        self.remove(item)

    def _restart_sum_attr(self):
        # sum attr
        self._sum = None

    def _restart_check_attr(self):
        # check x attrs
        self._length = None
        self._step = None
        self._x = None
        self._monotonicity = None

    def _restart_shift_attr(self):
        self._shift_calculated = {'values':np.array([np.NaN]*len(self)),
                                   'mode':None,
                                   'ref':None,
                                   'ref_value':None,
                                   'ranges':None,
                                  }

    def _restart_calib_attr(self):
        """Set relevant attributes to initial value."""
        self._calib_calculated = {'values':np.array([np.NaN]*len(self)),
                                  'centers':np.array([np.NaN]*len(self)),
                                   'mode':None,
                                   'popt':None,
                                   'value':None,
                                   'ranges':None,
                                   'func':None,
                                  }

    def _restart_factor_attr(self):
        """Set relevant attributes to initial value."""
        self._factor_calculated = {'values':np.array([np.NaN]*len(self)),
                                   'mode':None,
                                   'ref':None,
                                   'ref_value':None,
                                   'peak_idx':None,
                                  }

    def get_spectra_count(self):
        """Returns the number of spectra."""
        return len(self.data)

    def append(self, s=None, x=None, y=None, data=None):
        """Add spectrum to the spectrum list.

        Args:
            s (spectrum obj): :py:class:`Spectrum` object to be added.
            data (list or array, optional): two column list (or array).
            x (list or array, optional): x values (1D list/array). Overwrites `data`.
            y (list or array, optional): y values (1D list/array). Overwrites `data`.

        Returns:
            None

        See Also:
            :py:func:`Spectra.remove`.
        """
        if s is None:
            s=Spectrum(x=x, y=y, data=data)
            # print(s.x[0])
            self.append(s=s)
        elif s is not None:
            if isinstance(s, Iterable):
                # if isinstance(s, Spectrum) == False:
                #     raise ValueError('spectrum must be of type brixs.Spectrum.')
                self._data += copy.deepcopy(s)
            else:
                # for temp in s:
                #     if isinstance(temp, Spectrum) == False:
                #         raise ValueError('all entries must be of type brixs.spectrum.')
                self._data += [copy.deepcopy(s)]
        else:
            raise ValueError('No data to append.')
        # print(s.factor)
        self._restart_check_attr()
        self._restart_sum_attr()
        self._restart_shift_attr()
        self._restart_calib_attr()
        self._restart_factor_attr()

    def remove(self, idx):
        """Exclude spectrum from the spectrum list.

        Args:
            idx (int): index of the spectrum.

        Returns:
            None

        See Also:
            :py:func:`Spectra.append`.
        """
        del self.data[idx]
        self._restart_sum_attr()
        self._restart_check_attr()
        self._restart_shift_attr()
        self._restart_calib_attr()
        self._restart_factor_attr()


    def check_monotonicity(self):
        monotonicity = [None]*self.get_spectra_count()
        for i in range(self.get_spectra_count()):
            try:
                self[i].check_monotonicity()
                monotonicity[i] = self.data[i].monotonicity
            except ValueError:
                pass
        if all(x == 'increasing' for x in monotonicity):
            self._monotonicity = 'increasing'
        elif all(x == 'decreasing' for x in monotonicity):
            self._monotonicity = 'decreasing'
        else:
            text = ''
            for i in range(self.get_spectra_count()):
                text += f'spectrum: {i}, motonicity: {monotonicity[i]}\n'
            raise ValueError(f'some spectra have different monotonicity (increasing, decreasing) or no monotonicity at all (None): \n{text}')

    def fix_monotinicity(self, mode='increasing'):
        for s in self.data:
            s.fix_monotinicity(mode=mode)

    def check_length(self):
        """Checks if all spectra has the same length.

        If all spectra have the same length, it sets :py:attr:`Spectra.length` = length.
        Otherwise, it raises an error.

        Returns:
            None

        Raises:
            ValueError: spectra does not have the same length.

        See Also:
            :py:func:`Spectra.check_step_x`, :py:func:`Spectra.check_same_x`.
        """
        for idx, s in enumerate(self.data):
            try:
                if len(s.data) != len(self.data[idx+1].data):
                    self._length = None
                    raise ValueError(f"Spectrum {idx} and {idx+1} have the different length.")
            except IndexError:
                pass
        self._length = len(self.data[0].x)

    def check_step_x(self, max_error=0.1):
        """Check step between data points in the x-coordinates.

            If data has a well defined step size, it sets :py:attr:`Spectra.step` = step.
            Otherwise, it raises an error.

            Args:
                max_error (number, optional): percentage value (in terms of the
                    average x step) of the maximum allowed error.

            Three checks are performed:

            1. If not yet, checks if all spectra have same length by using
            :py:func:`Spectra.check_length`.

            2. Checks if the step between two data points is the same
            through out the x vector for each spectrum (vector uniformity).

            3. Checks if the step is the same between different spectra.

            Returns:
                None

            Raises:
                ValueError: If condition 1, 2, or 3 are not satisfied.

            See Also:
                :py:func:`Spectra.check_length`, :py:func:`Spectra.check_same_x`.
        """

        # 1) check length
        if self.length is None:
            self.check_length()

        # 2) check step uniformity
        steps = np.zeros(self.get_spectra_count())
        for idx, s in enumerate(self.data):
            if s.step is None:
                try:
                    s.check_step_x()
                except ValueError:
                    raise ValueError(f"Step in the x-coordinate of spectrum {idx} seems not to be uniform.")
            steps[idx] = s.step

        # 3) check step between spectra
        avg_step = np.mean(steps)
        if sum([abs(steps[i]-steps[i+1]) > abs(avg_step*max_error/100) for i in range(self.get_spectra_count()-1)]) > 0:
            self._step = None
            raise ValueError(f"Spectra seems to have different step size. Calculated step sizes = {steps}")
        self._step = avg_step

    def check_same_x(self, max_error=0.1):
        """Check if spectra have same x-coordinates.

        If data has same x-coordinates, it sets :py:attr:`Spectra.x` = x.
        Otherwise, it raises an error.

        Args:
            max_error (number, optional): percentage value (in terms of the
                x step) of the maximum allowed error.

        Returns:
            None

        Raises:
            ValueError: If any x-coodinate of any two spectrum is different.

        See Also:
            :py:func:`Spectra.check_length`, :py:func:`Spectra.check_step_x`.
        """
        if self.step is None:  # calculate step if it hasn't yet
            self.check_step_x(max_error=max_error)

        # check x between spectra
        for idx, s in enumerate(self.data):
            try:
                if max(abs(s.x - self.data[idx+1].x))*100/self.step > max_error:
                    self._x = None
                    raise ValueError(f"x axis of spectrum {idx} and {idx+1} seem to be different.\nTip: use brixs.Spectra.interp() to interpolate the data and make the x axis for different spectra match.")
            except IndexError:
                pass
        self._x = self.data[0].x

    def check_fit(self):
        temp = [None]*len(self)
        for i in range(len(self)):
                temp[i] = self[i].fit
        if None in temp:
            raise RuntimeError(f'Fit is not defined for some spectra. Use Spectra.fit_peaks().\nlist of fits: {temp}')

    def check_fitted_peaks(self):
        self.check_fit()

        n_peaks = [None]*len(self)
        for i in range(len(self)):
            if self[i].fit.peaks is not None:
                n_peaks[i] = len(self[i].fit.peaks)
        if None in n_peaks:
            raise RuntimeError(f'some fitted data have no peaks. Maybe try fitting again\nnumber of peaks for each spectrum: {n_peaks}')
        if any(x!=n_peaks[0] for x in n_peaks):
            return n_peaks
            # print(f'Different spectra have different number of fitted peaks. \nnumber of fitted peaks for each spectrum: {n_peaks}\n')


    def _gather_ys(self, ranges=None):
        """Return two lists, x and y's"""
        # check x
        if self.x is None:
            self.check_same_x()

        ys = np.zeros((self.length, self.get_spectra_count()))
        for i in range(self.get_spectra_count()):
            ys[:, i] = self.data[i].y
        if ranges is None:
            x = self.data[0].x
        else:
            x, ys = extract(self.data[0].x, ys, ranges=ranges)
        return x, ys


    def find_peaks(self, prominence=None, width=4, moving_average_window=8, ranges=None):
        for s in self.data:
            s.find_peaks(prominence=prominence, width=width, moving_average_window=moving_average_window, ranges=ranges)


    def fit_peaks(self, idx='all', asymmetry=False, fixed_m=0, ranges=None, offset=True):
        # n=self.get_spectra_count()
        # self.guess = Spectra(n=n)
        # self.fit = Spectra(n=n)
        # self.residue = Spectra(n=n)

        n = len(self)

        # check number of peaks is the same
        n_peaks = [None]*n
        for i in range(n):
            if self[i].peaks is not None:
                n_peaks[i] = len(self[i].peaks)
        if None in n_peaks:
            raise RuntimeError(f'peaks were not defined for some spectra. Use Spectra.find_peaks().\nnumber of peaks for each spectrum: {n_peaks}')
        if any(x!=n_peaks[0] for x in n_peaks):
            print(f'Different spectra have different number of peaks. \nnumber of peaks for each spectrum: {n_peaks}\nFitting anyway...')


        # fit peaks
        for i in range(n):
            self.data[i].fit_peaks(idx=idx, asymmetry=asymmetry, fixed_m=fixed_m, ranges=ranges, offset=offset)
            # self.guess[i] = self.data[i].guess
            # self.fit[i] = self.data[i].fit
            # self.residue[i] = self.data[i].residue


    def fit_peak(self, asymmetry=False, fixed_m=0, ranges=None, offset=True):
        for s in self:
            s.fit_peak(asymmetry=asymmetry, fixed_m=fixed_m, ranges=ranges, offset=offset)

    def calculate_shifts(self, ref_spectrum=0, mode='cross-correlation', ranges=None, verbose=False, idx=0, bypass=False, **kwargs):
        """Calculate how much shifted spectra are relatively to a reference spectrum.

        If :py:attr:`Spectra.step` is defined, it also sets :py:attr:`Spectra.shifts`.

        Args:
            ref (int, optional): index of reference spectrum. The shift of all
                other spectra is calculated based on the reference spectrum.
                Default is 0.
            mode (string, optional): method used to calculate the offsets.
                The current options are:

                1. 'cross-correlation' ('cc')
                2. 'max' ('m')
                3. 'elastic' ('e')

            ranges (list, optional): a pair of x values or a list of pairs. Each pair represents
                the start and stop of a data range. If None, the whole data set is used.
                If mode is peak, ranges is passed to Spectrum.fit_peak()
            verbose (bool,optional): turn verbose on/off.
            **kwargs: if ``mode = 'elastic'``, kwargs are passed to
                :py:func:`Spectra.guess_elastic_peak`.

        Returns:
            None

        Raises:
            ValueError: if ref_spectrum is not an integer or a non-valid mode is selected.

        Note:
            For ``mode = 'cc'``, spectra must have the same x-coordinates (this
            is checked before execution).

        See Also:
            :py:func:`Spectra.apply_shifts`
        """
        # check ref ============================================================
        if type(ref_spectrum) != int:
            raise ValueError('ref must be an integer.')
        else:
            ref = ref_spectrum

        # try checking if x is the same (execution is much faster if it is) ====
        if self.x is None:
            try:
                self.check_same_x()
            except ValueError:
                pass

        # common variables =====================================================
        shifts = np.array([0.0]*len(self))

        # CROSS-CORRELATION ====================================================
        if mode in cc:
            if bypass==False:
                # check if background is zero
                temp = [False]*len(self)
                for i in range(len(self)):
                    peaks, d = find_peaks(self[i].y, prominence=(max(self[i].y)-np.mean(self[i].y))*0.05, width=1)
                    bkg = list(self[i].y)
                    for j in range(len(peaks)):
                        del bkg[peaks[j]-int(d['widths'][j]*10):peaks[j]+int(d['widths'][j])*10]
                    temp[i] = np.mean(self[i].y)>max(self[i].y)*0.1
                if True in temp:
                    raise ValueError(f'some spectra have a background that seems to be too big for a reliable cross-corelation.\nSpectra with big bkg: {temp}\nTip: You can try using Spectra.floor()')

            # x must be the same for cc
            x, ys = self._gather_ys(ranges=ranges)
            # calculate cross-correlation
            for i in range(self.get_spectra_count()):
                if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                cross_correlation = np.correlate(ys[:, ref], ys[:, i], mode='full')
                shifts[i] = np.argmax(cross_correlation)
                # if self.data[i].shift_roll != 0: # fix shift in case there was a previous shift set
                #     shifts[i] += self.data[i].shift_roll
            shifts -= shifts[ref]
            # shifts *= -1

            # save values
            self._shift_calculated['ref_value'] = shifts[ref]
        # MAX ==================================================================
        elif mode == 'max':
            if self.x is not None: # if all data has the same x, execution is much faster
                x, ys = self._gather_ys(ranges=ranges)
                j_ref = np.argmax(ys[:, ref])
                self._shift_calculated['ref_value'] = x[j_ref]
                for i in range(self.get_spectra_count()):
                    if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                    j = np.argmax(ys[:, i])
                    shifts[i] = -(x[j] - x[j_ref])
                    # if self.data[i].shift != 0: # fix shift in case there was a previous shift set
                    #     if self.data[i].shift_mode in roll:
                    #         shifts[i] += self.data[i].shift*self.data[i].step
                    #     else:
                    #         shifts[i] += self.data[i].shift
            else:  # if x is not the same, execution is slower
                if ranges is None:   # if ranges is not defined, execution is faster
                    j_ref = np.argmax(self.data[ref].y)
                    self._shift_calculated['ref_value'] = self.data[ref].x[j_ref]
                    for i in range(self.get_spectra_count()):
                        if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                        j = np.argmax(self.data[i].y)
                        shifts[i] = -(self.data[i].x[j] - self.data[ref].x[j_ref])
                        # if self.data[i].shift != 0: # fix shift in case there was a previous shift set
                        #     if self.data[i].shift_mode in roll:
                        #         shifts[i] += self.data[i].shift*self.data[i].step
                        #     else:
                        #         shifts[i] += self.data[i].shift
                else: # if ranges is defined, execution is slower
                    x_ref, y_ref = extract(self.data[ref].x, self.data[ref].y, ranges=ranges)
                    j_ref = np.argmax(y_ref)
                    self._shift_calculated['ref_value'] = x_ref[j_ref]
                    for i in range(self.get_spectra_count()):
                        x, y = extract(self.data[i].x, self.data[i].y, ranges=ranges)
                        if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                        j = np.argmax(y)
                        shifts[i] = -(x[j] - x_ref[j_ref])
                        # if self.data[i].shift != 0: # fix shift in case there was a previous shift set
                        #     if self.data[i].shift_mode in roll:
                        #         shifts[i] += self.data[i].shift*self.data[i].step
                        #     else:
                        #         shifts[i] += self.data[i].shift
        # PEAK =================================================================
        elif mode == 'peak':
            n = len(self)

            # check if fit is defined for all spectra and if spectra has the same number of peaks
            temp = self.check_fitted_peaks()
            if temp is not None:
                print(f'Different spectra have different number of fitted peaks. \nnumber of fitted peaks for each spectrum: {temp}\n')

            # calculate
            self._shift_calculated['ref_value'] = self.data[ref].fit.peaks[idx]['c']
            for i in range(self.get_spectra_count()):
                if verbose:  print(f'({i}/{self.get_spectra_count()-1}) calculating...')
                shifts[i] = -(self.data[i].fit.peaks[idx]['c'] - self.data[ref].fit.peaks[idx]['c'])
                # if abs(self.data[ref].fit.peaks[idx]['c']) > abs(self.data[i].fit.peaks[idx]['c']):
                #     shifts[i] = -(self.data[i].fit.peaks[idx]['c'] - self.data[ref].fit.peaks[idx]['c'])
                # else:
                #     shifts[i] = (self.data[i].fit.peaks[idx]['c'] - self.data[ref].fit.peaks[idx]['c'])

                # if self.data[i].shift != 0:  # fix shift in case there was a previous shift set
                    # if self.data[i].shift_mode in roll:
                    #     shifts[i] -= self.data[i].shift*self.data[i].step
                    # else:
                    #     shifts[i] -= self.data[i].shift
        else:
            raise ValueError('mode not valid (cross-correlation, max, peak).')

        # save ranges ==========================================================
        if ranges is None:
            self._shift_calculated['ranges']  = ((self.get_min_x(), self.get_max_x()), )
        elif isinstance(ranges[0], Iterable):
            self._shift_calculated['ranges'] = ranges
        else:
            self._shift_calculated['ranges'] = (ranges, )

        # save stuff ===========================================================
        self._shift_calculated['values'] = shifts
        self._shift_calculated['ref']    = ref
        self._shift_calculated['mode']   = mode


        self._calculated_shifts = Spectrum(y=shifts)


        # finish ===============================================================
        if verbose:
            print('done!')

    def align(ref_spectrum=0, mode='cross-correlation', ranges=None, verbose=False, idx=0, bypass=False, **kwargs):
        self.calculate_shifts(ref_spectrum=0, mode='cross-correlation', ranges=None, verbose=False, idx=0, bypass=False, **kwargs)
        self.set_shift()

    def calculate_factors(self, ref_spectrum=0, mode='peak', idx=0):

        # if mode == 'peak':
        #     value = self[ref_spectrum].peaks[idx]['amp']
        #     for i in range(len(self)):
        #         self._factor_calculated['values'][i] = value/self[i].peaks[idx]['amp']
        if mode == 'peak':

            # check if fit is defined for all spectra and if spectra has the same number of peaks
            temp = self.check_fitted_peaks()
            if temp is not None:
                print(f'Different spectra have different number of fitted peaks. \nnumber of fitted peaks for each spectrum: {temp}\n')

            value = self[ref_spectrum].fit.peaks[idx]['amp']
            for i in range(len(self)):
                self._factor_calculated['values'][i] = value/self[i].fit.peaks[idx]['amp']
        elif mode == 'peak_area':
            value = self[ref_spectrum].peaks[idx]['amp']*self[ref_spectrum].fit.peaks[idx]['fwhm']
            for s in self:
                self._factor_calculated['values'][i] = value/(self[i].fit.peaks[idx]['amp']*self[i].fit.peaks[idx]['fwhm'])
        elif mode == 'value':
            raise NotImplementedError('sorry, not implemented yet.')
        elif mode == 'max':
            value = max(self[ref_spectrum].y)
            for i in range(len(self)):
                self._factor_calculated['values'][i] = value/max(self[i].y)
        else:
            raise ValueError('Invalid mode.')

        self._factor_calculated['mode'] = mode
        self._factor_calculated['ref'] = ref_spectrum
        self._factor_calculated['ref_value'] = value
        self._factor_calculated['peak_idx'] = idx

    def calculate_calib(self, start=None, stop=None, values=None, mode='peak', ranges=None, idx=None, **kwargs):
        """Calculate calibration factor (dispersion).

        It calculates the shift between spectra and infer the calibration factor.

        If mode is peak and idx is None, it will fit data with one peak
            using spectrum.fit_peak(). If idx is not None, it will look for
            defined FITTED peaks and use it to calculate calib.

        Args:
            start (number): value used for measuring the first spectrum.
            stop (number): value used for measuring the last spectrum.
            values (list): value list. It overwrites start and
                stop. Must be the same length as data.
            **kwargs: parameters to be passed to :py:func:`Spectra.calculate_shifts()`

        Returns:
            calibration factor (number)
        """
        # check number of values matches the numbre of spectra
        if values is None:
            values = np.linspace(start, stop, self.get_spectra_count())
        if self.get_spectra_count() != len(values):
            raise ValueError(f'number of values ({len(values)}) do not match the number of spectra ({self.get_spectra_count()})')

        if mode == 'peak':
            if idx is None:
                self.fit_peak(**kwargs)
                idx = 0
            centers = [0]*len(self)
            fwhms = [0]*len(self)
            for i in range(len(self)):
                centers[i] = self[i].fit.peaks[idx]['c']
                fwhms[i] = self[i].fit.peaks[idx]['fwhm']
        elif mode in cc:
            self.calculate_shifts(mode='cc', **kwargs)
            centers = (self.shift_calculated['values'] + self.shift_calculated['ref_value'])*self.step
            fwhms = [None]*len(self)
        elif mode == 'max':
            self.calculate_shifts(mode='max', **kwargs)
            centers = self.shift_calculated['values'] + self.shift_calculated['ref_value']
            fwhms = [None]*len(self)
        else:
            raise ValueError('mode not valid (cross-correlation, max, peak).')

        self._calib_calculated['mode'] = mode
        self._calib_calculated['values'] = values
        self._calib_calculated['centers'] = centers
        popt = np.polyfit(values, self._calib_calculated['centers'], deg=1)
        self._calib_calculated['popt'] = popt
        self._calib_calculated['value'] = 1/popt[0]
        self._calib_calculated['func'] = np.poly1d(popt)
        self._calib_calculated['ranges'] = ranges
        self._calib_calculated['fwhms'] = fwhms*1/popt[0]

        return 1/popt[0]





    def set_shift(self, value=None, mode=None):
        """Shift data recursively.

        Args:
            mode (string, optional): Three modes can be selected (default is 'roll'):

                1. 'x' or 'hard': y coordinates are fully preserved, while
                    x-coordinates are shifted.
                2. 'y', 'interp', or 'soft': x-coordinates are preserved, while
                    y-coordinates are interpolated with a shift.
                3. 'roll', 'rotate', or 'r': both x- and y-coordinates are preserved
                    as y-coordinates are simply rolled along the array. In this case,
                    shift values must be an integer and x-coordinates must be
                    evenly spaced.

        Returns:
            None

        Raises:
            ValueError: if a non-valid mode is selected, if shifts have not been
                calculated, or if x-coordinates are not evenly space for mode =
                'roll'.

        Warning:
            It is always better to use mode='hard' or 'roll' since y-coordinates are fully
            preserved (no interpolation). After applying a shift using the 'interp',
            one can apply a 'inverse' shift to retrieve the original data.
            The diference between the retrieved y-coordinates and the original
            data will give an ideia of the information loss caused by the interpolation.
            If statiscs are good enough, interpolating the data should not be a problem.

        See Also:
            :py:func:`Spectra.calculate_shifts`
        """
        if value is not None:
            if isinstance(value, Iterable):
                if len(value) == self.get_spectra_count():
                    for i in range(len(self)):
                        self[i].set_shift(value[i], mode=mode)
                else:
                    raise ValueError(f'Shifts must have the same length as number of spectra.\nNumber of spectra = {self.get_spectra_count()}\nNumber of shifts = {len(value)}.')

            else:
                for i in range(len(self)):
                    if mode in soft:
                        if self.monotonicity is None:
                            self.check_monotonicity()
                        if self.monotonicity != 'increasing':
                            raise ValueError('Arrays must be monotonicaly increasing.\nTip: use Spectra.fix_monotinicity()')
                    self[i].set_shift(value, mode=mode)
        else:

            # check if shifts were calculated
            if self.shift_calculated['mode'] is None:
                raise ValueError('Shift value not defined. Please, pass values or use Spectra.calculate_shifts().')

            # if all shifts are zero, does nothing
            if  all(x==0 for x in self.shift_calculated['values']):
                return

            # if shifts are different, reset x
            if all_equal(self.shift_calculated['values']) == False:
                self._x = None

            # if mode is not defined, auto pick the right one
            if mode is None:
                if self.shift_calculated['mode'] in cc:
                    mode = 'roll'
                elif self.shift_calculated['mode'] == 'max':
                    mode = 'x'
                elif self.shift_calculated['mode'] == 'peak':
                    mode = 'x'
                elif self.shift_calculated['mode'] == 'user-defined':
                    text = 'shift mode not define (please, select "roll", "soft", "hard").' +\
                            'Note that shifts were calculated elsewhere by the user, ' +\
                            'the script can has no way of knowing the most appropriate shift mode.'
                    raise ValueError(text)
                elif self.shift_calculated['mode'] in None:
                    raise ValueError('calculated shifts not defined. Use Spectra.calculate_shifts().')

            # if mode is roll, steps must be defined for all spectra
            if mode in roll:
                if self.step is None:
                    try:
                        self.check_step_x()
                    except ValueError:
                        temp = [self[i].step for i in  range(len(self))]
                        raise ValueError(f'Cannot apply roll because some spectra have no step defined.\nData must be homogeneous.\nSteps: {temp}\nUse Spectra.check_step_x()')

            # if shift were calculated by cc, one should use roll
            if self.shift_calculated['mode'] in cc and mode not in roll:
                # one could in principle do a roll even though the calculation mode
                # is cc by multiplying the shift value by the step, however,
                # the whole point of using cc is to be able to do a roll shift.
                # Therefore, when cc is used, roll shift is enforced.
                raise ValueError('shifts were calculated using cross-correlation. Only shift mode possible is "roll".')

            # if shift were calculated by peak or max, one can opt to shift via roll
            # as long as self.step is defined
            if self.shift_calculated['mode'] in ['peak', 'max'] and mode in roll:
                self.shift_calculated['values'] = np.array([int(round(x)) for x in self.shift_calculated['values']/self.step])

            # set shift
            for i in range(self.get_spectra_count()):

                # check for previous shifts
                if mode in roll:
                    value = self.shift_calculated['values'][i] + self[i].shift_roll
                elif mode in hard:
                    value = self.shift_calculated['values'][i] + self[i].shift
                elif mode in soft:
                    value = self.shift_calculated['values'][i] + self[i].shift_interp
                    # check monotonicity
                    if self.monotonicity is None:
                        self.check_monotonicity()
                    if self.monotonicity != 'increasing':
                        raise ValueError('Arrays must be monotonicaly increasing.\nTip: use Spectra.fix_monotinicity()')

                self.data[i].set_shift(value=value, mode=mode)

        self._restart_shift_attr()
        # self._restart_check_attr()
        self._restart_sum_attr()
        self._restart_calib_attr()

    def set_offset(self, value=None):
        if value is not None:
            if isinstance(value, Iterable):
                if len(value) == self.get_spectra_count():
                    for i in range(len(self)):
                        self[i].set_offset(value[i])
                else:
                    raise ValueError(f'value must have the same length as number of spectra.\nNumber of spectra = {self.get_spectra_count()}\nNumber of values = {len(value)}.')
            else:
                for i in range(len(self)):
                    self[i].set_offset(value)

            self._restart_sum_attr()

    def set_factor(self, value=None):
        if value is not None:
            if isinstance(value, Iterable):
                if len(value) == self.get_spectra_count():
                    for i in range(len(self)):
                        self[i].set_factor(value[i])
                else:
                    raise ValueError(f'value must have the same length as number of spectra.\nNumber of spectra = {self.get_spectra_count()}\nNumber of values = {len(value)}.')
            else:
                for i in range(len(self)):
                    self[i].set_factor(value)

            self._restart_sum_attr()
        else:
            # check if shifts were calculated
            if self.factor_calculated['mode'] is None:
                raise ValueError('Factor value not defined. Please, pass values or use Spectra.calculate_factor().')

            # if all shifts are zero, does nothing
            if  all(x==1 for x in self.factor_calculated['values']):
                return

            # set shift
            for i in range(self.get_spectra_count()):
                value = self.data[i].factor*self.factor_calculated['values'][i]
                self.data[i].set_factor(value=value)

            self._restart_sum_attr()

    def set_calib(self, value):
        """Calibrate data (from x-coordinates to energy).

        Args:
            value (number): dispersion of the diffraction grating in
                units of [energy/(unit of the x axis)].

        Returns:
            None
        """
        if value is not None:
            if isinstance(value, Iterable):
                if len(value) == self.get_spectra_count():
                    for i in range(len(self)):
                        self[i].set_calib(value[i])
                else:
                    raise ValueError(f'value must have the same length as number of spectra.\nNumber of spectra = {self.get_spectra_count()}\nNumber of values = {len(value)}.')
            else:
                for i in range(len(self)):
                    self[i].set_calib(value)

            self._restart_sum_attr()
            self._restart_check_attr()
            self._restart_shift_attr()
            self._restart_calib_attr()



    def calculate_sum(self):
        """Returns Spestrum object with the sum of all spectra.

        It also save it at :py:attr:`Spectra.sum`.

        Returns:
            :py:class:`Spectra` object.

        Note:
            All spectra have to have the same x-coordinates. This is verified
            before summing up the spectra.
        """

        if self.x is None:
            self.check_same_x()

        y = np.zeros(len(self.x))
        for i in range(self.get_spectra_count()):
            y += self.data[i].y

        self._sum = Spectrum(x=self.x, y=y)
        return self._sum

    def save(self, folderpath='./', prefix='spectrum_', suffix='.dat', zfill=None, delimiter=', ', header='', fmt='%.18e', verbose=False):
        r"""Saves all spectra in a folder.

        Args:
            folderpath (string or pathlib.Path, optional): path to directory.
                Default is the current directory.
            prefix (string, optional): prefix used for naming the files.
            suffix (string, optional): suffix used for naming the files.
            zfill (int, optional): number of digits for file numbering. If `None`,
                zfill will be determined.
            delimiter (str, optional): The string used to separate values.
                If whitespaces are used, consecutive whitespaces act as delimiter.
                Use ``\\t`` for tab. The default is comma (', ').
            header (string, optional): text to add at the beginning of each file.
                Use ``\n`` for new line. Comment flag ('# ') is added automatically.
            fmt (string, or list, optional): format for saving data.
                If string, the value is used for x- and y-coordinates. If tuple
                of strings, the first string is used for x-coordinates and the
                second for y-coordinates.

                    fmt = (%[flag]width[.precision]specifier)

                * `flag` can be:
                    1. '-' for left justify
                    2. '+', puts + or - in front of the numbers
                    3. '0' to Left pad the number with zeros instead of space (see width).

                * `width` is the minimum number of characters to be printed.

                * `precision` is the number of significant digits.

                * `specifier` is the type of notation. Tipically, either 'e' for
                scientific notation of 'f' for decimal floating point.

                * a common `fmt` strings is: '%.3f' for 3 decimal places.

                *  for more information see `np.savetxt <https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html?highlight=savetxt#numpy.savetxt>`_ documentation::

            verbose (bool, optional): turn verbose on and off. Default is `False`.

        Returns:
            None

        """
        folderpath = Path(folderpath)

        if zfill is None:
            zfill = n_digits(self.get_spectra_count()-1)[0]

        for idx, s in enumerate(self.data):
            filename = f'{prefix}' + f'{idx}'.zfill(zfill) + f'{suffix}'
            if verbose:  print(f'({i}/{self.get_spectra_count()-1}) saving {filename}')
            s.save(filepath=folderpath/filename, delimiter=delimiter, header=header, fmt='%.18e')
        if verbose: print('Done!')

    def interp(self, start=None, stop=None, num=None, step=None, x=None, ):
        """Interpolate data.

        Args:
            x (list or array, optional): The x-coordinates at which to
                evaluate the interpolated values. This overwrites all other arguments.
            start (number, optional): The starting value of the sequence. If `None`,
                the minium x value will be used.
            stop (number, optional): The end value of the sequence. If `None`,
                the maximum x value will be used.
            num (int, optional): Number of samples to generate.
            step (number, optional): Spacing between values. This overwrites ``num``.

        Returns:
            None
        """
        if x is None:
            if self.x is None:
                if start is None:
                    start = max([min(s.x) for s in self.data])
                if stop is None:
                    stop = min([max(s.x) for s in self.data])
                if num is None:
                    num = max([len(s.x) for s in self.data])
            else:
                if start is not None or stop is not None or num is not None or step is not None:
                    if start is None:
                        start = max(self.x)
                    if stop is None:
                        stop = min(self.x)
                    if num is None:
                        num = len(self.x)
                else:
                    return

        for s in self.data:
            s.interp(x=x, start=start, stop=stop, num=num, step=step)

        self._restart_sum_attr()
        self.check_same_x()

    def compress(self, ranges):
        """Slices data.

        Args:
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.

        Returns:
            None
        """
        if isinstance(ranges, Iterable):
            if isinstance(ranges[0], Iterable):
                for s in self.data:
                    s.compress(ranges=ranges)
                self._restart_sum_attr()
                self._restart_check_attr()
            else:
                self.crop(start=ranges[0], stop=ranges[0])

    def get_min_x(self):
        if self.x is None:
            return min([min(s.x) for s in self.data])
        else:
            return min(self.x)

    def get_max_x(self):
        if self.x is None:
            return max([max(s.x) for s in self.data])
        else:
            return max(self.x)

    def crop(self, start=None, stop=None):
        if self.x is None:
            # print('here')
            if start is None:
                start = min(self.data[0].x)
                for s in self.data:
                    temp = min(s.x)
                    if temp > start:
                        start = temp
            if stop is None:
                stop = max(self.data[0].x)
                for s in self.data:
                    temp = max(s.x)
                    if temp < stop:
                        stop = temp
            for s in self.data:
                s.crop(start, stop)
            self._restart_sum_attr()
            self._restart_check_attr()
        elif start is not None or stop is not None:
            for s in self.data:
                s.crop(start, stop)
            self._restart_sum_attr()
            self._restart_check_attr()

    def floor(self, value=None, n=8, ranges=None):
        """Sets zero value for y-coordinates (shifts data verticaly).

        Args:
            value (number, optional): x-coordinate to set y-coordinate value as
                zero. If `None`, the minium value of x is used.
            n (int, optional): number of data points to average and setting to
                zero.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                ranges overwrites value and n.

        Returns:
            None
        """
        for s in self.data:
            s.floor(value=value, n=n, ranges=ranges)

    def plot(self, ax=None, vertical_increment=0, shift=0, factor=1, ranges=None, **kwargs):
        """Plot spectra.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            vertical_increment (number, optional): defines
                the vertical offset between each curve in percentage of the
                y-coordinates variation (delta y). Default is 0.
            shift (float or int): horizontal shift value. Default is 0.
            factor (float or int): multiplicative factor. Default is 1.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            None
        """
        if ax is None:
            ax = plt

        # percentage wise increment ====================
        temp = [0]*self.get_spectra_count()
        for i in range(self.get_spectra_count()):
            temp[i] = self.data[i].get_deltay()
        vertical_increment = max(temp)*factor*vertical_increment/100

        for i in range(self.get_spectra_count()):
            self.data[i].plot(ax=ax, offset=vertical_increment*-i, shift=shift, factor=factor, ranges=ranges, label=i, **kwargs)

        plt.legend()

    def plot_shift_ranges(self, ax=None, **kwargs):

        if self.shift_ranges is None:
            raise ValueError('Shift range not defined. Use calculate_shifts().')

        if ax is None:
            ax = plt

        if 'lw' not in kwargs and 'linewidth' not in kwargs:
            kwargs['lw'] = 2
        if 'zorder' not in kwargs:
            kwargs['zorder'] = 10
        kwargs1 = copy.deepcopy(kwargs)
        kwargs2 = copy.deepcopy(kwargs)
        if 'color' not in kwargs1 and 'c' not in kwargs1:
            kwargs1['color'] = 'green'
        if 'color' not in kwargs2 and 'c' not in kwargs2:
            kwargs2['color'] = 'red'

        for r in self.shift_ranges:
            ax.axvline(r[0], **kwargs1)
            ax.axvline(r[1], **kwargs2)

    def plot_detected_peaks(self, ax=None, vertical_increment=0, shift=0, factor=1, **kwargs):

        if ax is None:
            ax = plt

        # check all spectra have peaks =================================
        for i in range(self.get_spectra_count()):
            if self.data[i].peaks is None:
                raise ValueError(f'find_peaks needs to be ran in Spectrum {i}.')

        # percentage wise increment ====================
        temp = [0]*self.get_spectra_count()
        for i in range(self.get_spectra_count()):
            temp[i] = self.data[i].get_deltay()
        vertical_increment = max(temp)*factor*vertical_increment/100

        for i in range(self.get_spectra_count()):
            self.data[i].plot_detected_peaks(ax=ax, offset=vertical_increment*-i, shift=shift, factor=factor, label=f'peaks for spectrun {i}', **kwargs)
        plt.legend()

    def plot_fit(self, ax=None, vertical_increment=0, shift=0, factor=1, ranges=None, **kwargs):
        if ax is None:
            ax = plt

        # percentage wise increment ====================
        temp = [0]*self.get_spectra_count()
        for i in range(self.get_spectra_count()):
            temp[i] = self.data[i].get_deltay()
        vertical_increment = max(temp)*factor*vertical_increment/100

        for i in range(self.get_spectra_count()):
            self.data[i].fit.plot(ax=ax, offset=vertical_increment*-i, shift=shift, factor=factor, ranges=ranges, **kwargs)

    def plot_guess(self, ax=None, vertical_increment=0, shift=0, factor=1, ranges=None, **kwargs):
        if ax is None:
            ax = plt

        # percentage wise increment ====================
        temp = [0]*self.get_spectra_count()
        for i in range(self.get_spectra_count()):
            temp[i] = self.data[i].get_deltay()
        vertical_increment = max(temp)*factor*vertical_increment/100

        for i in range(self.get_spectra_count()):
            self.data[i].guess.plot(ax=ax, offset=vertical_increment*-i, shift=shift, factor=factor, ranges=ranges, **kwargs)

    def plot_calib(self, ax=None, show_fit=True, **kwargs):
        """Plot position of elastic peak (or data shift) as function of energy.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            show_fit (bool, optional): show fit.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes
        """
        if ax is None:
            ax = plt

        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'
        if 'linewidth' not in kwargs and 'lw' not in kwargs:
            kwargs['linewidth'] = 0

        ax.plot(self.calib_calculated['values'], self.calib_calculated['centers'], **kwargs)
        if show_fit:
            x = np.linspace(self.calib_calculated['values'][0], self.calib_calculated['values'][-1], self.get_spectra_count()*10)
            y = self.calib_calculated['func'](x)
            ax.plot(x, y, color='black')

    def plot_map(self, start=None, stop=None, values=None, orientation='x', ranges=None, vlines=False, **kwargs):
        """Plot position of elastic peak (or data shift) as function of energy.

        x-coordinates must be the same for all spectra.

        Args:
            start (number): value (energy or momentum) used for measuring the first spectrum.
            stop (number): value (energy or momentum) used for measuring the last spectrum.
                Energies are them extrapolated from start to stop in
                steps of one.
            values (list): values list. It overwrites start and
                stop. Must be the same length as data.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                ranges overwrites value and n.
            vlines (bool, optional): show separation lines between spectra.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            matplotlib.axes
        """
        if self.x is None:
            self.check_same_x()

        rearrange_ys = False

        if values is None:
            if start is None or stop is None:
                values = np.arange(0, self.get_spectra_count())
            else:
                values = np.linspace(start, stop, self.get_spectra_count())
        else:
            if self.get_spectra_count() != len(values):
                raise ValueError(f'number of values ({len(values)}) do not match the number of spectra ({self.get_spectra_count()})')
            # check values monotonicity
            if np.all(np.diff(values) > 0) == True or np.all(np.diff(self.x) < 0) == True:
                pass
            else:
                values, ID, counts = np.unique(values, return_inverse=True, return_counts=True)
                rearrange_ys = True


        x, ys = self._gather_ys(ranges=ranges)

        # print(values)
        # print(counts)
        # print(ID)

        # if values is not monotonicaly
        if rearrange_ys:
            ys2 = copy.deepcopy(ys)
            for i, id in enumerate(ID):
                repeated = [k for k, x in enumerate(ID) if x == id]
                if len(repeated) > 1:
                    if min(repeated) == i:
                        ys2[:, id] = np.zeros(len(ys2[:, id]))
                        for j in repeated:
                            ys2[:, id] += ys[:, j]
                else:
                    ys2[:, id] = ys[:, i]
            n_del = sum(counts)-len(counts)
            ys = ys2[:, 0:len(ys2)-n_del]
            # print(n_del)
            # print(np.shape(ys))
            # print(len(counts))
            ys /= counts


        # print(len(values))
        # print(len(x))
        # print(len(ys[0, :]))
        # print(len(ys[:, 0]))


        if 'vmin' not in kwargs:
            kwargs['vmin'] = min([min(k) for k in ys])
        if 'vmax' not in kwargs:
            kwargs['vmax'] = max([max(k) for k in ys])
        if 'shading' not in kwargs:
            kwargs['shading'] = 'nearest'
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'

        # for i in range(self.get_spectra_count()):
        #     plt.plot(self.x, ys[:,i])
        fig = plt.figure()
        ax = fig.add_subplot(111)

        if orientation == 'x':
            ax.pcolormesh(values, x, ys, **kwargs)
        elif orientation == 'y':
            ax.pcolormesh(x, values,  ys.transpose(), **kwargs)
        else:
            raise ValueError('orientation should be either `x` or `y`.')

        if vlines:
            d = (np.diff(values)/2)[0]
            ax.vlines(values-d, min(self.x), max(self.x),)
            for i, x in enumerate(values):
                ax.text(x-d, max(self.x)-(max(self.x)-min(self.x))*0.1, i,)

        return ax


    def flip(self):
        for s in self:
            s.flip()
        self._restart_sum_attr()
        self._restart_check_attr()
        self._restart_shift_attr()
        self._restart_calib_attr()

    def concatenate(self):
        x = np.concatenate([s.x for s in self.data])
        y = np.concatenate([s.y for s in self.data])
        return Spectrum(x=x, y=y)

















# %%
