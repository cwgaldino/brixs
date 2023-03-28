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
import json

# backpack
from .backpack.filemanip import filelist
from .backpack.arraymanip import sort, all_equal
from .backpack.model_functions import voigt_fwhm, dirac_delta
from .backpack.figmanip import n_digits

# BRIXS
import brixs as br

# common definitions ===========================================================
relative = ['relative', 'r', 'rel']
absolute = ['a', 'abs', 'absolute']

# %% suport functions ==========================================================
def build_model_str(**kwargs):
    """Returns instructions for building peak functions.

    Args:
        asymmetry (bool): if True, the returned peak function will require two
            fwhm values, one for each half of the peak.
        fixed_m (False or number): m is the amount of lorentzian
            contribution for a peak. If False, m will be a required as an input
            argument on the returned peak function.
        idx (int, optional): number to be inprinted in the string

    Returns:
        function f(x), function f(x) as string, argument list as string
    """
    if 'asymmetry' in kwargs:
        asymmetry = kwargs.pop('asymmetry')
    else:
        asymmetry = False
    if 'idx' in kwargs:
        idx = kwargs.pop('idx')
    else:
        idx = 0

    # check kwargs
    delete = []
    for key in list(kwargs.keys()):
        if key not in ['amp', 'c', 'fwhm', 'm', 'fwhm1', 'm1', 'fwhm2', 'm2']:
            raise ValueError('Fixed parameter "' + str(key) + '" not valid.\nValid keys: "amp", "c", "fwhm", "m", "fwhm1", "m1", "fwhm2", "m2". ')
        else:
            if kwargs[key] is None:
                del kwargs[key]

    # check asymmetry
    # if 'fwhm1' in kwargs or 'fwhm2' in kwargs:
    #     if 'fwhm2' in kwargs and 'fwhm2' in kwargs:
    #         asymmetry = True
    #     else:
    #         raise ValueError('If fwhm1 or fwhm2 are defined, both must be defined.')
    # else:
    #     asymmetry = False


    # function string
    if asymmetry:
        # f_str = f'np.heaviside(x-c_{idx}, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w1_{idx}, m1_{idx}) + np.heaviside(c_{idx}-x, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w2_{idx}, m2_{idx}) + dirac_delta(x, amp_{idx}, c_{idx})'
        f_str = f'np.heaviside(c_{idx}-x, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w1_{idx}, m1_{idx}) + np.heaviside(x-c_{idx}, 0)*voigt_fwhm(x, amp_{idx}, c_{idx}, w2_{idx}, m2_{idx}) + dirac_delta(x, amp_{idx}, c_{idx})'
    else:
        f_str = f'voigt_fwhm(x, amp_{idx}, c_{idx}, w_{idx}, m_{idx})'

    # args string
    args_str = ''
    if 'amp' in kwargs:
        f_str = f_str.replace(f'amp_{idx}', str(kwargs['amp']))
    else:
        args_str += f'amp_{idx}, '
    if 'c' in kwargs:
        f_str = f_str.replace(f'c_{idx}', str(kwargs['c']))
    else:
        args_str += f'c_{idx}, '
    if asymmetry:
        if 'fwhm1' in kwargs:
            f_str = f_str.replace(f'w1_{idx}', str(kwargs['fwhm1']))
        else:
            args_str += f'w1_{idx}, '
        if 'm' in kwargs:
            if 'm1' not in kwargs:
                f_str = f_str.replace(f'm1_{idx}', str(kwargs['m']))
            if 'm1' in kwargs:
                f_str = f_str.replace(f'm1_{idx}', str(kwargs['m1']))
        elif 'm1' in kwargs:
            f_str = f_str.replace(f'm1_{idx}', str(kwargs['m1']))
        else:
            args_str += f'm1_{idx}, '
        if 'fwhm2' in kwargs:
            f_str = f_str.replace(f'w2_{idx}', str(kwargs['fwhm2']))
        else:
            args_str += f'w2_{idx}, '
        if 'm' in kwargs:
            if 'm2' not in kwargs:
                f_str = f_str.replace(f'm2_{idx}', str(kwargs['m']))
            if 'm2' in kwargs:
                f_str = f_str.replace(f'm2_{idx}', str(kwargs['m2']))
        elif 'm2' in kwargs:
            f_str = f_str.replace(f'm2_{idx}', str(kwargs['m2']))
    else:
        if 'fwhm' in kwargs:
            f_str = f_str.replace(f'w_{idx}', str(kwargs['fwhm']))
        else:
            args_str += f'w_{idx}, '
        if 'm' in kwargs:
            f_str = f_str.replace(f'm_{idx}', str(kwargs['m']))
        else:
            args_str += f'm_{idx}, '

    return f_str, args_str[:-2]

def build_model(**kwargs):
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
    f_str, args_str = build_model_str(**kwargs)
    # print(f_str)
    model_str = f'lambda x, {args_str}: {f_str}'
    return eval(model_str)


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
            asimetric peaks).
        m (number, optional): factor from 0 to 1 of the lorentzian amount.
            Default is 0 (fully gaussian peak). If peak is assymetric (m1 and m2
            are defined), m will be set the average m value.
        m1 (number, optional): m of the first half of the peak (for
            assimetric peaks). Default is 0.
        m2 (number, optional): m of the second half of the peak (for
            assimetric peaks). Default is 0.
        shift (number): initial shift value (does not affect the data).
        calib (number): initial calibration value (does not affect the data).
        offset( number): initial offset value (does not affect the data).
        factor (number): initial multiplicative factor (does not affect the data).


    Attributes:
        area (number): peak area.
        asymmetry (bool): True if peak is assymetric. This is defined when the
            object is created. If asymmetry=False, and fhwm1 (or fwhm2) is set,
            then asymmetry is set to True. If asymmetry=True, it cannot reverse
            to False, and another object must be created from scratch.
        store (dictionary):
        shift (number): shift value (value will be added to x-coordinates).
        calib (number): calibration value (x-coordinates will be multiplied by this value).
        offset (number): offset value (value will be added to y-coordinates).
        factor (number): multiplicative factor (y-coordinates will be multiplied by this value).
    """

    def __init__(self, *args, **kwargs):
        # core
        self._store = {'amp':None,
                      'area':None,
                      'c':None,
                      'fwhm':None,
                      'fwhm1':None,
                      'fwhm2':None,
                      'm':None,
                      'm1':None,
                      'm2':None,
                      }
        self.bounds = {'amp': [-np.inf, np.inf],
                      'c':    [-np.inf, np.inf],
                      'fwhm': [0, np.inf],
                      'fwhm1':[0, np.inf],
                      'fwhm2':[0, np.inf],
                      'm':    [0, 1],
                      'm1':   [0, 1],
                      'm2':   [0, 1],
                      }
        self.error = {'amp':  0,
                      'area': 0,
                      'c':    0,
                      'fwhm': 0,
                      'fwhm1':0,
                      'fwhm2':0,
                      'm':    0,
                      'm1':   0,
                      'm2':   0,
                      }

        # fit instructions
        if 'asymmetry' in kwargs:
            self._asymmetry = kwargs.pop('asymmetry')
        else:
            self._asymmetry = False
        if 'fixed' in kwargs:
            self._fixed = kwargs.pop('fixed')
        else:
            self._fixed = ['m']

        # modifiers
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
        if any([item not in self._store for item in kwargs.keys()]):
            for key in kwargs:
                if key not in self._store:
                    raise AttributeError(f'{key} is not a valid attribute.')

        # check kwargs
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
                kwargs['fwhm1'] = kwargs['fwhm']/2
                kwargs['fwhm2'] = kwargs['fwhm']/2
        if self.asymmetry:
            if 'm1' not in kwargs:
                kwargs['m1'] = 0
            else:
                self._check_m(kwargs['m1'])
            if 'm2' not in kwargs:
                kwargs['m2'] = 0
            else:
                self._check_m(kwargs['m1'])
            kwargs['m'] = (kwargs['m1']+kwargs['m2'])/2
        else:
            if 'm' not in kwargs:
                kwargs['m'] = 0
            else:
                self._check_m(kwargs['m'])

        # calculate area
        if self.asymmetry:
            c1 = kwargs['m1']/np.pi + (1-kwargs['m1'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            c2 = kwargs['m2']/np.pi + (1-kwargs['m2'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            kwargs['area'] = (kwargs['amp']*kwargs['fwhm1']*c1 + kwargs['amp']*kwargs['fwhm2']*c2)/2
        else:
            c = kwargs['m']/np.pi + (1-kwargs['m'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            kwargs['area'] = kwargs['amp']*kwargs['fwhm']*c

        # save kwargs
        self._store.update(dict(**kwargs))


    def __str__(self):
        # return str({name:self._store[name] for name in self._store if self._store[name] is not None})[1:-1].replace(', ', '\n')
        # return str({name:self._store[name] for name in self._store if self._store[name] is not None})[1:-1].replace(', ', '\n')
        return str(self._store)[1:-1].replace(', ', '\n')
        # return str(self._store).replace('}, ', '\n ')

    def __repr__(self):

        # return str({name:self._store[name] for name in self._store if self._store[name] is not None})
        return str(self._store)[1:-1].replace(', ', '\n')
        # return str({i:val for i, val in enumerate(self._store)})[1:-1].replace('}, ', '}\n')
        # return str(self._store).replace('}, ', '\n ')
        # return str(self._store)

    def __getitem__(self, name):
        return self._store[self._check_name(name)]

    def __setitem__(self, name, value):
        if name not in self._store:
            raise KeyError(f'{name} is not defined as a peak parameter and cannot be created.')
        if name == 'amp':
            self._check_amp(value)
            self._store['amp'] = value
        if name == 'c':
            self._store['c'] = value
        if name == 'area':
            raise ValueError('Area is read-only parameter. Change amp instead.')
        if name == 'fwhm':
            if self.asymmetry:
                raise ValueError('Cannot change fwhm value when asymmetry = True.\nChange fwhm1 and fhwm2 instead.\nOr change asymmetry to False.')
            else:
                self._check_fwhm(value)
                self._store['fwhm'] = value
        if name == 'fwhm1':
            if self.asymmetry:
                if value is None:
                    print('Peak asymmetry set to False.')
                    self.asymmetry = False
                else:
                    self._check_fwhm(value)
                    self._store['fwhm1'] = value
                    self._store['fwhm'] = self._store['fwhm1']+self._store['fwhm2']
            else:
                if value is None:
                    self._store['fwhm1'] = None
                else:
                    print('Peak asymmetry set to True.')
                    self.asymmetry = True

                    self._check_fwhm(value)
                    self._store['fwhm1'] = value
                    if self._store['fwhm']-self._store['fwhm1'] >= 0:
                        self._store['fwhm2'] = self._store['fwhm']-self._store['fwhm1']
                    else:
                        self._store['fwhm2'] = 0
                        self._store['fwhm'] = self._store['fwhm1']+self._store['fwhm2']
        if name == 'fwhm2':
            if self.asymmetry:
                if value is None:
                    print('Peak asymmetry set to False.')
                    self.asymmetry = False
                else:
                    self._check_fwhm(value)
                    self._store['fwhm2'] = value
                    self._store['fwhm'] = self._store['fwhm1']+self._store['fwhm2']
            else:
                if value is None:
                    self._store['fwhm2'] = None
                else:
                    print('Peak asymmetry set to True.')
                    self._asymmetry = True

                    self._check_fwhm(value)
                    self._store['fwhm2'] = value
                    if self._store['fwhm']-self._store['fwhm2'] >= 0:
                        self._store['fwhm1'] = self._store['fwhm']-self._store['fwhm2']
                    else:
                        self._store['fwhm1'] = 0
                        self._store['fwhm'] = self._store['fwhm1']+self._store['fwhm2']
        if name == 'm':
            if self.asymmetry:
                raise ValueError('cannot change m value when asymmetry = True.\nChange m1 and m2 instead.\nOr change asymmetry to False.')
            else:
                self._check_m(value)
                self._store['m'] = value
        if name == 'm2':
            if self.asymmetry:
                self._check_m(value)
                self._store['m2'] = value
                self._store['m'] = (self._store['m1']+self._store['m2'])/2
            else:
                raise ValueError('Cannot set m2 because asymmetry is False.')
        if name == 'm1':
            if self.asymmetry:
                self._check_m(value)
                self._store['m1'] = value
                self._store['m'] = (self._store['m1']+self._store['m2'])/2
            else:
                raise ValueError('Cannot set m1 because asymmetry is False.')

        # fix area
        self.calculate_area()

    def __delitem__(self, key):
        raise AttributeError('itens cannot be deleted')

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        # return len(self._store)
        return len({name:self._store[name] for name in self._store if self._store[name] is not None})


    def _check_name(self, name):
        if name not in self._store:
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

    def _check_fixed(self, value):
        for key in value:
            if key not in ['amp', 'c', 'fwhm', 'm']:
                raise ValueError('parameter "' + str(key) + '" not valid.\nValid keys: "amp", "c", "fwhm", "m". ')


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
        self.set_shifts(value)
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        return self._offset
    @offset.setter
    def offset(self, value):
        self.set_offsets(value)
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        return self._factor
    @factor.setter
    def factor(self, value):
        self.set_factors(value)
    @factor.deleter
    def factor(self):
            raise AttributeError('Cannot delete object.')

    @property
    def asymmetry(self):
        return self._asymmetry
    @asymmetry.setter
    def asymmetry(self, value):
        # set asymmetry
        if isinstance(value, bool):
            if self.asymmetry == value:
                return
            else:
                self._asymmetry = value
        else:
            raise TypeError('value must be bool (True or False).')
        # set everything else
        if self._asymmetry:
            self._store['fwhm1'] = self._store['fwhm']/2
            self._store['fwhm2'] = self._store['fwhm']/2
            self._store['m1']    = self._store['m']
            self._store['m2']    = self._store['m']
        else:
            self._store['fwhm1'] = None
            self._store['fwhm2'] = None
            self._store['m1']    = None
            self._store['m2']    = None
        # raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @asymmetry.deleter
    def asymmetry(self):
            raise AttributeError('Cannot delete object.')

    @property
    def fixed(self):
        return self._fixed
    @fixed.setter
    def fixed(self, value):
        self._check_fixed(value)
        self._fixed = value
    @fixed.deleter
    def fixed(self):
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


    def _string2save(self):
        """String used for saving peak to a file."""
        temp = {'_store': self._store,
                'error': self.error,
                'bounds': self.bounds,
                'asymmetry': self.asymmetry,
                'fixed': self.fixed,
                '_shift': self._shift,
                '_calib': self._calib,
                '_offset': self._offset,
                '_factor': self._factor}
        return json.dumps(temp, indent=4, sort_keys=False)

    def save(self, filepath='./Untitled.txt', check_overwrite=False):
        r"""Save peak to a text file. Wrapper for `json.dumps()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user want to overwrite file.

        Returns:
            None

        .. _json.dumps(): https://docs.python.org/3/library/json.html#json.dumps
        """
        filepath = Path(filepath)

        # check overwrite
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query('File already exists!! Do you wish to ovewrite it?', 'yes') == True:
                        pass
                    else:
                        warnings.warn('File not saved because user did not allow overwriting.')
                        return
                else:
                    warnings.warn('filepath is pointing to a folder. Saving file as Untitled.txt')
                    filepath = filepath/'Untitled.txt'

        with open(str(filepath), 'w') as file:
            file.write(self._string2save())

    def _obj_decode(self, obj):
        obj['_store'].pop('area')

        if obj['asymmetry']:
            obj['_store'].pop('fwhm')
            obj['_store'].pop('m')
        else:
            obj['_store'].pop('fwhm1')
            obj['_store'].pop('fwhm2')
            obj['_store'].pop('m1')
            obj['_store'].pop('m2')

        self.asymmetry = obj['asymmetry']
        for parameter in obj['_store']:
            # print(parameter)
            self[parameter] = obj['_store'][parameter]
        self.error = obj['error']
        self.bounds = obj['bounds']
        self.fixed = obj['fixed']
        self._shift = obj['_shift']
        self._calib = obj['_calib']
        self._offset = obj['_offset']
        self._factor = obj['_factor']

    def load(self, filepath):
        """Load peak from a text file. Wrapper for `json.load()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first decompressed.

        Returns:
            None

        .. _json.load(): https://docs.python.org/3/library/json.html#json.load
        """
        filepath = Path(filepath)

        with open(str(filepath), 'r') as file:
            obj = json.load(file)
        self._obj_decode(obj)


    def set_calib(self, value, type_='relative'):
        """Set calibration value.

        Args:
            value (number): calibration value (x-coordinates will be multiplied
                by this value).

        Returns:
            None
        """
        if type_ in relative:
            value = self.calib * value

        if self.calib != value:
            if self.calib != 1:
                self._store['c'] = self._store['c']*self.calib**-1
                self._store['fwhm'] = abs(self._store['fwhm']*self.calib**-1)
                if self.asymmetry:
                    self._store['fwhm1'] = abs(self._store['fwhm1']*self.calib**-1)
                    self._store['fwhm2'] = abs(self._store['fwhm2']*self.calib**-1)
            if value != 1:
                self._store['c'] = self._store['c']*value
                self._store['fwhm'] = abs(self._store['fwhm']*value)
                if self.asymmetry:
                    self._store['fwhm1'] = abs(self._store['fwhm1']*value)
                    self._store['fwhm2'] = abs(self._store['fwhm2']*value)
            self._calib = value

            # fix area
            self.calculate_area()

    def set_shifts(self, value, type_='relative'):
        """Set shift value.

        Args:
            value (float or int): shift value (value will be added to x-coordinates).

        Returns:
            None
        """
        if type_ in relative:
            value = self.shift + value

        if self.shift != value:
            if self.shift != 0:
                self._store['c'] = self._store['c']-self.shift
            if value != 0:
                self._store['c'] = self._store['c']+value
            self._shift = value

    def set_offsets(self, value, type_='relative'):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).

        Returns:
            None
        """
        if type_ in relative:
            value = self.offset + value

        if self.offset != value:
            if self.offset != 0:
                self._store['amp'] = self._store['amp']-self.offset
            if value != 0:
                self._store['amp'] = self._store['amp']+value
            self._offset = value

    def set_factors(self, value, type_='relative'):
        """Set y multiplicative factor.

        Args:
            value (number): multiplicative factor (y-coordinates will be
                multiplied by this value).

        Returns:
            None
        """
        if type_ in relative:
            value = self.factor * value

        if self.factor != value:
            if self.factor != 1:
                self._store['amp'] = self._store['amp']*self.factor**-1
            if value != 1:
                self._store['amp'] = self._store['amp']*value
            self._factor = value

            # fix area
            self.calculate_area()


    def calculate_area(self):
        """Calculate peak area.

        Returns:
            None
        """
        if self.asymmetry:
            c1 = self._store['m1']/np.pi + (1-self._store['m1'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            c2 = self._store['m2']/np.pi + (1-self._store['m2'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            self._store['area'] = (self._store['amp']*self._store['fwhm1']*c1 + self._store['amp']*self._store['fwhm2']*c2)/2
        else:
            c = self._store['m']/np.pi + (1-self._store['m'])*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            self._store['area'] = self._store['amp']*self._store['fwhm']*c


    def _find_suitable_x(self):
        return np.arange(self['c']-10*self['fwhm'], self['c']+10*self['fwhm'], self['fwhm']/20)

    def calculate_spectrum(self, x=None):
        """Return peak curve.

        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x, with at least 20 points within the peak,
                will be constructed.

        Returns:
            :py:class:`Spectrum`.
        """
        if x is None:
            x = self._find_suitable_x()
        # f = build_model(amp=self['amp'], c=self['c'],  fwhm=self['fwhm'], m=self['m'], fwhm1=self['fwhm1'], fwhm2=self['fwhm2'], m1=self['m1'], m2=self['m2'])
        f = self.build_model(fixed=['amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm1', 'm2', 'm'])
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
            if br.settings.FIGURE_FORCE_NEW_WINDOW:
                figure()

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


    def set_bounds(self, **kwargs):
        """Set percentage wise bounds.

        Args:


        Return:
            None
        """
        if 'type_' in kwargs:
            type_ = kwargs.pop('type_')
            if type_.startswith('a'):
                type_ = 'additive'
            elif type_.startswith('m'):
                type_ = 'multiplicative'
            elif type_.startswith('p'):
                type_ = 'percentage'
            else:
                raise ValueError(f'type_={type} not valid.\nValid type_s are: additive and multiplicative')

        if any([item not in self._store for item in kwargs.keys()]):
            for key in kwargs:
                if key not in self._store:
                    raise AttributeError(f'{key} is not a valid attribute.')

        for parameter in ['amp', 'c']:
            if parameter in kwargs:
                if kwargs[parameter] is None:
                    self.bounds[parameter] = [-np.inf, np.inf]
                elif isinstance(kwargs[parameter], Iterable):
                    if type_ == 'additive':
                        self.bounds[parameter] = [self[parameter]-kwargs[parameter][0], self[parameter]+kwargs[parameter][-1]]
                    elif type_ == 'multiplicative':
                        assert self[parameter] != 0,      f'bounds cannot be set via multiplicative factor for parameter "{parameter}", because its value is zero'
                        if self[parameter] < 0:
                            assert kwargs[parameter][1] <= 1, f'For negative values, the second bound multiplicative factor must be less than 1./nParameter: {parameter}/nParameter value: {self[parameter]}/nBound multiplicative factor: {kwargs[parameter]}'
                            assert kwargs[parameter][0] >= 1, f'For negative values, the first bound multiplicative factor must be higher than 1./nParameter: {parameter}/nParameter value: {self[parameter]}/nBound multiplicative factor: {kwargs[parameter]}'
                            self.bounds[parameter] = [self[parameter]*kwargs[parameter][0], self[parameter]*kwargs[parameter][-1]]
                        else:
                            assert kwargs[parameter][0] <= 1, f'For positive values, the first bound multiplicative factor must be less than 1./nParameter: {parameter}/nParameter value: {self[parameter]}/nBound multiplicative factor: {kwargs[parameter]}'
                            assert kwargs[parameter][1] >= 1, f'For positive values, the second bound multiplicative factor must be higher than 1./nParameter: {parameter}/nParameter value: {self[parameter]}/nBound multiplicative factor: {kwargs[parameter]}'
                            self.bounds[parameter] = [self[parameter]*kwargs[parameter][0], self[parameter]*kwargs[parameter][-1]]
                    else:
                        assert self[parameter] != 0,      f'bounds cannot be set via percentage wise factor because the parameter "{parameter}" is zero'
                        if self[parameter] < 0:
                            self.bounds[parameter] = [self[parameter]+self[parameter]*kwargs[parameter][0]/100, self[parameter]-self[parameter]*kwargs[parameter][-1]/100]
                        else:
                            self.bounds[parameter] = [self[parameter]-self[parameter]*kwargs[parameter][0]/100, self[parameter]+self[parameter]*kwargs[parameter][-1]/100]
                else:
                    assert kwargs[parameter] > 0, f'{parameter} cannot be negative or zero.'
                    if type_ == 'additive':
                        self.bounds[parameter] = [self[parameter]-kwargs[parameter], self[parameter]+kwargs[parameter]]
                    elif type_ == 'multiplicative':
                        raise ValueError(f'Value must be a tuple for type_ = multiplicative, not a number')
                    else:
                        if self[parameter] < 0:
                            self.bounds[parameter] = [self[parameter]+self[parameter]*kwargs[parameter]/100, self[parameter]-self[parameter]*kwargs[parameter]/100]
                        else:
                            self.bounds[parameter] = [self[parameter]-self[parameter]*kwargs[parameter]/100, self[parameter]+self[parameter]*kwargs[parameter]/100]
                assert self[parameter] >= self.bounds[parameter][0] and self[parameter] <= self.bounds[parameter][1], f'{parameter} value ('+ str(self[parameter]) +') is out of bounds.\nbounds = '+ str(self.bounds[parameter])

        for parameter in ['fwhm', 'fwhm1', 'fwhm2']:
            if parameter in kwargs:
                if kwargs[parameter] is None:
                    self.bounds[parameter] = [0, np.inf]
                elif isinstance(kwargs[parameter], Iterable):
                    if type_ == 'additive':
                        self.bounds[parameter] = [self[parameter]-kwargs[parameter][0], self[parameter]+kwargs[parameter][-1]]
                    elif type_ == 'multiplicative':
                        assert self[parameter] != 0,      f'bounds cannot be set via multiplicative factor for parameter "{parameter}", because its value is zero'
                        assert kwargs[parameter][0] <= 1, f'The first bound multiplicative factor must be less than 1./nParameter: {parameter}/nParameter value: {self[parameter]}/nBound multiplicative factor: {kwargs[parameter]}'
                        assert kwargs[parameter][1] >= 1, f'The second bound multiplicative factor must be higher than 1./nParameter: {parameter}/nParameter value: {self[parameter]}/nBound multiplicative factor: {kwargs[parameter]}'
                        self.bounds[parameter] = [self[parameter]*kwargs[parameter][0], self[parameter]*kwargs[parameter][-1]]
                    else:
                        self.bounds[parameter] = [self[parameter]-self[parameter]*kwargs[parameter][0]/100, self[parameter]+self[parameter]*kwargs[parameter][-1]/100]
                else:
                    assert kwargs[parameter] > 0, f'{parameter} cannot be negative or zero.'
                    if type_ == 'additive':
                        self.bounds[parameter] = [self[parameter]-kwargs[parameter], self[parameter]+kwargs[parameter]]
                    elif type_ == 'multiplicative':
                        raise ValueError(f'Value must be a tuple for type_ = multiplicative, not a number')
                    else:
                        self.bounds[parameter] = [self[parameter]-self[parameter]*kwargs[parameter]/100, self[parameter]+self[parameter]*kwargs[parameter]/100]
                assert self[parameter] >= self.bounds[parameter][0] and self[parameter] <= self.bounds[parameter][1], f'{parameter} value ('+ str(self[parameter]) +') is out of bounds.\nbounds = '+ str(self.bounds[parameter])

                if self.bounds[parameter][0] < 0: self.bounds[parameter][0] = 0

        for parameter in ['m', 'm1', 'm2']:
            if parameter in kwargs:
                if kwargs[parameter] is None:
                    self.bounds[parameter] = [0, 1]
                elif isinstance(kwargs[parameter], Iterable):
                    if type_ == 'additive':
                        self.bounds[parameter] = [self[parameter]-kwargs[parameter][0], self[parameter]+kwargs[parameter][-1]]
                    elif type_ == 'multiplicative':
                        assert self[parameter] != 0,      f'bounds cannot be set via multiplicative factor for parameter "{parameter}", because its value is zero'
                        assert kwargs[parameter][0] <= 1, f'The first bound multiplicative factor must be less than 1./nParameter: {parameter}/nParameter value: {self[parameter]}/nBound multiplicative factor: {kwargs[parameter]}'
                        assert kwargs[parameter][1] >= 1, f'The second bound multiplicative factor must be higher than 1./nParameter: {parameter}/nParameter value: {self[parameter]}/nBound multiplicative factor: {kwargs[parameter]}'
                        self.bounds[parameter] = [self[parameter]*kwargs[parameter][0], self[parameter]*kwargs[parameter][-1]]
                    else:
                        self.bounds[parameter] = [self[parameter]-self[parameter]*kwargs[parameter][0]/100, self[parameter]+self[parameter]*kwargs[parameter][-1]/100]
                else:
                    assert kwargs[parameter] > 0, f'{parameter} cannot be negative or zero.'
                    if type_ == 'additive':
                        self.bounds[parameter] = [self[parameter]-kwargs[parameter], self[parameter]+kwargs[parameter]]
                    elif type_ == 'multiplicative':
                        raise ValueError(f'Value must be a tuple for type_ = multiplicative, not a number')
                    else:
                        self.bounds[parameter] = [self[parameter]-self[parameter]*kwargs[parameter]/100, self[parameter]+self[parameter]*kwargs[parameter]/100]
                if self.bounds[parameter][0] < 0: self.bounds[parameter][0] = 0
                if self.bounds[parameter][1] > 1: self.bounds[parameter][1] = 1
                assert self[parameter] >= self.bounds[parameter][0] and self[parameter] <= self.bounds[parameter][1], f'{parameter} value ('+ str(self[parameter]) +') is out of bounds.\nbounds = '+ str(self.bounds[parameter])

        # amp_bounds, c_bounds, fwhm_bounds (tuple, optional): minimum and
        #     maximum multiplication factor for boundary values. For amp, the
        #     bounds are set between amp*amp_bounds[0] and amp*amp_bounds[1].
        #     For c, bounds are set c-fwhm*c_bound[0] and c+fwhm*c_bound[1].
        #     Finaly, for fwhm, the bounds are set between fwhm1*fwhm_bounds[0]
        #     and fwhm1*fwhm_bounds[0]. Note that if fwhm1*fwhm_bounds[0] is
        #     less than zero, the minimum fwhm boundary is set to zero as it
        #     cannot be negative.

    def build_guess(self, asymmetry=None, fixed=None):
        """Returns initial guess parameter for fitting.

        Args:
            asymmetry (bool, optional): Overwrites the Peak.asymmetry parameter.
                If True, the returned peak function will require two
                fwhm values, one for each half of the peak.
            fixed (bool, optional): Overwrites the Peak.fixed_m parameter.
                If False, m will be a required as an input
                argument on the returned peak function. Default is fixed m with
                value 0.

        Returns:
            three lists: p0, bounds_min, bounds_max
        """
        p0           = []
        bounds_min   = []
        bounds_max   = []

        # decode (for fitting purposes)
        decode_fixed = {}  # dictionary with parameter and value that will not be varied in a fit
        decode_free  = []  # parameters that shall be varied in a fit

        if asymmetry is None:
            asymmetry = self.asymmetry

        if fixed is None:
            fixed = self.fixed
        self._check_fixed(fixed)

        for parameter in ['amp', 'c']:
            if parameter not in fixed:
                assert self.bounds[parameter][1] != self.bounds[parameter][0], 'Minimum boundary cannot be equal to the maximum.\namp_bounds = ' + str(self.bounds[parameter])
                assert self[parameter] >= self.bounds[parameter][0] and self[parameter] <= self.bounds[parameter][1], f'Parameter ({parameter}) value ('+ str(self[parameter]) +') is out of bounds.\nbounds = '+ str(self.bounds[parameter])

                p0.append(self[parameter])
                bounds_min.append(self.bounds[parameter][0])
                bounds_max.append(self.bounds[parameter][1])
                decode_free.append(parameter)
            else:
                decode_fixed[parameter] = self[parameter]

        if asymmetry:
            if 'fwhm' not in fixed:
                for parameter in ['fwhm1', 'fwhm2']:
                    assert self.bounds[parameter][1] != self.bounds[parameter][0], 'Minimum boundary cannot be equal to the maximum.\namp_bounds = ' + str(self.bounds[parameter])
                    assert self[parameter] >= self.bounds[parameter][0] and self[parameter] <= self.bounds[parameter][1], f'Parameter ({parameter}) value ('+ str(self[parameter]) +') is out of bounds.\nbounds = '+ str(self.bounds[parameter])

                if 'm' not in fixed:
                    for parameter in ['m1', 'm2']:
                        if self.bounds[parameter][0] < 0: self.bounds[parameter][0] = 0
                        if self.bounds[parameter][1] > 1: self.bounds[parameter][1] = 1
                        assert self[parameter] >= self.bounds[parameter][0] and self[parameter] <= self.bounds[parameter][1], f'Parameter ({parameter}) value ('+ str(self[parameter]) +') is out of bounds.\nbounds = '+ str(self.bounds[parameter])

                    for parameter in ['fwhm1', 'm1', 'fwhm2', 'm2']:
                        p0.append(self[parameter])
                        bounds_min.append(self.bounds[parameter][0])
                        bounds_max.append(self.bounds[parameter][1])
                        decode_free.append(parameter)
                else:
                    for parameter in ['fwhm1', 'fwhm2']:
                        p0.append(self[parameter])
                        bounds_min.append(self.bounds[parameter][0])
                        bounds_max.append(self.bounds[parameter][1])
                        decode_free.append(parameter)
                    for parameter in ['m1', 'm2']:
                        decode_fixed[parameter] = self[parameter]
            else:
                if 'm' not in fixed:
                    for parameter in ['m1', 'm2']:
                        p0.append(self[parameter])
                        bounds_min.append(self.bounds[parameter][0])
                        bounds_max.append(self.bounds[parameter][1])
                        decode_free.append(parameter)
                    for parameter in ['fwhm1', 'fwhm2']:
                        decode_fixed[parameter] = self[parameter]
                else:
                    for parameter in ['fwhm1', 'm1', 'fwhm2', 'm2']:
                        decode_fixed[parameter] = self[parameter]
        else:
            if 'fwhm' not in fixed:
                if self.bounds['fwhm'][0] < 0: self.bounds['fwhm'][0] = 0
                assert self['fwhm'] >= self.bounds['fwhm'][0] and self['fwhm'] <= self.bounds['fwhm'][1], f'fwhm value ('+ str(self['fwhm']) +') is out of bounds.\nbounds = '+ str(self.bounds['fwhm'])

                if 'm' not in fixed:
                    if self.bounds['m'][0] < 0: self.bounds['m'][0] = 0
                    if self.bounds['m'][1] > 1: self.bounds['m'][1] = 1
                    assert self['m'] >= self.bounds['m'][0] and self['m'] <= self.bounds['m'][1], f'm value ('+ str(self['m']) +') is out of bounds.\nbounds = '+ str(self.bounds['m'])

                    for parameter in ['fwhm', 'm']:
                        p0.append(self[parameter])
                        bounds_min.append(self.bounds[parameter][0])
                        bounds_max.append(self.bounds[parameter][1])
                        decode_free.append(parameter)
                else:
                    p0.append(self['fwhm'])
                    bounds_min.append(self.bounds['fwhm'][0])
                    bounds_max.append(self.bounds['fwhm'][1])
                    decode_free.append('fwhm')
                    decode_fixed['m'] = self['m']
            else:
                if 'm' not in fixed:
                    if self.bounds['m'][0] < 0: self.bounds['m'][0] = 0
                    if self.bounds['m'][1] > 1: self.bounds['m'][1] = 1
                    assert self['m'] >= self.bounds['m'][0] and self['m'] <= self.bounds['m'][1], f'm value ('+ str(self['m']) +') is out of bounds.\nbounds = '+ str(self.bounds['m'])

                    p0.append(self['m'])
                    bounds_min.append(self.bounds['m'][0])
                    bounds_max.append(self.bounds['m'][1])
                    decode_free.append('m')
                else:
                    for parameter in ['fwhm', 'm']:
                        decode_fixed[parameter] = self[parameter]

        # decode function
        def decode(popt, psigma=None):
            peak  = {}
            peak.update(decode_fixed)

            # set initial error to zero
            if psigma is not None:
                error = {}
                for parameter in peak:
                    error[parameter] = 0

            # put parameters from p0 to a dictionary (peak)
            for i, parameter in enumerate(decode_free):
                peak[parameter]  = popt[decode_free.index(parameter)]
                if psigma is not None:
                    error[parameter] = psigma[decode_free.index(parameter)]

            # create peak object
            peak = Peak(**peak)
            if psigma is not None:
                peak.error.update(error)
            peak.asymmetry = copy.copy(self.asymmetry)
            peak.fixed     = copy.copy(self.fixed)
            peak.bounds    = copy.copy(self.bounds)
            return peak

        return p0, bounds_min, bounds_max, decode


    def build_model_str(self, fixed=None, idx=0):
        """Returns instructions for building peak functions.

        Args:
            asymmetry (bool): if True, the returned peak function will require two
                fwhm values, one for each half of the peak.
            fixed_m (False or number): m is the amount of lorentzian
                contribution for a peak. If False, m will be a required as an input
                argument on the returned peak function.
            idx (int, optional): number to be inprinted in the string

        Returns:
            function f(x), function f(x) as string, argument list as string
        """
        if fixed is None:
            fixed = self.fixed

        temp = {key:self[key] for key in fixed}
        temp['idx'] = idx
        temp['asymmetry'] = self.asymmetry
        f_str, args_str = build_model_str(**temp)

        return f_str, args_str

    def build_model(self, fixed=None, idx=0):
        if fixed is None:
            fixed = self.fixed

        temp = {key:self[key] for key in fixed}
        temp['idx'] = idx
        temp['asymmetry'] = self.asymmetry
        return build_model(**temp)


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

    Args:
        data (list): can be passed as positional arguments
        shift( number): initial shift value (does not have initial effect over the data).
        calib (number): initial calibration value (does not have initial effect over the data).
        offset( number): initial offset value (does not have initial effect over the data).
        factor (number): initial multiplicative factor (does not have initial effect over the data).
        *args: Peak objects
    """

    def __init__(self, *args, **kwargs):
        # core
        self._store = []

        # modifiers
        self._factor       = 1
        self._offset       = 0
        self._calib        = 1
        self._shift        = 0

        data, filepath = self._sort_args(args, kwargs)

        # data
        if data is not None:
            if isinstance(data, list):
                for peak in data:
                    self.append(peak)
            else:
                raise TypeError('data must be a list')
        elif filepath is not None:
            self.load(filepath)

    def __str__(self):
        # return str({i:val for i, val in enumerate(self._store)})[1:-1].replace('}, ', '}\n')
        return str({i:val for i, val in enumerate(self._store)})[1:-1].replace(', ', '\n\n')

    def __repr__(self):
        # print('f')
        # return str({i:val for i, val in enumerate(self._store)})[1:-1].replace('}, ', '}\n')
        return str({i:val for i, val in enumerate(self._store)})[1:-1].replace(', ', '\n\n')
        # return str(self._store)


    def __getitem__(self, key):
        if isinstance(key, int):
            return self._store[key]
        elif isinstance(key, slice):
            return Peaks(self._store[key])
        else:
            raise TypeError('Index must be int, not {}'.format(type(key).__name__))

    def __setitem__(self, key, value):
        if isinstance(value, Peak):
            self._store[key] = value
        elif isinstance(value, dict):
            self._store[key] = Peak(**value)
        else:
            raise ValueError('valuea must be a dict or a peak object')
        self._fix_order()

    def __delitem__(self, key):
        del self._store[key]
        self._fix_order()

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return len(self._store)

    def _sort_args(self, args, kwargs):
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
        if any([item not in ['data', 'filepath'] for item in kwargs.keys()]):
            raise AttributeError(f'invalid attributes.\nValid atributes are `data`, `dirpath`, and `filepaths`\nInput attributes: {kwargs.keys()}')

        data     = None
        filepath = None
        if 'data' in kwargs:
            data = kwargs['data']
        elif 'filepath' in kwargs:
            filepath = kwargs['filepath']
        elif len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepath = Path(args[0])
            elif isinstance(args[0], Peak):
                    data = [args[0], ]
            elif isinstance(args[0], Iterable):
                data = args[0]
        elif len(args) > 1:
            data = args
        return data, filepath

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
                self._store = sort(c, self._store)


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
    def shift(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].shift
        return temp
    @shift.setter
    def shift(self, value):
        self.set_shifts(value)
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].offset
        return temp
    @offset.setter
    def offset(self, value):
        self.set_offsets(value)
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        # return self._factor
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].factor
        return temp
    @factor.setter
    def factor(self, value):
        self.set_factors(value)
    @factor.deleter
    def factor(self):
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

    @property
    def spectra(self):
        return self.calculate_spectra()
    @spectra.setter
    def spectra(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @spectra.deleter
    def spectra(self):
            raise AttributeError('Cannot delete object.')


    def save(self, filepath='./Untitled.txt', check_overwrite=False):
        r"""Save peak to a text file. Wrapper for `json.dumps()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user want to overwrite file.

        Returns:
            None

        .. _json.dumps(): https://docs.python.org/3/library/json.html#json.dumps
        """
        filepath = Path(filepath)
        # check overwrite
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query('File already exists!! Do you wish to ovewrite it?', 'yes') == True:
                        pass
                    else:
                        warnings.warn('File not saved because user did not allow overwriting.')
                        return
                else:
                    warnings.warn('filepath is pointing to a folder. Saving file as Untitled.txt')
                    filepath = filepath/'Untitled.txt'

        string2save = ''
        for peak in self:
            string2save += peak._string2save()
            string2save += '\nend_of_peak\n'

        with open(str(filepath), 'w') as file:
            file.write(string2save)

    def load(self, filepath):
        """Load peak from a text file. Wrapper for `json.load()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first decompressed.

        Returns:
            None

        .. _json.load(): https://docs.python.org/3/library/json.html#json.load
        """
        filepath = Path(filepath)

        # remove all peaks
        for i in range(len(self)):
            self.remove(0)

        with open(str(filepath), 'r') as file:
            try:
                peaks_str = file.read().split('end_of_peak')
                for peak_str in peaks_str:
                    temp = br.Peak(amp=0, c=0, fwhm=1)
                    obj = json.loads(peak_str)
                    temp._obj_decode(obj)
                    self.append(temp)
            except json.JSONDecodeError as e:
                if peak_str == '\n':
                    pass
                else:
                    raise e


    def append(self, value):
        """Append peak. Peak order is reassigned.

        Args:
            value (Peak or dict): peak to be appended.

        Returns:
            None
        """
        if isinstance(value, Peak):
            self._store.append(value)
        elif isinstance(value, dict):
            peak = Peak(**value)
            # peak._shift  = self.shift
            # peak._calib  = self.calib
            # peak._offset = self.offset
            # peak._factor = self.factor
            self._store.append(peak)
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
                del self._store[k]
        else:
            del self._store[key]

    def copy(self, value):
        if isinstance(value, Peaks):
            calib  = self.calib
            shift  = self.shift
            offset = self.offset
            factor = self.factor

            # modifiers must be the same 
            # in case you have a different number of peaks between the current obj
            # and the one you are coping from.
            if all_equal(calib) == False:
                raise RuntimeError('calib have different values for different peaks. All values must be the same.')
            if all_equal(shift) == False:
                raise RuntimeError('shift have different values for different peaks. All values must be the same.')
            if all_equal(offset) == False:
                raise RuntimeError('offset have different values for different peaks. All values must be the same.')
            if all_equal(factor) == False:
                raise RuntimeError('factor have different values for different peaks. All values must be the same.')

            self._store = copy.deepcopy(value._store)
            for peak in self:
                peak._calib  = calib
                peak._shift  = shift
                peak._offset = offset
                peak._factor = factor
        else:
            raise ValueError('obj to copy must be of type br.Peaks.')

    def set_calib(self, value, type_='relative'):
        """Set calibration value.

        Args:
            value (number): calibration value (x-coordinates will be multiplied
                by this value).

        Returns:
            None
        """
        # check if value is a number
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        for i in range(len(self)):
            self[i].set_calib(value=value[i], type_=type_)
        self._fix_order()

    def set_shifts(self, value, type_='relative'):
        """Set shift value.

        Args:
            value (float or int): shift value (value will be added to x-coordinates).

        Returns:
            None
        """

        # check if value is a number
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        for i in range(len(self)):
            # peak.calib = value
            self[i].set_shifts(value=value[i], type_=type_)
        self._fix_order()

        # for peak in self._store:
        #     peak.shift = value
        # self._shift = value

    def set_offsets(self, value, type_='relative'):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).

        Returns:
            None
        """
        # check if value is a number
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        for i in range(len(self)):
            self[i].set_offsets(value=value[i], type_=type_)

        # for peak in self._store:
        #     peak.offset = value
        # self._offset = value

    def set_factors(self, value, type_='relative'):
        """Set y multiplicative factor.

        if value is a list, type_ is applied. if value is a number, type_ is set to relative.

        Args:
            value (number): multiplicative factor (y-coordinates will be
                multiplied by this value).

        Returns:
            None
        """
        # check if value is a number
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        for i in range(len(self)):
            self[i].set_factors(value=value[i], type_=type_)

        # if isinstance(value, Iterable):
        #     # value must be the right length
        #     assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        #     for i in range(len(self)):
        #         self[i].set_factors(value=value[i], type_=type_)
        # else:
        #     type_='relative'
        #     if self.factor != value:
        #         if self.factor != 0:
        #             for i in range(len(self)):
        #                 self[i].set_factors(value=self.factor, type_=type_)
        #         if value != 0:
        #             for i in range(len(self)):
        #                 self[i].set_factors(value=value, type_=type_)
        #     self._factor = value


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
            for i in range(len(self._store)):
                count = len([k2 for k2 in key if k2 == i])   # counts same key
                if count > 0:
                    for n in range(count):
                        temp = copy.deepcopy(self._store[key])
                        temp['c']+=temp['fwhm']/4
                        self.append(temp)
            for value in temp:
                self.append(value)
        else:
            temp = copy.deepcopy(self._store[key])
            self._store[key]['c'] -= self._store[key]['fwhm']/4
            temp['c'] += temp['fwhm']/4
            self.append(temp)


    def _find_suitable_x(self):
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

        return np.arange(vmin, vmax, fwhm_min/20)

    def calculate_spectra(self, x=None):
        """Return each peak spectrum separetely.

        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be constructed.

        Returns:
            :py:class:`Spectra`.
        """
        if len(self) > 0:
            if x is None:
                x = self._find_suitable_x()

            ss = br.Spectra(n=len(self))
            for i in range(len(self)):
                ss[i] = self[i].calculate_spectrum(x=x)
            return ss
        else:
            raise ValueError('No peaks defined.')

    def calculate_spectrum(self, x=None):
        """Return the spectrum with all peaks summed up.

        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be constructed.

        Returns:
            :py:class:`Spectrum`.
        """
        ss = self.calculate_spectra(x=x)
        s = br.Spectrum(x=ss[0].x, y=np.zeros(len(ss[0].x)))
        for s1 in ss:
            s += s1
        return s


    def set_bounds(self, **kwargs):
        """Set percentage wise bounds.

        Args:


        Return:
            None
        """
        for peak in self:
            peak.set_bounds(**kwargs)

    def build_guess(self):
        p0         = []
        bounds_min = []
        bounds_max = []

        decode_func = []
        decode_len  = []

        for p in self:
            p0_temp, bounds_min_temp, bounds_max_temp, decode_temp = p.build_guess()
            p0         = p0 + p0_temp
            bounds_min = bounds_min + bounds_min_temp
            bounds_max = bounds_max + bounds_max_temp
            decode_len.append(len(p0_temp))
            decode_func.append(decode_temp)


        def decode(popt, psigma=None):
            peaks = Peaks()

            cumulative_len_list = [sum(decode_len[0:x:1]) for x in range(0, len(decode_len)+1)]
            popt = [popt[cumulative_len_list[i]:cumulative_len_list[i+1]] for i in range(0, len(cumulative_len_list)-1)]
            if psigma is not None:
                psigma = [psigma[cumulative_len_list[i]:cumulative_len_list[i+1]] for i in range(0, len(cumulative_len_list)-1)]

            for i in range(len(decode_func)):
                if psigma is not None:
                    peaks.append(decode_func[i](popt[i], psigma[i]))
                else:
                    peaks.append(decode_func[i](popt[i]))
            return peaks

        return p0, bounds_min, bounds_max, decode

    def build_model_str(self):
        f_str = ''
        args_str = ''

        for idx, p in enumerate(self):
            f_str_temp, args_str_temp = p.build_model_str(idx=idx)
            f_str = f_str + f_str_temp + ' + '
            args_str  = args_str  + args_str_temp  + ', '

        return f_str[:-3], args_str[:-2]

    def build_model(self):
        f_str, args_str = self.build_model_str()
        # print(f_str)
        # print(args_str)
        model_str = f'lambda x, {args_str}: {f_str}'
        return eval(model_str)


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
            list with `ErrorbarContainer`_

        .. matplotlib.pyplot.errorbar(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
        .. ErrorbarContainer: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
        """
        if ax is None:
            ax = plt
            if br.settings.FIGURE_FORCE_NEW_WINDOW:
                figure()
               
        elif type(ax) == str:
            raise ValueError(f'ax parameter cannot be type str ("{ax}").')
        # elif type(ax) == module:
        #     ax = plt

        r = [None]*len(self)
        for i in range(len(self)):
            r[i] = self[i].plot(ax=ax, offset=offset, shift=shift, factor=factor, **kwargs)
            plt.text(self[i]['c']*calib+shift, self[i]['amp']*factor+offset, i, fontsize=14)
        return r


class Collection(MutableMapping):

    def __init__(self, *args, **kwargs):
        # core
        self._store = []

        # argument parsing
        data, dirpath, filepaths = self._sort_args(args, kwargs)

        # data
        if data is not None:
            if isinstance(data, Peaks):
                self.append(data)
            elif isinstance(data, Iterable):
                for peaks in data:
                    self.append(peaks)
            else:
                raise ValueError('data must be a list of type brixs.Peaks.')
        elif dirpath is not None:
            self.load(dirpath)
        elif filepaths is not None:
            self.load(filepaths)

    def __str__(self):
        # temp = str([val for i, val in enumerate(self._store)])[1:-1].replace('], [', ']\n=======\n[')
        # return temp.replace('}, {', '}\n{')
        return str([val for i, val in enumerate(self._store)])[1:-1].replace(', ', '\n=======\n')

    def __repr__(self):
        # temp = str([val for i, val in enumerate(self._store)])[1:-1].replace('], [', ']\n=======\n[')
        # return temp.replace('}, {', '}\n{')
        # return str(self._store)
        # return str([val for i, val in enumerate(self._store)])[1:-1].replace('}, ', '}\n=====\n')
        return str([val for i, val in enumerate(self._store)])[1:-1].replace(', ', '\n=======\n')

    def __getitem__(self, key):
        if type(key) == int:
            return self._store[key]
        elif type(key) == str:
            peaks_attr = ['amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2', 'area']
            if key not in peaks_attr:
                raise KeyError("Collection indices must be integers or a peak attribute ('amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2', 'area')")
            n_peaks = [len(peaks) for peaks in self]
            if all_equal(n_peaks) == False:
                print('WARNING: number of peaks is different between spectra')
            if max(n_peaks) == 0:
                return {}
            else:
                return [[peaks[j][key] for peaks in self] for j in range(max(n_peaks))]
        else:
            raise TypeError("Collection indices must be integers or a peak attribute ('amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2', 'area')")

    def __setitem__(self, key, value):

        if isinstance(value, Peaks):
            self._store[key] = value
        else:
            raise ValueError('value must be a brixs.Peaks object')

    def __delitem__(self, key):
        del self._store[key]

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return len(self._store)

    def _sort_args(self, args, kwargs):
        """checks initial arguments.

         Keyword arguments (kwargs) cannot be mixed with positional arguments.

        The hierarchy for Keyword arguments is: 1) data, 2) y (and x), and finally
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
        if any([item not in ['data', 'dirpath'] for item in kwargs.keys()]):
            raise AttributeError(f'invalid attributes.\nValid atributes are `data`, `dirpath`, and `filepaths`\nInput attributes: {kwargs.keys()}')

        data      = None
        dirpath   = None
        filepaths = None
        if 'data' in kwargs:
            data = kwargs['data']
        elif 'dirpath' in kwargs:
            dirpath = kwargs['dirpath']
        elif 'filepaths' in kwargs:
            filepaths = kwargs['filepaths']
        elif len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                temp = Path(args[0])
                if temp.is_file():
                    filepaths = [temp, ]
                elif temp.is_dir():
                    dirpath = temp
                else:
                    raise ValueError(f'cannot read dirpath or filepath.\nError: {dirpath}')
            elif isinstance(args[0], Iterable):
                if isinstance(args[0][0], str) or isinstance(args[0][0], Path):
                    filepaths = [Path(x) for x in args[0]]
                else:
                    data = args[0]
        elif len(args) > 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                filepaths = [Path(x) for x in args]
            elif isinstance(args[0], Iterable):
                data = args
        return data, dirpath, filepaths



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
    def shift(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].shift
        return temp
    @shift.setter
    def shift(self, value):
        self.set_shifts(value)
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].offset
        return temp
    @offset.setter
    def offset(self, value):
        self.set_offsets(value)
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].factor
        return temp
    @factor.setter
    def factor(self, value):
        self.set_factors(value)
    @factor.deleter
    def factor(self):
            raise AttributeError('Cannot delete object.')


    def save(self, dirpath='./', prefix='peaks_', suffix='.dat', zfill=None, verbose=False, **kwargs):
        r"""Save peak to a text file. Wrapper for `json.dumps()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user want to overwrite file.

        Returns:
            None

        .. _json.dumps(): https://docs.python.org/3/library/json.html#json.dumps
        """
        dirpath = Path(dirpath)

        # check if dirpath is a directory
        assert dirpath.exists(), f'dirpath does not exists.\ndirpath: {dirpath}'
        assert dirpath.is_dir(), f'dirpath is not a directory.\ndirpath: {dirpath}'

        # set filenames
        if zfill is None:
            zfill = n_digits(len(self)-1)[0]

        # saving
        if verbose: print('saving files...')
        for i, peaks in enumerate(self):
            filename = f'{prefix}' + f'{i}'.zfill(zfill) + f'{suffix}'
            if verbose:  print(f':{i}/{len(self)-1}: {filename}')
            peaks.save(filepath=dirpath/filename, **kwargs)
        if verbose: print('Done!')

    def _load(self, dirpath, string='*', verbose=True):
        """THIS FUNCTION WORKS, BUT IT DOES NOT WORK WELL ON A SPECTRA OBJECT.
        
        Load peak from a text file. Wrapper for `json.load()`_.

        Args:
            dirpath (string or path object, optional): filepath, list of 
                filepaths (or folderpaths), or folderpath.
                If filepath, peak is loaded from file and appended. 
                If list, each filepath within the list is loaded.
                If folderpath,
                All files inside folder are loaded and peaks are appended.
                If the filename extension is .gz or .bz2, the file is first decompressed.
            string (str, optional): file names without this string will be ignored.
                Use '*' for matching anything. Default is '*'.
            verbose (bool, optional): verbose. Default is True.


        Returns:
            None

        .. _json.load(): https://docs.python.org/3/library/json.html#json.load
        """
        # reset data
        self._store     = []

        # if dirpath is str
        if isinstance(dirpath, str) or isinstance(dirpath, Path):
            dirpath = Path(dirpath)
            if dirpath.is_file():
                if verbose: print('dirpath is a file')
                if verbose: print('Loading...')
                self.append(Peaks(filepath=dirpath))
                if verbose: print('Done!')
                return
            elif dirpath.is_dir():
                dirpath = [dirpath, ]
            else:
                raise ValueError(f'cannot read dirpath.\ndirpath: {dirpath}')

        # if dirpath is iterable
        if isinstance(dirpath, Iterable):
            if verbose: print('Loading...')
            for j, filepath in enumerate(dirpath):
                if verbose: print(f'{j+1}/{len(dirpath)}: {filepath}')
                
                if Path(filepath).is_dir():
                    fl = filelist(dirpath=filepath, string=string)
                    for i, f in enumerate(fl):
                        if verbose: print(f'    {i+1}/{len(fl)}: {f}')
                        self.append(Peaks(filepath=f))

                elif Path(filepath).is_file():
                    self.append(Peaks(filepath=filepath))
                else:
                    raise ValueError(f'cannot read filepath.\nfilepath: {dirpath}')
        if verbose: print('Done!')
        # print(self._store)


    def append(self, value):
        """Append peak. Peak order is reassigned.

        Args:
            value (Peak or dict): peak to be appended.

        Returns:
            None
        """
        if isinstance(value, Peaks):
            # print('here')
            self._store.append(value)
        else:
            raise ValueError('value must be a brixs.Peaks object')

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
                del self._store[k]
        else:
            del self._store[key]


    def set_calib(self, value, type_='relative'):
        """Set calibration value.

        Args:
            value (number): calibration value (x-coordinates will be multiplied
                by this value).

        Returns:
            None
        """
        # check if value is a number
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        for peaks in self._store:
            peaks.set_calib(value=value, type_=type_)

    def set_shifts(self, value, type_='relative'):
        """Set shift value.

        Args:
            value (float or int): shift value (value will be added to x-coordinates).

        Returns:
            None
        """
        # check if value is a number
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        for peaks in self._store:
            peaks.set_shifts(value=value, type_=type_)

    def set_offsets(self, value, type_='relative'):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).

        Returns:
            None
        """
        # check if value is a number
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        for peaks in self._store:
            peaks.set_offsets(value=value, type_=type_)

    def set_factors(self, value, type_='relative'):
        """Set y multiplicative factor.

        Args:
            value (number): multiplicative factor (y-coordinates will be
                multiplied by this value).

        Returns:
            None
        """
        # check if value is a number
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        for peaks in self._store:
            peaks.set_factors(value=value, type_=type_)


    def get_errors(self, key):
        if type(key) == str:
            peaks_attr = ['amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2', 'area']
            if key not in peaks_attr:
                raise KeyError("Collection indices must be integers or a peak attribute ('amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2', 'area')")
            n_peaks = [len(peaks) for peaks in self]
            if all_equal(n_peaks) == False:
                print('WARNING: number of peaks is different between spectra')
            if max(n_peaks) == 0:
                return {}
            else:
                return [[peaks[j].error[key] for peaks in self] for j in range(max(n_peaks))]
        else:
            raise TypeError("Collection indices must be integers or a peak attribute ('amp', 'c', 'fwhm', 'fwhm1', 'fwhm2', 'm', 'm1', 'm2', 'area')")

    def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, **kwargs):
        """Place a marker at the maximum of every peak position. Wrapper for `matplotlib.pyplot.errorbar()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number, optionakl): defines a vertical offset. Default is 0.
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
            if br.settings.FIGURE_FORCE_NEW_WINDOW:
                figure()
        elif type(ax) == str:
            raise ValueError(f'ax parameter cannot be type str ("{ax}").')
        # elif type(ax) == module:
        #     ax = plt

        # percentage wise increment ====================
        if 'vi' in kwargs and 'vertical_increment' in kwargs:
            raise SyntaxError('keyword argument repeated: vertical increment/vi')
        elif 'vi' in kwargs or 'vertical_increment' in kwargs:
            if 'vi' in kwargs:
                vi = kwargs['vi']
                del kwargs['vi']
            if 'vertical_increment' in kwargs:
                vi = kwargs['vertical_increment']
                del kwargs['vertical_increment']
            temp = [0]*len(self)
            for i in range(len(self)):
                temp[i] = max(self.data[i].y) - min(self.data[i].y)
            vi = max(temp)*factor*vi/100
        else:
            vi = 0

        # percentage wise horizontal increment ====================
        if 'hi' in kwargs and 'horizontal_increment' in kwargs:
            raise SyntaxError('keyword argument repeated: horizontal increment and hi')
        elif 'hi' in kwargs or 'horizontal_increment' in kwargs:
            if 'hi' in kwargs:
                hi = kwargs['hi']
                del kwargs['hi']
            if 'horizontal_increment' in kwargs:
                hi = kwargs['horizontal_increment']
                del kwargs['horizontal_increment']
            temp = [0]*len(self)
            for i in range(len(self)):
                temp[i] = max(self.data[i].x) - min(self.data[i].x)
            hi = max(temp)*factor*hi/100
        else:
            hi = 0

        # offset
        if isinstance(offset, Iterable):
            if len(offset) != len(self):
                raise ValueError(f'offset must be a number of a list with length compatible with the number of spectra.\nnumber of offsets: {len(offset)}\nnumber of spectra: {len(self)}')
        else:
            offset = [offset]*len(self)
            for i in range(len(self)):
                offset[i] = offset[i]+(vi*i)

        # shift
        if isinstance(shift, Iterable):
            if len(shift) != len(self):
                raise ValueError(f'shift must be a number of a list with length compatible with the number of spectra.\nnumber of shift: {len(shift)}\nnumber of spectra: {len(self)}')
        else:
            shift = [shift]*len(self)
            for i in range(len(self)):
                shift[i] = shift[i]+(hi*i)

        # calib
        if isinstance(calib, Iterable):
            if len(calib) != len(self):
                raise ValueError(f'calib must be a number of a list with length compatible with the number of spectra.\nnumber of calib: {len(calib)}\nnumber of spectra: {len(self)}')
        else:
            calib = [calib]*len(self)

        # factor
        if isinstance(factor, Iterable):
            if len(factor) != len(self):
                raise ValueError(f'factor must be a number of a list with length compatible with the number of spectra.\nnumber of factor: {len(factor)}\nnumber of spectra: {len(self)}')
        else:
            factor = [factor]*len(self)


        # plot
        r = [0]*len(self)
        for i in range(len(self)):
            r[i] = self[i].plot(ax=ax, offset=offset[i], shift=shift[i], factor=factor[i], calib=calib[i], **kwargs)
        return r
