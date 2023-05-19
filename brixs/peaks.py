#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Peak objects"""

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
import lmfit
from scipy.signal import find_peaks
from collections.abc import Iterable, MutableMapping
# import json

# backpack
from .backpack.filemanip import filelist
from .backpack.arraymanip import all_equal
from .backpack.model_functions import voigt_fwhm, voigt_area_fwhm, dirac_delta
from .backpack.figmanip import n_digits

# BRIXS
import brixs as br

# common definitions ===========================================================
names = ['amp', 'area', 'c', 'w', 'w1', 'w2', 'm', 'm1', 'm2']

cc   = ['cross-correlation', 'cc']
roll = ['roll', 'rotate', 'r', 'rot']
hard = ['hard', 'x', 'h', 'Hard']
soft = ['soft', 'Soft', 'interp', 'y', 's']
relative = ['relative', 'r', 'rel']
absolute = ['a', 'abs', 'absolute']
increasing = ['inc', 'i', 'up', 'increasing', 'increasingly']
decreasing = ['dec', 'd', 'down', 'decreasing', 'decreasingly']


# %% Peaks =====================================================================
class Peaks(lmfit.Parameters):
    def __init__(self, *args, **kwargs):
        super().__init__()

        self.next_index = 0

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

    def append(self, amp=None, c=None, w=None, w1=None, w2=None, m=0, m1=0, m2=0, area=None):
        index = self.next_index

        if area is None and amp is not None:
            self.add(f'amp_{index}', value=amp, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)
        elif amp is None and area is not None:
            self.add(f'area_{index}', value=area, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)
        else:
            raise ValueError('amp or area not defined')
        
        if c is not None:
            self.add(f'c_{index}',   value=c,   vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)
        else:
            raise ValueError('c not defined')

        if w is not None:
            self.add(f'w_{index}', value=w, vary=True, min=0, max=np.inf, expr=None, brute_step=None)
            self.add(f'm_{index}', value=m, vary=False, min=0, max=1, expr=None, brute_step=None)
        elif w1 is not None and w2 is not None:
            self.add(f'w1_{index}', value=w1, vary=True, min=0, max=np.inf, expr=None, brute_step=None)
            self.add(f'w2_{index}', value=w2, vary=True, min=0, max=np.inf, expr=None, brute_step=None)
            self.add(f'm1_{index}', value=m1, vary=False, min=0, max=1, expr=None, brute_step=None)
            self.add(f'm2_{index}', value=m1, vary=False, min=0, max=1, expr=None, brute_step=None)
        else:
            raise ValueError('Either w or (w1 and w2) must be defined')

        self.next_index += 1

    def add(self, name, value=None, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None):
        # # check if name is valid
        # if not name.startswith(tuple(names)):
        #     raise ValueError(f'name= {name} is not valid.\nValid names:{names}')

        # add
        super().add(name, value=value, vary=vary, min=min, max=max, expr=expr, brute_step=brute_step)

    def __setitem__(self, name, value):
        # # check if name is valid
        # if not name.startswith(tuple(names)):
        #     raise ValueError(f'name= {name} is not valid.\nValid names:{names}')

        # # area and amp
        # if self._enable_update_area_amp:
        #     if self.use_area:
        #         if name == 'amp':
        #             raise ValueError('amp cannot be edited if Peak.use_area is True.\nUse Peak.update_area_amp()\nOr change Peak.use_area to False.')
        #     else:
        #         if name == 'area':
        #             raise ValueError('area cannot be edited if Peak.use_area is False.\nUse Peak.update_area_amp().\nOr change Peak.use_area to True.')
        
        # set values
        super().__setitem__(name, value)
    
        # if name == 'w':
        #     self['w1'].value = self['w']/2
        #     self['w2'].value = self['w']/2

        # # check
        # if self._enable_update_area_amp:
        #     self.update_area_amp()

    def __getitem__(self, name):
        if isinstance(name, int):
            return self.get_params_with_index(i=name)
        else:
            if name in ['amp', 'area', 'c', 'w', 'w1', 'w2', 'm', 'm1', 'm2']:
                if self.length() == 1:
                    index = self.get_indexes()[0]
                    name = name + f'_{index}'
                else:
                    raise ValueError(f'This list has multiple peaks, please, identify the peak number')
            if name in self:
                return super().__getitem__(name)
            elif name.startswith('area'):
                index = name.split('_')[1]
                return lmfit.Parameter(name=name, value=self.calculate_area(i=index), vary=False)
            elif name.startswith('amp'):
                index = name.split('_')[1]
                return lmfit.Parameter(name=name, value=self.calculate_amp(i=index), vary=False)
            elif name.startswith('w'):
                index = name.split('_')[1]
                return lmfit.Parameter(name=name, value=self.calculate_w(i=index), vary=False)
            else:
                raise ValueError(f'Cannot find parameter: {name}')

    def __delitem__(self, i):
        names = self._get_with_index(i=i)
        for name in names:
            del self[name]
    
    def length(self):
        return len(self.get_indexes())
    
    # support
    def copy(self, *args, **kwargs):
        i = None
        if 'i' in kwargs:
            i = kwargs['i']
        elif len(args) > 0:
            if isinstance(args[0], int):
                i = args[0]
            elif isinstance(args[0], Peaks):
                self.clear()
                for name in args[0]:
                    self[name] = copy.deepcopy(args[0][name])
                return
        elif len(args) == 0:
            pass
        else:
            raise ValueError('input not recognized')
        
        peaks = Peaks()
        if i is not None:
            names = self._get_with_index(i=i)
        else:
            names = list(self.keys())
        for name in names:
            peaks[name] = copy.deepcopy(self[name])
        return peaks
    
    def _get_with_index(self, i):
        """return list of names with index i"""
        final = []
        for name in self:
            index = int(name.split('_')[1])
            if index == i:
                final.append(name)
        return final
    
    def get_params_with_index(self, i):
        names = self._get_with_index(i=i)
        final = {}
        for name in names:
            final[name.split('_')[0]] = self[name]
        return final

    def split_params(self):
        indexes = self.get_indexes()
        peaks = {}
        for i in indexes:
            peaks[i] = self.get_params_with_index(i)
        return peaks

    def get_indexes(self):
        """return list of all indexes."""
        final = []
        for name in self:
            index = int(name.split('_')[1])
            if index not in final:
                final.append(index) 
        return final

    def is_asymmetric(self, i):
        if f'w1_{i}' in self and f'w2_{i}' in self:
            return True
        elif f'w_{i}' in self:
            return False
        else:
            raise ValueError(f'w not found for i={i}')
    
    def use_area(self, i):
        if f'area_{i}' in self:
            return True
        elif f'amp_{i}' in self:
            return False
        else:
            raise ValueError(f'area/amp not found for i={i}')
    
    def _find_suitable_x(self, i=None):
        if i is None:
            vmin  = []
            vmax  = []
            step  = []

            for i in self.get_indexes():
                temp = self._find_suitable_x(i=i)
                vmin.append(min(temp))
                vmax.append(max(temp))
                step.append(abs(temp[1]-temp[0]))
            return np.arange(min(vmin), max(vmax), min(step))
        else:
            # get width
            if self.is_asymmetric(i):
                w = self[f'w1_{i}'].value + self[f'w2_{i}'].value
                m = max([self[f'm1_{i}'].value, self[f'm2_{i}'].value])
            else:
                w = self[f'w_{i}'].value
                m = self[f'm_{i}'].value

            if w == 0:
                w = self[f'c_{i}'].value*0.1
            if w == 0:
                w = 1

            vmin = self[f'c_{i}'].value - w * (m*4 + 4)
            vmax = self[f'c_{i}'].value + w * (m*4 + 4)

            if vmin == 0 and vmax == 0:
                vmin = -10
                vmax = 10

            return np.arange(vmin, vmax, w/20)

    # save and load (TODO)*****
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

    # basic
    def remove(self, i):
        """Remove peak.

        Args:
            key (int or list): peaks to be removed.

        Returns:
            None
        """
        for name in self._get_with_index(i=i):
            del self[name]

    def clear(self):
        super().clear()
        self.next_index = 0

    # plot and visualization
    def plot(self, i=None, ax=None, offset=0, shift=0, factor=1, calib=1, **kwargs):
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
                br.figure()

        if i is None:
            for i in self.get_indexes():
                self.plot(i=i, ax=ax, offset=offset, shift=shift, factor=factor, calib=calib, **kwargs)
        else:
            # data
            c    = self[f'c_{i}'].value
            amp  = self.calculate_amp(i)
            xerr = self.calculate_w(i)/2

            if 'lw' not in kwargs and 'linewidth' not in kwargs:
                kwargs['lw'] = 0
            if 'elinewidth' not in kwargs :
                kwargs['elinewidth'] = 2
            if 'marker' not in kwargs :
                kwargs['marker'] = 'o'
            if 'markersize' not in kwargs and 'ms' not in kwargs:
                kwargs['markersize'] = 5

            return ax.errorbar((c*calib)+shift, amp*factor+offset, xerr=xerr*calib, **kwargs)

    # extractors
    def split_peaks(self):
        """return a dict with type Peaks."""
        indexes = self.get_indexes()
        peaks = {}
        for i in indexes:
            peaks[i] = Peaks()
            temp = self.get_params_with_index(i)
            for param in temp:
                temp[param] = {'name':param + '_0',  
                               'value':temp[param].value, 
                               'min':temp[param].min, 
                               'max':temp[param].max,}
                peaks[i].add(**temp[param])
        return peaks
    
    def calculate_spectrum(self, i=None, x=None):
        """Return peak curve.

        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x, with at least 20 points within the peak,
                will be constructed.

        Returns:
            :py:class:`Spectrum`.
        """
        if x is None:
            x = self._find_suitable_x(i=i)
        f = self.model(i=i)

        s = br.Spectrum(x=x, y=f(x))

        # # copy modifiers
        # s._shift  = self.shift
        # # s._shift_roll    = self.shift_roll  # cannot copy shift_roll because step might be different
        # s._shift_interp  = self.shift_interp
        # s._factor = self.factor
        # s._calib  = self.calib
        # s._offset = self.offset
        return s

    def calculate_spectra(self, x=None):
        """Return each peak spectrum separately.

        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be constructed.

        Returns:
            :py:class:`Spectra`.
        """
        if x is None:
            x = self._find_suitable_x(i=i)

        indexes = self.get_indexes()
        n_peaks = len(indexes)

        if n_peaks > 0:
            ss = br.Spectra(n_peaks)
            for j, i in enumerate(indexes):
                ss[i] = self.calculate_spectrum(i=i, x=x)
            return ss
        else:
            raise ValueError('No peaks defined.')
    
    # calculation and info
    def _model_str(self, i=None):
        """Returns string for building peak function.

        Returns:
            function f(x) as string
        """
        if i is None:
            final = ''
            for i in self.get_indexes():
                final += self._model_str(i) + ' + '
            
            return final[:-3]

        if self.is_asymmetric(i):
            if self.use_area(i):
                return f"np.heaviside(c_{i}-x, 0)*br.voigt_area_fwhm(x, area_{i}, c_{i}, w1_{i},  m1_{i}) + np.heaviside(x-c_{i}, 0)*br.voigt_area_fwhm(x, area_{i}, c_{i},  w2_{i},  m2_{i}) + dirac_delta(x, amp_{i}, c_{i})" 
            else:
                return f"np.heaviside(c_{i}-x, 0)*br.voigt_fwhm(x, amp_{i}, c_{i}, w1_{i},  m1_{i}) + np.heaviside(x-c_{i}, 0)*br.voigt_fwhm(x, amp_{i}, c_{i},  w2_{i},  m2_{i}) + dirac_delta(x, amp_{i}, c_{i})" 
        else:
            if self.is_asymmetric(i):
                return f"br.voigt_area_fwhm(x, area_{i}, c_{i}, w_{i}, m_{i})"
            else:
                return f"br.voigt_fwhm(x, amp_{i}, c_{i}, w_{i}, m_{i})"

    def _generate_residual_function(self, x, y, yerr=None, i=None):
        # params = self
        if yerr is None:
            def residual(params, x, y):
                # for name in params:
                #     eval(f'{name} = self[{name}]')
                model = params.model(i=i)
                # pvals = params.valuesdict()
                # model = eval(model_str)
                return model(x)-y
            return residual
        
        else:
            def residual(params, x, y, yerr):
                model = params.model(i=i)
                # pdict = params.valuesdict()
                # model = eval(model_str)
                return (model-y)/yerr
            return residual
        
    def model(self, i=None):
        """Returns a function f(x) for the peak."""
        var_str = ''
        if i is None:
            for name in self:
                var_str += f'{name} = self["{name}"].value' + ', '
        else:
            for name in self._get_with_index(i=i):
                var_str += f'{name} = self["{name}"].value' + ', '

        model_str = f'lambda x, {var_str[:-2]}: {self._model_str(i=i)}'
        return eval(model_str)

    def calculate_area(self, i):
        """Updates the area amp values.

        if use_area = True, amp is updated based on the area.
        if use_area = False, area is updated based on the amp.

        Only necessary if amp and area values are changed by hand.        
        """
        if self.use_area(i=i):
            return self[f'area_{i}']
        else:
            if self.is_asymmetric(i=i):
                c1 = self[f'm1_{i}'].value/np.pi + (1-self[f'm1_{i}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                c2 = self[f'm2_{i}'].value/np.pi + (1-self[f'm2_{i}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return self[f'amp_{i}'].value * (self[f'w1_{i}'].value*c1 + self[f'w2_{i}'].value*c2)/2
            else:
                c = self[f'm_{i}'].value/np.pi + (1-self[f'm_{i}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return self[f'amp_{i}'].value * self[f'w_{i}'].value * c

    def calculate_amp(self, i):
        """Updates the area amp values.

        if use_area = True, amp is updated based on the area.
        if use_area = False, area is updated based on the amp.

        Only necessary if amp and area values are changed by hand.        
        """
        if self.use_area(i=i) == False:
            return self[f'amp_{i}']
        else:
            if self.is_asymmetric(i=i):
                c1 = self[f'm1_{i}'].value/np.pi + (1-self[f'm1_{i}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                c2 = self[f'm2_{i}'].value/np.pi + (1-self[f'm2_{i}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return self[f'area_{i}'].value * 2 /(self[f'w1_{i}'].value*c1 + self[f'w2_{i}'].value*c2)
            else:
                c = self[f'm_{i}'].value/np.pi + (1-self[f'm_{i}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return self[f'area_{i}'].value/self[f'w_{i}'].value/c

    def calculate_w(self, i):
        if self.is_asymmetric(i):
            return self[f'w1_{i}'].value + self[f'w2_{i}'].value
        else:
            return self[f'w_{i}'].value

    def fit(self, x, y, yerr=None, i=None, method='least_squares', return_fitted_peaks=False):
        """output, peaks = fit()
        
        I think x and y must be monotonic.
        """
        residual = self._generate_residual_function(x=x, y=y, yerr=yerr, i=i)
        if yerr is None:
            out = lmfit.minimize(residual, method=method, params=self, args=(x, y))
        else:
            out = lmfit.minimize(residual, method=method, params=self, args=(x, y, yerr))

        # return peaks2
        if return_fitted_peaks:
            return out, out.params
        else:
            for name in self:
                self[name] = out.params[name]
            return out
    
    def find(self, x, y, prominence=5, width=4, moving_average_window=8):
        """I think x and y must be monotonic and uniform."""
        # check width and moving_average_window
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
            
        # data smoothing
        if moving_average_window > 1:
            y2 = br.moving_average(y, moving_average_window)
            x2 = br.moving_average(x, moving_average_window)
        else:
            x2 = x
            y2 = y

        # parameters
        if prominence is None:
            prominence = (max(y2)-min(y2))*0.1
        else:
            prominence = (max(y2)-min(y2))*prominence/100

        try:
            peaks, d = find_peaks(y2, prominence=prominence, width=width)
            assert len(peaks) > 0, 'No peaks found.'
            self.clear()
            for i in range(len(peaks)):
                amp = d['prominences'][i]+max([y2[d['right_bases'][i]], y2[d['left_bases'][i]]])
                c = x2[peaks[i]]
                w = abs(d['widths'][i]*np.mean(np.diff(x)))

                self.append(amp=amp, c=c, w=w)
        except IndexError:
            pass

    # %% attributes and properties  (TODO)*****
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
    def shift_roll(self):
        return self._shift_roll
    @shift_roll.setter
    def shift_roll(self, value):
        # raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shifts(value, mode='roll').")
        self.set_shift(value, mode='roll')
    @shift_roll.deleter
    def shift_roll(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift_interp(self):
        return self._shift_interp
    @shift_interp.setter
    def shift_interp(self, value):
        # raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shifts(value, mode='interp').")
        self.set_shift(value, mode='interp')
    @shift_interp.deleter
    def shift_interp(self):
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
    def step(self):
        return self._step
    @step.setter
    def step(self, value):
        raise AttributeError('Cannot edit object.')
    @step.deleter
    def step(self):
        raise AttributeError('Cannot delete object.')

    # modifiers  (TODO)*****
    def set_calib(self, value, type_='relative'):
        """Set calibration value.

        Args:
            value (number): calibration value (x-coordinates will be multiplied
                by this value).

        Returns:
            None
        """
        # is relative?
        if type_ in relative:
            value = self.calib * value

        # apply
        if self.calib != value:
            if self.calib != 1:
                self['c'].value = self['c'].value*self.calib**-1
                self['w'].value = abs(self['w'].value*self.calib**-1)
                if self.asymmetry:
                    self['w1'].value = abs(self['w1'].value*self.calib**-1)
                    self['w2'].value = abs(self['w2'].value*self.calib**-1)
            if value != 1:
                self['c'].value = self['c'].value*value
                self['w'].value = abs(self['w'].value*value)
                if self.asymmetry:
                    self['w1'].value = abs(self['w1'].value*value)
                    self['w2'].value = abs(self['w2'].value*value)
            self._calib = value

            # check
            self.update_area_amp()

    def set_shift(self, value, mode, type_='relative'):
        """Set shift value.

        Args:
            value (float or int): shift value (value will be added to x-coordinates).

        Returns:
            None
        """
        # is relative?
        if type_ in relative:
            value = self.shift + value
        
        # apply
        if mode in roll: 
            value = self.step*value

            if self.shift_roll != value:
                if self.shift_roll != 0:
                    self['c'].value = self['c'].value - self.shift_roll
                if value != 0:
                    self['c'].value = self['c'].value + value
                self._shift_roll = value
        elif mode in hard:
            if self.shift != value:
                if self.shift != 0:
                    self['c'].value = self['c'].value - self.shift
                if value != 0:
                    self['c'].value = self['c'].value + value
                self._shift = value
        elif mode in soft:
            if self.shift_interp != value:
                if self.shift_interp != 0:
                    self['c'].value = self['c'].value - self.shift_interp
                if value != 0:
                    self['c'].value = self['c'].value + value
                self._shift_interp = value
        else:
            raise ValueError(f'Invalid mode. Valid options are `roll`, `x`, `interp`.') 

    def set_offset(self, value, type_='relative'):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).

        Returns:
            None
        """
        # is relative?
        if type_ in relative:
            value = self.offset + value

        # apply
        if self.offset != value:
            if self.offset != 0:
                self['amp'].value = self['amp'].value - self.offset
            if value != 0:
                self['amp'].value = self['amp'].value + value
            self._offset = value

        # check
        self.update_area_amp()

    def set_factor(self, value, type_='relative'):
        """Set y multiplicative factor.

        Args:
            value (number): multiplicative factor (y-coordinates will be
                multiplied by this value).

        Returns:
            None
        """
        # is relative?
        if type_ in relative:
            value = self.factor * value

        if self.factor != value:
            if self.factor != 1:
                self['amp'].value  = self['amp'].value * self.factor**-1
                self['area'].value = self['area'].value * self.factor**-1
            if value != 1:
                self['amp'].value  = self['amp'].value * value
                self['area'].value = self['area'].value * value
            self._factor = value

        # check
        self.update_area_amp()

    # extractors  (TODO)*****
    def _split_params(self):
        # index list
        indexes = []
        for key in self.keys():
            index = int(key.split('_')[-1])
            if index not in indexes:
                indexes.append(index)
        indexes = br.sort(indexes, indexes)

        # allocation
        peaks2 = Peaks()

        # sort
        for index in indexes:
            peak = Peak()
            peak.index = index
            # peak.pretty_print()
            for key in self.keys():
                # print(index)
                if index == int(key.split('_')[-1]):
                    peak[key] = self[key]
            # peak.pretty_print()
            peaks2.append(peak)
        
        return peaks2

    # #############  (TODO)*****

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



# %%
class Collection(lmfit.Parameters):
    def __init__(self, *args, **kwargs):
        super().__init__()

        self.next_index = 0

        if 'data' in kwargs:
            data = kwargs['data']
        elif len(args) == 1:
            data = args[0]
        else:
            data = args

        # append peaks
        for peaks in data:
            self.append(peaks)        
    
    def __getitem__(self, name):
        if type(name) == int:
            return self.split_peaks_per_spectrum()[name]
        
        elif type(name) == str:
            if len(name.split('_')) == 2:
                return self.get_values(attr=name)
            if len(name.split('_')) == 3:
                if name in self:
                    return super().__getitem__(name)
                else:
                    raise ValueError(f'Cannot find parameter: {name}')
    
    def length(self):
        return len(self._get_indexes_i2())
    
    # basic
    # def append(self, peaks):
    #     temp = peaks.copy()
    #     for name in temp:
    #         temp[name].set(expr='')
    #         self[name + f'_{i}'] = temp[name]

    #     for name in peaks:
    #         self[name + f'_{self.next_index}'] = peaks[name]
    #     self.next_index += 1 


    # support   
    def copy(self, *args):
        if len(args) == 0:
            return copy.deepcopy(self)
        elif len(args) == 1:
            if isinstance(args[0], Collection):
                raise NotImplementedError('sorry, not implemented')
        else:
            raise ValueError('input not recognized')
    
    def copy_from_spectra(self, ss):
        for i, s in enumerate(ss):
            temp = s.peaks.copy()
            for name in temp:
                temp[name].set(expr='')
                self[name + f'_{i}'] = temp[name]

    def copy_to_spectra(self, ss):
        temp = self.copy()
        for par in temp:
            temp[par].set(expr='')

        for i2, s in enumerate(ss):
            s.peaks.clear()
            names = temp._get_with_index(i2=i2)
            for name in names:
                _name = '_'.join(name.split('_')[:2])
                s.peaks[_name] = temp[name]

    def _get_indexes_i2(self):
        """return list of all indexes i2."""
        final = []
        for name in self:
            index = int(name.split('_')[2])
            if index not in final:
                final.append(index) 
        return final

    def _get_indexes_i1(self, i2):
        """return list of all indexes i1 for a specific i2."""
        final = []
        for name in self:
            split = name.split('_')
            _i1 = int(split[1])
            _i2 = int(split[2])
            if _i1 not in final and _i2 == i2:
                final.append(_i1) 
        return final

    def _get_with_index(self, i1=None, i2=None):
        """Return list of names with i2."""
        assert i1 is not None or i2 is not None, 'i1 and i2 cannot both be None'
        
        final = []
        if i1 is None:
            for name in self:
                split = name.split('_')
                # _i1 = int(split[1])
                _i2 = int(split[2])
                if _i2 == i2:
                    final.append(name)    
        elif i2 is None:
            for name in self:
                split = name.split('_')
                _i1 = int(split[1])
                # _i2 = int(split[2])
                if _i1 == i1:
                    final.append(name)
        else:
            for name in self:
                split = name.split('_')
                _i1 = int(split[1])
                _i2 = int(split[2])
                if _i1 == i1 and _i2 == i2:
                    final.append(name)
        return final

    def get_params(self, attr, i1=None):
        final = []
        if i1 is None:
            split = attr.split('_')
            i1 = split[1]
            attr = split[0]

        for i2 in self._get_indexes_i2():
            name = attr + f"_{i1}" + f"_{i2}"
            if name in self:
                final.append(self[name])
            else:
                raise ValueError(f'Cannot find parameter: {name}')
        return final

    def split_peaks_per_spectrum(self):
        temp = self.copy()
        for par in temp:
            temp[par].set(expr='')

        final = {}
        for i2 in self._get_indexes_i2():
            final[i2] = Peaks()
            for i1 in self._get_indexes_i1(i2=i2):
                names = self._get_with_index(i1=i1, i2=i2)
                for name in names:
                    _name = '_'.join(name.split('_')[:2])
                    final[i2][_name] = temp[name] 
        return final
    
        

    def _get_params_with_index(self, i1, i2):
        """Return peaks with all names with i2."""
        names = self._get_with_index(i1=i1, i2=i2)
        final = {}
        for name in names:
            final[name.split('_')[0]] = self[name]
        return final

    def is_asymmetric(self, i1, i2):
        if f'w1_{i1}_{i2}' in self and f'w2_{i1}_{i2}' in self:
            return True
        elif f'w_{i1}_{i2}' in self:
            return False
        else:
            raise ValueError(f'w not found for i1={i1}, i2={i2}')
    
    def use_area(self, i1, i2):
        if f'area_{i1}_{i2}' in self:
            return True
        elif f'amp_{i1}_{i2}' in self:
            return False
        else:
            raise ValueError(f'amp/area not found for i1={i1}, i2={i2}')
   
    def get_values(self, attr, i1=None):
        final = []
        if i1 is None:
            i1 = attr.split('_')[1]
            for i2 in self._get_indexes_i2():
                name = attr + f"_{i2}"
                if name in self:
                    final.append(self[name].value)
                elif name.startswith('area'):
                    final.append(self.calculate_area(i1=i1, i2=i2))
                elif name.startswith('amp'):
                    final.append(self.calculate_amp(i1=i1, i2=i2))
                elif name.startswith('w'):
                    final.append(self.calculate_w(i1=i1, i2=i2))
                else:
                    raise ValueError(f'Cannot find parameter: {name}')
        else:
            for i2 in self._get_indexes_i2():
                name = attr + f"_{i1}" + f"_{i2}"
                if name in self:
                    final.append(self[name].value)
                elif name.startswith('area'):
                    final.append(self.calculate_area(i1=i1, i2=i2))
                elif name.startswith('amp'):
                    final.append(self.calculate_amp(i1=i1, i2=i2))
                elif name.startswith('w'):
                    final.append(self.calculate_w(i1=i1, i2=i2))
                else:
                    raise ValueError(f'Cannot find parameter: {name}')
        return final
    
    def _find_suitable_x(self, i2, i1=None):
        if i1 is None:
            vmin  = []
            vmax  = []
            step  = []

            for i1 in self._get_indexes_i1(i2=i2):
                temp = self._find_suitable_x(i2=i2, i1=i1)
                vmin.append(min(temp))
                vmax.append(max(temp))
                step.append(abs(temp[1]-temp[0]))
            return np.arange(min(vmin), max(vmax), min(step))
        else:
            # get width
            if self.is_asymmetric(i1=i1, i2=i2):
                w = self[f'w1_{i1}_{i2}'].value + self[f'w2_{i1}_{i2}'].value
                m = max([self[f'm1_{i1}_{i2}'].value, self[f'm2_{i1}_{i2}'].value])
            else:
                w = self[f'w_{i1}_{i2}'].value
                m = self[f'm_{i1}_{i2}'].value

            if w == 0:
                w = self[f'c_{i1}_{i2}'].value*0.1
            if w == 0:
                w = 1

            vmin = self[f'c_{i1}_{i2}'].value - w * (m*4 + 4)
            vmax = self[f'c_{i1}_{i2}'].value + w * (m*4 + 4)

            if vmin == 0 and vmax == 0:
                vmin = -10
                vmax = 10

            return np.arange(vmin, vmax, w/20)
        
    # calculation and info
    def calculate_spectrum(self, i2, i1=None, x=None):
        """Return peak curve.

        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x, with at least 20 points within the peak,
                will be constructed.

        Returns:
            :py:class:`Spectrum`.
        """
        if x is None:
            x = self._find_suitable_x(i2=i2, i1=i1)
        f = self.model(i2=i2, i1=i1)

        s = br.Spectrum(x=x, y=f(x))
        return s
    
    def calculate_area(self, i1, i2):
        """Updates the area amp values.

        if use_area = True, amp is updated based on the area.
        if use_area = False, area is updated based on the amp.

        Only necessary if amp and area values are changed by hand.        
        """
        if self.use_area(i1=i1, i2=i2):
            return self[f'area_{i1}_{i2}']
        else:
            if self.is_asymmetric(i1=i1, i2=i2):
                c1 = self[f'm1_{i1}_{i2}'].value/np.pi + (1-self[f'm1_{i1}_{i2}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                c2 = self[f'm2_{i1}_{i2}'].value/np.pi + (1-self[f'm2_{i1}_{i2}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return self[f'amp_{i1}_{i2}'].value * (self[f'w1_{i1}_{i2}'].value*c1 + self[f'w2_{i1}_{i2}'].value*c2)/2
            else:
                c = self[f'm_{i1}_{i2}'].value/np.pi + (1-self[f'm_{i1}_{i2}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return self[f'amp_{i1}_{i2}'].value * self[f'w_{i1}_{i2}'].value * c

    def calculate_amp(self, i1, i2):
        """Updates the area amp values.

        if use_area = True, amp is updated based on the area.
        if use_area = False, area is updated based on the amp.

        Only necessary if amp and area values are changed by hand.        
        """
        if self.use_area(i1=i1, i2=i2) == False:
            return self[f'amp_{i1}_{i2}']
        else:
            if self.is_asymmetric(i1=i1, i2=i2):
                c1 = self[f'm1_{i1}_{i2}'].value/np.pi + (1-self[f'm1_{i1}_{i2}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                c2 = self[f'm2_{i1}_{i2}'].value/np.pi + (1-self[f'm2_{i1}_{i2}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return self[f'area_{i1}_{i2}'].value * 2 /(self[f'w1_{i1}_{i2}'].value*c1 + self[f'w2_{i1}_{i2}'].value*c2)
            else:
                c = self[f'm_{i1}_{i2}'].value/np.pi + (1-self[f'm_{i1}_{i2}'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return self[f'area_{i1}_{i2}'].value/self[f'w_{i1}_{i2}'].value/c

    def calculate_w(self, i1, i2):
        if self.is_asymmetric(i1=i1, i2=i2):
            return self[f'w1_{i1}_{i2}'].value + self[f'w2_{i1}_{i2}'].value
        else:
            return self[f'w_{i1}_{i2}'].value

    def _model_str(self, i2, i1=None):
        """Returns string for building peak function.

        i1

        Returns:
            function f(x) as string
        """
        if i1 is None:
            final = ''
            for i1 in self._get_indexes_i1(i2=i2):
                final += self._model_str(i1=i1, i2=i2) + ' + '
            return final[:-3]

        if self.is_asymmetric(i1=i1, i2=i2):
            if self.use_area(i1=i1, i2=i2):
                return f"np.heaviside(c_{i1}_{i2}-x, 0)*br.voigt_area_fwhm(x, area_{i1}_{i2}, c_{i1}_{i2}, w1_{i1}_{i2},  m1_{i1}_{i2}) + np.heaviside(x-c_{i1}_{i2}, 0)*br.voigt_area_fwhm(x, area_{i1}_{i2}, c_{i1}_{i2},  w2_{i1}_{i2},  m2_{i1}_{i2}) + dirac_delta(x, amp_{i1}_{i2}, c_{i1}_{i2})" 
            else:
                return f"np.heaviside(c_{i1}_{i2}-x, 0)*br.voigt_fwhm(x, amp_{i1}_{i2}, c_{i1}_{i2}, w1_{i1}_{i2},  m1_{i1}_{i2}) + np.heaviside(x-c_{i1}_{i2}, 0)*br.voigt_fwhm(x, amp_{i1}_{i2}, c_{i1}_{i2},  w2_{i1}_{i2},  m2_{i1}_{i2}) + dirac_delta(x, amp_{i1}_{i2}, c_{i1}_{i2})" 
        else:
            if self.is_asymmetric(i1=i1, i2=i2):
                return f"br.voigt_area_fwhm(x, area_{i1}_{i2}, c_{i1}_{i2}, w_{i1}_{i2}, m_{i1}_{i2})"
            else:
                return f"br.voigt_fwhm(x, amp_{i1}_{i2}, c_{i1}_{i2}, w_{i1}_{i2}, m_{i1}_{i2})"

    def model(self, i2, i1=None):
        """Returns a function f(x) for the peak."""
        var_str = ''
        for name in self._get_with_index(i1=i1, i2=i2):
            var_str += f'{name} = self["{name}"].value' + ', '

        model_str = f'lambda x, {var_str[:-2]}: {self._model_str(i1=i1, i2=i2)}'
        return eval(model_str)
    
    def _generate_residual_function(self, yerr=None):

        def residual(params, xs, ys):
            residual = []

            # make residual per data set
            for i2, x in enumerate(xs):
                residual.append(ys[i2] - params.model(i2=i2)(x))

            residual = np.concatenate(residual).ravel()
            return residual

        return residual
    
    def fit(self, xs, ys, yerr=None, method='least_squares', return_fitted_peaks=False):
        """output, peaks = fit()
        
        I think x and y must be monotonic.
        """
        residual = self._generate_residual_function(yerr=yerr)
        if yerr is None:
            out = lmfit.minimize(residual, method=method, params=self, kws={'xs':xs, 'ys':ys})
        else:
            raise NotImplementedError('sorry, not implemented yet')
            # out = lmfit.minimize(residual, method=method, params=self, args=(x, y, yerr))

        # return peaks2
        if return_fitted_peaks:
            return out, out.params
        else:
            for name in self:
                self[name] = out.params[name]
            return out
    


# %%


# %%
class Peak(MutableMapping):
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

        # sort
        data, filepath = self._sort_args(args, kwargs)

        # data
        if data is not None:
            if isinstance(data, Iterable):
                for peak in data:
                    self.append(peak)
            else:
                raise TypeError('data must be a list or tuple')
        elif filepath is not None:
            self.load(filepath)
        
    # basic
    def __str__(self):
        # return str({i:val for i, val in enumerate(self._store)})[1:-1].replace('}, ', '}\n')
        # return str({i:val for i, val in enumerate(self._store)})[1:-1].replace(', ', '\n\n')
        return str({i:peak for i, peak in enumerate(self._store)})#[1:-1].replace(', ', '\n\n')
 

    def __repr__(self):
        return str({i:peak for i, peak in enumerate(self._store)})#[1:-1].replace(', ', '\n\n')
        # return str({i:val for i, val in enumerate(self._store)})[1:-1].replace('}, ', '}\n')
        # return str({i:val for i, val in enumerate(self._store)})[1:-1].replace(', ', '\n\n')
        # return str(self._store)

    def __getitem__(self, key):
        if isinstance(key, int):
            return self._store[key]
        elif isinstance(key, slice):
            return Peaks(self._store[key])
        else:
            raise TypeError('Index must be int or a slice, not {}'.format(type(key).__name__))

    def __setitem__(self, key, value):
        if isinstance(key, int):
            if isinstance(value, Peak):
                self._store[key] = value
            else:
                raise TypeError('value must be br.Peak, not {}'.format(type(key).__name__))
        else:
            raise TypeError('Index must be int, not {}'.format(type(key).__name__))

    def __delitem__(self, key):
        del self._store[key]

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return len(self._store)

    # support
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

    def reorder(self):
        """Change the order of the peaks in the list according to its position c."""
        if self == []:
            return
        else:
            l = len(self)
            c = [0]*l
            for i, peak in enumerate(self):
                c[i] = peak['c'].value
            if sum(1 for test in np.diff(c) if test < 0) > 0:
                self._store = br.sort(c, self._store)
        
        # update index
        for i in range(len(self)):
            self[i].index = i

class _Collection(MutableMapping):

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

    # basic
    def __str__(self):
        # temp = str([val for i, val in enumerate(self._store)])[1:-1].replace('], [', ']\n=======\n[')
        # return temp.replace('}, {', '}\n{')
        # return str([val for i, val in enumerate(self._store)])[1:-1].replace(', ', '\n=======\n')
        return str(self._store)
    
    def __repr__(self):
        return str(self._store)
        # temp = str([val for i, val in enumerate(self._store)])[1:-1].replace('], [', ']\n=======\n[')
        # return temp.replace('}, {', '}\n{')
        # return str(self._store)
        # return str([val for i, val in enumerate(self._store)])[1:-1].replace('}, ', '}\n=====\n')
        # return str([val for i, val in enumerate(self._store)])[1:-1].replace(', ', '\n=======\n')

    def __getitem__(self, key):
        if type(key) == int:
            return self._store[key]
        
        elif type(key) == str:
            if key not in names:
                raise KeyError(f"Collection indices must be integers or a peak attribute ({names})")
            n_peaks = [len(peaks) for peaks in self]
            if all_equal(n_peaks) == False:
                raise ValueError(f'Number of peaks is different between spectra.\nNumber of peaks: {n_peaks}')
            if max(n_peaks) == 0:
                return {}
            else:
                # final = []
                # for j in range(max(n_peaks)):
                #     if j <= 
                #         temp = [peaks[j][key].value for peaks in self]
                #     else:
                #         pass
                # return final
                return [[peaks[j][key].value for peaks in self] for j in range(max(n_peaks))]
        else:
            raise KeyError(f"Collection indices must be integers or a peak attribute ({names})")

    def __setitem__(self, key, value):
        print('nothing to set')
        # if isinstance(value, Peaks):
        #     self._store[key] = value
        # else:
        #     raise ValueError('value must be a brixs.Peaks object')

    def __delitem__(self, key):
        print('nothing to del')
        # del self._store[key]

    def __iter__(self):
        return iter(self._store)

    def __len__(self):
        return len(self._store)

    # support
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

    # save and load
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

    # calculation and info
    def calculate_spectra(self, x=None):
        ss = br.Spectra()
        for peaks in self:
            s = peaks.calculate_spectrum(x=x)
            ss.append(s)
        return ss

    def errors(self, key):
        if type(key) == int:
            return self._store[key]
        
        elif type(key) == str:
            if key not in names:
                raise KeyError(f"Collection indices must be integers or a peak attribute ({names})")
            n_peaks = [len(peaks) for peaks in self]
            if all_equal(n_peaks) == False:
                print('WARNING: number of peaks is different between spectra')
            if max(n_peaks) == 0:
                return {}
            else:
                return [[peaks[j][key].stderr for peaks in self] for j in range(max(n_peaks))]
        else:
            raise KeyError(f"Collection indices must be integers or a peak attribute ({names})")
        
    # plotting and visualization
    def pretty_print(self):
        for i in range(len(self)):
            print('='*5 + ' Spectrum ' + str(i) + ' ' + '='*5)
            self[i].pretty_print()
            print('\n')

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