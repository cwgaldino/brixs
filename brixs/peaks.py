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
import numbers

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
    """m vary is false by default
    """

    def __init__(self, *args, **kwargs):
        super().__init__()
        self._reset_modifiers()

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
        # define index
        indexes = self.get_indexes()
        if indexes == []:
            index = 0
        else:
            index = max(indexes) + 1

        # check parameters
        assert area is not None or amp is not None, 'amp or area not defined'
        assert c is not None, 'c not defined'
        assert w is not None or (w1 is not None and w2 is not None), 'Either w or (w1 and w2) must be defined'

        assert area is None or isinstance(area, numbers.Number), 'area must be a number'
        assert amp is None or isinstance(amp, numbers.Number), 'amp must be a number'
        assert c is None or isinstance(c, numbers.Number), 'c must be a number'
        assert w is None or isinstance(w, numbers.Number), 'w must be a number'
        assert w1 is None or isinstance(w1, numbers.Number), 'w1 must be a number'
        assert w2 is None or isinstance(w2, numbers.Number), 'w2 must be a number'
        assert m is None or isinstance(m, numbers.Number), 'm must be a number'
        assert m1 is None or isinstance(m1, numbers.Number), 'm1 must be a number'
        assert m2 is None or isinstance(m2, numbers.Number), 'm2 must be a number'

        if isinstance(w, numbers.Number):
            assert w >= 0, 'w must be positive'
        if isinstance(w1, numbers.Number):
            assert w1 >= 0, 'w1 must be positive'
        if isinstance(w2, numbers.Number):
            assert w2 >= 0, 'w2 must be positive'

        if isinstance(m, numbers.Number):
            assert m <= 1 and m >= 0, 'm must be positive less than 1'
        if isinstance(m1, numbers.Number):
            assert m1 <= 1 and m1 >= 0, 'm1 must be positive less than 1'
        if isinstance(m2, numbers.Number):
            assert m2 <= 1 and m2 >= 0, 'm2 must be positive less than 1'

        # add peak
        if area is None:
            self.add(f'amp_{index}', value=amp, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)
        else:
            self.add(f'area_{index}', value=area, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)

        self.add(f'c_{index}',   value=c,   vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)

        if w1 is not None and w2 is not None:
            self.add(f'w1_{index}', value=w1, vary=True, min=0, max=np.inf, expr=None, brute_step=None)
            self.add(f'w2_{index}', value=w2, vary=True, min=0, max=np.inf, expr=None, brute_step=None)
            self.add(f'm1_{index}', value=m1, vary=False, min=0, max=1, expr=None, brute_step=None)
            self.add(f'm2_{index}', value=m2, vary=False, min=0, max=1, expr=None, brute_step=None)
        else:
            self.add(f'w_{index}', value=w, vary=True, min=0, max=np.inf, expr=None, brute_step=None)
            self.add(f'm_{index}', value=m, vary=False, min=0, max=1, expr=None, brute_step=None)

    def add(self, name, value=None, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None):
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
            names = self._get_names_with_index(i=i)
        else:
            names = list(self.keys())
        for name in names:
            peaks[name] = copy.deepcopy(self[name])
        return peaks
    
    def _get_names_with_index(self, i):
        """return list of names with index i
        if i is negative, gets indexes from high to low
        """
        if i < 0:
            try: 
                i = np.sort(self.get_indexes())[i]
            except IndexError:
                raise AssertionError(f'index {i} does not exist\nPeak indexes available: {self.get_indexes()}')
            
        assert i in self.get_indexes(), f'index {i} does not exist\nPeak indexes available: {self.get_indexes()}'

        final = []
        for name in self:
            index = int(name.split('_')[1])
            if index == i:
                final.append(name)
        return final
    
    def get_params_with_index(self, i):
        names = self._get_names_with_index(i=i)
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
            key (int): peaks to be removed.

        Returns:
            None
        """
        names = self._get_names_with_index(i=i)
        for name in names:
            del self[name]

    def clear(self):
        """Erase all peaks."""
        super().clear()

    def replace_index(self, i, i_new):
        """Change the index of a peak."""
        # check i exists

        # check i_new is int

        # get params
        params = self.get_params_with_index(i)

        # create new
        for p in params:
            name_new = params[p].name.split('_')[0] + '_' + str(i_new)
            self.add(name=name_new, 
                     value=params[p].value, 
                     vary=params[p].vary, 
                     min=params[p].min, 
                     max=params[p].max, 
                     expr=params[p].expr, 
                     brute_step=params[p].brute_step)
        
        # delete old
        self.remove(i)

    def fix_indexes(self):
        """Fix indexes so they start from zero."""
        indexes = self.get_indexes()
        final = len(indexes)
        
        # replace indexes
        for i_new, index in enumerate(np.sort(indexes)):
            if i_new != index:
                self.replace_index(index, i_new)

        # remove
        for index in self.get_indexes():
            if index >= final:
                self.remove(index)

    def reorder(self, attr='c'):
        """Change peak indexes based on the position of the peak or other attr."""
        raise NotImplementedError('sorry, not implemented yet.')

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
            for name in self._get_names_with_index(i=i):
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
        self.check_feasibility()
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

    def check_feasibility(self):
        """Check if start values are within boundaries."""
        for p in self:
            value = self[p].value
            vmax  = self[p].max
            vmin  = self[p].min

            if value > vmax or value < vmin:
                raise ValueError(f'Start value for parameter {p} is out of bounds\nbound max = {vmax}\nstart value = {value}\nbound min = {vmin}')


    # %% attributes and properties
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
    def roll(self):
        return self._roll
    @roll.setter
    def roll(self, value):
        # raise AttributeError("Attribute is 'read only'. Cannot set attribute.\nPlease, use Spectrum.set_shifts(value, mode='roll').")
        self.set_roll(value)
    @roll.deleter
    def roll(self):
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

    # modifiers
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
                for name in self:
                    if name.split('_')[0] == 'c':
                        self[name].value = self[name].value*self.calib**-1
                    elif name.split('_')[0] == 'w':
                        self[name].value = abs(self[name].value*self.calib**-1)
            if value != 1:
                for name in self:
                    if name.split('_')[0] == 'c':
                        self[name].value = self[name].value*value
                    elif name.split('_')[0] == 'w':
                        self[name].value = abs(self[name].value*value)
            self._calib = value

    def set_roll(self, value, type_='relative'):
        """Set roll value (compatible with Spectrum roll).

        Args:
            value (float or int): shift value (value will be added to x-coordinates).

        Returns:
            None
        """
        # is relative?
        if type_ in relative:
            value = self.roll + value
        
        # Fix step
        value = self.step*value

        # apply
        if self.roll != value:
            if self.roll != 0:
                for name in self:
                    if name.split('_')[0] == 'c':
                        self[name].value = self[name].value - self.roll
            if value != 0:
                for name in self:
                    if name.split('_')[0] == 'c':
                        self[name].value = self[name].value + value
            self._roll = value
            
    def set_shift(self, value, type_='relative'):
            # is relative?
            if type_ in relative:
                value = self.shift + value

            if self.shift != value:
                if self.shift != 0:
                    for name in self:
                        if name.split('_')[0] == 'c':
                            self[name].value = self[name].value - self.shift
                if value != 0:
                    for name in self:
                        if name.split('_')[0] == 'c':
                            self[name].value = self[name].value + value
                self._shift = value

    def set_offset(self, value, type_='relative'):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).

        Returns:
            None
        """
        # to fix offset we have to add a background option
        pass

    def set_factor(self, value, type_='relative'):
        """Set y multiplicative factor.

        Args:
            value (number): multiplicative factor (y-coordinates will be
                multiplied by this value).

        Returns:
            None
        """
        pass
        # is relative?
        if type_ in relative:
            value = self.factor * value

        if self.factor != value:
            if self.factor != 1:
                for name in self:
                        if name.split('_')[0] == 'amp' or name.split('_') == 'area':
                            self[name].value  = self[name].value * self.factor**-1
            if value != 1:
                for name in self:
                    if name.split('_')[0] == 'amp' or name.split('_') == 'area':
                        self[name].value  = self[name].value * value
            self._factor = value

    def _reset_modifiers(self):
        self._factor = 1
        self._offset = 0
        self._calib  = 1
        self._shift  = 0
        self._roll   = 0

    # save and load (TODO)*****

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
            if len(name.split('_')) == 1:
                raise NotImplementedError('sorry, calling __getitem___ like this is not implemented yet.')
            if len(name.split('_')) == 2:
                return self.get_values(attr=name)
            if len(name.split('_')) == 3:
                if name in self:
                    return super().__getitem__(name)
                else:
                    raise ValueError(f'Cannot find parameter: {name}')
    
    def length(self):
        return len(self._get_indexes_i2())

    def append(self, i2='all', amp=None, c=None, w=None, w1=None, w2=None, m=0, m1=0, m2=0, area=None):

        # check parameters
        assert area is not None or amp is not None, 'amp or area not defined'
        assert c is not None, 'c not defined'
        assert w is not None or (w1 is not None and w2 is not None), 'Either w or (w1 and w2) must be defined'

        assert area is None or isinstance(area, numbers.Number), 'area must be a number'
        assert amp is None or isinstance(amp, numbers.Number), 'amp must be a number'
        assert c is None or isinstance(c, numbers.Number), 'c must be a number'
        assert w is None or isinstance(w, numbers.Number), 'w must be a number'
        assert w1 is None or isinstance(w1, numbers.Number), 'w1 must be a number'
        assert w2 is None or isinstance(w2, numbers.Number), 'w2 must be a number'
        assert m is None or isinstance(m, numbers.Number), 'm must be a number'
        assert m1 is None or isinstance(m1, numbers.Number), 'm1 must be a number'
        assert m2 is None or isinstance(m2, numbers.Number), 'm2 must be a number'

        if isinstance(w, numbers.Number):
            assert w >= 0, 'w must be positive'
        if isinstance(w1, numbers.Number):
            assert w1 >= 0, 'w1 must be positive'
        if isinstance(w2, numbers.Number):
            assert w2 >= 0, 'w2 must be positive'

        if isinstance(m, numbers.Number):
            assert m <= 1 and m >= 0, 'm must be positive less than 1'
        if isinstance(m1, numbers.Number):
            assert m1 <= 1 and m1 >= 0, 'm1 must be positive less than 1'
        if isinstance(m2, numbers.Number):
            assert m2 <= 1 and m2 >= 0, 'm2 must be positive less than 1'

        # run
        if i2 == 'all':
            for i2 in self._get_indexes_i2():
                self.append(i2=i2, amp=amp, c=c, w=w, w1=w1, w2=w2, m=m, m1=m1, m2=m2, area=area)
            return
        else:
            # define index
            indexes_i1 = self._get_indexes_i1(i2=i2)
            if indexes_i1 == []:
                i1 = 0
            else:
                i1 = max(indexes_i1) + 1

            # add parameters
            if area is None:
                self.add(f'amp_{i1}_{i2}', value=amp, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)
            else:
                self.add(f'area_{i1}_{i2}', value=area, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)
            
            self.add(f'c_{i1}_{i2}',   value=c,   vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)

            if w1 is not None and w2 is not None:
                self.add(f'w1_{i1}_{i2}', value=w1, vary=True, min=0, max=np.inf, expr=None, brute_step=None)
                self.add(f'w2_{i1}_{i2}', value=w2, vary=True, min=0, max=np.inf, expr=None, brute_step=None)
                self.add(f'm1_{i1}_{i2}', value=m1, vary=False, min=0, max=1, expr=None, brute_step=None)
                self.add(f'm2_{i1}_{i2}', value=m1, vary=False, min=0, max=1, expr=None, brute_step=None)
            else:
                self.add(f'w_{i1}_{i2}', value=w, vary=True, min=0, max=np.inf, expr=None, brute_step=None)
                self.add(f'm_{i1}_{i2}', value=m, vary=False, min=0, max=1, expr=None, brute_step=None)

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
            names = temp._get_names_with_index(i2=i2)
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

    def _get_names_with_index(self, i1=None, i2=None):
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
                names = self._get_names_with_index(i1=i1, i2=i2)
                for name in names:
                    _name = '_'.join(name.split('_')[:2])
                    final[i2][_name] = temp[name] 
        return final
    
        

    def _get_params_with_index(self, i1, i2):
        """Return peaks with all names with i2."""
        names = self._get_names_with_index(i1=i1, i2=i2)
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
        """returns a list with values of a peak for all spectra"""
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
        for name in self._get_names_with_index(i1=i1, i2=i2):
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
    
    def check_feasibility(self):
        """Check if start values are within boundaries."""
        for p in self:
            value = self[p].value
            vmax  = self[p].max
            vmin  = self[p].min

            if value > vmax or value < vmin:
                raise ValueError(f'Start value for parameter {p} is out of bounds\nbound max = {vmax}\nstart value = {value}\nbound min = {vmin}')

    def fit(self, xs, ys, yerr=None, method='least_squares', return_fitted_peaks=False):
        """output, peaks = fit()
        
        I think x and y must be monotonic.
        """
        self.check_feasibility()

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

