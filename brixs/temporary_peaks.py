# %%

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
import lmfit
from collections.abc import Iterable, MutableMapping
import json

# BRIXS
import brixs as br

# %% common definitions ===========================================================
names = ['amp', 'area', 'c', 'w', 'w1', 'w2', 'm', 'm1', 'm2']

br.cc   = ['cross-correlation', 'cc']
br.roll = ['roll', 'rotate', 'r', 'rot']
br.hard = ['hard', 'x', 'h', 'Hard']
br.soft = ['soft', 'Soft', 'interp', 'y', 's']
br.relative = ['relative', 'r', 'rel']
br.absolute = ['a', 'abs', 'absolute']
br.increasing = ['inc', 'i', 'up', 'increasing', 'increasingly']
br.decreasing = ['dec', 'd', 'down', 'decreasing', 'decreasingly']

# %% Peaks =====================================================================
class Peak(lmfit.Parameters):
    def __init__(self, *args, **kwargs):
        super().__init__()

        # when init, avoid run update_area_amp()
        self._enable_update_area_amp = False

        # index
        if 'index' in kwargs:
            self._index = kwargs.pop('index')
        else:
            self._index = 0

        # initial definitions
        self._use_area  = False
        self._step      = None
        self._asymmetry = False
        self._use_peak  = True
        # self.bkg        = ''

        # parameters
        self.add(f'amp_{self.index}',  value=0)
        self.add(f'area_{self.index}', value=0, vary=False)
        self.add(f'c_{self.index}',    value=0)
        self.add(f'w_{self.index}',    value=0, min=0)
        self.add(f'w1_{self.index}',   value=0, min=0, vary=False)
        self.add(f'w2_{self.index}',   value=0, min=0, vary=False)
        self.add(f'm_{self.index}',    value=0, min=0, max=1)
        self.add(f'm1_{self.index}',   value=0, min=0, max=1, vary=False)
        self.add(f'm2_{self.index}',   value=0, min=0, max=1, vary=False)

        # use area
        if 'use_area' in kwargs:
            self.use_area = kwargs.pop('use_area')
        # asymmetry
        if 'asymmetry' in kwargs:
            self.asymmetry = kwargs.pop('asymmetry')
        # step
        if 'step' in kwargs:
            self._step = kwargs.pop('step')


        # values
        if 'amp' in kwargs:
            self['amp'].value  = kwargs.pop('amp')
        if 'area' in kwargs:
            self['area'].value = kwargs.pop('area')
        if 'c' in kwargs:
            self['c'].value   = kwargs.pop('c')
        if 'w' in kwargs:
            self['w'].value   = kwargs.pop('w')
        if 'w1' in kwargs:
            self['w1'].value  = kwargs.pop('w1')
        if 'w2' in kwargs:
            self['w2'].value  = kwargs.pop('w2')
        if 'm' in kwargs:
            self['m'].value   = kwargs.pop('m')
        if 'm1' in kwargs:
            self['m1']  = kwargs.pop('m1')
        if 'm2' in kwargs:
            self['m2']  = kwargs.pop('m2')

        # modifiers
        if 'shift' in kwargs:
            self._shift = kwargs.pop('shift')
        else:
            self._shift = 0
        if 'shift_interp' in kwargs:
            self._shift_interp = kwargs.pop('shift_interp')
        else:
            self._shift_interp = 0
        if 'shift_roll' in kwargs:
            self._shift_roll = kwargs.pop('shishift_rollft')
        else:
            self._shift_roll = 0
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

        # check
        self._enable_update_area_amp = True
        self.update_area_amp()

    def add(self, name, value=None, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None):
        # check if name is valid
        if not name.startswith(tuple(names)):
            raise ValueError(f'name= {name} is not valid.\nValid names:{names}')
        
        # add
        super().add(name, value=value, vary=vary, min=min, max=max, expr=expr, brute_step=brute_step)

    def __setitem__(self, name, value):
        # check if name is valid
        if not name.startswith(tuple(names)):
            raise ValueError(f'name= {name} is not valid.\nValid names:{names}')

        # area and amp
        if self._enable_update_area_amp:
            if self.use_area:
                if name == 'amp':
                    raise ValueError('amp cannot be edited if Peak.use_area is True.\nUse Peak.update_area_amp()\nOr change Peak.use_area to False.')
            else:
                if name == 'area':
                    raise ValueError('area cannot be edited if Peak.use_area is False.\nUse Peak.update_area_amp().\nOr change Peak.use_area to True.')
        
        # set values
        super().__setitem__(name, value)
    
        # if name == 'w':
        #     self['w1'].value = self['w']/2
        #     self['w2'].value = self['w']/2

        # check
        if self._enable_update_area_amp:
            self.update_area_amp()

    def __getitem__(self, name):
        if name in names:
            name = name + f'_{self.index}'
        return super().__getitem__(name)

    # support
    def _find_suitable_x(self):
        if self.asymmetry:
            w = self['w1'].value + self['w2'].value
        else:
            w = self['w'].value

        if w == 0:
            w = self['c'].value*0.1
        if w == 0:
            w = 1

        vmin = self['c'].value-10*w
        vmax = self['c'].value+10*w

        if vmin == 0 and vmax == 0:
            vmin = -10
            vmax = 10
        num = int((vmax-vmin)/w*20)+1
        return np.linspace(vmin, vmax, num)

    # attributes and properties
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

    @property
    def use_area(self):
        return self._use_area
    @use_area.setter
    def use_area(self, value):
        assert isinstance(value, bool), 'value must be True or False.'
        self._use_area = value
        if value:
            for item in ['area']:
                self[item+f'_{self.index}'].vary = True 
            for item in ['amp']:
                self[item+f'_{self.index}'].vary = False
        else:
            for item in ['area']:
                self[item+f'_{self.index}'].vary = False 
            for item in ['amp']:
                self[item+f'_{self.index}'].vary = True
    @use_area.deleter
    def use_area(self):
        raise AttributeError('Cannot delete object.')

    @property
    def asymmetry(self):
        return self._asymmetry
    @asymmetry.setter
    def asymmetry(self, value):
        assert isinstance(value, bool), 'value must be True or False.'
        self._asymmetry = value
        if value:
            for item in ['w1', 'w2', 'm1', 'm1']:
                self[item+f'_{self.index}'].vary = True 
            for item in ['w', 'm']:
                self[item+f'_{self.index}'].vary = False
        else:
            for item in ['w1', 'w2', 'm1', 'm1']:
                self[item+f'_{self.index}'].vary = False 
            for item in ['w', 'm']:
                self[item+f'_{self.index}'].vary = True
    @asymmetry.deleter
    def asymmetry(self):
        raise AttributeError('Cannot delete object.')

    @property
    def use_peak(self):
        return self._use_peak
    @use_peak.setter
    def use_peak(self, value):
        assert isinstance(value, bool), 'value must be True or False.'
        self._use_peak = value
    @use_peak.deleter
    def use_peak(self):
        raise AttributeError('Cannot delete object.')

    @property
    def index(self):
        return self._index
    @index.setter
    def index(self, value):
        assert isinstance(value, int), 'value must be int.'
        if self._index != value:
            self._enable_update_area_amp = False
            for name in names:
                # self[name].name = name+f'_{value}'
                old_name = name+f'_{self._index}'
                self.add(name+f'_{value}', value=self[old_name].value, 
                                        min=self[old_name].min, 
                                        max=self[old_name].max, 
                                        vary=self[old_name].vary, 
                                        expr=self[old_name].expr)
                del self[old_name]
            self._index = value
            self._enable_update_area_amp = True
    @index.deleter
    def index(self):
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
        if type_ in br.relative:
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
        if type_ in br.relative:
            value = self.shift + value
        
        # apply
        if mode in br.roll: 
            value = self.step*value

            if self.shift_roll != value:
                if self.shift_roll != 0:
                    self['c'].value = self['c'].value - self.shift_roll
                if value != 0:
                    self['c'].value = self['c'].value + value
                self._shift_roll = value
        elif mode in br.hard:
            if self.shift != value:
                if self.shift != 0:
                    self['c'].value = self['c'].value - self.shift
                if value != 0:
                    self['c'].value = self['c'].value + value
                self._shift = value
        elif mode in br.soft:
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
        if type_ in br.relative:
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
        if type_ in br.relative:
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

    # calculation and info
    def update_area_amp(self):
        """Updates the area amp values.

        if use_area = True, amp is updated based on the area.
        if use_area = False, area is updated based on the amp.

        Only necessary if amp and area values are changed by hand.        
        """
        if self.asymmetry:
            c1 = self['m1'].value/np.pi + (1-self['m1'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            c2 = self['m2'].value/np.pi + (1-self['m2'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            if self.use_area:
                self['amp'].value = self['area'].value * 2 /(self['w1'].value*c1 + self['w2'].value*c2)
            else:
                self['area'].value = self['amp'].value * (self['w1'].value*c1 + self['w2'].value*c2)/2
        else:
            c = self['m'].value/np.pi + (1-self['m'].value)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
            if self.use_area:
                self['amp'].value = self['area'].value/self['w'].value/c
            else:
                self['area'].value = self['amp'].value * self['w'].value * c

    # plot and visualization
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

        self.update_area_amp()

        # data
        c = self['c'].value
        amp = self['amp'].value
        xerr = self['fwhm'].value/2

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
        f = self.model()
        s = br.Spectrum(x=x, y=f(x))

        # copy modifiers
        s._shift  = self.shift
        # s._shift_roll    = self.shift_roll  # cannot copy shift_roll because step might be different
        s._shift_interp  = self.shift_interp
        s._factor = self.factor
        s._calib  = self.calib
        s._offset = self.offset
        return s

    def _split_peaks(self):
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

    # calculation and info
    def _model_str(self):
        """Returns string for building peak function.

        Returns:
            function f(x) as string
        """
        i = self.index

        if self.asymmetry:
            if self.use_area:
                return f"np.heaviside(parvals['c_{i}']-x, 0)*br.voigt_area_fwhm(x, parvals['area_{i}'], parvals['c_{i}'], parvals['w1_{i}'],  parvals['m1_{i}']) + np.heaviside(x-parvals['c_{i}'], 0)*br.voigt_area_fwhm(x, parvals['area_{i}'], parvals['c_{i}'],  parvals['w2_{i}'],  parvals['m2_{i}']) + dirac_delta(x, parvals['amp_{i}'], parvals['c_{i}'])" 
            else:
                return f"np.heaviside(parvals['c_{i}']-x, 0)*br.voigt_fwhm(x, parvals['amp_{i}'], parvals['c_{i}'], parvals['w1_{i}'],  parvals['m1_{i}']) + np.heaviside(x-parvals['c_{i}'], 0)*br.voigt_fwhm(x, parvals['amp_{i}'], parvals['c_{i}'],  parvals['w2_{i}'],  parvals['m2_{i}']) + dirac_delta(x, parvals['amp_{i}'], parvals['c_{i}'])" 
        else:
            if self.use_area:
                return f"br.voigt_area_fwhm(x, parvals['area_{i}'], parvals['c_{i}'], parvals['w_{i}'], parvals['m_{i}'])"
            else:
                return f"br.voigt_fwhm(x, parvals['amp_{i}'], parvals['c_{i}'], parvals['w_{i}'], parvals['m_{i}'])"

    def model(self):
        """Returns a function f(x) for the peak."""
        parvals = self.valuesdict()
        model_str = f'lambda x, parvals: {self._model_str()}'
        final = eval(model_str)
        return lambda x: final(x, parvals)

    #############

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

    # OBSOLETE

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

# %%
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

    def _find_suitable_x(self):
        vmin  = []
        vmax  = []
        step  = []
        for peak in self:
            temp = peak._find_suitable_x()
            vmin.append(min(temp))
            vmax.append(max(temp))
            step.append(abs(temp[1]-temp[0]))
        num = int((max(vmax)-min(vmin))/max(step))+1
        return np.linspace(min(vmin), max(vmax), num)

    # attributes and properties
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
            # print(type(self[i]))
            temp[i] = self[i].shift
        return temp
    @shift.setter
    def shift(self, value):
        self.set_shifts(value, mode='x')
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
        self.set_shifts(value, mode='roll')
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
        self.set_shifts(value, mode='interp')
    @shift_interp.deleter
    def shift_interp(self):
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

    # save and load
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
    def append(self, value):
        """Append peak.

        Args:
            value (Peak or dict): peak to be appended.

        Returns:
            None
        """
        # value.pretty_print()
        assert isinstance(value, Peak), 'Only type br.Peak can be appended.'
        next_index = len(self)
        value.index = next_index
        self._store.append(value)
        

        # update modifiers
        self._store[-1]._shift         = self.shift[0]
        self._store[-1]._shift_roll    = self.shift_roll[0]
        self._store[-1]._shift_interp  = self.shift_interp[0]
        self._store[-1]._calib   = self.calib[0]
        self._store[-1]._offset  = self.offset[0]
        self._store[-1]._factor  = self.factor[0]

        # # check
        # self.reorder()

    def remove(self, key):
        """Remove peak.

        Args:
            key (int or list): peaks to be removed.

        Returns:
            None
        """
        del self._store[key]

    # modifiers
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

    def set_shifts(self, value, mode, type_='relative'):
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
            self[i].set_shifts(value=value[i], mode=mode, type_=type_)

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

    # plot and visualization
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

        r = [None]*len(self)
        for i in range(len(self)):
            r[i] = self[i].plot(ax=ax, offset=offset, shift=shift, factor=factor, **kwargs)
            plt.text(self[i]['c'].value*calib+shift, self[i]['amp'].value*factor+offset, i, fontsize=14)
        return r

    # calculation and info
    def _model_str(self):
        for i in range(len(self)):
            if i == 0:
                model_str = self[i]._model_str() + '+'
                p = self[i]
            else:
                model_str += self[i]._model_str() + '+'
                p += self[i]
        return model_str[:-1]

    def _add_peaks(self):
        for i in range(len(self)):
            if i == 0:
                p = self[i]
            else:
                p += self[i]
        return p
   
    def model(self):
        """Returns a function f(x) for the peak."""
        parvals = self._add_peaks().valuesdict()        
        model_str = f'lambda x, parvals: {self._model_str()}'
        final = eval(model_str)
        return lambda x: final(x, parvals)

    def generate_residual_function(self, x, y):
        model_str = self._model_str()  
        p = self._add_peaks()

        def residual(p, x, y):
            parvals = p.valuesdict()
            model = eval(model_str)
            return model-y
        
        return residual
    
    def fit(self, x, y, method='least_squares'):
        """output, peaks = fit()"""
        residual = self.generate_residual_function(x, y)
        out = lmfit.minimize(residual, method='least_squares', params=self._add_peaks(), args=(x, y))
        return out, out.params._split_peaks()
    
    # extractors
    def calculate_spectra(self, x=None):
        """Return each peak spectrum separately.

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
    
    # OBSOLETE  
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


# %%

# %% data to plot
x = np.linspace(-2, 6, 100)
y1 = br.voigt_fwhm(x=x, amp=10, c=2, w=1, m=0)#+br.voigt_fwhm(x=x, amp=12, c=4, w=1, m=0)
y2 = br.voigt_fwhm(x=x, amp=10, c=2, w=1, m=0)+br.voigt_fwhm(x=x, amp=12, c=4, w=1, m=0)


# %% peak 1 ++++++++++++
p1 = Peak()
p1.index = 0
# p1 = lmfit.Parameters()
# p1.add('amp')
# p1.add('c')
# p1.add('w')
# p1.add('m')
# p1.add('w1')
# p1.add('w2')
# p1.add('m1')
# p1.add('m2')
# p1.add('area')
# p1['w1'].vary = False
# p1['w2'].vary = False
# p1['m'].vary  = False
# p1['m1'].vary = False
# p1['m2'].vary = False
# p1['area'].vary = False

p1['amp'].value = 8
p1['c'].value = 1.5
p1['w'].value = 0.3

# p1._find_suitable_x()

p1.pretty_print()
# s = p1.calculate_spectrum()
# s.plot()
# plt.show()

# %% peak 2 ===================
p2 = Peak()
p2.index = 1

p2['amp'].value = 6
# p2['amp'].expr = 'amp_0'
# p2['amp'].expr = ''
p2['c'].value = 5
p2['w'].value = 0.3

# p2.pretty_print()
# s = p2.calculate_spectrum()
# s.plot()
# plt.show()

# %% fit
peaks = Peaks(p1, p2)

out, peaks2 = peaks.fit(x, y2)

yfit = peaks2.calculate_spectrum()

plt.figure()
plt.plot(x, y2, marker='o')
yfit.plot()
plt.show()

# %%