#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Model object for fitting

`Model` is a lmfit.Parameters object (note the `s` in Parameters). 

The parent of `Model` is always a br.Spectrum or br.Spectra object.


brixs     Model   Component  Entries(i1)   Name    Tag   lmfit-parameter

                                        ┌── amp ── amp1 ── amp1_0_(i2)
                            ┌── entry 0 ├──  c  ──  c1  ──  c1_0_(i2)
                            |           └──  w  ──  w1  ──  w1_0_(i2)
                            |
                            |           ┌── amp ── amp1 ── amp1_0_(i2)
                 ┌──  peak  ├── entry 1 ├──  c  ──  c1  ──  c1_0_(i2)
                 |          |           └──  w  ──  w1  ──  w1_0_(i2)
                 |          └── ...  
                 |
Spectra ┐        |          ┌── entry 0
Spectum ├─ Model ├── arctan ├── entry 1
                 |          └── ...
                 |          ┌── entry 0
                 ├──  step  ├── entry 1
                 |          └── ...
                 └──  ...


#####################
# indexes i2 and i1 #
#####################

i2 is the spectrum number. When using br.Spectra object we have multiple spectra.
The index i2 indicates which entries are associated with which spectrum in br.Spectra.

i1 is the entry number. For example, one spectrum may have multiple peaks. Each 
peak will have a number i1.  

Typically, there will be two combinations of i1 and i2 that make sense when 
setting the arguments of a function/method:

(1) i2=None, i1='all'
    If i2 is None, i2 is assumed to be unique. Otherwise, raise an error
    If i1 is defined, check if i1 exists. Otherwise, raise an error
    Use self._check_i2_i1()

(2) i2='all', i1=None
    If i1 is None, i1 will be replaced by the max(i1) + 1
    If i1 is defined, check if i1 exits. If yes, raise an error.
    If i2 is defined, check if i2 exists. Otherwise, raise an error.
    Use self._check_i2_i1()
    I expect that we only use this one inside the `add()` function.


############################
# Creating a new Component #
############################

1. Every time a new component is created. It needs to be include inside Model.components.

2. A new component must have _ComponentTemplate as a parent class 
(i.e. every `Component` must have at least the methods and attributes defined 
in this template class.

#### setting up the __init__ function ####

__init__(self, parent) -> must receive the brixs object as an argument (either 
br.Spectrum or br.Spectra)

One must run _ComponentTemplate._initialize() inside __init__()

_ComponentTemplate._initialize() requires 6 things to run:

parent: Model object
nametags: a dictionary for converting from names to tags. Names are parameter
    names that will be used by the user. Tags are also parameters names which are
    going to be used by lmfit. Names must be simple and easy to remember. Tags
    can be complicated as long the do not have the underscore ('_') character. 
    Also, Tags must be unique even between different `Components`.
function: a function f(i1, i2) that returns a string. This string must be a 
    mathematical expression y = f(x, parameters) which will return a value based on x
    and the parameters.
    Tags must be accompanied of the indexes i1 and i2. See Peaks() and Linear() 
     Components below.
min, max, step: a function f(i1, i2) that returns a list of three strings. The 
    first string must be an expression that returns the minimum x value necessary
    to consider when 'plotting' the y = f(x, parameters). Likewise, the second
    string is the maximum value necessary and the third value is a reasonable x
    step. For example (See Peaks() and Linear() Components), for a `Component` where y = f(x, amp, c, w)) = gaussian_peak(x, amp, c, w)
    where amp is the amplitude of the peak, c is the center, and w is the width,
    we can define the minimum value the c-10*w (the center of the peak
    minus 10 times the width of the peak) and the maximum value as c+10*w (the 
    center of the peak plus 10 time the width of the peak. 
    This way, the most 'interesting' part 
    of the function f(x), i.e. the part where the peak is, will be contained inside
    min < x < max. Finally, the step value is the x step size necessary so we 
    have enough x point to describe/visualize y=f(x). In the example (See Peaks() and Linear() Components), a reasonable 
    step size would be w/20, i.e., the 'peak' will have 20 x points within width.

#### setting up the add() method ####

The add() methods inside a component is responsible for saving new parameters 
inside the parent model. See Peaks() and Linear() Components as an example. 


#### setting up the base-modifier functions ####

When one use set_shift() in br.Spectrum.set_shift(), the spectrum shifts by a value.
If Components have a set_shift() function defined, then whenever a spectrum is shifted,
then the parameters relating to that Component also shift. See Peaks() and Linear() 
Components as an example. This is also valid for set_offset(), set_factor(), and
set_calib(). Note that these functions are not required for a component to work.



########        
# TODO #
########        
s.model.peaks.find_peaks()
save and load
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from collections.abc import Iterable
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import copy

# %% -------------------------- Special Imports --------------------------- %% #
from scipy.signal import find_peaks
import lmfit

# %% --------------------------- brixs Imports ---------------------------- %% #
import brixs as br
# %%

# %% ------------------------- support functions -------------------------- %% #
def _name_parser(name):
    temp = name.split('_')
    name = temp[0]
    i1   = int(temp[1])
    i2   = int(temp[2])
    return name, i1, i2

# %% -------------------------------- Model ------------------------------- %% #
class Model_(lmfit.Parameters):
    
    def __init__(self, parent=None):  # talk about why parent has to be None (only from within the fit - minimizer - function)
        super().__init__()
        self.parent = parent

        ##############
        # components #
        ##############
        self.components = {'peaks':  Peaks,
                           'linear': Linear,
                           'asymmetricpeaks': AsymmetricPeaks,
                           'arctan': Arctan,
                           }

        # check repeated tags
        unique_tags = []
        for component in self.components:
            for tag in list(self.components[component](None).nametags.values()):
                if tag in unique_tags:
                    raise KeyError(f'tag `{tag}` in `{component}` already exist. Tags must be unique. Please, change this tag')
                else:
                    unique_tags.append(tag)
        del unique_tags

        # start components
        for name, obj in self.components.items():
            self.__setattr__(name, obj(self))

    #################
    # magic methods #
    #################
    def __setitem__(self, value, _value):
        super().__setitem__(value, _value)
    
    def __getitem__(self, value):
        #######
        # int #
        #######
        # if isinstance(value, int):
        #     return self._get_parameters_for_given_i2(i2=value) 
        ##########
        # string #
        ##########
        if isinstance(value, str):
            return super().__getitem__(value)
        ##########
        # error #
        ##########
        else:
            raise ValueError('key must be str')# or int')
        
    def __delitem__(self, value):
        """delete"""
        super().__delitem__(value)
        # self.remove(i2=value)

    #########
    # basic #
    #########
    def remove_parameters_for_given_i2(self, i2):
        """remove all components with index i2"""
        # valid i2
        i2s = self._get_list_of_indexes_i2()
        assert i2 in i2s, f'i2={i2} not found\navailable i2: {i2s}'

        # remove
        for name in self:
            if name.split('_')[2] == str(i2):
                del self[name]

    def get_available_components(self):
        """returns list of available components"""
        return list(self.components.keys())

    ###########
    # support #
    ###########
    def _check_if_model_have_only_one_i2(self):
        """return True if model has only one i2, i.e., parameters for only one spectrum"""
        if len(self._get_list_of_indexes_i2()) < 2: return True
        else: return False
    
    def _get_list_of_indexes_i2(self):
        """return list of all secondary (spectrum number) indexes i2
        
        Returns:
            list
        """
        return list(np.unique([int(name.split('_')[2]) for name in self]))

    # _get_indexes_i1_for_each_component_and_each_i2
    def _get_indexes_dict(self):
        """get indexes

        Returns:
            dict of dicts of list
                {component name: {i2: [i1]}}
        """
        return {name: self.__getattribute__(name)._get_indexes_dict() for name in self.components}
    
    # def _get_parameters_for_given_i2(self, i2):
    #     """returns parameters for a given i2

    #     Args:
    #         i2 (int or None, optional): i2 index

    #     Returns:
    #         model object
    #     """
    #     # valid i2
    #     i2s = self._get_list_of_indexes_i2()
    #     assert i2 in i2s, f'i2={i2} not found\navailable i2: {i2s}'
        
    #     # get item
    #     final = Model(self.parent)
    #     for tag in self:
    #         if tag.split('_')[2] == str(i2):
    #             final[tag] = self[tag]
    #     return final 
    
    #########
    # model #
    #########
    def check_feasibility(self):
        """Check if start values are within boundaries."""
        final = []
        for name in self:
            value = self[name].value
            vmax  = self[name].max
            vmin  = self[name].min

            if value > vmax or value < vmin:
                final.append(f'parameter `{name}`\n    bound max = {vmax}\n    start value = {value}\n    bound min = {vmin}\n')
        if final != []:
            raise ValueError(f'Feasibility error. Start values outside boundaries:\n\n' + '\n'.join(final))
    
    def get_model_str(self, i2=None):
        """Returns string for building model.

        Args:
            i2 (int, optional): secondary (spectrum) index i2. If None, i2 must 
                be unique. Default is None.

        Returns:
            two strings: variables and function f(x) as string
        """
        # check i2
        if i2 is None:
            if self._check_if_model_have_only_one_i2():
                i2 = self._get_list_of_indexes_i2()[0]
            else:
                raise ValueError('i2 not defined')
        else:
            assert i2 in self._get_list_of_indexes_i2(), f'i2={i2} not found'
        
        # Initialization
        indexes   = self._get_indexes_dict()
        final     = ''
        variables = ''

        # model
        for component in indexes:
            if i2 in indexes[component]:
                for _i1 in indexes[component][i2]:
                    _variables, _final = self.__getattribute__(component).get_model_str(i1=_i1, i2=i2)
                    final     += _final     + ' + '
                    variables += _variables + ', '
        variables = variables.replace('self.parent[', 'self[')
        return variables[:-2], final[:-3]

    def get_model(self, i2=None):
        """Returns a function f(x) for the peak.
        
        Args:
            i2 (int, optional): secondary (spectrum) index i2. If None, i2 must 
                be unique. Default is None.

        Returns:
            function f(x)
        """
        variables, _model = self.get_model_str(i2=i2)
        model = f'lambda x, {variables}: {_model}'
        return eval(model)

    def fit(self, method='least_squares', limits=None, update=True):
        """fit spectrum with model function and return lmfit.minimizer output

        Warning: 

        Warning:
            I am unsure if the data to be fitted (x, y) needs to be monotonic. 
            After some testing, it does not seem to be a requirement for lmfit.
            However, I would keep an eye on it since this is a requirement for
            scipy.curve_fit().

        Args:
            method (str, optional): Name of the fitting method to use. See methods
                available on `lmfit.minimize()`_ documentation.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            update (bool, optional): if True, Model will be updated with 
                fitted values.
        
        Returns:
            `lmfit.minimize()`_ output object  

        .. _lmfit.minimize(): https://lmfit.github.io/lmfit-py/fitting.html      
        """
        ###############
        # feasibility #
        ###############
        self.check_feasibility()

        ##########
        # ranges #
        ##########
        temp = self.parent.copy(limits=limits)

        ################
        # x and y data #
        ################
        if isinstance(temp, br.Spectrum):
            xs = [temp.x, ]
            ys = [temp.y, ]
        elif isinstance(temp, br.Spectra):
            xs = [s.x for s in temp]
            ys = [s.y for s in temp]
        else:
            raise('there is no parent (br.Spectrum/br.Spectra) object associated with this fitting model')

        #####################
        # residual function #
        #####################
        def residual(params, xs, ys):
            """residual function for lmfit.minimize function

            Args:
                params (lmfit.parameters)
                xs, ys (list): list with x and y arrays

            Returns:
                residual
            """
            residual = []

            # make residual per data set
            for i2, x in enumerate(xs):
                residual.append(ys[i2] - params.get_model(i2=i2)(x))

            residual = np.concatenate(residual).ravel()
            return residual
        
        ##################
        # lmfit minimize #
        ##################
        out = lmfit.minimize(residual, method=method, params=self, kws={'xs':xs, 'ys':ys})        
        
        #####################
        # update parameters #
        #####################
        if update:
            for name in self:
                self[name] = out.params[name]
        return lmfit.fit_report(out), out
    
    ########
    # copy #
    ########
    def copy(self, parent=None):
        model = copy.deepcopy(self)
        model.parent = parent
        return model

    #################
    # save and load #
    #################
    # TODO

    ##################
    # base modifiers #
    ##################
    def set_shift(self, value):
        """Shift peaks to components

        Args:
            value (float or int): shift value (value will be added to center c).

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        # if isinstance(value, Iterable) == False:
        #     value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        # assert len(value) == len(self._get_list_of_indexes_i2()), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self._get_list_of_indexes_i2())}'

        ##################
        # applying value #
        ##################
        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:  # checking if there are any components 
                obj.set_shift(value)

    def set_calib(self, value):
        """Set calibration value to components

        Args:
            value (number): calibration value (centers and widths will be multiplied
                by this value).

        Returns:
            None
        """      
        ###################################
        # asserting validity of the input #
        ################################### 
        if value == 0:  # calib cannot be zero
            raise ValueError('cannot set calib = 0')
        elif value == 1:   # if calib is 1, do nothing
            return
        
        ##################
        # applying value #
        ##################
        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:
                obj.set_calib(value)
                
    def set_offset(self, value):
        """Set offset value to component

        Args:
            value (value): offset value.

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if value == 0:
            return

        ##################
        # applying value #
        ##################
        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:           
                obj.set_offset(value)

    def set_factor(self, value):
        """Set multiplicative factor to component.

        Args:
            value (number): multiplicative factor (amplitudes and areas will be
                multiplied by this value).

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if value == 0:
            raise ValueError('cannot set factor = 0.')
        elif value == 1:
            return

        ##################
        # applying value #
        ##################
        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:                      
                obj.set_factor(value)

    ######################
    # calculate spectrum #
    ######################
    def get_suitable_x(self, i2=None, i1='all'):
        """estimates a suitable array with x values for plotting.

        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.

        Returns:
            array
        """
        # empty model
        assert len(self) > 0, 'cannot estimate suitable x in a empty model'

        vmin = None
        vmax = None
        step = None

        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:
                x = obj.get_suitable_x(i2=i2, i1=i1, verbose=False)
                if vmin is None or vmin > min(x):
                    vmin = min(x)
                if vmax is None or vmax < max(x):
                    vmax = max(x)
                if step is None or step > x[1]-x[0]:
                    step = x[1]-x[0]

            # _vmin, _vmax, _step = obj._get_min_max_step(i2=i2, i1=i1)
            # _vmin = _vmin
            # _vmax = _vmax
            # _step = _step
            # if vmin is None or vmin > _vmin:
            #     vmin = _vmin
            # if vmax is None or vmax < _vmax:
            #     vmax = _vmax
            # if step is None or step < _step:
            #     step = _step

        return np.arange(vmin, vmax, step)
    
    def calculate_spectrum(self, i2=None, x=None):
        """Return spectrum

        Args:
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be calculated

        Returns:
            :py:class:`Spectrum`.
        """        
        # get model
        model = self.get_model(i2=i2)

        # get x
        if x is None:
            x = self.get_suitable_x(i2=i2)

        return br.Spectrum(x=x, y=model(x))

    def calculate_spectra(self, x=None):
        """Return a Spectra object with spectra for each i2 index 

        Note:
            If model is absent for a certain i2, this will return an empty spectrum
        
        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be constructed.

        Returns:
            :py:class:`Spectra`
        """
        ss = br.Spectra()
        for i2 in range(len(self.parent)):
            if i2 not in self._get_list_of_indexes_i2():
                ss.append(br.Spectrum())
            else:
                if x is None:
                    x = self.get_suitable_x(i2=i2)
                ss.append(self.calculate_spectrum(i2=i2, x=x))
        return ss

    def check_values_close_to_bounds(self, max_error=0.1):
        """Return dict with parameter values close (or equal to boundaries)
        
            Args:
                max_error (number, optional): percentage value 
                    of the max distance allowed between parameters value and 
                    boundary conditions. Default is 0.01 %
            
            Returns:
                Dict
        """
        final = {}
        for parameter in self:
            if self[parameter].vary:
                value = self[parameter].value
                _min  = self[parameter].min
                _max  = self[parameter].max
                if value + value*max_error/100 >= _max:
                    final[parameter] = 'max'
                    if value - value*max_error/100 <= _min:
                        final[parameter] = 'narrow'
                elif value - value*max_error/100 <= _min:
                    final[parameter] = 'min'

        return final

    ##########################        
    # plot and visualization #
    ##########################
    def pretty_print(self, *args, **kwargs):
        """Wrapper for `lmfit.parameter.Parameters.pretty_print()`_.

        .. _lmfit.parameter.Parameters.pretty_print(): https://lmfit.github.io/lmfit-py/parameters.html#lmfit.parameter.Parameters.pretty_print
        """
        if len(self) == 0:
            return None
        return super().pretty_print(*args, **kwargs)

    def plot(self, i2=None, ax=None, offset=0, shift=0, roll=0, factor=1, calib=1, smooth=1, label=None, limits=None, switch_xy=False, **kwargs):
        """Plot spectrum. Wrapper for `matplotlib.pyplot.plot()`_.

        Note:
            If `label` is `None` and if spectrum have attr 
            `label`, this attr will be used as label, e.g., 
            `plt.plot(s.x, s.y, label=s.label)`.  

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            calib, shift (number, optional): multiplicative and additive factor
                 on the x-coordinates. calib is applied first.
            factor, offset (number, optional): multiplicative and additive factor
                 on the y-coordinates. Factor is applied first.
            roll (int, optional): Roll value of array elements of the x-coordinates
            smooth (int, optional): number of points to average data. Default is 1.
            label (str, number, optional): if str or number, this label will be 
                applied. If None and if spectrum have attr `label`, 
                this attr will be used as label, e.g., `plt.plot(s.x, s.y, label=s.label)`.
                Default is None. 
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            switch_xy (bool, optional): Switch x and y axis.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            `Line2D`_

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        ################
        # get spectrum #
        ################
        x = self.get_suitable_x(i2=i2)
        f = self.get_model(i2=i2)
        y = f(x)
        s = br.Spectrum(x=x, y=y)

        ########
        # plot #
        ########
        return s.plot(ax=ax, offset=offset, shift=shift, roll=roll, factor=factor, calib=calib, smooth=smooth, label=label, limits=limits, switch_xy=switch_xy, **kwargs)


class Model(lmfit.Parameters):

    def __init__(self, parent=None):  # talk about why parent has to be None (only from within the fit - minimizer - function)
        super().__init__()
        self.parent = parent

        ##############
        # components #
        ##############
        self.components = {'peaks':  Peaks,
                           'linear': Linear,
                           'asymmetricpeaks': AsymmetricPeaks,
                           'arctan': Arctan,
                           'erf': Erf,
                           }

        # check repeated tags
        unique_tags = []
        for component in self.components:
            for tag in list(self.components[component](None).nametags.values()):
                if tag in unique_tags:
                    raise KeyError(f'tag `{tag}` in `{component}` already exist. Tags must be unique. Please, change this tag')
                else:
                    unique_tags.append(tag)
        del unique_tags

        # start components
        for name, obj in self.components.items():
            self.__setattr__(name, obj(self))

    #################
    # magic methods #
    #################
    def __setitem__(self, value, _value):
        super().__setitem__(value, _value)
    
    def __getitem__(self, value):
        #######
        # int #
        #######
        if isinstance(value, int):
            return self._get_dict_for_given_i2(i2=value) 
        ##########
        # string #
        ##########
        elif isinstance(value, str):
            return super().__getitem__(value)
        ##########
        # error #
        ##########
        else:
            raise ValueError('key must be str or int')
        
    def __delitem__(self, value):
        """delete"""
        super().__delitem__(value)
        # self.remove(i2=value)

    #########
    # basic #
    #########
    def remove_parameters_for_given_i2(self, i2):
        """remove all components with index i2"""
        # valid i2
        i2s = self._get_list_of_indexes_i2()
        assert i2 in i2s, f'i2={i2} not found\navailable i2: {i2s}'

        # remove
        for name in self:
            if name.split('_')[2] == str(i2):
                del self[name]

    def get_available_components(self):
        """returns list of available components"""
        return list(self.components.keys())

    def copy_from_spectra(self):
        assert isinstance(self.parent, br.Spectra), 'parent object must be of type br.Spectra'
        print('still not sure this works 100%')

        for i2, _s in enumerate(self.parent):
            for name in self.parent[i2].model:
                attr = self.parent[i2].model[name]
                new_name = '_'.join(name.split('_')[:2]) + '_' + str(i2)
                self[new_name] = attr
                self[new_name].name = new_name
                self[new_name].expr = ''
        return

    def split_parameters_per_spectra(self):
        assert isinstance(self.parent, br.Spectra), 'parent object must be of type br.Spectra'
        print('still not sure this works 100%')

        # remove expr from all parameters
        temp = copy.deepcopy(self)
        for attr in temp:
            temp[attr].expr = ''

        for i2, _s in enumerate(self.parent):
            for name in temp:
                _name, _i1, _i2 = _name_parser(name)
                if _i2 == i2:
                    new_name = _name + '_' + str(_i1) + '_0'
                    _s.model[new_name] = temp[name]
        return

    ###########
    # support #
    ###########
    def _check_if_model_have_only_one_i2(self):
        """return True if model has only one i2, i.e., parameters for only one spectrum"""
        if len(self._get_list_of_indexes_i2()) < 2: return True
        else: return False
    
    def _get_list_of_indexes_i2(self):
        """return list of all secondary (spectrum number) indexes i2
        
        Returns:
            list
        """
        return list(np.unique([int(name.split('_')[2]) for name in self]))

    # _get_indexes_i1_for_each_component_and_each_i2
    def _get_indexes_dict(self):
        """get indexes

        Returns:
            dict of dicts of list
                {component name: {i2: [i1]}}
        """
        return {name: self.__getattribute__(name)._get_indexes_dict() for name in self.components}
    
    def _get_dict_for_given_i2(self, i2):
        i2s = self._get_list_of_indexes_i2()
        assert i2 in i2s, f'i2={i2} not found\navailable i2: {i2s}'
        final = {}
        for tag in self:
            if tag.split('_')[2] == str(i2):
                final[tag] = self[tag].value
        return final
    
    # def _get_parameters_for_given_i2(self, i2):
    #     """returns parameters for a given i2

    #     Args:
    #         i2 (int or None, optional): i2 index

    #     Returns:
    #         model object
    #     """
    #     # valid i2
    #     i2s = self._get_list_of_indexes_i2()
    #     assert i2 in i2s, f'i2={i2} not found\navailable i2: {i2s}'
        
    #     # get item
    #     final = Model(self.parent)
    #     for tag in self:
    #         if tag.split('_')[2] == str(i2):
    #             temp = self[tag]
    #             temp.expr = ''
    #             final[tag] = temp
    #     return final
    
    #########
    # model #
    #########
    def check_feasibility(self):
        """Check if start values are within boundaries."""
        final = []
        for name in self:
            value = self[name].value
            vmax  = self[name].max
            vmin  = self[name].min

            if value > vmax or value < vmin:
                final.append(f'parameter `{name}`\n    bound max = {vmax}\n    start value = {value}\n    bound min = {vmin}\n')
        if final != []:
            raise ValueError(f'Feasibility error. Start values outside boundaries:\n\n' + '\n'.join(final))
    
    def get_model_str(self, i2=None):
        """Returns string for building model.

        Args:
            i2 (int, optional): secondary (spectrum) index i2. If None, i2 must 
                be unique. Default is None.

        Returns:
            two strings: variables and function f(x) as string
        """
        # check i2
        if i2 is None:
            if self._check_if_model_have_only_one_i2():
                i2 = self._get_list_of_indexes_i2()[0]
            else:
                raise ValueError('i2 not defined')
        else:
            assert i2 in self._get_list_of_indexes_i2(), f'i2={i2} not found'
        
        # Initialization
        indexes   = self._get_indexes_dict()
        final     = ''
        variables = ''

        # model
        for component in indexes:
            if i2 in indexes[component]:
                for _i1 in indexes[component][i2]:
                    _variables, _final = self.__getattribute__(component).get_model_str(i1=_i1, i2=i2)
                    final     += _final     + ' + '
                    variables += _variables + ', '
        variables = variables.replace('self.parent[', 'self[')
        return variables[:-2], final[:-3]

    def get_model(self, i2=None):
        """Returns a function f(x) for the peak.
        
        Args:
            i2 (int, optional): secondary (spectrum) index i2. If None, i2 must 
                be unique. Default is None.

        Returns:
            function f(x)
        """
        variables, _model = self.get_model_str(i2=i2)
        model = f'lambda x, {variables}: {_model}'
        return eval(model)

    def fit(self, method='least_squares', limits=None, update=True):
        """fit spectrum with model function and return lmfit.minimizer output

        Warning: 

        Warning:
            I am unsure if the data to be fitted (x, y) needs to be monotonic. 
            After some testing, it does not seem to be a requirement for lmfit.
            However, I would keep an eye on it since this is a requirement for
            scipy.curve_fit().

        Args:
            method (str, optional): Name of the fitting method to use. See methods
                available on `lmfit.minimize()`_ documentation.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            update (bool, optional): if True, Model will be updated with 
                fitted values.
        
        Returns:
            `lmfit.minimize()`_ output object  

        .. _lmfit.minimize(): https://lmfit.github.io/lmfit-py/fitting.html      
        """
        ###############
        # feasibility #
        ###############
        self.check_feasibility()

        ##########
        # ranges #
        ##########
        temp = self.parent.copy(limits=limits)

        ################
        # x and y data #
        ################
        if isinstance(temp, br.Spectrum):
            xs = [temp.x, ]
            ys = [temp.y, ]
        elif isinstance(temp, br.Spectra):
            xs = [s.x for s in temp]
            ys = [s.y for s in temp]
        else:
            raise('there is no parent (br.Spectrum/br.Spectra) object associated with this fitting model')

        #####################
        # residual function #
        #####################
        def residual(params, xs, ys):
            """residual function for lmfit.minimize function

            Args:
                params (lmfit.parameters)
                xs, ys (list): list with x and y arrays

            Returns:
                residual
            """
            residual = []

            # make residual per data set
            for i2, x in enumerate(xs):
                residual.append(ys[i2] - params.get_model(i2=i2)(x))

            residual = np.concatenate(residual).ravel()
            return residual
        
        ##################
        # lmfit minimize #
        ##################
        out = lmfit.minimize(residual, method=method, params=self, kws={'xs':xs, 'ys':ys})        
        
        #####################
        # update parameters #
        #####################
        if update:
            for name in self:
                self[name] = out.params[name]
        return lmfit.fit_report(out), out
    
    ########
    # copy #
    ########
    def copy(self, parent=None):
        model = copy.deepcopy(self)
        model.parent = parent
        # if isinstance(parent, br.Spectra):
        #     for _s in parent:
        #         _s.model
        return model
    
    def copy_from(self, object):
        # assert isinstance(self, object), 'cannot copy model from an'
        if isinstance(object, br.Spectrum) or isinstance(object, br.Spectra):
            model = object.model
        elif isinstance(object, Model):
            model = object
        parent = self.parent 
        self = copy.deepcopy(model)
        self.parent = parent

    #################
    # save and load #
    #################
    def save(self, filepath):
        """Writes JSON representation of Parameters to a file"""
        filepath = Path(filepath)
        with filepath.open('w') as f:
            final = self.dump(f, indent=4)
        return final

    def load(self, filepath):
        """Loads JSON representation of Parameters from a file"""
        filepath = Path(filepath)
        with filepath.open('r') as f:
            _ = super().load(f)
        return
    
    ##################
    # base modifiers #
    ##################
    def set_shift(self, value):
        """Shift peaks to components

        Args:
            value (float or int): shift value (value will be added to center c).

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        # if isinstance(value, Iterable) == False:
        #     value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        # assert len(value) == len(self._get_list_of_indexes_i2()), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self._get_list_of_indexes_i2())}'

        ##################
        # applying value #
        ##################
        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:  # checking if there are any components 
                obj.set_shift(value)

    def set_calib(self, value):
        """Set calibration value to components

        Args:
            value (number): calibration value (centers and widths will be multiplied
                by this value).

        Returns:
            None
        """      
        ###################################
        # asserting validity of the input #
        ################################### 
        if value == 0:  # calib cannot be zero
            raise ValueError('cannot set calib = 0')
        elif value == 1:   # if calib is 1, do nothing
            return
        
        ##################
        # applying value #
        ##################
        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:
                obj.set_calib(value)
                
    def set_offset(self, value):
        """Set offset value to component

        Args:
            value (value): offset value.

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if value == 0:
            return

        ##################
        # applying value #
        ##################
        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:           
                obj.set_offset(value)

    def set_factor(self, value):
        """Set multiplicative factor to component.

        Args:
            value (number): multiplicative factor (amplitudes and areas will be
                multiplied by this value).

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if value == 0:
            raise ValueError('cannot set factor = 0.')
        elif value == 1:
            return

        ##################
        # applying value #
        ##################
        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:                      
                obj.set_factor(value)

    ######################
    # calculate spectrum #
    ######################
    def get_suitable_x(self, i2=None, i1='all'):
        """estimates a suitable array with x values for plotting.

        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.

        Returns:
            array
        """
        # empty model
        assert len(self) > 0, 'cannot estimate suitable x in a empty model'

        vmin = None
        vmax = None
        step = None

        for component in self.components:
            obj = self.__getattribute__(component)
            if len(obj._get_all_tags()) > 0:
                x = obj.get_suitable_x(i2=i2, i1=i1, verbose=False)
                if vmin is None or vmin > min(x):
                    vmin = min(x)
                if vmax is None or vmax < max(x):
                    vmax = max(x)
                if step is None or step > x[1]-x[0]:
                    step = x[1]-x[0]

            # _vmin, _vmax, _step = obj._get_min_max_step(i2=i2, i1=i1)
            # _vmin = _vmin
            # _vmax = _vmax
            # _step = _step
            # if vmin is None or vmin > _vmin:
            #     vmin = _vmin
            # if vmax is None or vmax < _vmax:
            #     vmax = _vmax
            # if step is None or step < _step:
            #     step = _step

        return np.arange(vmin, vmax, step)
    
    def calculate_spectrum(self, i2=None, x=None):
        """Return spectrum

        Args:
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be calculated

        Returns:
            :py:class:`Spectrum`.
        """        
        # get model
        model = self.get_model(i2=i2)

        # get x
        if x is None:
            x = self.get_suitable_x(i2=i2)

        return br.Spectrum(x=x, y=model(x))

    def calculate_spectra(self, x=None):
        """Return a Spectra object with spectra for each i2 index 

        Note:
            If model is absent for a certain i2, this will return an empty spectrum
        
        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be constructed.

        Returns:
            :py:class:`Spectra`
        """
        ss = br.Spectra()
        for i2 in range(len(self.parent)):
            if i2 not in self._get_list_of_indexes_i2():
                ss.append(br.Spectrum())
            else:
                if x is None:
                    x = self.get_suitable_x(i2=i2)
                ss.append(self.calculate_spectrum(i2=i2, x=x))
        return ss

    def check_values_close_to_bounds(self, max_error=0.1):
        """Return dict with parameter values close (or equal to boundaries)
        
            Args:
                max_error (number, optional): percentage value 
                    of the max distance allowed between parameters value and 
                    boundary conditions. Default is 0.01 %
            
            Returns:
                Dict
        """
        final = {}
        for parameter in self:
            if self[parameter].vary:
                value = self[parameter].value
                _min  = self[parameter].min
                _max  = self[parameter].max
                if value + value*max_error/100 >= _max:
                    final[parameter] = 'max'
                    if value - value*max_error/100 <= _min:
                        final[parameter] = 'narrow'
                elif value - value*max_error/100 <= _min:
                    final[parameter] = 'min'

        return final

    ##########################        
    # plot and visualization #
    ##########################
    def pretty_print(self, *args, **kwargs):
        """Wrapper for `lmfit.parameter.Parameters.pretty_print()`_.

        .. _lmfit.parameter.Parameters.pretty_print(): https://lmfit.github.io/lmfit-py/parameters.html#lmfit.parameter.Parameters.pretty_print
        """
        if len(self) == 0:
            return None
        return super().pretty_print(*args, **kwargs)

    def plot(self, i2=None, ax=None, offset=0, shift=0, roll=0, factor=1, calib=1, smooth=1, label=None, limits=None, switch_xy=False, **kwargs):
        """Plot spectrum. Wrapper for `matplotlib.pyplot.plot()`_.

        Note:
            If `label` is `None` and if spectrum have attr 
            `label`, this attr will be used as label, e.g., 
            `plt.plot(s.x, s.y, label=s.label)`.  

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            calib, shift (number, optional): multiplicative and additive factor
                 on the x-coordinates. calib is applied first.
            factor, offset (number, optional): multiplicative and additive factor
                 on the y-coordinates. Factor is applied first.
            roll (int, optional): Roll value of array elements of the x-coordinates
            smooth (int, optional): number of points to average data. Default is 1.
            label (str, number, optional): if str or number, this label will be 
                applied. If None and if spectrum have attr `label`, 
                this attr will be used as label, e.g., `plt.plot(s.x, s.y, label=s.label)`.
                Default is None. 
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            switch_xy (bool, optional): Switch x and y axis.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            `Line2D`_

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        ################
        # get spectrum #
        ################
        x = self.get_suitable_x(i2=i2)
        f = self.get_model(i2=i2)
        y = f(x)
        s = br.Spectrum(x=x, y=y)

        ########
        # plot #
        ########
        return s.plot(ax=ax, offset=offset, shift=shift, roll=roll, factor=factor, calib=calib, smooth=smooth, label=label, limits=limits, switch_xy=switch_xy, **kwargs)
# %%

# %% ------------------------- Component Template ------------------------- %% #
class _ComponentTemplate(object):

    ##################
    # initialization #
    ##################
    def _initialize(self, parent, nametags, function, _min, _max, _step):
        """
        """
        self.parent   = parent
        self.nametags = nametags
        self.tagnames = dict((v, k) for k, v in nametags.items())
        self.function = function
        self._min     = _min
        self._max     = _max
        self._step    = _step

    #########
    # basic #
    #########
    def add(self, *args, **kwargs):
        raise NotImplementedError(f'add function is not defined for model component of type {type(self)}')

    def remove(self, i2=None, i1='all'):
        """delete all parameters associated with an entry

        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.

        Returns:
            None
        """
        # get tags
        tags = self._get_tags_dict_for_given_i2_i1(i2=i2, i1=i1)

        # delete tags
        for _i1 in tags:
            for tag in tags[_i1]:
                del self.parent[tag]
    
    def clear(self):
        """remove all entries for this component"""
        tags = self._get_all_tags()
        for tag in tags:
            del self.parent[tag]

    #################
    # magic methods #
    #################
    def __repr__(self):
        return '\n'.join(["'" + self.tagnames[tag.split('_')[0]] + '_' + '_'.join(tag.split('_')[1:]) + "'" + ': ' + str(self.parent[tag])  for tag in self._get_all_tags()])

    def __setitem__(self, value, _value):
        ##########
        # string #
        ##########
        if isinstance(value, str):
            split = value.split('_')
            assert len(split) == 3, f'`{name}` is not a recognizable attr for component of type `{type(self)}`'
            name = split[0]
            i1   = split[1]
            i2   = split[2]

            # check if name is a valid attribute
            assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

            # check if i2 is valid
            assert br.numanip.is_integer(i2), f'i2 must be a integer, not {type(i2)}'
            i2 = int(i2)
            indexes = self._get_indexes_dict()
            assert i2 in indexes, f'i2={i2} is not valid\nvalid i2: {list(indexes.keys())}'

            # check i1 is valid
            assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
            i1 = int(i1)
            assert i1 in indexes[i2], f'i1={i1} is not valid\nvalid i1: {indexes[i2]}'

            # set parameter value
            self.parent[self.nametags[name] + '_' + str(i1) + '_' + str(i2)].value = _value
        ##########
        # error #
        ##########
        else:
            raise ValueError('key must be str or int')
        
    def __getitem__(self, value):
        """getitem

        for value int: value will be considered i1 (only works if i2 is unique)
        for value str: value must match the parameter name
        """
        #######
        # int #
        #######
        if isinstance(value, int):
            # return self._get_parameters_for_given_i2_i1(i1=value, i2=None) 
            return self._get_dict_for_given_i2_i1(i1=value, i2=None) 
        ##########
        # string #
        ##########
        elif isinstance(value, str):
            split = value.split('_')
        if len(split) == 1:
            name = split[0]

            # check if name is a valid attribute
            assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'
            
            # check i2 is unique
            assert self._check_if_component_have_only_one_i2(), 'for now, this syntax requires i2 to be unique'
            i2 = self._get_list_of_indexes_i2()[0]

            # check i1 is unique
            i1s = self._get_indexes_dict()[i2]
            assert len(i1s) == 1, 'for now, this syntax requires i1 to be unique'
            i1 = int(i1s[0])

            return self.parent[self.nametags[name] + '_' + str(i1) + '_' + str(i2)]

        elif len(split) == 3:
            # assert len(split) == 3, f'`{name}` is not a recognizable attr for component of type `{type(self)}`'
            name = split[0]
            i1   = split[1]
            i2   = split[2]

            # check if name is a valid attribute
            assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

            # check if i2 is valid
            assert br.numanip.is_integer(i2), f'i2 must be a integer, not {type(i2)}'
            i2      = int(i2)
            indexes = self._get_indexes_dict()
            assert i2 in indexes, f'i2={i2} is not valid\nvalid i2: {list(indexes.keys())}'

            # check i1 is valid
            assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
            i1 = int(i1)
            assert i1 in indexes[i2], f'i1={i1} is not valid\nvalid i1: {indexes[i2]}'

            # get parameter
            return self.parent[self.nametags[name] + '_' + str(i1) + '_' + str(i2)]
        ##########
        # error #
        ##########
        elif len(split) != 1 and len(split) != 3:
            raise KeyError(f'`{name}` is not a recognizable attr for component of type `{type(self)}`')
        else:
            raise KeyError('key must be str')# or int')
        
    def __delitem__(self, value):
        """cannot be used"""
        raise ValueError('Please, use Peaks.remove() to delete Component entries')

    ##################
    # base modifiers #
    ##################
    def set_shift(self, value):
        """Shift peaks.

        Adds a value to the center parameter.

        Args:
            value (float or int): shift value (value will be added to center c).

        Returns:
            None
        """
        #print(f'Warning: Cannot shift Component `{type(self).__name__}` as it does not have set_shift() implemented')
        ###################################
        # asserting validity of the input #
        ###################################
        # if value == 0:
        #     return

        # ##################
        # # applying value #
        # ##################
        # for name in self._get_all_tags():
        #     if name.split('_')[0] in self.modifier_instructions['shift']:
        #         self.parent[name].value += value
        return 
    
    def set_offset(self, value):
        """Set offset value (in fact, it does nothing)

        Applying an offset, should not change the amp or area of peaks.

        Args:
            value (value): offset value.

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        # if value == 0:
        #     return

        # ##################
        # # applying value #
        # ##################
        # for name in self._get_all_tags():
        #     if name.split('_')[0] in self.modifier_instructions['offset']:
        #         self.parent[name].value += value
        return
    
    def set_calib(self, value):
        """Set calibration value.

        Args:
            value (number): calibration value (centers and widths will be multiplied
                by this value).

        Returns:
            None
        """      
        ###################################
        # asserting validity of the input #
        ################################### 
        # if value == 0:  # calib cannot be zero
        #     raise ValueError('cannot set calib = 0')
        # elif value == 1:   # if calib is 1, do nothing
        #     return
        
        # ##################
        # # applying value #
        # ##################
        # for name in self._get_all_tags():
        #     raw = name.split('_')[0]
        #     if raw in self.modifier_instructions['calib']:
        #         self.parent[name].value = self[name].value*value
        #     elif raw in self.modifier_instructions['calib_abs']:
        #         self.parent[name].value = abs(self[name].value*value)
        return
                  
    def set_factor(self, value):
        """Set multiplicative factor.

        Args:
            value (number): multiplicative factor (amplitudes and areas will be
                multiplied by this value).

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        # if value == 0:
        #     raise ValueError('cannot set factor = 0.')
        # elif value == 1:
        #     return

        # ##################
        # # applying value #
        # ##################
        # for name in self._get_all_tags():
        #     if name.split('_')[0] in self.modifier_instructions['factor']:
        #         self.parent[name].value += value
        return
    
    ###########
    # support #
    ###########
    def _check_if_component_have_only_one_i2(self):
        """return True if component has only one i2, i.e., parameters for only one spectrum"""
        if len(list(self._get_indexes_dict().keys())) < 2: return True
        else: return False

    def _get_list_of_all_available_i2(self):
        """return list of all secondary (spectrum number) indexes i2 available
        
        Returns:
            list
        """
        if isinstance(self.parent.parent, br.Spectra):
            return np.arange(len(self.parent.parent))
        else:
            return [0, ]
    
    def _get_list_of_indexes_i2(self):
        """return list of all secondary (spectrum number) indexes i2 with peaks defined
        
        Returns:
            list
        """
        return list(np.unique([int(tag.split('_')[2]) for tag in self._get_all_tags()]))

    def _get_indexes_dict(self):
        """get indexes for this component

        Returns:
            dict of lists
                {i2: [i1]}
        """
        tags = self._get_all_tags()
        _i2  = self._get_list_of_indexes_i2()
        return {i2:list(np.unique([int(tag.split('_')[1]) for tag in tags if int(tag.split('_')[2]) == i2])) for i2 in _i2}
      
    def _get_dict_for_given_i2_i1(self, i2=None, i1='all'):
        """returns dict with parameters name and value for a given i1 and i2

        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.

        Returns:
            dict
        """
        # get tags
        tags = self._get_tags_dict_for_given_i2_i1(i2=i2, i1=i1)

        # get parameters
        final = {}
        for _i1 in tags:
            final.update({tag: self.parent[tag].value for tag in tags[_i1]})
        
        return final 

    # def _get_parameters_for_given_i2_i1(self, i2=None, i1='all'):
    #     """returns parameters for a given i1 and i2

    #     Args:
    #         i1 (int or str, optional): i1 index. Use `all` to get data from all
    #             valid i1 indexes. Default is `all`.
    #         i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
    #             unique. Default is None.

    #     Returns:
    #         a component type object
    #     """
    #     # get tags
    #     tags = self._get_tags_dict_for_given_i2_i1(i2=i2, i1=i1)

    #     # get parameters
    #     parent = {}
    #     for _i1 in tags:
    #         parent.update({tag: self.parent[tag] for tag in tags[_i1]})
        
    #     # final
    #     final = self.__class__(parent)
    #     return final 

    def _get_all_tags(self):
        """returns a list of all tags associated with this component"""
        return [tag for tag in self.parent if tag.split('_')[0] in self.tagnames]

    def _get_tags_dict_for_given_i2_i1(self, i2=None, i1='all'):
        """get tags for a given i1 and i2

        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.

        Returns:
            dict {i1: tags}
        """
        # check i1, i2
        i1, i2 = self._check_i2_i1(i2=i2, i1=i1)

        # get names
        tags = self._get_all_tags()

        # filter names
        if isinstance(i1, Iterable) == False:
            i1 = [i1, ]
        return {_i1: [tag for tag in tags if int(tag.split('_')[2]) == i2 and int(tag.split('_')[1]) == _i1] for _i1 in i1}

    def _check_i2_i1(self, i2=None, i1='all'):
        """check/return valid i1's based on i2. Raises error if i1 and i2 are not valid

        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.

        Note:
            usually used in function where the args are of type `i2=None, i1='all'`
            i.e. i2 must be defined or unique.

        Returns:
            valid i1 (int or list) and i2 (int)
        """
        indexes = self._get_indexes_dict()

        # check i2
        if i2 is None:
            if self._check_if_component_have_only_one_i2():
                i2 = list(indexes.keys())[0]
            else:
                raise ValueError('i2 not defined')
        else:
            assert br.numanip.is_integer(i2), f'i2 must be a integer, not {type(i2)}'
            assert i2 in indexes, f'i2={i2} not found'

        # check i1
        if i1 == 'all':
            i1 = indexes[i2]
        else:
            assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
            assert i1 in indexes[i2], f'i1 not found for i2={i2}'

        return i1, i2

    def _check_i1_i2(self, i1=None, i2='all'):
        """Check/return valid i2's based on i1. Raises error if i1 and i2 are not valid

        Args:
            i1 (int or str, optional): i1 index. If None, i1 will be defined as
                the max(i1)+1. Default is None.
            i2 (int, list, None, optional): i2 index. Use 'all' to add components to all i2 
                available. Can also be a list with i2 values. Default is 'all'

        Returns:
            valid i1 (int) and i2 (list)
        """
        # initial
        i2s     = self._get_list_of_all_available_i2()
        indexes = self._get_indexes_dict()

        # empty Spectrum().model
        if isinstance(self.parent.parent, br.Spectrum):
            if len(indexes) == 0:
                if i1 is None:
                    i1 = 0
                else:
                    assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
                return i1, [0, ]

        # check i2
        if i2 == 'all':
            i2 = i2s
        else:
            if isinstance(i2, Iterable):
                for _i2 in i2:
                    assert br.numanip.is_integer(_i2), f'i2={_i2} must be a integer, not {type(i2)}'
                    assert _i2 in i2s, f'Spectra does not have i2={_i2} not found. available i2={i2s}'
            else:
                assert br.numanip.is_integer(i2), f'i2 must be a integer, not {type(i2)}'
                assert i2 in i2s, f'Spectra does not have i2={i2} not found. available i2={i2s}'
                i2 = [i2, ]

        # check i1
        if i1 is None:
            i1 = -1
            for _i2 in i2:
                if _i2 in indexes:
                    if len(indexes[_i2]) > 0:
                        if i1 < max(indexes[_i2]):
                            i1 = max(indexes[_i2])
            i1 += 1
            # i1 = max([max(indexes[_i2]) for _i2 in i2 if _i2 in indexes]) + 1
        else:
            assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
            
            check = {}
            for _i2 in i2:
                check[_i2] = False
                if _i2 in indexes:
                    if i1 in indexes[_i2]:
                        check[_i2] = True
            assert True not in list(check.values()), f'i1={i1} already exists for some i2\n{check}'
        return i1, i2
    
    #########
    # model #
    #########
    def check_feasibility(self):
        """Check if start values are within boundaries."""
        final = []
        for tag in self._get_all_tags():
            value = self.parent[tag].value
            vmax  = self.parent[tag].max
            vmin  = self.parent[tag].min

            if value > vmax or value < vmin:
                final.append(f'type `{type(self)}`\npar `{self.tagnames[tag]}`\ntag `{tag}`\n    max = {vmax}\n    start = {value}\n    min = {vmin}\n')
        if final != []:
            raise ValueError(f'Feasibility error. Start values outside boundaries:\n\n' + '\n'.join(final))
    
    def get_model_str(self, i2=None, i1='all'):
        """Returns string for building model.

        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.

        Returns:
            two strings: variables and function f(x) as string
        """
        # Initialization
        final     = ''
        variables = ''

        # check i1, i2
        i1, i2 = self._check_i2_i1(i2=i2, i1=i1)
        
        # model
        if isinstance(i1, Iterable):
            for _i1 in i1:
                _variables, _final  = self.get_model_str(i1=_i1, i2=i2)
                final     +=  _final     + ' + '
                variables +=  _variables + ', '
        else:
            # final = f"br.voigt_fwhm(x, amp_{i1}_{i2}, c_{i1}_{i2}, w_{i1}_{i2}, m_{i1}_{i2}) + "
            final = self.function(i1=i1, i2=i2) + ' + '
            for tag in self.tagnames:
                variables += f'{tag}_{i1}_{i2} = self.parent["{tag}_{i1}_{i2}"].value, '
        return variables[:-2], final[:-3]

    def get_model(self, i2=None, i1='all'):
        """Returns a function f(x) for the peak.
        
        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.

        Returns:
            function f(x)
        """
        variables, _model = self.get_model_str(i1=i1, i2=i2)
        model = f'lambda x, {variables}: {_model}'
        return eval(model)
   
    #############
    # modifiers #
    #############
    def reorder_i1_by_parameter_value(self, name, i2='all', increasing=True):
        """Change component i1 indexes based on parameter value
        
        args:
            name (str): parameter name
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.
            increasing (bool, optional): Default is True.

        Returns:
            None        
        """
        raise NotImplementedError('sorry, not implemented yet.')

    #################
    # save and load #
    #################
    # TODO

    ######################
    # calculate spectrum #
    ######################
    def _get_suitable_min_max_step_str(self, i2=None, i1='all'):
        """return string function to get min, max, and step values
        
        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.

        Returns:
            variables, vmin, vmax, step
        """
        # # check if empty
        # assert len(self._get_all_tags()) > 0, 'cannot calculate suitable min max step for empty Component'
    
        # check i1, i2
        i1, i2 = self._check_i2_i1(i2=i2, i1=i1)

        # Initialization
        vmin      = 'np.min([_ for _ in ['
        vmax      = 'np.max([_ for _ in ['
        step      = 'np.min([_ for _ in ['
        variables = ''

        # min_max_step
        if isinstance(i1, Iterable):
            for _i1 in i1:
                _variables, _vmin, _vmax, _step = self._get_suitable_min_max_step_str(i1=_i1, i2=i2)
                if _vmin is not None:
                    vmin  +=  _vmin + ', '
                if _vmax is not None:
                    vmax  +=  _vmax + ', '
                if _step is not None:
                    step  +=  _step + ', '
                variables +=  _variables
            return variables[:-2], vmin[:-2] + '] if _ is not None])', vmax[:-2] + '] if _ is not None])', step[:-2] + '] if _ is not None])'
        else:
            vmin = self._min(i1=i1, i2=i2)
            vmax = self._max(i1=i1, i2=i2)
            step = self._step(i1=i1, i2=i2)
            for tag in self.tagnames:
                variables += f'{tag}_{i1}_{i2} = self.parent["{tag}_{i1}_{i2}"].value, '
            return variables, vmin, vmax, step

    def _get_min_max_step(self, i2=None, i1='all', verbose=True):
        """return min, max, and step values
        
        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.
            verbose(bool, optional): if True, warns when min, max, and step values 
                cannot be estimated unambiguously. Default is True.

        Returns:
            vmin, vmax, step
        """
        variables, _vmin, _vmax, _step = self._get_suitable_min_max_step_str(i1=i1, i2=i2)
        vmin = f'lambda {variables}: {_vmin}'
        vmax = f'lambda {variables}: {_vmax}'
        step = f'lambda {variables}: {_step}'

        flag = False
        try:
            vmin = eval(vmin)()
        except ValueError:
            vmin = -1
            flag = True
        try:
            vmax = eval(vmax)()
        except ValueError:
            vmax = 1
            flag = True
        try:
            step = eval(step)()
        except ValueError:
            step = 2/100
            flag = True

        if flag and verbose: print('warning: suitable x cannot be found unambiguously. Defaulting to x = np.arange(-1, 1, 200)')
        # return eval(vmin), eval(vmax), eval(step)
        return vmin, vmax, step
    
    def get_suitable_x(self, i2=None, i1='all', verbose=True):
        """returns a array with x values for plotting. Wrapper for `np.arange(vmin, vmax, step)`

        Note:
            if min, max, and step values used in np.arange(min, max, step) cannot 
                be estimated unambiguously, it defaults to np.arange(-1, 1, 0.2)

        Args:
            i1 (int or str, optional): i1 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
            i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None.
            verbose(bool, optional): if True, warns when min, max, and step values 
                used in np.arange(min, max, step) cannot be estimated unambiguously. 
                Default is True.

        Returns:
            x array
        """
        # check if empty
        assert len(self._get_all_tags()) > 0, 'cannot calculate suitable x for empty Component'
    
        # return suitable x
        vmin, vmax, step = self._get_min_max_step(i2=i2, i1=i1, verbose=verbose)
        return np.arange(vmin, vmax, step)
    
    def calculate_spectrum(self, i2=None, i1='all', x=None):
        """Return spectrum

        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be calculated

        Returns:
            :py:class:`Spectrum`.
        """        
        # get model
        model = self.get_model(i2=i2, i1=i1)

        # get x
        if x is None:
            x = self.get_suitable_x(i2=i2, i1=i1)

        return br.Spectrum(x=x, y=model(x))

    def calculate_spectrum_for_each_i1(self, i2=None, x=None):
        """return spectra with spectra for each component for a certain i2
        
        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be constructed.

        Returns:
            :py:class:`Spectra`
        """
        # check i1, i2
        i1, i2 = self._check_i2_i1(i2=i2, i1='all')
        _x = x

        # calculate spectra
        ss = br.Spectra()
        for _i1 in i1:
            if x is None:
                _x = self.get_suitable_x(i2=i2, i1=_i1)
            ss.append(self.calculate_spectrum(i2=i2, i1=_i1, x=_x))
        return ss

    def calculate_spectra(self, x=None):
        """Return a Spectra object with spectrum for each i2 index 

        If component is absent for a certain i2, this will return an empty spectrum
        
        Args:
            x (list, optional): x values to which the curve will be calculated.
                If None, a suitable x will be constructed.

        Returns:
            :py:class:`Spectra`
        """
        ss = br.Spectra()
        for i2 in self.parent._get_list_of_indexes_i2():
            if i2 not in self._get_list_of_indexes_i2():
                ss.append(br.Spectrum())
            else:
                if x is None:
                    x = self.get_suitable_x(i2=i2)
                ss.append(self.calculate_spectrum(i2=i2, x=x))
        return ss

    ##############
    # get values #
    ##############
    # def get_values_for_each_i2(self, name, i1=None):
    #     """Returns a list with values of a component parameter for each i2.

    #     Args:
    #         name (str): parameter name.
    #         i1   (int): i1 index (component number). If i1 is None, i1 is assumed
    #             to be unique for all i2. Default is None.
        
    #     Returns:
    #         list of values
    #     """
    #     # check if name is a valid attribute
    #     assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

    #     # check i1 is valid
    #     indexes = self._get_indexes_dict()
    #     if i1 is None:
    #         check = {i2: len(indexes[i2]) for i2 in indexes}
    #         assert sum(list(check.values())) == len(indexes), f'i1 is not unique for all i2\ni1 must be defined\nindexes = {indexes}'
    #     else:
    #         check = {i2: i1 in indexes[i2] for i2 in indexes}
    #         assert False not in list(check.values()), f'some i2 indexes do not have i1 component\ni2 error = {check}\nindexes = {indexes}'

    #     # get values
    #     return [self._get_parameters_for_given_i2_i1(i2=i2, i1=i1)[name].value for i2 in indexes]

    # def get_errors_for_each_i2(self, name, i1=None):
    #     """Returns a list with errors of a component parameter for each i2.

    #     Args:
    #         name (str): parameter name.
    #         i1   (int): i1 index (component number). If i1 is None, i1 is assumed
    #             to be unique for all i2. Default is None.
        
    #     Returns:
    #         list of values
    #     """
    #     # check if name is a valid attribute
    #     assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

    #     # check i1 is valid
    #     indexes = self._get_indexes_dict()
    #     if i1 is None:
    #         check = {i2: len(indexes[i2]) for i2 in indexes}
    #         assert sum(list(check.values())) == len(indexes), f'i1 is not unique for all i2\ni1 must be defined\nindexes = {indexes}'
    #     else:
    #         check = {i2: i1 in indexes[i2] for i2 in indexes}
    #         assert False not in list(check.values()), f'some i2 indexes do not have i1 component\ni2 error = {check}\nindexes = {indexes}'

    #     # get values
    #     return [self._get_parameters_for_given_i2_i1(i2=i2, i1=i1)[name].err for i2 in indexes]

    ##########################        
    # plot and visualization #
    ##########################
    def pretty_print(self):
        """print parameters"""
        print(self.__repr__())
        return

    def plot(self, i2=None, i1='all', ax=None, offset=0, shift=0, roll=0, factor=1, calib=1, smooth=1, label=None, limits=None, switch_xy=False, **kwargs):
        """Plot spectrum. Wrapper for `matplotlib.pyplot.plot()`_.

        Note:
            If `label` is `None` and if spectrum have attr 
            `label`, this attr will be used as label, e.g., 
            `plt.plot(s.x, s.y, label=s.label)`.  

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            calib, shift (number, optional): multiplicative and additive factor
                 on the x-coordinates. calib is applied first.
            factor, offset (number, optional): multiplicative and additive factor
                 on the y-coordinates. Factor is applied first.
            roll (int, optional): Roll value of array elements of the x-coordinates
            smooth (int, optional): number of points to average data. Default is 1.
            label (str, number, optional): if str or number, this label will be 
                applied. If None and if spectrum have attr `label`, 
                this attr will be used as label, e.g., `plt.plot(s.x, s.y, label=s.label)`.
                Default is None. 
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            switch_xy (bool, optional): Switch x and y axis.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            `Line2D`_

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        ################
        # get spectrum #
        ################
        x = self.get_suitable_x(i2=i2, i1=i1)
        f = self.get_model(i2=i2, i1=i1)
        y = f(x)
        s = br.Spectrum(x=x, y=y)

        ########
        # plot #
        ########
        return s.plot(ax=ax, offset=offset, shift=shift, roll=roll, factor=factor, calib=calib, smooth=smooth, label=label, limits=limits, switch_xy=switch_xy, **kwargs)
# %%

# %% ----------------------------- Components ----------------------------- %% #
# %% peaks
class Peaks(_ComponentTemplate):
    def __init__(self, parent):
        super().__init__()

        # parameter names
        # name = name used by the user
        # tag  = name used by lmfit
        # name: tag
        nametags = {'amp':'amp1', 'c':'c1', 'w':'w1', 'm':'m1'}

        # functions
        def function(i1, i2):
            """function that returns f(x)"""
            return f"br.model.voigt_fwhm(x, amp1_{i1}_{i2}, c1_{i1}_{i2}, w1_{i1}_{i2}, m1_{i1}_{i2})"
        
        def _min_suitable_x_str(i1, i2):
            """return min x value (as string) to be used for quickly plotting this component"""
            return f'c1_{i1}_{i2} - w1_{i1}_{i2}*10'
        
        def _max_suitable_x_str(i1, i2):
            """return max x value (as string) to be used for quickly plotting this component"""
            return f'c1_{i1}_{i2} + w1_{i1}_{i2}*10'
        
        def _suitable_x_step_str(i1, i2):
            """return x step value (as string) to be used for quickly plotting this component"""
            return f'w1_{i1}_{i2}/20'
    
        # initialization 
        self._initialize(parent=parent, nametags=nametags, function=function, _min=_min_suitable_x_str, _max=_max_suitable_x_str, _step=_suitable_x_step_str)

    #################
    # add parameter #
    #################
    def add(self, amp=None, c=None, w=None, m=0, 
            amp_min=-np.inf, amp_max=np.inf, amp_vary=True, amp_expr=None,
            c_min=-np.inf, c_max=np.inf, c_vary=True, c_expr=None,
            w_min=0, w_max=np.inf, w_vary=True, w_expr=None,
            m_min=0, m_max=1, m_vary=False, m_expr=None,
            i1=None, i2='all'):  
        """add parameters
        
        Args:
            i1 (int or str, optional): i1 index. If None, i1 will be defined as
                the max(i1)+1. Default is None.
            i2 (int, list, None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None. Use 'all' to add component to all i2 
                available.
        
        Return:
            None
        """        
        ###################
        # check i1 and i2 #
        ###################
        i1, i2 = self._check_i1_i2(i1=i1, i2=i2)

        #######################################
        # assert validity of parameter values #
        #######################################
        assert amp >= 0,          'amp cannot be negative'

        assert w >= 0,            'w cannot be negative'
        assert w_min >= 0,        'w_min cannot be negative'

        assert m >= 0 and m <= 1, 'm must be between 0 and 1'
        assert m_min >= 0 and m_min <= 1, 'm_min must be between 0 and 1'
        assert m_max >= 0 and m_max <= 1, 'm_max must be between 0 and 1'


        ##################
        # add parameters #
        ##################
        for _i2 in i2:
            tag = self.nametags['amp']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=amp, vary=amp_vary,  min=amp_min, max=amp_max, expr=amp_expr, brute_step=None)
            
            tag = self.nametags['c']
            self.parent.add(f'{tag}_{i1}_{_i2}',   value=c, vary=c_vary,  min=c_min, max=c_max, expr=c_expr, brute_step=None)
            
            tag = self.nametags['w']
            self.parent.add(f'{tag}_{i1}_{_i2}',   value=w, vary=w_vary,  min=w_min, max=w_max, expr=w_expr, brute_step=None)
            
            tag = self.nametags['m']
            self.parent.add(f'{tag}_{i1}_{_i2}',   value=m, vary=m_vary, min=m_min, max=m_max, expr=m_expr, brute_step=None)
        return

    ##################
    # base modifiers #
    ##################
    def set_shift(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'c1':
                # bounds
                self.parent[_name].min   += value[i2]
                self.parent[_name].max   += value[i2]
                # value
                self.parent[_name].value += value[i2]
                # error
                # self.parent[_name].stderr += abs(value[i2])

        return 
    
    def set_offset(self, value):
        pass
        return 
    
    def set_factor(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'amp1':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])
        return 
    
    def set_calib(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'c1':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])
                
            elif name == 'w1':
                # bounds
                self.parent[_name].min   *= abs(value[i2])
                self.parent[_name].max   *= abs(value[i2])
                # value
                self.parent[_name].value = abs(self.parent[_name].value*value[i2])
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr   *= abs(value[i2])
        return 

    #########
    # EXTRA #
    #########
    def find(self, prominence=5, width=4, moving_average_window=8):
        """uses scipy function find_peaks() to estimate number of peaks and positions

        Peak parameters are saved inside model.peaks()

        Args:
            prominence (number, optional): difference between baseline and 
                amplitude of the peaks
            width (int, optional): width of the peaks in terms of number of datapoints
            moving_average_window (int, optional): moving average before searching for peaks

        Returns:
            None
        """
        if isinstance(self.parent.parent, br.Spectra):
            raise NotImplementedError('find peaks is not implemented for br.Spectra.model.peaks yet. Only for br.Spectrum.model.peaks')

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
            y2 = br.arraymanip.moving_average(self.parent.parent.y, moving_average_window)
            x2 = br.arraymanip.moving_average(self.parent.parent.x, moving_average_window)
        else:
            x2 = self.parent.parent.x
            y2 = self.parent.parent.y

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
                amp = d['prominences'][i] + max([y2[d['right_bases'][i]], y2[d['left_bases'][i]]])
                c   = x2[peaks[i]]
                w   = abs(d['widths'][i]*np.mean(np.diff(self.parent.parent.x)))
                try:
                    self.add(amp=amp, c=c, w=w)
                except AssertionError:
                    pass
        except IndexError:
            pass
# %%

# %% peaks
class AsymmetricPeaks(_ComponentTemplate):
    def __init__(self, parent):
        super().__init__()

        # parameter names
        # name = name used by the user
        # tag  = name used by lmfit (cannot contain '_')
        # name: tag
        nametags = {'amp':'ampasy', 'c':'casy', 'w1':'w1asy', 'w2':'w2asy', 'm1':'m1asy', 'm2':'m2asy'}

        # functions
        def function(i1, i2):
            """function that returns f(x)"""
            # return f"br.model.voigt_fwhm(x, amp1_{i1}_{i2}, c1_{i1}_{i2}, w1_{i1}_{i2}, m1_{i1}_{i2})"
            return f"np.heaviside(casy_{i1}_{i2}-x, 0)*br.model.voigt_fwhm(x, ampasy_{i1}_{i2}, casy_{i1}_{i2}, w1asy_{i1}_{i2},  m1asy_{i1}_{i2}) + np.heaviside(x-casy_{i1}_{i2}, 0)*br.model.voigt_fwhm(x, ampasy_{i1}_{i2}, casy_{i1}_{i2},  w2asy_{i1}_{i2},  m2asy_{i1}_{i2}) + br.model.dirac_delta(x, ampasy_{i1}_{i2}, casy_{i1}_{i2})" 
        
        def _min_suitable_x_str(i1, i2):
            """return min x value (as string) to be used for quickly plotting this component"""
            return f'casy_{i1}_{i2} - w1asy_{i1}_{i2}*10'
        
        def _max_suitable_x_str(i1, i2):
            """return max x value (as string) to be used for quickly plotting this component"""
            return f'casy_{i1}_{i2} + w2asy_{i1}_{i2}*10'
        
        def _suitable_x_step_str(i1, i2):
            """return x step value (as string) to be used for quickly plotting this component"""
            return f'(w1asy_{i1}_{i2} + w2asy_{i1}_{i2})/(2*20)'
    
        # initialization 
        self._initialize(parent=parent, nametags=nametags, function=function, _min=_min_suitable_x_str, _max=_max_suitable_x_str, _step=_suitable_x_step_str)

    #################
    # add parameter #
    #################
    def add(self, amp=None, c=None, w1=None, w2=None, m1=0, m2=0, 
            amp_min=-np.inf, amp_max=np.inf, amp_vary=True, amp_expr=None,
            c_min=-np.inf, c_max=np.inf, c_vary=True, c_expr=None,
            w1_min=0, w1_max=np.inf, w1_vary=True, w1_expr=None,
            w2_min=0, w2_max=np.inf, w2_vary=True, w2_expr=None,
            m1_min=0, m1_max=1, m1_vary=False, m1_expr=None,
            m2_min=0, m2_max=1, m2_vary=False, m2_expr=None,
            i1=None, i2='all'):  
        """add parameters
        
        Args:
            i1 (int or str, optional): i1 index. If None, i1 will be defined as
                the max(i1)+1. Default is None.
            i2 (int, list, None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None. Use 'all' to add component to all i2 
                available.
        
        Return:
            None
        """        
        ###################
        # check i1 and i2 #
        ###################
        i1, i2 = self._check_i1_i2(i1=i1, i2=i2)

        #######################################
        # assert validity of parameter values #
        #######################################
        assert amp >= 0,          'amp cannot be negative'

        assert w1 >= 0,            'w1 cannot be negative'
        assert w2 >= 0,            'w2 cannot be negative'
        assert w1_min >= 0,        'w1_min cannot be negative'
        assert w2_min >= 0,        'w2_min cannot be negative'

        assert m1 >= 0 and m1 <= 1, 'm1 must be between 0 and 1'
        assert m2 >= 0 and m2 <= 1, 'm2 must be between 0 and 1'
        assert m1_min >= 0 and m1_min <= 1, 'm1_min must be between 0 and 1'
        assert m1_max >= 0 and m1_max <= 1, 'm1_max must be between 0 and 1'
        assert m2_min >= 0 and m2_min <= 1, 'm2_min must be between 0 and 1'
        assert m2_max >= 0 and m2_max <= 1, 'm2_max must be between 0 and 1'

        ##################
        # add parameters #
        ##################
        for _i2 in i2:
            tag = self.nametags['amp']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=amp, vary=amp_vary,  min=amp_min, max=amp_max, expr=amp_expr, brute_step=None)
            
            tag = self.nametags['c']
            self.parent.add(f'{tag}_{i1}_{_i2}',   value=c, vary=c_vary,  min=c_min, max=c_max, expr=c_expr, brute_step=None)
            
            tag = self.nametags['w1']
            self.parent.add(f'{tag}_{i1}_{_i2}',   value=w1, vary=w1_vary,  min=w1_min, max=w1_max, expr=w1_expr, brute_step=None)
            
            tag = self.nametags['w2']
            self.parent.add(f'{tag}_{i1}_{_i2}',   value=w2, vary=w2_vary,  min=w2_min, max=w2_max, expr=w2_expr, brute_step=None)
            
            tag = self.nametags['m1']
            self.parent.add(f'{tag}_{i1}_{_i2}',   value=m1, vary=m1_vary, min=m1_min, max=m1_max, expr=m1_expr, brute_step=None)
            
            tag = self.nametags['m2']
            self.parent.add(f'{tag}_{i1}_{_i2}',   value=m2, vary=m2_vary, min=m2_min, max=m2_max, expr=m2_expr, brute_step=None)
        return

    ##################
    # base modifiers #
    ##################
    def set_shift(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)
        
        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'casy':
                # bounds
                self.parent[_name].min   += value[i2]
                self.parent[_name].max   += value[i2]
                # value
                self.parent[_name].value += value[i2]
                # error
                # self.parent[_name].stderr += abs(value[i2])
        return 
    
    def set_offset(self, value):
        # super().set_offset(value)
        return 
    
    def set_factor(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'ampasy':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])

        return 
    
    def set_calib(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'casy':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])

            elif name == 'w1asy':
                # swap w1 with w2
                if value[i2] < 0:
                    _temp = self.parent[f'w1asy_{i1}_{i2}']
                    self.parent[f'w1asy_{i1}_{i2}'] = self.parent[f'w2asy_{i1}_{i2}']
                    self.parent[f'w2asy_{i1}_{i2}'] = _temp

                # bounds
                self.parent[_name].min   *= abs(value[i2])
                self.parent[_name].max   *= abs(value[i2])

                self.parent[f'w2asy_{i1}_{i2}'].min   *= abs(value[i2])
                self.parent[f'w2asy_{i1}_{i2}'].max   *= abs(value[i2])
  
                # value
                self.parent[_name].value = abs(self.parent[f'w1asy_{i1}_{i2}'].value*value[i2])
                self.parent[f'w2asy_{i1}_{i2}'].value = abs(self.parent[f'w2asy_{i1}_{i2}'].value*value[i2])
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr   *= abs(value[i2])
                if self.parent[f'w2asy_{i1}_{i2}'].stderr is not None:
                    self.parent[f'w2asy_{i1}_{i2}'].stderr   *= abs(value[i2])

        return 

    #########
    # EXTRA #
    #########
    pass
# %%

# %% linear polynomial
class Linear(_ComponentTemplate):
    def __init__(self, parent):
        super().__init__()

        # parameter names
        # name = name used by the user
        # tag  = name used by lmfit
        # name: tag
        nametags = {'linear':'a1', 'const':'a0'}  # make it immutable

        def function(i1, i2):
            """function that returns f(x)"""
            return f"a1_{i1}_{i2}*x + a0_{i1}_{i2}"
        
        def _min_suitable_x_str(i1, i2):
            """return min x value (as string) to be used for quickly plotting this component"""
            return 'None'
        
        def _max_suitable_x_str(i1, i2):
            """return max x value (as string) to be used for quickly plotting this component"""
            return 'None'
        
        def _suitable_x_step_str(i1, i2):
            """return x step value (as string) to be used for quickly plotting this component"""
            return 'None'
    
        # initialization 
        self._initialize(parent=parent, nametags=nametags, function=function, _min=_min_suitable_x_str, _max=_max_suitable_x_str, _step=_suitable_x_step_str)

    #################
    # add parameter #
    #################
    def add(self, linear, const, i1=None, i2='all'):  
        """add parameters
        
        Args:
            i1 (int or str, optional): i1 index. If None, i1 will be defined as
                the max(i1)+1. Default is None.
            i2 (int or None, optional): i2 index. Use `all` to get data from all
                valid i1 indexes. Default is `all`.
        
        Return:
            None
        """
        ###################
        # check i1 and i2 #
        ###################
        i1, i2 = self._check_i1_i2(i1=i1, i2=i2)

        #######################################
        # assert validity of parameter values #
        #######################################
        pass

        ##################
        # add parameters #
        ##################
        for _i2 in i2:
            tag = self.nametags['linear']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=linear, vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
            
            tag = self.nametags['const']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=const,   vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)    
        return

    ##################
    # base modifiers #
    ##################
    def set_shift(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)
        
        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'a0':
                # bounds
                # self.parent[_name].min   += value[i2]
                # self.parent[_name].max   += value[i2]
                self.parent[_name].min   -= self.parent[f'a1_{i1}_{i2}'].value*value[i2]
                self.parent[_name].max   -= self.parent[f'a1_{i1}_{i2}'].value*value[i2]
                # value
                self.parent[_name].value -= self.parent[f'a1_{i1}_{i2}'].value*value[i2]
        return 
    
    def set_offset(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)
        
        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'a0':
                # bounds
                self.parent[_name].min   += value[i2]
                self.parent[_name].max   += value[i2]
                # value
                self.parent[_name].value += value[i2]
        return 
    
    def set_factor(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'a1':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])

            if name == 'a0':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])
        return 
    
    def set_calib(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'a1':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= 1/value[i2]
                    self.parent[_name].max   *= 1/value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max/value[i2]
                    self.parent[_name].max   = _min/value[i2]
                
                # value
                self.parent[_name].value *= 1/(value[i2])
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= 1/abs(value[i2])
        return 

# %% Arctangent
class Arctan(_ComponentTemplate):
    def __init__(self, parent):
        super().__init__()

        # parameter names
        # name = name used by the user
        # tag  = name used by lmfit
        # name: tag
        nametags = {'amp':'amparc', 'c':'carc', 'w':'warc'}  # make it immutable

        def function(i1, i2):
            """function that returns f(x)"""
            return f"br.model.arctan_fwhm(x, amparc_{i1}_{i2}, carc_{i1}_{i2}, warc_{i1}_{i2})"
        
        def _min_suitable_x_str(i1, i2):
            """return min x value (as string) to be used for quickly plotting this component"""
            return f'carc_{i1}_{i2} - warc_{i1}_{i2}*10'
        
        def _max_suitable_x_str(i1, i2):
            """return max x value (as string) to be used for quickly plotting this component"""
            return f'carc_{i1}_{i2} + warc_{i1}_{i2}*10'
        
        def _suitable_x_step_str(i1, i2):
            """return x step value (as string) to be used for quickly plotting this component"""
            return f'warc_{i1}_{i2}/20'
    
        # initialization 
        self._initialize(parent=parent, nametags=nametags, function=function, _min=_min_suitable_x_str, _max=_max_suitable_x_str, _step=_suitable_x_step_str)

    #################
    # add parameter #
    #################
    def add(self, amp=None, c=None, w=None,
            amp_min=-np.inf, amp_max=np.inf, amp_vary=True, amp_expr=None,
            c_min=-np.inf, c_max=np.inf, c_vary=True, c_expr=None,
            w_min=0, w_max=np.inf, w_vary=True, w_expr=None,
            i1=None, i2='all'):  
        """add parameters
        
        Args:
            i1 (int or str, optional): i1 index. If None, i1 will be defined as
                the max(i1)+1. Default is None.
            i2 (int, list, None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None. Use 'all' to add component to all i2 
                available.
        
        Return:
            None
        """        
        ###################
        # check i1 and i2 #
        ###################
        i1, i2 = self._check_i1_i2(i1=i1, i2=i2)

        #######################################
        # assert validity of parameter values #
        #######################################
        assert amp >= 0,          'amp cannot be negative'

        assert w >= 0,            'w cannot be negative'
        assert w_min >= 0,        'w_min cannot be negative'

        ##################
        # add parameters #
        ##################
        for _i2 in i2:
            tag = self.nametags['amp']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=amp, vary=amp_vary,  min=amp_min, max=amp_max, expr=amp_expr, brute_step=None)
            
            tag = self.nametags['c']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=c, vary=c_vary,  min=c_min, max=c_max, expr=c_expr, brute_step=None)
            
            tag = self.nametags['w']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=w, vary=w_vary,  min=w_min, max=w_max, expr=w_expr, brute_step=None)
        return

    ##################
    # base modifiers #
    ##################
    def set_shift(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'carc':
                # bounds
                self.parent[_name].min   += value[i2]
                self.parent[_name].max   += value[i2]
                # value
                self.parent[_name].value += value[i2]

        return 
    
    def set_offset(self, value):
        pass
        return 
    
    def set_factor(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'amparc':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])
        return 
    
    def set_calib(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'carc':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])
            elif name == 'warc':
                # bounds
                self.parent[_name].min   *= abs(value[i2])
                self.parent[_name].max   *= abs(value[i2])
                # value
                self.parent[_name].value = abs(self.parent[_name].value*value[i2])
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr   *= abs(value[i2])
        return 

# %% Arctangent
class Erf(_ComponentTemplate):
    def __init__(self, parent):
        super().__init__()

        # parameter names
        # name = name used by the user
        # tag  = name used by lmfit
        # name: tag
        nametags = {'amp':'amperf', 'c':'cerf', 'w':'werf'}  # make it immutable

        def function(i1, i2):
            """function that returns f(x)"""
            return f"br.model.erf_fwhm(x, amperf_{i1}_{i2}, cerf_{i1}_{i2}, werf_{i1}_{i2})"
        
        def _min_suitable_x_str(i1, i2):
            """return min x value (as string) to be used for quickly plotting this component"""
            return f'cerf_{i1}_{i2} - werf_{i1}_{i2}*5'
        
        def _max_suitable_x_str(i1, i2):
            """return max x value (as string) to be used for quickly plotting this component"""
            return f'cerf_{i1}_{i2} + werf_{i1}_{i2}*5'
        
        def _suitable_x_step_str(i1, i2):
            """return x step value (as string) to be used for quickly plotting this component"""
            return f'werf_{i1}_{i2}/20'
    
        # initialization 
        self._initialize(parent=parent, nametags=nametags, function=function, _min=_min_suitable_x_str, _max=_max_suitable_x_str, _step=_suitable_x_step_str)

    #################
    # add parameter #
    #################
    def add(self, amp=None, c=None, w=None,
            amp_min=-np.inf, amp_max=np.inf, amp_vary=True, amp_expr=None,
            c_min=-np.inf, c_max=np.inf, c_vary=True, c_expr=None,
            w_min=0, w_max=np.inf, w_vary=True, w_expr=None,
            i1=None, i2='all'):  
        """add parameters
        
        Args:
            i1 (int or str, optional): i1 index. If None, i1 will be defined as
                the max(i1)+1. Default is None.
            i2 (int, list, None, optional): i2 index. If None, i2 is assumed to be 
                unique. Default is None. Use 'all' to add component to all i2 
                available.
        
        Return:
            None
        """        
        ###################
        # check i1 and i2 #
        ###################
        i1, i2 = self._check_i1_i2(i1=i1, i2=i2)

        #######################################
        # assert validity of parameter values #
        #######################################
        assert amp >= 0,          'amp cannot be negative'

        assert w >= 0,            'w cannot be negative'
        assert w_min >= 0,        'w_min cannot be negative'

        ##################
        # add parameters #
        ##################
        for _i2 in i2:
            tag = self.nametags['amp']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=amp, vary=amp_vary,  min=amp_min, max=amp_max, expr=amp_expr, brute_step=None)
            
            tag = self.nametags['c']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=c, vary=c_vary,  min=c_min, max=c_max, expr=c_expr, brute_step=None)
            
            tag = self.nametags['w']
            self.parent.add(f'{tag}_{i1}_{_i2}', value=w, vary=w_vary,  min=w_min, max=w_max, expr=w_expr, brute_step=None)
        return

    ##################
    # base modifiers #
    ##################
    def set_shift(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'cerf':
                # bounds
                self.parent[_name].min   += value[i2]
                self.parent[_name].max   += value[i2]
                # value
                self.parent[_name].value += value[i2]

        return 
    
    def set_offset(self, value):
        pass
        return 
    
    def set_factor(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'amperf':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])
        return 
    
    def set_calib(self, value):
        ##############################
        # check if value is a number #
        ##############################
        if self.parent is not None: 
            if isinstance(self.parent, br.Spectrum):
                value = [value, ]
            else:
                if isinstance(value, Iterable) == False:
                    value = [value]*len(self.parent)

        ##################
        # applying value #
        ##################
        for _name in self._get_all_tags():
            name, i1, i2 = _name_parser(_name)
            if name == 'cerf':
                # bounds
                if value[i2] > 0:
                    self.parent[_name].min   *= value[i2]
                    self.parent[_name].max   *= value[i2]
                else:
                    _min = self.parent[_name].min
                    self.parent[_name].min   = self.parent[_name].max*value[i2]
                    self.parent[_name].max   = _min*value[i2]
                # value
                self.parent[_name].value *= value[i2]
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr *= abs(value[i2])
            elif name == 'werf':
                # bounds
                self.parent[_name].min   *= abs(value[i2])
                self.parent[_name].max   *= abs(value[i2])
                # value
                self.parent[_name].value = abs(self.parent[_name].value*value[i2])
                # error
                if self.parent[_name].stderr is not None:
                    self.parent[_name].stderr   *= abs(value[i2])
        return 
# %%




    



    
    


