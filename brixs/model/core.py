#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Model object for fitting

###################
# Developer notes #
###################

1) Two types of default i1, i2 combination for functions:

i2=None, i1='all'
Description: 
    If i2 is None, i2 is assumed to be unique. Otherwise, raise an error
    If i1 is defined, check if i1 exists. Otherwise, raise an error
    Use self._check_i2_i1()

i2='all', i1=None
Description: 
    If i1 is None, i1 will be replaced by the max(i1)+1
    If i1 is defined, check if i1 exits. If yes, raise an error.
    If i2 is defined, check if i2 exists. Otherwise, raise an error.
    Use self._check_i2_i1()

"""

# %% ------------------------- Standard Imports --------------------------- %% #
import copy
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
from collections.abc import Iterable
from scipy.signal import find_peaks

lmfitok = False
try:
    import lmfit
    lmfitok = True
except:
    pass    

# BRIXS
import brixs as br

# TODO ================================
# s.model.peaks.find_peaks()
# s.model.peaks.pretty_print()
# s.model.pretty_print()
# __str__() and __repr__()

# copy

# save and load

# go to brixs and fix spectrum and spectra
# on spectra, fix __repr__ for empty spectra

# done
if lmfitok:
    # %% Peaks ========================================================================
    class ComponentTemplate(object):

        ##################
        # initialization #
        ##################
        def _initialize(self, parent, nametags, modifier_instructions, function, min_max_step):
            self.parent                = parent
            self.nametags              = nametags
            self.modifier_instructions = modifier_instructions
            self.function              = function
            self.min_max_step          = min_max_step
            self.tagnames              = dict((v, k) for k, v in nametags.items())

        #########
        # basic #
        #########
        def add(self, *args, **kwargs):
            raise NotImplementedError(f'add function is not defined for model component of type {type(self)}')

        def remove(self, i2=None, i1='all'):
            """delete all parameters associated with this component

            Args:
                i1 (int or str, optional): i1 index. Use `all` to get data from all
                    valid i1 indexes. Default is `all`.
                i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                    unique. Default is None.

            Returns:
                None
            """
            # get tags
            tags = self._get_tags(i2=i2, i1=i1)

            # delete tags
            for tag in tags:
                del self.parent[tag]
        
        def clear(self):
            """remove all components"""
            tags = self._get_all_tags()
            for tag in tags:
                del self.parent[tag]

        #################
        # magic methods #
        #################
        # def __str__(self):
        #     return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

        # def __repr__(self):
        #     return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

        def __setitem__(self, value, _value):
            ##########
            # string #
            ##########
            if isinstance(value, str):
                # (1) Object only one component: user just needs the attr name without i1 and i2
                # (2) Object have multiple components for one i2: user must indicate i1 only
                # (3) Object have multiple components for multiple i2: i1 and i2 must be indicated
                
                ##############
                # split name #
                ##############
                split = value.split('_')

                ##############
                # name alone #
                ##############
                if len(split) == 1:
                    name = split[0]

                    # check if name is a valid attribute
                    assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

                    # check if i2 is unique
                    assert self._unique_i2(), f'i2 is not unique and must be defined'
                    i2      = self._get_indexes_i2()[0]
                    indexes = self._get_indexes()

                    # check if i1 is unique
                    assert len(indexes[i2])==1, f'i1 is not unique and must be defined'
                    i1 = indexes[i2][0]

                    # set parameter value
                    self.parent[self.nametags[name] + '_' + str(i1) + '_' + str(i2)].value = _value
                ###############
                # name and i1 #
                ###############
                elif len(split) == 2:
                    name = split[0]
                    i1   = split[1]

                    # check if name is a valid attribute
                    assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

                    # check if i2 is unique
                    assert self._unique_i2(), f'i2 is not unique and must be defined'
                    i2      = self._get_indexes_i2()[0]
                    indexes = self._get_indexes()

                    # check i1 is valid
                    assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
                    i1 = int(i1)
                    assert i1 in indexes[i2], f'i1={i1} is not valid\nvalid i1: {indexes[i2]}'

                    # set parameter value
                    self.parent[self.nametags[name] + '_' + str(i1) + '_' + str(i2)].value = _value
                ####################
                # name, i1, and i2 #
                ####################
                elif len(split) == 3:
                    name = split[0]
                    i1   = split[1]
                    i2   = split[2]

                    # check if name is a valid attribute
                    assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

                    # check if i2 is valid
                    assert br.numanip.is_integer(i2), f'i2 must be a integer, not {type(i2)}'
                    i2 = int(i2)
                    indexes = self._get_indexes()
                    assert i2 in indexes, f'i2={i2} is not valid\nvalid i2: {list(indexes.keys())}'

                    # check i1 is valid
                    assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
                    i1 = int(i1)
                    assert i1 in indexes[i2], f'i1={i1} is not valid\nvalid i1: {indexes[i2]}'

                    # set parameter value
                    self.parent[self.nametags[name] + '_' + str(i1) + '_' + str(i2)].value = _value
                else:
                    raise KeyError(f'`{name}` is not a recognizable attr for component of type `{type(self)}`')
            ##########
            # error #
            ##########
            else:
                raise ValueError('key must be str or int')
            
        def __getitem__(self, value):
            #######
            # int #
            #######
            if isinstance(value, int):
                # (1) Object have multiple components for one i2: int will be considered i1
                # (2) Object have multiple components for multiple i2: int will be considered as i2
                
                if self._unique_i2(): 
                    i2 = self._get_indexes_i2()[0]  
                    i1 = value
                else:
                    i2 = value      
                    i1 = 'all'
                return self._get_parameters(i1=i1, i2=i2) 
            ##########
            # string #
            ##########
            elif isinstance(value, str):
                # (1) Object only one component: user just needs the attr name without i1 and i2
                # (2) Object have multiple components for one i2: user must indicate i1 only
                # (3) Object have multiple components for multiple i2: i1 and i2 must be indicated
                
                ##############
                # split name #
                ##############
                split = value.split('_')

                ##############
                # name alone #
                ##############
                if len(split) == 1:
                    name = split[0]

                    # check if name is a valid attribute
                    assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

                    # check if i2 is unique
                    assert self._unique_i2(), f'i2 is not unique and must be defined'
                    i2      = self._get_indexes_i2()[0]
                    indexes = self._get_indexes()

                    # check if i1 is unique
                    assert len(indexes[i2])==1, f'i1 is not unique and must be defined'
                    i1 = indexes[i2][0]

                    # get parameter
                    return self.parent[self.nametags[name] + '_' + str(i1) + '_' + str(i2)]
                ###############
                # name and i1 #
                ###############
                if len(split) == 2:
                    name = split[0]
                    i1   = int(split[1])

                    # check if name is a valid attribute
                    assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

                    # check if i2 is unique
                    assert self._unique_i2(), f'i2 is not unique and must be defined'
                    i2      = self._get_indexes_i2()[0]
                    indexes = self._get_indexes()

                    # check i1 is valid
                    assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
                    i1 = int(i1)
                    assert i1 in indexes[i2], f'i1={i1} is not valid\nvalid i1: {indexes[i2]}'

                    # get parameter
                    return self.parent[self.nametags[name] + '_' + str(i1) + '_' + str(i2)]
                ####################
                # name, i1, and i2 #
                ####################
                elif len(split) == 3:
                    name = split[0]
                    i1   = split[1]
                    i2   = split[2]

                    # check if name is a valid attribute
                    assert name in self.nametags, f'`{name}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

                    # check if i2 is valid
                    assert br.numanip.is_integer(i2), f'i2 must be a integer, not {type(i2)}'
                    i2 = int(i2)
                    indexes = self._get_indexes()
                    assert i2 in indexes, f'i2={i2} is not valid\nvalid i2: {list(indexes.keys())}'

                    # check i1 is valid
                    assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
                    i1 = int(i1)
                    assert i1 in indexes[i2], f'i1={i1} is not valid\nvalid i1: {indexes[i2]}'

                    # get parameter
                    return self.parent[self.nametags[name] + '_' + str(i1) + '_' + str(i2)]
                else:
                    raise KeyError(f'`{name}` is not a recognizable attr for component of type `{type(self)}`')
            ##########
            # error #
            ##########
            else:
                raise ValueError('key must be str or int')
            
        def __delitem__(self, value):
            """same as remove"""
            # (1) Object have multiple components for one i2: int will be considered i1
            # (2) Object have multiple components for multiple i2: int will be considered as i2
                
            if self._unique_i2(): 
                i2 = self._get_indexes_i2()[0]  
                i1 = value
            else:
                i2 = value      
                i1 = 'all'
            self.remove(i2=i2, i1=i1)

        ##################
        # base modifiers #
        ##################
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
            if value == 0:  # calib cannot be zero
                raise ValueError('cannot set calib = 0')
            elif value == 1:   # if calib is 1, do nothing
                return
            
            ##################
            # applying value #
            ##################
            for name in self._get_all_tags():
                raw = name.split('_')[0]
                if raw in self.base_modifiers['calib'] or raw in self.base_modifiers['calib_abs']:
                    if isinstance(value, Iterable):
                        f = lambda x: np.polyval(value, x)
                        new_value = np.float(f(self[name].value))
                    elif callable(value):
                        new_value = np.float(f(self[name].value))
                    else:
                        new_value = self[name].value*value
                    
                    self[name].value = new_value

                if raw in self.base_modifiers['calib_abs']:
                    self[name].value = abs(new_value)
                    
        def set_shift(self, value):
            """Shift peaks.

            Adds a value to the center parameter.

            Args:
                value (float or int): shift value (value will be added to center c).

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
            for name in self._get_all_tags():
                if name.split('_')[0] in self.base_modifiers['shift']:
                    self[name].value += value

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
            if value == 0:
                return

            ##################
            # applying value #
            ##################
            for name in self._get_all_tags():
                if name.split('_')[0] in self.base_modifiers['offset']:
                    self[name].value += value

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
            if value == 0:
                raise ValueError('cannot set factor = 0.')
            elif value == 1:
                return

            ##################
            # applying value #
            ##################
            for name in self._get_all_tags():
                if name.split('_')[0] in self.base_modifiers['factor']:
                    self[name].value += value

        ###########
        # support #
        ###########
        def _get_all_tags(self):
            """ok"""
            return [tag for tag in self.parent if tag.split('_')[0] in self.tagnames]

        def _unique_i2(self):
            if len(self._get_indexes_i2()) < 2: return True
            else: return False
        
        def _check_i2_i1(self, i2=None, i1='all'):
            """raises error if i2 and i1 are not valid

            Args:
                i1 (int or str, optional): i1 index. Use `all` to get data from all
                    valid i1 indexes. Default is `all`.
                i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                    unique. Default is None.

            Returns:
                valid i1 and i2
                i1 --> list
                i2 --> int
            """
            indexes = self._get_indexes()

            # check i2
            if i2 is None:
                if self._unique_i2():
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
            """raises error if i2 and i1 are not valid

            Args:
                i1 (int or str, optional): i1 index. If None, i1 will be defined as
                    the max(i1)+1. Default is None.
                i2 (int or None, optional): i2 index. Use `all` to get data from all
                    valid i1 indexes. Default is `all`.

            Returns:
                valid i1 and i2
                i1 --> int
                i2 --> list
            """
            indexes = self._get_indexes()

            # empty object
            if indexes == {}:
                i1 = 0
                i2 = 0
                return i1, [i2, ]

            # check i2
            if i2 == 'all':
                i2 = list(indexes.keys())
            else:
                assert br.numanip.is_integer(i2), f'i2 must be a integer, not {type(i2)}'
                assert i2 in indexes, f'i2={i2} not found'
                i2 = [i2, ]

            # check i1
            if i1 == None:
                i1 = max([max(indexes[_i2]) for _i2 in i2]) + 1
            else:
                assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
                check = {_i2: i1 in indexes[_i2] for _i2 in i2}
                assert False not in check, f'i1={i1} already exists for some i2\ni2 has i1={check}'
            return i1, i2

        def _get_indexes_i2(self):
            """return list of all secondary (spectrum number) indexes i2
            
            Returns:
                list
            """
            return list(np.unique([int(tag.split('_')[2]) for tag in self._get_all_tags()]))

        def _get_indexes(self):
            tags = self._get_all_tags()
            _i2 = np.unique([int(tag.split('_')[2]) for tag in tags])
            return {i2:list(np.unique([int(tag.split('_')[1]) for tag in tags if int(tag.split('_')[2]) == i2])) for i2 in _i2}
            
        def _get_tags(self, i2=None, i1='all'):
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

        def _get_parameters(self, i2=None, i1='all'):
            """returns parameters for a given i1 and i2

            Args:
                i1 (int or str, optional): i1 index. Use `all` to get data from all
                    valid i1 indexes. Default is `all`.
                i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                    unique. Default is None.

            Returns:
                a component type object
            """
            # get tags
            tags = self._get_tags(i2=i2, i1=i1)

            # obsolete
            # if len(tags) > 1:
            #     return {_i1: {self.tagnames[tag.split('_')[0]]: self.parent[tag] for tag in tags[_i1]} for _i1 in tags}
            # else:
            #     return {self.tagnames[tag.split('_')[0]]: self.parent[tag] for tag in tags[i1]}
            
            # get parameters
            parent = {}
            for _i1 in tags:
                parent.update({tag: self.parent[tag] for tag in tags[_i1]})

            # obsolete
            # if len(tags) > 1:
            #     parent = {}
            #     for _i1 in tags:
            #         parent.update({tag: self.parent[tag] for tag in tags[_i1]})
            # else:
            #     parent = {tag: self.parent[tag] for tag in tags[i1]}
            
            # final
            final = self.__class__(parent)
            return final 
        
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
                function f(x) as string
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

        def get_min_max_step_str(self, i2=None, i1='all'):
            """return string function to get min, max, and step values
            
            Args:
                i1 (int or str, optional): i1 index. Use `all` to get data from all
                    valid i1 indexes. Default is `all`.
                i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                    unique. Default is None.

            Returns:
                variables, vmin, vmax, step
            """
            # Initialization
            vmin      = 'np.min(['
            vmax      = 'np.max(['
            step      = 'np.max(['
            variables = ''

            # check i1, i2
            i1, i2 = self._check_i2_i1(i2=i2, i1=i1)

            # min_max_step
            if isinstance(i1, Iterable):
                for _i1 in i1:
                    _variables, _vmin, _vmax, _step = self.get_min_max_step_str(i1=_i1, i2=i2)
                    vmin      +=  _vmin     + ', '
                    vmax      +=  _vmax     + ', '
                    step      +=  _step     + ', '
                    variables +=  _variables
                return variables[:-2], vmin[:-2] + '])', vmax[:-2] + '])', step[:-2] + '])'
            else:
                vmin, vmax, step = self.min_max_step(i1=i1, i2=i2)
                for tag in self.tagnames:
                    variables += f'{tag}_{i1}_{i2} = self.parent["{tag}_{i1}_{i2}"].value, '
                return variables, vmin, vmax, step

        def get_min_max_step(self, i2=None, i1='all'):
            """return min, max, and step values
            
            Args:
                i1 (int or str, optional): i1 index. Use `all` to get data from all
                    valid i1 indexes. Default is `all`.
                i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                    unique. Default is None.

            Returns:
                vmin, vmax, step
            """
            variables, _vmin, _vmax, _step = self.get_min_max_step_str(i1=i1, i2=i2)
            vmin = f'lambda {variables}: {_vmin}'
            vmax = f'lambda {variables}: {_vmax}'
            step = f'lambda {variables}: {_step}'

            return eval(vmin), eval(vmax), eval(step)
        
        #############
        # modifiers #
        #############
        def reorder_by_attr(self, attr, increasing=True):
            """Change component i1 indexes based on parameter value
            
            args:
                attr (str): attr name
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
        def get_suitable_x(self, i2=None, i1='all'):
            """returns a array with x values for plotting.

            Args:
                i1 (int or str, optional): i1 index. Use `all` to get data from all
                    valid i1 indexes. Default is `all`.
                i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                    unique. Default is None.

            Returns:
                array
            """
            vmin, vmax, step = self.get_min_max_step(i2=i2, i1=i1)
            return np.arange(vmin(), vmax(), step())
        
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

            # calculate spectra
            ss = br.Spectra()
            for _i1 in i1:
                if x is None:
                    x = self.get_suitable_x(i2=i2, i1=_i1)
                ss.append(self.calculate_spectrum(i2=i2, i1=_i1, x=x))
            return ss

        def calculate_spectra(self, x=None):
            """Return a Spectra object with spectra for each i2 index 

            If component is absent for a certain i2, this will return an empty spectrum
            
            Args:
                x (list, optional): x values to which the curve will be calculated.
                    If None, a suitable x will be constructed.

            Returns:
                :py:class:`Spectra`
            """
            ss = br.Spectra()
            for i2 in self.parent._get_indexes_i2():
                if i2 not in self._get_indexes_i2():
                    ss.append(br.Spectrum())
                else:
                    if x is None:
                        x = self.get_suitable_x(i2=i2)
                    ss.append(self.calculate_spectrum(i2=i2, x=x))
            return ss

        ##############
        # get values #
        ##############
        def get_values_for_each_i2(self, attr, i1=None):
            """Returns a list with values of a component parameter for each i2.

            Args:
                attr (str): attribute name.
                i1   (int): i1 index (component number). If i1 is None, i1 is assumed
                    to be unique for all i2. Default is None.
            
            Returns:
                list of values
            """
            # check if name is a valid attribute
            assert attr in self.nametags, f'`{attr}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

            # check i1 is valid
            indexes = self._get_indexes()
            if i1 is None:
                check = {i2: len(indexes[i2]) for i2 in indexes}
                assert sum(list(check.values())) == len(indexes), f'i1 is not unique for all i2\ni1 must be defined\nindexes = {indexes}'
            else:
                check = {i2: i1 in indexes[i2] for i2 in indexes}
                assert False not in list(check.values()), f'some i2 indexes do not have i1 component\ni2 error = {check}\nindexes = {indexes}'

            # get values
            return [self._get_parameters(i2=i2, i1=i1)[attr].value for i2 in indexes]

        def get_errors_for_each_i2(self, attr, i1=None):
            """Returns a list with errors of a component parameter for each i2.

            Args:
                attr (str): attribute name.
                i1   (int): i1 index (component number). If i1 is None, i1 is assumed
                    to be unique for all i2. Default is None.
            
            Returns:
                list of values
            """
            # check if name is a valid attribute
            assert attr in self.nametags, f'`{attr}` is not a valid parameter for type `{type(self)}`\nvalid names: {list(self.nametags.keys())}'

            # check i1 is valid
            indexes = self._get_indexes()
            if i1 is None:
                check = {i2: len(indexes[i2]) for i2 in indexes}
                assert sum(list(check.values())) == len(indexes), f'i1 is not unique for all i2\ni1 must be defined\nindexes = {indexes}'
            else:
                check = {i2: i1 in indexes[i2] for i2 in indexes}
                assert False not in list(check.values()), f'some i2 indexes do not have i1 component\ni2 error = {check}\nindexes = {indexes}'

            # get values
            return [self._get_parameters(i2=i2, i1=i1)[attr].err for i2 in indexes]

        ##########################        
        # plot and visualization #
        ##########################
        def _check_modifiers_format(self, offset, shift, factor, calib):
            number_of_i2 = len(self._get_indexes_i2())

            # offset
            if isinstance(offset, Iterable):
                if len(offset) != number_of_i2:
                    raise ValueError(f'offset must be a number of a list with length compatible with the number of spectra.\nnumber of offsets: {len(offset)}\nnumber of spectra: {number_of_i2}')
            else:
                offset = [offset]*number_of_i2

            # shift
            if isinstance(shift, Iterable):
                if len(shift) != number_of_i2:
                    raise ValueError(f'shift must be a number of a list with length compatible with the number of spectra.\nnumber of shift: {len(shift)}\nnumber of spectra: {number_of_i2}')
            else:
                shift = [shift]*number_of_i2

            # calib
            if isinstance(calib, Iterable):
                if len(calib) != number_of_i2:
                    raise ValueError(f'calib must be a number of a list with length compatible with the number of spectra.\nnumber of calib: {len(calib)}\nnumber of spectra: {number_of_i2}')
            else:
                calib = [calib]*number_of_i2

            # factor
            if isinstance(factor, Iterable):
                if len(factor) != number_of_i2:
                    raise ValueError(f'factor must be a number of a list with length compatible with the number of spectra.\nnumber of factor: {len(factor)}\nnumber of spectra: {number_of_i2}')
            else:
                factor = [factor]*number_of_i2
            
            return offset, shift, factor, calib
        
        def plot(self, *args, **kwargs):
            raise NotImplementedError(f'plot function is not defined for model component of type {type(self)}')

    class Model(lmfit.Parameters):

        def __init__(self, parent=None):
            """
            talk about why parent has to be None (only from within the fit - minimizer - function)
            """
            super().__init__()
            self.parent = parent

            ##############
            # components #
            ##############
            self.components = {'peaks':      Peaks,
                            #    'asympeaks' : AsymmetricPeaks,
                            #    'areapeaks' : AreaPeaks,
                            #    'asym_area_peaks' : AsymmetricAreaPeaks
                            }

            # TODO
            # check repeated tags

            # start components
            for name, obj in self.components.items():
                self.__setattr__(name, obj(self))

        #################
        # magic methods #
        #################
        # TODO
        # def __str__(self):
        #     return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

        # TODO
        # def __repr__(self):
        #     return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

        def __setitem__(self, value, _value):
            super().__setitem__(value, _value)
        
        def __getitem__(self, value):
            #######
            # int #
            #######
            if isinstance(value, int):
                return self._get_parameters(i2=value) 
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
        def remove(self, i2):
            """remove all components with index i2"""
            # valid i2
            i2s = self._get_indexes_i2()
            assert i2 in i2s, f'i2={i2} not found\navailable i2: {i2s}'

            # remove
            for name in self:
                if name.split('_')[2] == i2:
                    del self[name]

        ###########
        # support #
        ###########
        def _unique_i2(self):
            if len(self._get_indexes_i2()) < 2: return True
            else: return False
        
        def _get_indexes_i2(self):
            """return list of all secondary (spectrum number) indexes i2
            
            Returns:
                list
            """
            return list(np.unique([int(name.split('_')[2]) for name in self]))

        def _get_indexes(self):
            """get indexes

            Returns:
                dict of dicts
                    {component name: {i2: [i1]}}
            """
            return {name: self.__getattribute__(name)._get_indexes() for name in self.components}
        
        def _get_parameters(self, i2):
            """returns parameters for a given i2

            Args:
                i2 (int or None, optional): i2 index

            Returns:
                a model object
            """
            # valid i2
            i2s = self._get_indexes_i2()
            assert i2 in i2s, f'i2={i2} not found\navailable i2: {i2s}'
            
            # get item
            final = Model(self.parent)
            for name in self:
                if name.split('_')[2] == i2:
                    final[name] = self[name]

            # # get tags
            # tags = self._get_tags(i2=i2)

            # # get parameters
            # if len(tags) > 1:
            #     parent = {}
            #     for _i1 in tags:
            #         parent.update({tag: self.parent[tag] for tag in tags[_i1]})
            # else:
            #     parent = {tag: self.parent[tag] for tag in tags[i1]}
            
            # # final
            # final = self.__class__(parent)
            return final 
        
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
                function f(x) as string
            """
            # check i2
            if i2 is None:
                if self._unique_i2():
                    i2 = self._get_indexes_i2()[0]
                else:
                    raise ValueError('i2 not defined')
            else:
                assert i2 in self._get_indexes_i2(), f'i2={i2} not found'
            
            # Initialization
            indexes   = self._get_indexes()
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

        def fit(self, ranges=None, method='least_squares', update=True):
            """fit spectrum with peaks model

            ############################# TODO
            Developers note: I think x and y must be monotonic, but I am not sure.
                If this is the case, include a check_monotonic in the future.

            Args:
                x, y (array or list of arrays): x and y coordinates to fit. If 
                    peaks has secondary (spectrum) indexes i2, then x and y must be
                    a list os x's and y's in the respective order.
                method (str, optional): Name of the fitting method to use. See methods
                    available on `lmfit.minimize()`_ documentation.
                update_peaks (bool, optional): if True, peaks will be updated with 
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
            if ranges is None:
                temp = self.parent
            else:
                temp = self.parent._extract(ranges=ranges)

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
                raise('there is not parent br.Spectrum/br.Spectra object associated with this fitting model')

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
            self
            out = lmfit.minimize(residual, method=method, params=self, kws={'xs':xs, 'ys':ys})        

            #####################
            # update parameters #
            #####################
            if update:
                for name in self:
                    self[name] = out.params[name]
            return out
        
        #################
        # save and load #
        #################
        # TODO

        ##################
        # base modifiers #
        ##################
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
                self.__getattribute__(component).set_calib(value)
                    
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
            if value == 0:
                return
            
            ##################
            # applying value #
            ##################
            for component in self.components:
                self.__getattribute__(component).set_shift(value)

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
                self.__getattribute__(component).set_offset(value)

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
                self.__getattribute__(component).set_factor(value)

        ##############
        # extractors #
        ##############
        # TODO
        def copy(self, *args, **kwargs):
            """Return a copy of the data.

            Usage:
                p2.copy(p1)                 # s2 is now a copy of s1
                p2.copy(peaks=p1)
                p2 = p1.copy()              

                peak = p1.copy(i1=2)        # peak is a copy of peak with index i1, assuming p1 does not have multiple spectra
                peak = p1.copy(2)

                peak = p1.copy(i1=2, i2=3)  # peak is a copy of peak with index i1 from spectrum i2
                peak = p1.copy(2, 3)

                peak = p1.copy(i2=3)        # peak ia a copy of all peaks from spectrum i2.

            Args:
                peaks (Peaks, optional): Peaks object is copied. See usage.
                i1 (int): primary index i1. Use i1='A' for additional parameters.
                i2 (int, optional): secondary (spectrum) index i2. Only used if 
                    peaks has secondary indexes.

            Returns:
                :py:attr:`Peaks`
            """
            ###################################
            # asserting validity of the input #
            ###################################
            error_message = 'Wrong input. Peaks object could not be copied. Please, use one ' +\
                            'of the examples below:\n' +\
                            '\n' +\
                            'p2.copy(peaks=p1)\n' +\
                            'p2.copy(p1)\n' +\
                            'p2 = p1.copy()' +\
                            '\n' +\
                            'peak = p1.copy(i1=2)\n' +\
                            'peak = p1.copy(2)\n' +\
                            '\n' +\
                            'peak = p1.copy(i1=2, i2=3)\n' +\
                            'peak = p1.copy(2, 3)\n' +\
                            '\n' +\
                            'peak = p1.copy(i2=3)'
            if kwargs != {} and args != ():
                raise AttributeError(error_message)
            if any([item not in ['peaks', 'i1', 'i2'] for item in kwargs.keys()]):
                raise AttributeError(error_message)
            if len(args) > 2 or len(kwargs) > 2:
                raise AttributeError(error_message)
            if ('i1' in kwargs or 'i2' in kwargs) and 'peaks' in kwargs:
                raise AttributeError(error_message)

            #################
            # sorting input #
            #################
            peaks = None
            i1    = None
            i2    = None
            # keyword arguments
            if 'peaks' in kwargs:
                peaks = kwargs['peaks']
            elif 'i1' in kwargs:
                i1 = kwargs['i1']
                if 'i2' in kwargs:
                    i2 = kwargs['i2']
            elif 'i2' in kwargs:
                i2 = kwargs['i2']

            # positional arguments
            if len(args) == 1:
                if isinstance(args[0], Peaks):
                    peaks = args[0]
                else:
                    i1 = args[0]
            elif len(args) == 2:
                i1 = args[0]
                i2 = args[1]

            # validating types
            if peaks is not None:
                assert isinstance(peaks, Peaks), 'Only type br.Peaks can be copied to type br.Peaks'
            if i1 is not None:
                assert br.numanip.is_integer(i1), f'i1 must be a integer, not {type(i1)}'
            if i2 is not None:
                assert br.numanip.is_integer(i2), f'i1 must be a integer, not {type(i2)}'

            ##################################
            # if peaks is passed, copy peaks #
            ##################################
            if peaks is not None:
                self.clear()
                for name in peaks:
                    self[name] = copy.deepcopy(peaks[name])

                # user defined attrs
                for attr in self._get_user_attrs():
                    self.__delattr__(attr)
                for attr in peaks._get_user_attrs():
                    self.__setattr__(attr, peaks.__dict__[attr])
                # additional attr
                peaks._additional = self._additional
                return
            
            ###############################
            # Otherwise, self is returned #
            ###############################
            else:
                peaks = self._get_peaks_by_index(i1=i1, i2=i2)
                # user defined attrs
                for attr in self._get_user_attrs():
                    peaks.__setattr__(attr, self.__dict__[attr])
                # additional attr
                if 'A' in peaks._get_indexes_i1():
                    peaks._additional = self._additional
                return copy.deepcopy(peaks)
            
        ######################
        # calculate spectrum #
        ######################
        def get_suitable_x(self, i2=None, i1='all'):
            """returns a array with x values for plotting.

            Args:
                i1 (int or str, optional): i1 index. Use `all` to get data from all
                    valid i1 indexes. Default is `all`.
                i2 (int or None, optional): i2 index. If None, i2 is assumed to be 
                    unique. Default is None.

            Returns:
                array
            """
            vmin = None
            vmax = None
            step = None

            for component in self.components:
                obj = self.__getattribute__(component)
                _vmin, _vmax, _step = obj.get_min_max_step(i2=i2, i1=i1)
                _vmin = _vmin()
                _vmax = _vmax()
                _step = _step()
                if vmin is None or vmin > _vmin:
                    vmin = _vmin
                if vmax is None or vmax < _vmax:
                    vmax = _vmax
                if step is None or step < _step:
                    step = _step

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

            If model is absent for a certain i2, this will return an empty spectrum
            
            Args:
                x (list, optional): x values to which the curve will be calculated.
                    If None, a suitable x will be constructed.

            Returns:
                :py:class:`Spectra`
            """
            ss = br.Spectra()
            for i2 in range(len(self.parent)):
                if i2 not in self._get_indexes_i2():
                    ss.append(br.Spectrum())
                else:
                    if x is None:
                        x = self.get_suitable_x(i2=i2)
                    ss.append(self.calculate_spectrum(i2=i2, x=x))
            return ss

    class Peaks(ComponentTemplate):
        def __init__(self, parent):
            super().__init__()

            ##############
            # EDIT THESE #
            ##############

            # parameter names
            # name = name used by the user
            # tag  = name used by lmfit
            # name: tag
            nametags = {'amp':'amp1', 'c':'c1', 'w':'w1', 'm':'m1'}  # make it immutable

            # base modifier instructions (use tags, not names)
            modifier_instructions = {'shift':     ['c1', ],    # additive factor to x axis
                                    'offset':    [],          # additive factor to y axis
                                    'calib':     ['c1'],      # multiplicative factor to x axis
                                    'calib_abs': ['w1'],      # same as calib, but avoid negative values (absolute)
                                    'factor':    ['amp1', ]}  # multiplicative factor to y axis
            
            def function(i1, i2):
                """function that returns f(x)"""
                return f"br.model.voigt_fwhm(x, amp1_{i1}_{i2}, c1_{i1}_{i2}, w1_{i1}_{i2}, m1_{i1}_{i2})"
            def min_max_step(i1, i2):
                """function that returns the min, max and step x value used when creating spectrum"""
                return f'c1_{i1}_{i2}-w1_{i1}_{i2}*10',  f'c1_{i1}_{i2}+w1_{i1}_{i2}*10', f'w1_{i1}_{i2}/20'
        
            # initialization 
            self._initialize(parent=parent, nametags=nametags, modifier_instructions=modifier_instructions, function=function, min_max_step=min_max_step)

        ####################
        # must be modified #
        ####################
        def add(self, amp=None, c=None, w=None, m=0, i1=None, i2='all'):  
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
            assert amp >= 0,          'amp cannot be negative'
            assert w >= 0,            'w cannot be negative'
            assert m >= 0 and m <= 1, 'm must be between 0 and 1'

            ##################
            # add parameters #
            ##################
            for _i2 in i2:
                tag = self.nametags['amp']
                self.parent.add(f'{tag}_{i1}_{_i2}', value=amp, vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
                
                tag = self.nametags['c']
                self.parent.add(f'{tag}_{i1}_{_i2}',   value=c,   vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
                
                tag = self.nametags['w']
                self.parent.add(f'{tag}_{i1}_{_i2}',   value=w,   vary=True,  min=0,       max=np.inf, expr=None, brute_step=None)
                
                tag = self.nametags['m']
                self.parent.add(f'{tag}_{i1}_{_i2}',   value=m,   vary=False, min=0,       max=1,      expr=None, brute_step=None)
            return

        ##########################        
        # plot and visualization #
        ########################## 
        def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, **kwargs):
            """Place a marker at the maximum of every peak position. Wrapper for `matplotlib.pyplot.errorbar()`_.

            If `Peaks` have secondary indexes i2, then of offset, shift, factor, and calib
            can be list, where each value will be used for the

            Args:
                ax (matplotlib.axes, optional): axes for plotting on.
                offset (number or list, optional): defines a vertical offset. Default is 0.
                shift (number or list, optional): horizontal shift value. Default is 0.
                factor (number or list, optional): multiplicative factor on the y axis.
                    Default is 1.
                calib (number or list, optional): multiplicative factor on the x axis.
                    Default is 1.
                **kwargs: kwargs are passed to `matplotlib.pyplot.errorbar()`_ that plots the data.

            Returns:
                `ErrorbarContainer`_

            .. matplotlib.pyplot.errorbar(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
            .. ErrorbarContainer: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.errorbar.html
            """
            ##########
            # figure #
            ##########
            if ax is None:
                ax = plt
                if br.settings.FIGURE_FORCE_NEW_WINDOW:
                    br.figure()

            ########################
            # check base modifiers #
            ########################
            offset, shift, factor, calib = self._check_modifiers_format(offset=offset, shift=shift, factor=factor, calib=calib)

            ##################
            # Initialization #
            ##################
            indexes = self._get_indexes()
            final = []

            ##################
            # default kwargs #
            ##################
            if 'lw' not in kwargs and 'linewidth' not in kwargs:
                kwargs['lw'] = 0
            if 'elinewidth' not in kwargs :
                kwargs['elinewidth'] = 2
            if 'marker' not in kwargs :
                kwargs['marker'] = 'o'
            if 'markersize' not in kwargs and 'ms' not in kwargs:
                kwargs['markersize'] = 5

            ########
            # plot #
            ########
            for _i, i2 in enumerate(indexes):
                # gather parameters
                c    = []
                amp  = []
                w    = []
                for i1 in indexes[i2]:
                    component = self._get_parameters(i2=i2, i1=i1)

                    c    += [component['c'].value]
                    amp  += [component['amp'].value]
                    w    += [component['w'].value/2]

                # apply modifiers
                c    = (np.array(c)*calib[_i])  + shift[_i]
                amp  = np.array(amp)*factor[_i] + offset[_i]
                w    = abs(np.array(w)*calib[_i])

                # plot
                final += [ax.errorbar(c, amp, xerr=w, **kwargs), ]

                # save modifiers
                final[0].offset = offset
                final[0].shift  = shift
                final[0].calib  = calib
                final[0].factor = factor
            return final

        #########
        # EXTRA #
        #########
        # TODO
        def find(self, prominence=5, width=4, moving_average_window=8):
            """
            
            ############################# TODO
            Developers note: I think x and y must be monotonic and uniform, but I am not sure.
                If this is the case, include a check_monotonic in the future.
                
            """

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
                y2 = br.arraymanip.moving_average(y, moving_average_window)
                x2 = br.arraymanip.moving_average(x, moving_average_window)
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
                    amp = d['prominences'][i] + max([y2[d['right_bases'][i]], y2[d['left_bases'][i]]])
                    c   = x2[peaks[i]]
                    w   = abs(d['widths'][i]*np.mean(np.diff(x)))

                    self.append(amp=amp, c=c, w=w)
            except IndexError:
                pass
else:
    class Model():
        pass

# %%
if __name__ == "__main__":
    import lmfit
    import numpy as np
    import matplotlib.pyplot as plt
    import brixs as br
    plt.ion()
    get_ipython().run_line_magic('load_ext', 'autoreload')
    get_ipython().run_line_magic('autoreload', '2')

    # TODO
    # for peak in self.model.peaks

    # unique i2
    s = br.Spectrum()
    s.model.peaks.add(amp=10, c=-4, w=4)
    s.model.peaks.add(amp=5, c=6, w=5)

    # check tags
    s.model.pretty_print()

    # test component
    s.model.peaks._unique_i2()
    s.model.peaks._get_indexes()

    s.model.peaks.check_feasibility()

    var, vmin, vmax, step = s.model.peaks.get_min_max_step_str()
    vmin, vmax, step = s.model.peaks.get_min_max_step()
    x = s.model.peaks.get_suitable_x()

    s.model.peaks._get_all_tags()
    s.model.peaks._get_tags()
    s.model.peaks._get_parameters()

    s.model.peaks[0]
    s.model.peaks[0]['c']

    s.model.peaks.get_values_for_each_i2('amp')

    # get model
    var, model = s.model.peaks.get_model_str(0, 1)
    model = s.model.peaks.get_model(0)
    y = model(x)

    s2 = s.model.peaks.calculate_spectrum()




    # multiple i2
    s = br.Spectrum()
    s.model.peaks.add(    1, 0, 0.1, 0, 0, 0)
    s.model.peaks.add(    2, 0, 0.1, 0, 1, 0)
    s.model.peaks.add(    1, 0, 0.1, 0, 0, 1)

    s.model
    s.model.check_feasibility()
    var, model = s.model.get_model_str(0)
    var, model = s.model.get_model_str(1)

    var, model = s.model.get_model_str()


    s.model.peaks.check_feasibility()
    var, model = s.model.peaks.get_model_str()
    model = s.model.peaks.get_model()
    
    s.model.fit()


# %%    

    # def _get_tags_by_index(self, i1='all', i2='all'):
    #     """return list of names with index i1 and i2.

    #     Args:
    #         i1 (int, optional): primary index i1
    #         i2 (int, optional): secondary (spectrum) index i2
        
    #     Returns:
    #         list
    #     """           

    #     # transform in list
    #     if i2 is None:
    #         i2 = self._get_indexes_i2()
    #         if i1 is None:
    #             i1 = [self._get_indexes_i1(i2=_i2) for _i2 in i2]
    #         else:
    #             # certify that all i2 have i1
    #             i1_exists_for_all_i2 = {f'i2={_i2}': i1 in self._get_indexes_i1(i2=_i2) for _i2 in i2}
    #             for _i2 in i1_exists_for_all_i2:
    #                 if i1_exists_for_all_i2[_i2] == False: 
    #                     raise ValueError(f'i1 does not exists for all i2\n{i1_exists_for_all_i2}')
    #             # if passes the test, then:
    #             i1 = [[i1, ] for _ in i2]
    #     else:
    #         if i1 is None:
    #             i1 = [self._get_indexes_i1(i2=i2), ]
    #         else:
    #             i1 = [[i1], ]
    #         i2 = [i2, ]
      
    #     # gather attrs
    #     final = []
    #     for i, _i2 in enumerate(i2):
    #         for _i1 in i1[i]:
    #             for name in self:
    #                 split = name.split('_')
    #                 raw  = split[0]
    #                 if raw != 'A':
    #                     __i1 = int(split[1])
    #                     __i2 = int(split[2])
    #                     if __i1 == _i1 and __i2 == _i2:
    #                         final.append(name)
    #     return final
    


    



    
    


