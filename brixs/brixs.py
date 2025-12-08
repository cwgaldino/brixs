#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core brixs module. Defines the main objects: 

    Spectrum, Spectra, PhotonEvents, Image, Dummy

Every BRIXS object has 5 types of attributes: 

    `Core`: 
        stores the main data of the object

    `Check`: 
        asserts data quality

    `Modifiers` (or `base modifiers`): 
        absolute values of 'shift', 'factor', 'offset', 'calib'

    `Labels`: 
        labels for datapoints

    `User`: 
        user-defined attributes

########################        
# NOTES FOR DEVELOPERS #  
########################       

Initializing objects and empty objects:
    Empty objects should have their core attrs set to None, e.g., Spectrum 
    has x and y data. Initially, s.x and s.y should be None. If the user sets
    s.x or s.y to an empty state, i.e., s.x = [], then s.x should be set to None
    (initial state).

Check attrs:
    
    1. `check` attrs cannot be user modifiable and shall only 
    be modified via 'check methods'

    2. check attrs exist so one does not waste time running checks all the time.
    Once a check is done, it's result is saved in a check attr.

    3. one should avoid using `check` methods inside functions/methods. It's
    preferable to raise an error and let the user deal with it or at least check 
    if check attrs are defined already before calling the check methods. 


Writing new methods:
    
    1. Methods shall avoid changing the `core` attrs inside 
    `self`. Instead, a copy of `self` shall be created, modified, and 
    returned to the user. If `core` attrs are modified directly on `self`, 
    this should be explicitly clear in the docstring.

    2. Methods that return a copy of `self` must (as much as possible) copy all 
    suitable attrs (`Modifiers`, `Labels`, and `User`). As for `Check` attrs, 
    methods must ensure that `check` attrs are still valid after operating
    on `self` or on the copy of `self`. If `check` attrs are not valid 
    anymore, reset them to None.  

    3. Methods should take into consideration what happens in case the object is 
    'empty' (object with no data).  

Copy methods:
    copy() must copy all attrs including all 4 attrs types
    (data, check, modifiers, labels) within _copy() and user attrs included in copy().

    
"""

# %% ------------------------- Standard Imports -------------------------- %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import warnings
import bisect
import copy

# %% -------------------------- Special Imports -------------------------- %% #
from collections.abc import Iterable, MutableMapping
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy.lib.stride_tricks import sliding_window_view

# %% ----------------------------- backpack ------------------------------ %% #
from .backpack import filemanip, arraymanip, figmanip, numanip, vectormanip, other

# %% ------------------------------ settings ----------------------------- %% #
from .config import settings
# %%

# %% =========================== support class =========================== %% #
class _Meta(type):
    """Metaclass to facilitate creation of read-only and non-removable attributes."""
    def __new__(self, class_name, bases, attrs):

        ###################
        # read only attrs #
        ###################
        def lazy_read_only(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    raise AttributeError('Attribute is "read only". Cannot set attribute.')
                def deleter(self):
                    raise AttributeError('Attribute is "read only". Cannot delete object.')
            return getter, setter, deleter#, 'read only attribute'

        #######################
        # non removable attrs #
        #######################
        def lazy_non_removable(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    return getattr(self, variable)
                def setter(self, value):
                    return setattr(self, variable, value)
                def deleter(self):
                    raise AttributeError('Attribute cannot be deleted.')
            return getter, setter, deleter#, 'non-removable attribute'

        #################
        # reserved word #
        #################
        def lazy_reserved(_attr):
            variable = '_' + _attr
            if not hasattr(self, variable):
                def getter(self):
                    raise AttributeError(f'Cannot get attribute. `{_attr}` is a reserved word.')
                def setter(self, value):
                    raise AttributeError(f'Cannot set attribute. `{_attr}` is a reserved word.')
                def deleter(self):
                    raise AttributeError('Attribute cannot be deleted.')
            return getter, setter, deleter

        ################
        # create attrs #
        ################
        new_attrs = {}

        # if '_reserved' in attrs:
        #     print('reserved')
        #     for attr in attrs['_reserved']:
        #         print(attr)
        #         _property = property(*lazy_reserved(attr))
        #         new_attrs[attr] = _property
        if '_non_removable' in attrs:
            for attr in attrs['_non_removable']:
                _property = property(*lazy_non_removable(attr))
                new_attrs[attr] = _property
        if '_read_only' in attrs:
            for attr in attrs['_read_only']:
                _property = property(*lazy_read_only(attr))
                new_attrs[attr] = _property
        for name, value in attrs.items():
            if name not in ('_non_removable', '_read_only', '_reserved'):
                new_attrs[name] = value

        return type(class_name, bases, new_attrs)

class _BrixsObject(object): 
    """Parent class defining common attrs and methods"""   
    #########
    # attrs #
    #########       
    def get_core_attrs(self): 
        """return a list of core attrs"""
        if isinstance(self, Spectrum):
            name = 'Spectrum'
        elif isinstance(self, Spectra):
            name = 'Spectra'
        elif isinstance(self, Image):
            name = 'Image'
        elif isinstance(self, PhotonEvents):
            name = 'PhotonEvents'
        elif isinstance(self, Dummy):
            name = 'Dummy'
        else:
            raise ValueError('What?')
        return settings._reserved_words[name]['pseudovars']
        
    def get_attrs(self):
        """return attrs that are user defined.""" 
        # return [key for key in self.__dict__.keys() if key.startswith('_') == False and key not in settings._reserved_words]
        if isinstance(self, Spectrum):
            name = 'Spectrum'
        elif isinstance(self, Spectra):
            name = 'Spectra'
        elif isinstance(self, Image):
            name = 'Image'
        elif isinstance(self, PhotonEvents):
            name = 'PhotonEvents'
        elif isinstance(self, Dummy):
            name = 'Dummy'
        else:
            raise ValueError('What?')

        return [_ for _ in dir(self)
                if _ not in self._get_methods() 
                and (_.startswith('__')==False and _.endswith('__')==False)
                and _ not in settings._reserved_words[name]['pseudovars']
                and _ not in settings._reserved_words[name]['vars']]

    def get_attrs_dict(self):
        """return a dict of user defined attrs and their values"""
        return {key: self.__getattribute__(key) for key in self.get_attrs()}

    def _get_methods(self):
        """return a list of methods available including hidden ones"""
        if isinstance(self, Spectrum):
            name = 'Spectrum'
        elif isinstance(self, Spectra):
            name = 'Spectra'
        elif isinstance(self, Image):
            name = 'Image'
        elif isinstance(self, PhotonEvents):
            name = 'PhotonEvents'
        elif isinstance(self, Dummy):
            name = 'Dummy'
        else:
            raise ValueError('What?')
        
        methodList = []
        for method_name in dir(self):
            if method_name in settings._reserved_words[name]['pseudovars'] + settings._reserved_words[name]['vars']:
                pass
            elif method_name in settings._reserved_words[name]['methods']:
                    methodList.append(str(method_name))
            else:
                try:
                    if callable(getattr(self, method_name)):
                        methodList.append(str(method_name))
                except Exception:
                    methodList.append(str(method_name))
        return methodList
    
    def get_methods(self):
        """return a list of methods available"""
        return [_ for _ in self._get_methods() 
                if _.startswith('_') == False]
                #and _.endswith('_') == False]
    
    ############################
    # attrs that return a copy #
    ############################
    def remove_attrs(self):
        """Delete all user defined attrs.
        
        Returns:
            brixs object with removed attrs
        """

        s2 = self.copy()
        for attr in self.get_attrs():
            s2.__delattr__(attr)

        return s2
    
    def copy_attrs_from(self, s):
        """Copy user defined attributes from another brixs object.

        Usage:
            s3 = s2.copy_attrs_from(s)

        Args:
            s (brixs object): Either a Spectrum, Spectra, Image, or PhotonEvents
                to copy user defined attributes from.
        
        Returns:
            brixs object with copied attrs
        """
        # check type
        if isinstance(s, Spectrum) or isinstance(s, Spectra) or isinstance(s, Image) or isinstance(s, PhotonEvents) or isinstance(s, Dummy):
            pass
        else:
            raise TypeError(f'type {type(s)} not valid\nCan only copy user attrs from type br.Spectrum, br.Spectra, br.Image, or br.PhotonEvents')

        # transfer attrs
        s2 = self.copy()
        for attr in s.get_attrs():
            value = copy.deepcopy(s.__getattribute__(attr))
            s2.__setattr__(attr, value)

        return s2
    
    def _copy_attrs_from(self, s):
        """Copy user defined attributes from another brixs object.

        Usage:
            s1._copy_attrs_from(s)

        Args:
            s (brixs object): Either a Spectrum, Spectra, Image, or PhotonEvents
                to copy user defined attributes from.
        
        Returns:
            None
        """
        # check type
        if isinstance(s, Spectrum) or isinstance(s, Spectra) or isinstance(s, Image) or isinstance(s, PhotonEvents) or isinstance(s, Dummy):
            pass
        else:
            raise TypeError(f'type {type(s)} not valid\nCan only copy user attrs from type br.Spectrum, br.Spectra, br.Image, or br.PhotonEvents')

        # transfer attrs
        for attr in s.get_attrs():
            value = copy.deepcopy(s.__getattribute__(attr))
            self.__setattr__(attr, value)

        return 
# %%

# %% ====================== common support functions ===================== %% #
def get_functions():
    """return a list of methods available including hidden ones"""        
    methodList = []

    for _module in (filemanip, arraymanip, figmanip, numanip, vectormanip, other):
        methodList.append(f'')
        methodList.append(f'===== {_module.__name__} =====')
        for _method_name in dir(_module):
            try:
                if callable(getattr(_module, _method_name)) and _method_name.startswith('_') == False:
                    methodList.append(str(_method_name))
            except Exception:
                methodList.append(str(_method_name))
    return methodList

def _attr2str(attrs_dict, verbose):
    """returns a list with strings ("name: value") for each attr and attr value

        Note:
            This function is similar to json.dump(), but formatting is more 
            appropriate and gives us more flexibility.

        Note: 
            nested dictionaries will be indented with 4 spaces.

        Warning:
            numpy arrays will return a `list`.
        
        Args:
            attrs_dict (dict): a dictionary with attr names and values, e.g., 
                {'attr1': 10, 'attr2': [1, 2, 3]}. Values can be numbers, lists,
                arrays, dict, None, datetime.
            verbose (bool): if True, message is printed when attr cannot be 
                converted to string.
        
        Returns:
            list
    """
    final = []

    #######################
    # collect attr values #
    #######################
    for name in attrs_dict:
        try:
            ##############
            # type: None #
            ##############
            if attrs_dict[name] is None:
                final.append(f'{name}: None')
            #############
            # type: str #
            #############
            elif isinstance(attrs_dict[name], str):
                temp2 = str(attrs_dict[name]).replace('\n','\\n')
                final.append(f'{name}: \"{temp2}\"')
            ##############
            # type: dict #
            ##############
            elif isinstance(attrs_dict[name], dict) or isinstance(attrs_dict[name], MutableMapping):
                final.append(f'{name}: ' + '{')
                for line in _attr2str(attrs_dict[name], verbose):
                    final.append('    ' + line)
                final.append('}')
            ##################
            # type: Spectrum #
            ##################
            elif isinstance(attrs_dict[name], Spectrum):
                x = list(attrs_dict[name].x)
                y = list(attrs_dict[name].y)
                final.append(f'{name}: ' + f'Spectrum(x={x}, y={y})')
            #################
            # type: Spectra #
            #################
            elif isinstance(attrs_dict[name], Spectra):
                _internal = ''
                for _s in attrs_dict[name]:
                    x = list(_s.x)
                    y = list(_s.y)
                    _internal += f'Spectrum(x={x}, y={y}), '
                final.append(f'{name}: ' + f'Spectra(data=[' + _internal + '])')
            ######################
            # type: PhotonEvents #
            ######################
            elif isinstance(attrs_dict[name], PhotonEvents):
                x = list(attrs_dict[name].x)
                y = list(attrs_dict[name].y)
                final.append(f'{name}: ' + f'PhotonEvents(x={x}, y={y})')
            ######################
            # type: PhotonEvents #
            ######################
            elif isinstance(attrs_dict[name], Image):
                _data = [list(_) for _ in attrs_dict[name].data]
                final.append(f'{name}: ' + f'Image(data={_data})')
            ########################
            # type: list and tuple #
            ########################
            elif isinstance(attrs_dict[name], Iterable):
                final.append(f'{name}: {arraymanip.array2list(attrs_dict[name])}')
            ##################
            # type: filepath #
            ##################
            elif isinstance(attrs_dict[name], Path):
                final.append(f'{name}: {attrs_dict[name]}')           
            ################
            # type: number #
            ################
            elif numanip.is_number(attrs_dict[name]):
                tosave = str(attrs_dict[name])
                if tosave[-1] == '\n':
                    tosave = tosave[:-1]
                final.append(f'{name}: {tosave}')
            else:
                temp2 = str(attrs_dict[name]).replace('\n','\\n')
                final.append(f'{name}: \"{temp2}\"')
        except:
            if verbose:
                type_ = str(type(attrs_dict[name]))
                print('\nWarning: Cannot save attr of type: ' + type_ + '\nattr name: '+ name + '\nTo turn off this warning, set verbose to False.')
    return final

def _str2attr(header, indentation=0, verbose=False):
    """returns a dictionary with (names: values) given a list of strings "name: value"

    Args:
        header (list): list of strings "name: value". This should be the output
            of the _attr2str() function.
        indentation (int, optional): number of spaces before characters start
        verbose (bool): if True, message is printed when name/value cannot be 
                converted from string to a python object.
    """
    output = {}
    for i, line in enumerate(header):
        line = line[indentation:]
        if line[0] == ' ' or line[0] == '}':
            pass
        elif ':' not in line:
            pass
        else:
            # get name and value
            _split = line.split(':')
            name   = _split[0].strip()
            value  = ':'.join(_split[1:]).strip()

            # parse dictionaries
            if value == '{':
                for j, line2 in enumerate(header[i+1:]):
                    line2 = line2[indentation:]
                    if line2[0] == '}':
                        break
                indentation2 = indentation + 4
                value = _str2attr(header[i+1:i+j+1], indentation=indentation2, verbose=False)
            else:
                try:
                    value = eval(str(value).strip())
                except NameError:
                    if value.startswith('Spectrum'):
                        value = eval(str(value).replace('nan', 'np.nan').strip())
                except:# Exception as e:
                    value = str(value).strip()
                
            try:
                output[name] = value

            except Exception as e:
                if verbose:
                    print(f'cannot read attr ({name}: {value})\nAttribute not set\n{e}\n')
    return output
# %%

# %% ====================== modified figure function ===================== %% #
def figure(*args, **kwargs):
    """same as br.figmanip.figure(), but the following br.settings affect figure:

        br.settings.FIGURE_POSITION
        br.settings.FIGURE_FORCE_ON_TOP
        br.settings.FIGURE_DPI
        br.settings.FIGURE_SIZE
        br.settings.FIGURE_GRID
    """
    fig = figmanip.figure(*args, **kwargs)

    ############
    # position #
    ############
    if settings.FIGURE_POSITION is not None:
        figmanip.set_window_position(settings.FIGURE_POSITION)

    ###############
    # set onclick #
    ###############
    if settings.FIGURE_ONCLICK is not None:
        cid1 = fig.canvas.mpl_connect('button_press_event', settings.FIGURE_ONCLICK)

    # #############
    # # force top #
    # #############
    # if br.settings.FIGURE_FORCE_ON_TOP:
    #     print('aa')
    #     bring2top()

    ##############
    # figure DPI #
    ##############
    if 'dpi' not in kwargs:
        if settings.FIGURE_DPI is not None:
            fig.set_dpi(settings.FIGURE_DPI)

    ########################
    # figure size and grid #
    ########################
    if 'figsize' not in kwargs:
        if settings.FIGURE_SIZE is not None:
            figmanip.set_window_size(settings.FIGURE_SIZE)

        # grid
        if settings.FIGURE_GRID:
            rows    = settings.FIGURE_GRID[0]
            columns = settings.FIGURE_GRID[1]

            if (rows > 1 and columns > 0) or (columns > 1 and rows > 0):
                count  = settings._figure_count - 1
                row    = int((count/columns)%rows)
                column = count%columns

                if settings.FIGURE_SIZE is None:
                    height, width = figmanip.get_window_size()
                else:
                    height = settings.FIGURE_SIZE[0]
                    width  = settings.FIGURE_SIZE[1]

                position = (settings.FIGURE_POSITION[0]+row*(height+settings.FIGURE_GRID_OFFSET[0]), settings.FIGURE_POSITION[1]+column*(width+settings.FIGURE_GRID_OFFSET[1]))
                figmanip.set_window_position(position)

                settings._figure_count += 1
        else:
            # set_window_position()
            settings._figure_count = 0
    return fig

def subplots(*args, **kwargs):
    """same as br.figmanip.subplots(), but the following br.settings affect figure:

        br.settings.FIGURE_POSITION
        br.settings.FIGURE_FORCE_ON_TOP
        br.settings.FIGURE_DPI
        br.settings.FIGURE_SIZE
        br.settings.FIGURE_GRID

    """
    fig, axes = figmanip.subplots(*args, **kwargs)

    ############
    # position #
    ############
    if settings.FIGURE_POSITION is not None:
        figmanip.set_window_position(settings.FIGURE_POSITION)

    ###############
    # set onclick #
    ###############
    if settings.FIGURE_ONCLICK is not None:
        cid1 = fig.canvas.mpl_connect('button_press_event', settings.FIGURE_ONCLICK)

    # #############
    # # force top #
    # #############
    # if br.settings.FIGURE_FORCE_ON_TOP:
    #     print('aa')
    #     bring2top()

    ##############
    # figure DPI #
    ##############
    if 'dpi' not in kwargs:
        if settings.FIGURE_DPI is not None:
            fig.set_dpi(settings.FIGURE_DPI)

    ########################
    # figure size and grid #
    ########################
    if 'figsize' not in kwargs:
        if settings.FIGURE_SIZE is not None:
            figmanip.set_window_size(settings.FIGURE_SIZE)

        # grid
        if settings.FIGURE_GRID:
            rows    = settings.FIGURE_GRID[0]
            columns = settings.FIGURE_GRID[1]

            if (rows > 1 and columns > 0) or (columns > 1 and rows > 0):
                count  = settings._figure_count - 1
                row    = int((count/columns)%rows)
                column = count%columns

                if settings.FIGURE_SIZE is None:
                    height, width = figmanip.get_window_size()
                else:
                    height = settings.FIGURE_SIZE[0]
                    width  = settings.FIGURE_SIZE[1]

                position = (settings.FIGURE_POSITION[0]+row*(height+settings.FIGURE_GRID_OFFSET[0]), settings.FIGURE_POSITION[1]+column*(width+settings.FIGURE_GRID_OFFSET[1]))
                figmanip.set_window_position(position)

                settings._figure_count += 1
        else:
            # set_window_position()
            settings._figure_count = 0
    return fig, axes
# %%

# %% ============================= Spectrum ============================== %% #
class Spectrum(_BrixsObject, metaclass=_Meta):
    """Returns a ``Spectrum`` object. 

    Args:
        x (list or array, optional): x values (1D list/array).
        y (list or array, optional): y values (1D list/array).
        filepath (str or Path, optional): filepath.
        **kwargs: kwargs are passed to :py:func:`Spectrum.load` function

    Usage:        
            >>> s = br.Spectrum()
            >>> s = br.Spectrum(x, y)
            >>> s = br.Spectrum(x=x, y=y)
            >>> s = br.Spectrum(y=y)
            >>> s = br.Spectrum().load(filepath=<filepath>)
            >>> s = br.Spectrum().load(filepath=<filepath>, delimiter=',')
            >>>
            >>> print(s.get_core_attrs()) # print list of core attrs
            >>> print(s.get_attrs())      # print list of attrs
            >>> print(s.get_methods())    # print list of methods available

    Attributes:
        Every BRIXS object has 5 types of attributes: 
        `Core` , `Check`, `Modifiers`, `Labels`, `User`.

        *1. Core*
            x, y (array): 1D arrays.
        
        *2. Check*
            has_nan (bool): True if x or y has non-numeric values (NaN).
                Can only be modified by s.check_nan().
            step (number): None or a number if the step between two data points 
                is the same through out the x vector. Can only be modified by 
                s.check_step() method.
            monotonicity (string): None if data is not monotonic or 'increasing'
                 or 'decreasing'. Can only be modified by s.check_monotonicity()
                 method.

        *3. Modifiers*
            calib, factor, shift, offset (number): absolute values of 
            modifications made to the data. 
        
        *4. Labels*
            None

        *5. User*
            anything that the user defined on the fly.        
   """
    _read_only     = ['step', 'monotonicity', 'has_nan']
    _non_removable = []

    def __init__(self, x=None, y=None, filepath=None, **kwargs):
        """Initialize the object instance"""
        ###########################
        # Initializing attributes #
        ###########################
        # core
        self._x = None
        self._y = None

        # check
        self._step         = None
        self._monotonicity = None
        self._has_nan      = None

        # modifiers
        self._calib  = 1
        self._factor = 1
        self._shift  = 0
        self._offset = 0

        # labels
        pass

        # extra
        for extra in settings._init['Spectrum']:
            self.__setattr__(extra, settings._init['Spectrum'][extra](self))

        ################
        # loading data #
        ################
        if y is not None:
            self.y = y
        if x is not None:
            self.x = x
        elif filepath is not None and y is None:
            self.load(filepath, **kwargs)
        return

    ###################
    # core attributes #
    ###################
    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, value):
        # check None
        if value is None:
            if self.y is not None:
                self._x = np.arange(0, len(self.y))
            else:
                self._x = None
            self._step         = None
            self._monotonicity = None
            self._has_nan      = None
            return
        
        ###################################
        # asserting validity of the input #
        ###################################
        # check type
        if not isinstance(value, Iterable):
            raise TypeError(f'The x-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a x-array which is not an Iterable.\nThe type of the variable you passed is: {type(value)}\nAccepted types are: list, array, ...')
    
        # check if iterable is made of numbers
        try:
            _ = sum(value)
        except:
            raise TypeError(f'The x-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a x-array which is Iterable, but is NOT ENTIRELY made of numbers because we fail to sum all the elements of the array you are passing.')
        
        # check length
        if self.y is not None:
            assert len(value) == len(self.y), f'Length of x array (len={len(value)}) you are trying to set is not compatible with current length of the y array (len={len(self.y)}).'
        else:
            if len(value) == 0:
                self._x = None
                self._step         = None
                self._monotonicity = None
                self._has_nan      = None
                return
            else:
                self._y = np.arange(0, len(value))
              
        #################
        # set attribute #
        #################
        self._x = np.array(value, dtype='float')

        ##########################
        # reset check attributes #
        ##########################
        self._step         = None
        self._monotonicity = None
        self._has_nan      = None
    @x.deleter
    def x(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def y(self):
        return self._y
    @y.setter
    def y(self, value):
        ###################################
        # asserting validity of the input #
        ###################################
        # check None
        if value is None:
            if self.x is not None:
                self._y = np.arange(0, len(self.x))
            else:
                self._y = None
            return
            
        # check type
        if not isinstance(value, Iterable):
            raise TypeError(f'The y-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a y-array which is not an Iterable.\nThe type of the variable you passed is: {type(value)}\nAccepted types are: list, array, ...')
    
        # check if iterable is made of numbers
        try:
            _ = sum(value)
        except:
            raise TypeError(f'The y-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a y-array which is Iterable, but is NOT ENTIRELY made of numbers because we fail to sum all the elements of the array you are passing.')
        
        # check length
        if self.x is not None:
            assert len(value) == len(self.x), f'Length of y-array (len={len(value)}) you are trying to set is not compatible with current length of the x-array (len={len(self.x)}).'
        else:
            if len(value) == 0:
                self._y = None
                return
            else:
                self._x = np.arange(0, len(value))
                self._step         = 1
                self._monotonicity = 'increasing'

        #################
        # set attribute #
        #################
        self._y = np.array(value, dtype='float')

        ##########################
        # reset check attributes #
        ##########################
        # self._step         = None
        # self._monotonicity = None
        self._has_nan      = None
    @y.deleter
    def y(self):
        raise AttributeError('Cannot delete object.')

    ###################################
    # computed (read-only) attributes #
    ###################################
    @property
    def data(self):
        if self.x is not None:
            return np.vstack((self.x, self.y)).transpose()
        else:
            return
    @data.setter
    def data(self, value):
        raise AttributeError('Attribute is "read only".')
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    #########################
    # write-only attributes #
    #########################
    pass

    #######################
    # modifier attributes #
    #######################
    @property
    def calib(self):
        return self._calib
    @calib.setter
    def calib(self, value):
        _s = self.set_calib(1/self._calib).set_calib(value)
        self.x      = _s.x
        self._calib = value
        self._shift = _s.shift
    @calib.deleter
    def calib(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift(self):
        return self._shift
    @shift.setter
    def shift(self, value):
        _s = self.set_shift(-self._shift).set_shift(value)
        self.x      = _s.x
        self._shift = value
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        return self._offset
    @offset.setter
    def offset(self, value):
        _s = self.set_offset(-self._offset).set_offset(value)
        self.y      = _s.y
        self._offset = value
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        return self._factor
    @factor.setter
    def factor(self, value):
        _s = self.set_factor(1/self._factor).set_factor(value)
        self.y      = _s.y
        self._factor= value
        self._offset = _s.offset
    @factor.deleter
    def factor(self):
        raise AttributeError('Cannot delete object.')

    #################
    # magic methods #
    #################
    def __setattr__(self, name, value):
        if name in settings._forbidden_words['Spectrum']:
            raise AttributeError(f'`{name}` is a reserved word and cannot be set as an attribute')
        super().__setattr__(name, value)

    def __len__(self):
        if self.x is None:
            return 0
        else:
            return len(self.x)

    def __add__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot operate on spectra. x axis is different.\nMaybe try interpolating the x axis.')
            s    = self.copy()
            s._y = self.y + object.y
            return s
        elif isinstance(object, (np.floating, float, int)):
            return self.set_offset(object)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')

    def __sub__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot operate on spectra. x axis is different.\nMaybe try interpolating the x axis.')
            s    = self.copy()
            s._y = self.y - object.y
            return s
        elif isinstance(object, (np.floating, float, int)):
            return self.set_offset(-object)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')

    def __mul__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot operate on spectra. x axis is different.\nMaybe try interpolating the x axis.')
            s    = self.copy()
            s._y = self.y*object.y
            return s
        elif isinstance(object, (np.floating, float, int)):
            return self.set_factor(object)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')

    def __div__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot divide spectra. x axis is different.\nMaybe try interpolating the x axis.')
            if 0 in object.y:
                raise ZeroDivisionError(f'y axis contain zeros. Cannot divide by zero.')
            else:
                s    = self.copy()
                s._y = self.y/object.y
                return s
        elif isinstance(object, (np.floating, float, int)):
            if object == 0:
                raise ZeroDivisionError(f'Cannot divide by zero.')
            else:
                return self.set_factor(1/object)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')

    def __truediv__(self, object):
        return self.__div__(object)

    # this is what makes possible max(s) and min(s)
    def __iter__(self):
        for y in self.y:
            yield y
        # for x, y in zip(self.x, self.y):
        #     yield x, y

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.x[item], self.y[item]
        elif isinstance(item, slice):
            x = self.x[item]
            y = self.y[item]

            s = Spectrum(x=x, y=y)

            # transfer attrs
            for attr in self.get_attrs():
                value = copy.deepcopy(self.__dict__[attr])
                s.__setattr__(attr, value)
            return s
        else:
            raise TypeError('Index must be int or a slice, not {}'.format(type(item).__name__))

    def __delitem__(self, item):
        if isinstance(item, int) or isinstance(item, slice):
            self._x = np.delete(self.x, item)
            self._y = np.delete(self.y, item)

            ##########################
            # reset check attributes #
            ##########################
            self._step         = None
            self._monotonicity = None
            # self._has_nan      = None
        else:
            raise TypeError('Index must be int or a slice, not {}'.format(type(item).__name__))
        
    # This is what makes possible s1 > s2 ?
    # Also makes possible max(ss)
    def __gt__(self, object):
        return max(self.y) > max(object.y)
    
    # This is what makes possible s1 < s2 ?
    # Also makes possible min(ss)
    def __lt__(self, object):
        return min(self.y) < min(object.y)
    
    #########
    # attrs #
    #########
    pass

    ###########
    # support #
    ###########
    def _check_limits(self, limits):
        """returns limits in the right format.

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            None or limits in the following format:
                ((xi_1, xf_1), (xi_2, xf_2), ...)               
        """
        #####################
        # if limits is None #
        #####################
        if limits is None:
            return None
        
        ##################################
        # assert that limits is Iterable #
        ##################################
        assert isinstance(limits, Iterable), f'`limits` must be an Iterable, not {type(limits)}'
        
        ################
        # empty object #
        ################
        if self.x is None:
            raise ValueError('cannot operate on empty spectrum')

        ################################
        # get min and max range values #
        ################################
        vmin = min(self.x)
        vmax = max(self.x)

        ################
        # empty limits #
        ################
        if len(limits) == 0:
            return ((vmin, vmax),)

        ##############
        # fix format #
        ##############
        # one pair
        if len(limits) == 1: # ((xi, xf), )
            assert isinstance(limits[0], Iterable), f'wrong format for limits={limits}'
        # two pairs 
        elif len(limits) == 2: # (xi, xf), or ((xi1, xf1), (xi2, xf2))
            if isinstance(limits[0], Iterable) == False:
                if isinstance(limits[1], Iterable) == False:
                    limits = [limits, ]
                else:
                    raise ValueError(f'wrong format for limits={limits}')
        final = []
        # three or more pairs 
        for lim in limits:
            assert isinstance(lim, Iterable), f'wrong format for limits={limits}'
            temp = [lim[0], lim[1]]
            if temp[0] == None: temp[0] = vmin
            if temp[1] == None: temp[1] = vmax
            final.append((temp[0], temp[1]))   

        return final

    ################
    # core methods #
    ################
    def append(self, x, y):
        """return spectrum with added new data points
        
        Args:
            x, y (number): number to be appended
            
        Returns:
            Spectrum
        """
        assert isinstance(x, Iterable) == False, 'x must be a number'
        try:
            x + 3.2, 
        except:
            raise ValueError('x must be a number')
        assert isinstance(y, Iterable) == False, 'y must be a number'
        try:
            y + 3.2, 
        except:
            raise ValueError('y must be a number')

        # append
        s = self.copy()
        s._x = np.append(self._x, x)
        s._y = np.append(self._y, y)

        ##########################
        # reset check attributes #
        ##########################
        # self._step         = None
        # self._monotonicity = None
        # self._has_nan      = None
        return s
    
    ########
    # copy #
    ########
    def _copy(self, limits=None):
        """Same as self.copy(), but attrs are NOT copied."""
        ################
        # copy x and y #
        ################
        x = copy.deepcopy(self.x)
        y = copy.deepcopy(self.y)

        #####################
        # if limits is None #
        #####################
        limits = self._check_limits(limits=limits)
        if limits is None:
            s = Spectrum(x=x, y=y)
            s._step         = self.step
            s._monotonicity = self.monotonicity
            s._shift        = self.shift
            s._calib        = self.calib
            s._offset       = self.offset
            s._factor       = self.factor
            return s
        
        #####################
        # if empty spectrum #
        #####################
        if self.x is None:
            # print(f'Warning: No datapoints within limits={limits}')
            s = Spectrum(x=x, y=y)
            s._step         = self.step
            s._monotonicity = self.monotonicity
            s._shift        = self.shift
            s._calib        = self.calib
            s._offset       = self.offset
            s._factor       = self.factor
            return s
        
        ########################################
        # check if extract is really necessary #
        ########################################
        if len(limits) == 1:
            if min(limits[0]) <= min(self.x) and max(limits[0]) >= max(self.x):
                s = Spectrum(x=x, y=y)
                s._step         = self.step
                s._monotonicity = self.monotonicity
                s._shift        = self.shift
                s._calib        = self.calib
                s._offset       = self.offset
                s._factor       = self.factor
                return s
        
        ###########
        # extract #
        ###########
        # try:
        x, y = arraymanip.extract(x, y, limits)
        s = Spectrum(x=x, y=y)
        s._shift  = self.shift
        s._calib  = self.calib
        s._offset = self.offset
        s._factor = self.factor
        return s
        # except RuntimeError:
        #     print(f'Warning: No datapoints within limits={limits}')
        #     return Spectrum()

    def copy(self, limits=None):
        """Return a copy of the data within a range limits.

        Usage:
            >>> # full copy
            >>> s2 = s1.copy()  # s2 is now a copy of s2
            >>>
            >>> # s2 will have only data between x=0 and 10 and between x=90 and 100
            >>> s2 = s1.copy(limits=((0, 10), (90, 100)))  

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Note:
            spectrum is not required to be monotonic. 

        Returns:
            :py:class:`Spectrum`
        """
        s = self._copy(limits=limits)
        s._copy_attrs_from(self)

        # extra
        for extra in settings._copy['Spectrum']:
            if hasattr(s, extra):
                s.__setattr__(extra, self.__getattribute__(extra).copy(s))

        return s

    #################
    # save and load #
    #################     
    def save(self, filepath, number_of_decimal_places=None, only_data=False, check_overwrite=False, verbose=False, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Warning:
            Attrs are saved as comments if only_data is False. Saving attrs to file
            is not always reliable because requires converting variables to string. 
            Attrs that are of type: strings, numbers, arrays, lists, strings, 
             None, dicts should work fine. More exotic types must be tested. 

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            number_of_decimal_places (int or list, optional): if not None, this argument
                defines the number of decimal places to save the data (it does
                not affect attrs, only x and y arrays). If list, it must have two 
                values corresponding to the number_of_decimal_places to be used for 
                to x and y separately. The 'fmt' argument overwrites
                number_of_decimal_places. Default is None. If None, the number
                of decimal places will be such that data will be saved with the 
                best precision necessary.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user wants to overwrite file.
            verbose (bool, optional): turn verbose on and off. Default is `False`.
            *kwargs: kwargs are passed to ``np.savetxt()`` that saves the data.

        If not specified, the following parameters are passed to `numpy.savetxt()`_:

        Args:
            fmt (string, or list, optional): A single format (like ``%10.5f``), or a
                sequence of formats. See numpy's documentation for more information.
                If not specified, best fmt is calculated based on the
                number of decimal places of the data. Specifying fmt makes the code
                runs a little faster (not much tough, but it might make a difference
                if saving a lot of files).
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (", ").
            newline (str, optional): String or character separating lines.
                Default is ``\n``.
            header (bool, optional): String that will be written at the beginning of the file.
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        See Also:
            :py:func:`Spectrum.load`

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        ############
        # filepath #
        ############
        filepath = Path(filepath)

        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'
        
        ###################
        # check overwrite #
        ###################
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if other.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        ##########
        # kwargs #
        ##########
        if 'fmt' in kwargs: # pick best format
            pass
        elif number_of_decimal_places is not None:
            if isinstance(number_of_decimal_places, Iterable):
                kwargs['fmt'] = (f'%.{number_of_decimal_places[0]}f', f'%.{number_of_decimal_places[1]}f')
            else:
                kwargs['fmt'] = f'%.{number_of_decimal_places}f'
        else:
            if self.has_nan is None:
                self.check_nan()
            if self.has_nan:
                temp = self.remove_nan()
            else:
                temp = self.copy()
            if self.x is not None:
                number_of_decimal_places_x = max([numanip.n_decimal_places(x) for x in temp.x])
                number_of_decimal_places_y = max([numanip.n_decimal_places(y) for y in temp.y])
                kwargs['fmt'] = (f'%.{number_of_decimal_places_x}f', f'%.{number_of_decimal_places_y}f')
        kwargs.setdefault('delimiter', ', ')
        kwargs.setdefault('newline', '\n')
        kwargs.setdefault('comments', '# ')

        #####################
        # header and footer #
        #####################
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
        else:
            attrs_dict = {_:self.__getattribute__(_) for _ in settings._reserved_words['Spectrum']['vars'] if _ not in ['_x', '_y']}
            attrs_dict.update(self.get_attrs_dict())
            header = '\n'.join(_attr2str(attrs_dict, verbose)) + '\n'

            if 'header' not in kwargs:
                kwargs['header'] = header
            else:
                if kwargs['header'] == '':
                    kwargs['header'] = header
                elif kwargs['header'][-1] != '\n':
                    kwargs['header'] += '\n'
        
        ########
        # save #
        ########
        if self.x is None:
            np.savetxt(Path(filepath), [], **kwargs)
        else:
            np.savetxt(Path(filepath), self.data, **kwargs)
        return

    def load(self, filepath, only_data=False, verbose=False, **kwargs):
        """Return spectrum with loaded data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Warning:
            This a very simple loading function that works well with column text
            files (xy files). If file has more columns than two columns, the first two columns
            will be loaded by default. Use `usecols` to select columns with x and y data.

        Note:
            If file was saved by br.Spectrum.save(), then the metadata (comments) can be 
            recovered. If not, one might get better results by setting `only_data = True`.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first 
                decompressed. Last used filepath is saved to an attr s.filepath.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is loaded. Default is False.
            verbose (book, optional): Default is False. If True, it will print
                an warning when attributes cannot be loaded from the file.
            **kwargs: kwargs are passed to ``numpy.genfromtxt()`` that loads the data.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``#``. Attributes picked up
                from the header will be loaded too.
            usecols (tuple, optional): Default is (0, 1).

        Returns:
            Spe trum

        See Also:
            :py:func:`Spectrum.save`

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        ############
        # filepath #
        ############
        filepath = Path(filepath)
        assert filepath.is_file(), f'filepath must point to a file, {filepath}'
        assert filepath.exists(),  f'filepath does not exist, {filepath}'

        ##########
        # kwargs #
        ##########
        kwargs.setdefault('delimiter', ', ')
        kwargs.setdefault('comments', '# ')
        kwargs.setdefault('usecols', (0, 1))

        ########
        # read #
        ########
        data = np.genfromtxt(Path(filepath), **kwargs)

        s = Spectrum()
        if len(data) != 0:
            ##############
            # check data #
            ##############
            x = data[:, 0]
            y = data[:, 1]
            assert len(x) == len(y), f'Length of x array (len={len(x)}) is not compatible with y array (len={len(y)}).'
            
            ##########
            # assign #
            ##########
            s._x = x
            s._y = y
        else:
            s._x = None
            s._y = None

        ##########################
        # reset check attributes #
        ##########################
        # self._step         = None
        # self._monotonicity = None
        # self._has_nan      = None

        ###################
        # reset modifiers #
        ###################
        # self._calib  = 1
        # self._factor = 1
        # self._shift  = 0
        # self._offset = 0

        ###############
        # read header #
        ###############
        if only_data is False:
            # get header
            header = filemanip.load_comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            
            # remove comment flag (#)
            comment_flag_len = len(kwargs['comments'])
            if header:
                for i, line in enumerate(header):
                    header[i] = line[comment_flag_len:]

                # attrs dict
                attrs_dict = _str2attr(header[:-1], verbose=verbose)

                # set attrs
                for attr in attrs_dict:
                    s.__setattr__(attr, attrs_dict[attr])
        return s

    #########
    # check #
    #########
    def check_nan(self):
        """Check if x or y have non-numeric (NaN) values

        Result (True or False) is stored in s.has_nan attribute
        
        Returns:
            None
        """
        ################
        # empty object #
        ################
        if self.x is None:
            self._has_nan = False
            return
        
        #############
        # check nan #
        #############
        if np.isnan(self.x).any():
            self._has_nan = True
        elif np.isnan(self.y).any():
            self._has_nan = True
        else:
            self._has_nan = False
        return None
    
    def find_nan(self):
        """Return a list with indexes where non-numeric (NaN) values were found
        
        Returns:
            list
        """
        if self.has_nan is None:
            self.check_nan()

        if self.has_nan:
            # xi = list(np.argwhere(np.isnan([0, 1, np.nan, 2, 3]))[:, 0])
            # yi = list(np.argwhere(np.isnan([0, 1, np.nan, np.nan, 3]))[:, 0])
            xi = list(np.argwhere(np.isnan(self.x))[:, 0])
            yi = list(np.argwhere(np.isnan(self.y))[:, 0])
            return np.unique(np.concatenate((xi, yi), 0))
        return []
    
    def remove_nan(self):
        """remove data points (x, y) that contain non-numeric (NaN) values
        
        Returns:
            Spectrum
        """
        if self.has_nan is None:
            self.check_nan()

        if self.has_nan == False:
            return self.copy()
        else:
            index2remove = list(np.argwhere(np.isnan(self.x))) + list(np.argwhere(np.isnan(self.y)))
            x = np.delete(self.x, index2remove)
            y = np.delete(self.y, index2remove)

            # copy
            s = self.copy()
            s._x = x
            s._y = y
            # reset checks
            s._step = None
            s._has_nan = False
            return s
    
    def check_step(self, max_error=0.1):
        """Checks vector uniformity of the x-coordinates.

            If the step between two data points is the same through out the
            x vector, it sets :py:attr:`s.step` with the value of the step size.

            Args:
                max_error (number, optional): percentage (of the x step) value 
                of the max error. Default is 0.1 %
                
            Note:
                Step uniformity is verified by the following equation:

                    (max(steps) - min(steps))/np.mean(steps) * 100 < max_error
            
            Note:
                If s.step is well defined, then s.monotonicity is also well defined. 
                If spectrum has only one datapoint, step will be set to None and 
                monotonicity set to 'increasing'.

            Returns:
                None

            Raises:
                ValueError: If x-coordinates are not uniform.
        
            See Also:
                :py:func:`Spectrum.check_monotonicity`
        """
        ########################
        # check empty spectrum #
        ########################
        if self.x is None:
            raise ValueError('cannot check step for empty spectrum')
        
        #####################################
        # check spectrum with one datapoint #
        #####################################
        if len(self.x) == 1:
            self._step = None
            self._monotonicity = 'increasing'
            return
        
        # if data is not monotonic, than it is not uniform
        if self.monotonicity is None:
            try:
                self.check_monotonicity()
            except ValueError:
                raise ValueError(f"Step in the x-coordinate seems not to be uniform. In fact, it is not even monotonic. Use Spectrum.fix_monotonicity().")

        # check step uniformity
        d = np.diff(self.x)
        if abs((max(d) - min(d))*100/np.mean(np.diff(self.x))) > max_error:
            self._step = None
            raise ValueError(f"Step in the x-coordinate seems not to be uniform.")

        # set step
        self._step = np.mean(d)

        return

    def check_monotonicity(self):
        """Sets monotonicity attribute to 'increasing' or 'decreasing'.

        Raises:
            ValueError if data is not monotonic.

        Note:
            If spectrum has only one datapoint, step will be set to None and 
            monotonicity set to 'increasing'.

        Returns:
            None

        See Also:
                :py:func:`Spectrum.check_step`, :py:func:`Spectrum.fix_monotonicity`
        """
        ########################
        # check empty spectrum #
        ########################
        if self.x is None:
            raise ValueError('cannot check monotonicity for empty spectrum')
        
        if np.all(np.diff(self.x) > 0) == True:
            self._monotonicity = 'increasing'
        elif np.all(np.diff(self.x) < 0) == True:
            self._monotonicity = 'decreasing'
        else:
            self._monotonicity = None
            raise ValueError('x array is not monotonic. Use Spectrum.fix_monotonicity()')

    def fix_monotonicity(self, mode='increasing'):
        """Rearrange (x, y) such as x array is monotonically increasing or decreasing.

            Args:
                mode (str, optional): increasing or decreasing.

            Note:
                duplicated datapoints are averaged.

            Returns:
                Spectrum
            
            See Also:
                :py:func:`Spectrum.check_monotonicity`
        """
        ########################
        # check empty spectrum #
        ########################
        if self.x is None:
            raise ValueError('cannot operate on empty spectrum')
        
        # check mode
        increasing = ['inc', 'i', 'up', 'increasing', 'increasingly']
        decreasing = ['dec', 'd', 'down', 'decreasing', 'decreasingly']
        if mode not in increasing and mode not in decreasing:
            raise ValueError('mode should be "decreasing" or "increasing".')
        
        if mode in increasing: mode = 'increasing'
        else: mode = 'decreasing'

        # copy
        s = self.copy()

        # turn array into monotonic
        if s.monotonicity is None:
            try:
                s.check_monotonicity()
            except ValueError:
                s._step = None
                unqa, ID, counts = np.unique(s.x, return_inverse=True, return_counts=True)
                data        = np.column_stack((unqa, np.bincount(ID, s.y)/counts))
                
                s._x = data[:, 0]
                s._y = data[:, 1]
                # check
                s.check_monotonicity()
            
        # make it decreasing or increasing
        if s.monotonicity != mode:
            s._step = None
            s._x = self.x[::-1]
            s._y = self.y[::-1]
        s.check_monotonicity()

        return s

    ##################
    # BASE modifiers #
    ##################
    def set_shift(self, value):
        """Add value to x-coordinates.

        Args:
            value (float or int): shift value.

        Returns:
            :py:class:`Spectrum`
        """
        s         = self.copy()
        s._x     += value
        s._shift += value

        # extra
        for extra in settings._shift['Spectrum']:
            if hasattr(s, extra):
                s.__getattribute__(extra).set_shift(value)
            
        return s

    def set_offset(self, value):
        """Add value to y-coordinates.

        Args:
            value (value): offset value.

        Returns:
            :py:class:`Spectrum`
        """
        s          = self.copy()
        s._y      += value
        s._offset += value

        # extra
        for extra in settings._offset['Spectrum']:            
            if hasattr(s, extra):
                s.__getattribute__(extra).set_offset(value)
            
        return s

    def set_factor(self, value):
        """Multiply y-coordinates by value.

        Args:
            value (number): multiplicative factor

        Raises:
            AttributeError: if value is 0.

        Returns:
            :py:class:`Spectrum`
        """
        if value == 0:
            raise AttributeError('cannot set factor = 0.')
        s          = self.copy()
        s._y      *= value
        s._factor *= value
        s._offset *= value

        # extra
        for extra in settings._factor['Spectrum']:            
            if hasattr(s, extra):
                s.__getattribute__(extra).set_factor(value)
            
        return s

    def set_calib(self, value):
        """Multiply x-coordinates by value.

        Args:
            value (number): calibration value

        Raises:
            AttributeError: if value is 0.
                
        Returns:
            :py:class:`Spectrum`
        """
        if value == 0:
            raise AttributeError('cannot set calib = 0.')
        s         = self.copy()
        s._x     *= value
        s._calib *= value
        s._shift *= value

        # extra
        for extra in settings._calib['Spectrum']:            
            if hasattr(s, extra):
                s.__getattribute__(extra).set_calib(value)
            
        return s

    #############
    # modifiers #
    #############
    def set_x_via_polyval(self, p):
        """Set x to np.polyval(p, x).

        Args:
            p (array): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term.

        Returns:
            :py:class:`Spectrum`
        """
        s    = self.copy()
        f    = lambda x: np.polyval(p, x)
        s._x = np.array([f(x) for x in self.x])
        s._step         = None
        s._monotonicity = None
        return s

    def set_x_via_function(self, f):
        """Set x to f(x).

        Args:
            f (function): function where argument is x-coordinate elements

        Returns:
            :py:class:`Spectrum`
        """
        s    = self.copy()
        s._x = np.array([f(x) for x in s.x])
        s._step         = None
        s._monotonicity = None
        return s
    
    def set_y_via_function(self, f):
        """Set y to f(x, y).

        Args:
            f (function): function where argument is x- and y-coordinate elements

        Returns:
            :py:class:`Spectrum`
        """
        s    = self.copy()
        s._y = np.array([f(x, y) for x, y in zip(s.x, s.y)])
        s._step         = None
        s._monotonicity = None
        return s

    def set_roll(self, value):
        """Roll array elements of the y-coordinate [same as s.set_y_roll()]

        Note:
            Elements that roll beyond the last position are re-introduced at the
            first.

        Args:
            value (float or int): roll value.

        Returns:
            :py:class:`Spectrum`
        """    
        return self.set_y_roll(value)

    def set_y_roll(self, value):
        """Roll array elements of the y-coordinate

        Note:
            Elements that roll beyond the last position are re-introduced at the
            first.

        Args:
            value (float or int): roll value.

        Returns:
            :py:class:`Spectrum`
        """    
        ###################################
        # asserting validity of the input #
        ###################################
        # x axis must be uniform if mode is roll
        if self.step is None:
            try:
                self.check_step()
            except ValueError:
                raise ValueError(f'Cannot roll data, because x-coordinates are not uniform. Use s.interp() to make data uniform, or shift data using s.set_shift()')
        # value must be an integer
        if numanip.is_integer(value) == False:
            raise ValueError('roll value must be an integer.')
        
        ########
        # roll #
        ########
        s         = self.copy()
        s._y      = np.roll(self.y, int(value))
        s._shift += value*s.step
        return s

    def set_x_roll(self, value):
        """Roll array elements of the x-coordinates

        Note:
            Elements that roll beyond the last position are re-introduced at the
             first.

        Args:
            value (float or int): roll value.

        Returns:
            :py:class:`Spectrum`
        """    
        ###################################
        # asserting validity of the input #
        ###################################
        # x axis must be uniform if mode is roll
        if self.step is None:
            try:
                self.check_step()
            except ValueError:
                raise ValueError(f'Cannot roll data, because x-coordinates are not uniform. Use s.interp() to make data uniform, or shift data using s.set_shift()')
        # value must be an integer
        if numanip.is_integer(value) == False:
            raise ValueError('roll value must be an integer.')
        
        ########
        # roll #
        ########
        s         = self.copy()
        s._x      = np.roll(self.x, int(value))
        s._shift -= value*s.step
        return s

    def flip_x(self):
        """Flip sign of x axis.

        Returns:
            :py:class:`Spectrum`
        """
        return self.set_calib(value=-1)

    def flip_y(self):
        """Flip sign of y axis.

        Returns:
            :py:class:`Spectrum`
        """
        return self.set_factor(value=-1)

    def floor(self, limits=None):
        """Sets zero value for y-coordinates.

        Usage:
            >>> # brings the avg of all data points to zero.
            >>> s.floor()
            >>>
            >>> # Brings the avg between x=0 and 10 and between x=90 and 100 to zero
            >>> s.floor(limits=((0, 10), (90, 100)))  

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectrum`
        """
        if limits is None: limits = (None, None)
        temp  = self._copy(limits=limits)
        assert temp.y is not None, 'no data points within limits = {limits}'
        value = -np.mean(temp.y)
        return self.set_offset(value)
    
    def normalize(self, value=1, limits=None):
        """Set a factor such as the average y between limits is equal to value.

        Args:
            value (number): value. Default is 1.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectrum`
        """
        if limits is None: limits = (None, None)
        temp = self._copy(limits=limits)
        assert temp.y is not None, 'no data points within limits = {limits}'
        return self.set_factor(value/np.mean(temp.y))

    ############
    # advanced #
    ############
    def interp(self, start=None, stop=None, num=None, step=None, x=None):
        """Interpolate data.

        Args:
            start (number, optional): The starting value of the sequence. If None,
                the minium x value will be used.
            stop (number, optional): The end value of the sequence. If None,
                the maximum x value will be used.
            num (int, optional): Number of samples to generate.
            step (number, optional): Spacing between values. This overwrites ``num``.
            x (list or array, optional): The x-coordinates at which to
                evaluate the interpolated values. This overwrites all other arguments.
            
        None:
            If none arguments is given, this function will adjust x to have regular spacing.

        Returns:
            :py:class:`Spectrum`
        """
        ################################
        # assert spectrum is not empty #
        ################################
        assert len(self) > 0, 'Cannot operate on an empty Spectrum'

        ##############################################
        # np.interp requires that array is monotonic #
        ##############################################
        self.check_monotonicity()
        if self.monotonicity != 'increasing':
            raise ValueError('monotonicity of data must be strictly increasing. Use Spectrum.fix_monotonicity().')

        ##############
        # new x axis #
        ##############
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
                temp = self._copy(limits=(start, stop))
                x = np.linspace(start, stop, num=len(temp))
        else:
            x = np.array(x)
            
        #################
        # interpolation #
        #################
        y = np.interp(x, self.x, self.y)

        ########
        # copy #
        ########
        s        = self.copy()
        s._x    = x
        s._y    = y
        s._step = None
        s._monotonicity = None
        s._has_nan = None

        return s

    def derivative(self, order=1):
        """Returns the derivative of y-coordinates as a function of x-coordinates.

        Args:
            order (number, optional): derivative order. Default is 1.

        Returns:
            :py:class:`Spectrum`
        """
        # assert spectrum is not empty
        assert len(self) > 0, 'Cannot calculate derivative of an empty Spectrum'

        # copy x and y and calculate derivative
        x, y = arraymanip.derivative(self.x, self.y, order=order)
        
        # copy
        s       = self.copy()
        s._x    = x
        s._y    = y
        s._step = None
        s._monotonicity = None
        s._has_nan = None

        return s
    
    def moving_average(self, n):
        """Returns the moving average of the data.

        Args:
            n (int): number of points to average.

        Returns:
            :py:class:`Spectrum` of length given by (len(x)-n+1)
        """
        # copy x and y and smooth
        x = arraymanip.moving_average(self.x, n=n)
        y = arraymanip.moving_average(self.y, n=n)

        # copy
        s       = self.copy()
        s._x    = x
        s._y    = y
        s._step = None
        s._monotonicity = None
        s._has_nan = None

        return s

    def smooth(self, n, force_divisible=False):
        """Returns Spectrum which is the average over every n elements

        Args:
            n (int): number of points to average.
            force_divisible (bool, optional): if True, the length of the data
                must be divisible by n. Default is False.

        Warning:        
            If force_divisible=False, the last data point might be averaged by 
            less points then n.

        Returns:
            :py:class:`Spectrum` of length given by (len(x)/n).
        """
        # check if data is divisible by n
        if force_divisible:
            assert len(self) % n == 0, f"The length of data ({len(self)}) is not evenly divisible by n={n}\nPlease, set force_divisible=False or pick one of the following numbers: {np.sort(list(numanip.factors(len(self))))}"

        if len(self) % n == 0:
            x = np.mean(self.x.reshape(-1, n), axis=1)
            y = np.mean(self.y.reshape(-1, n), axis=1)
        else:
            ids = np.arange(len(self.x))//n
            x = np.bincount(ids, self.x)/np.bincount(ids)

            ids = np.arange(len(self.y))//n
            y = np.bincount(ids, self.y)/np.bincount(ids)

        # copy
        s       = self.copy()
        s._x    = x
        s._y    = y
        s._step = None
        s._monotonicity = None
        s._has_nan = None
        
        return s

    def smooth_with_multiple_n(self, intervals):
        """Returns Spectrum which is the average over n elements where n can be vary over different data ranges.

        Args:
            intervals (list of 3-element lists): a list of the type 
            [[xstart_0, xstop_0, n_0], ..., [xstart_m, xstop_m, n_m]]

        Returns:
            smoothed out spectrum    

        See Also:
            :py:func:`Spectrum.smooth`
        """
        ss_smooth_intervals  = Spectra()
        for xstart, xstop, n in intervals:
            _s = self._copy(limits=(xstart, xstop))
            ss_smooth_intervals.append(_s.smooth(n))
        s_smooth = ss_smooth_intervals.concatenate()
        s_smooth._copy_attrs_from(self)
        return s_smooth

    def crop(self, start=None, stop=None):
        """Crop edges of the dataset.

        Args:
            start (number, optional): start x value. If None, the minimum value of
                x will be used.
            stop (number, optional): final x value. If None, the maximum value of
                x will be used.
        
        Note:
            spectrum is not required to be monotonic. 

        Returns:
            :py:class:`Spectrum`
        """
        # assert spectrum is not empty
        assert len(self) > 0, 'Cannot crop an empty Spectrum'
  
        # sort start and stop
        if start is None and stop is None:
            return self.copy()
        
        if start is None:
            start = min(self.x)
        if stop is None:
            stop = max(self.x)
        
        if start <= min(self.x) and stop >= max(self.x):
            return self.copy()
        
        # copy
        return self.copy(limits=(start, stop))

    def switch_xy(self):
        """Switch x and y axis.
        
        Returns:
            :py:class:`Spectrum`
        """
        s       = self.copy()
        s._x    = copy.deepcopy(self.y)
        s._y    = copy.deepcopy(self.x)
        s._step = None
        s._monotonicity = None
        
        return s

    def remove(self, limits):
        """Remove datapoints within a range limits.

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectrum`
        """
        limits = self._check_limits(limits=limits)
        x, y   = arraymanip.extract(self.x, self.y, limits, invert=True)
        
        s       = self.copy()
        s._x    = x
        s._y    = y
        s._step = None
        s._monotonicity = None
        s._has_nan = None
        
        return s

    ########################
    # calculation and info #
    ########################
    def calculate_area(self, limits=None):
        """Returns the calculated area under the curve. Wrapper for `numpy.trapz()`_.

        Usage:
            >>> s.calculate_area()                # returns the area for the whole dataset
            >>> s.calculate_area(limits=(0, 10))  # returns the area between x=0 and 10
            
        Warning:
            Calculate_area for broken limits, i.e., `s.calculate_area(limits=((0, 10), (90, 100)))`
             works, but will 
            typically not return a desirable value, because the ``trapz()``
            algorithm will consider the area between 10 and 90 as a rectangle.
            The correct approach in this case would be to calculated the area
            between 0 and 10 and between 90 and 100 separately.

            >>> area = s.calculate_area(0, 10) + s.calculate_area(90, 100)

            this behavior might change in the future.

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            number

        .. _numpy.trapz(): https://numpy.org/doc/stable/reference/generated/numpy.trapz.html
        """
        ################################
        # assert spectrum is not empty #
        ################################
        assert len(self) > 0, 'Cannot operate on an empty Spectrum'

        ##############################################
        # np.interp requires that array is monotonic #
        ##############################################
        self.check_monotonicity()
        if self.monotonicity != 'increasing':
            raise ValueError('monotonicity of data must be strictly increasing. Use Spectrum.fix_monotonicity().')

        s = self._copy(limits=limits)
        return np.trapz(y=s.y, x=s.x)

    def calculate_x_sum(self, limits=None):
        """Returns sum of x elements within a range.
        
        Usage:
            >>> s.calculate_x_sum()  # returns the x sum for the whole dataset
            >>> s.calculate_x_sum((0, 10), (90, 100))  # returns the x sum from data between x=0 and 10 and between x=90 and 100

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            number
        """
        s = self._copy(limits=limits)
        return sum(s.x)

    def calculate_y_sum(self, limits=None):
        """Returns sum of y elements within a range.
        
        Usage:
            >>> s.calculate_y_sum()  # returns the y sum for the whole dataset
            >>> s.calculate_y_sum((0, 10), (90, 100))  # returns the y sum from data between x=0 and 10 and between x=90 and 100

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            number
        """
        s = self._copy(limits=limits)
        return sum(s.y)

    def calculate_x_average(self, limits=None):
        """Returns the average value of x elements within a range.
        
        Usage:
            >>> s.calculate_x_average()  # returns the y average for the whole dataset
            >>> s.calculate_x_average((0, 10), (90, 100))  # returns the y average from data between x=0 and 10 and between x=90 and 100

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            number
        """
        s = self._copy(limits=limits)
        return np.mean(s.x)
    
    def calculate_y_average(self, limits=None):
        """Returns the average value of y elements within a range.
        
        Usage:
            >>> s.calculate_y_average()  # returns the y average for the whole dataset
            >>> s.calculate_y_average((0, 10), (90, 100))  # returns the y average from data between x=0 and 10 and between x=90 and 100

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            number
        """
        s = self._copy(limits=limits)
        return np.mean(s.y)

    def calculate_calib(self, values, mode='peaks', deg=1, limits=None, **kwargs):
        """return a calibration factor (`spectral feature vs value`).

        Note:
            The calibration factor is the shift value (x-coord) as a function of
             the `values` array. Assuming each spectral feature (i.e., peak) was collected 
             based on a certain experimental condition (i.e, with a photon energy E), 
             than one can find a calibration factor to multiply the x axis so 
             that the x-coordinates reflect this experimental condition.

        Args:
            values (list): value list for each spectral feature
            mode (string, optional): method used. Options are: 

                 1) 'peaks': Fit multiple peaks and calculate the distance between them 
                 (requires that `brixs.addons.fitting` is imported)
            
            deg (int, optional): a polynomial degree order used to fit the curve
                `shift vs value`.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            :py:class:`Spectrum` -> s

            s.x     = values
            s.y     = calculated shifts
            s.popt  = polynomial coeff. (highest power first.) that fit the calibration curve.
            s.model = function f(x) of the calibration curve.
            s.R2    = R2 factor of the fitting.
        """
        raise NotImplementedError('not implemented yet')

    def polyfit(self, deg, limits=None):
        """Fit data with a polynomial. Wrapper for `numpy.polyfit()`_.

        Usage:
            >>> out = polyfit(deg)
            >>> out['fit']    # --> Spectrum
            >>> out['popt']   # --> optimized parameters
            >>> out['R2']     # --> R2 error [1 - (sum( (y-fit(x) )**2)  / sum( (y - mean_y )**2) )]
            >>> out['model']  # --> function f(x) of the fitted spectrum

        Args:
            deg (int): degree of the fitting polynomial.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
        
        Returns:
            dictionary {fit, popt, R2, f(x)}

                fit (spectrum): polynomial fit spectrum with 100x more intepolated points

                popt (np.array): 1D array of polynomial coefficients 
                    (including coefficients equal to zero) from highest degree to 
                    the constant term.

                R2 (number): R2 error

                model (function): funcion f(x_centers)

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        s = self._copy(limits=limits)
        x = s.x
        y = s.y

        # fit
        popt  = np.polyfit(x, y, deg=deg)
        model = lambda x: np.polyval(popt, x)
        R2 =  1 - (sum((self.y-model(self.x))**2)/sum((self.y-np.mean(self.y))**2))
    
        start = min(x)-abs(max(x)-min(x))*0.1
        stop  = max(x)+abs(max(x)-min(x))*0.1
        _x    = np.linspace(start, stop, len(x)*100)
        arr100 = Spectrum(x=_x, y=model(_x))

        return {'fit': arr100, 'popt': popt, 'R2': R2, 'model': model}
   
    def index(self, x, closest=True):
        """Return the index value for a given x.
        
        Args:
            x (number): x value.
            closest (bool, optional): if True, x does not have to be exact and
                the index which is closest to the x value is returned. 
            
        Returns:
            int
        """
        return arraymanip.index(self.x, x, closest=closest)

    def x2y(self, x, closest=True):
        """Return the y value associated with a given x.

        Args:
            x (number): x value.
            closest (bool, optional): if True, x does not have to be exact and
                the index which is closest to the x value is returned. 
            
        Returns:
            y number
        """
        return self.y[self.index(x=x, closest=closest)]

    def get_x_where_y_is_max(self):
        """return x value where y is max"""
        return self.x[np.argmax(self.y)]
    
    def get_x_where_y_is_min(self):
        """return x value where y is min"""
        return self.x[np.argmin(self.y)]
    
    ############
    # composed #
    ############
    pass

    ##########################        
    # plot and visualization #
    ##########################        
    def plot(self, ax=None, smooth=1, label=None, limits=None, switch_xy=False, verbose=True, **kwargs):
        """Plot spectrum. Wrapper for `matplotlib.pyplot.plot()`_.

        Note:
            If `label` is `None` and if spectrum have attr 
            `label`, this attr will be used as label, e.g., 
            `plt.plot(s.x, s.y, label=s.label)`.  

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
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
            verbose (bool, optional): if True, prints warning if ploted data has
                nun-numeric values (NaN). Default is True.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            `Line2D`_

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        ##########
        # limits #
        ##########
        if limits is not None:
            temp = self._copy(limits=limits)
            if hasattr(self, 'label'):
                temp.label = self.label
            self = temp
        if verbose:
            if self.has_nan is None:
                self.check_nan()
            if self.has_nan:
                print('Warning: ploting spectrum with NaN values') 
        x = self.x
        y = self.y

        #############
        # switch xy #
        #############
        if switch_xy:
            _x = x
            x = y
            y = _x

        ##########
        # smooth #
        ##########
        if smooth > 1:
            ids = np.arange(len(x))//int(smooth)
            x = np.bincount(ids, x)/np.bincount(ids)

            ids = np.arange(len(y))//int(smooth)
            y = np.bincount(ids, y)/np.bincount(ids)

        #########
        # label #
        #########
        if label is None:
            if hasattr(self, 'label'):
                kwargs['label'] = self.label
        else:
            kwargs['label'] = label    

        ###################
        # figure and axes #
        ###################
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW or len(plt.get_fignums()) == 0:
                figure()

        ########
        # plot #
        ########
        line = ax.plot(x, y, **kwargs)
        line[0].smooth = smooth

        return line[0]

# %% ============================== Spectra ============================== %% #
class Spectra(_BrixsObject, metaclass=_Meta):
    """Returns a ``spectra`` object.

    Args:
        data (list or array, optional): list of :py:class:`spectrum` objects.

    Usage:        
            >>> ss = br.Spectra()
            >>> ss = br.Spectra([s1, s2, ...])
            >>> ss = br.Spectra(data=[s1, s2, ...])
            >>> ss = br.Spectrum().load(folderpath=<folderpath>)
            >>> ss = br.Spectrum().load(folderpath=<folderpath>, delimiter=',')
            >>> ss = br.Spectrum().load_from_single_file(filepath=<filepath>)
            >>> ss = br.Spectrum().load_from_single_file(filepath=<filepath>, delimiter=',')
            >>>
            >>> print(ss.get_core_attrs()) # print list of core attrs
            >>> print(ss.get_attrs())      # print list of attrs
            >>> print(ss.get_methods())    # print list of methods available

        Attributes:
            Every BRIXS object has 5 types of attributes: 
            `Core` , `Check`, `Modifiers`, `Labels`, `User`.

            *1. Core*
                data (list): list of Spectrum objects
            
            *2. Check*
                None

            *3. Modifiers*
                None
            
            *4. Labels*
                None

            *5. User*
                anything that the user defined on the fly.        
    """
    _read_only     = []
    _non_removable = []

    def __init__(self, data=None):
        """Initialize the object instance"""
        ###########################
        # Initializing attributes #
        ###########################
        # core
        self._data = []

        # check
        pass

        # modifiers
        pass

        # labels
        pass

        # extra
        for extra in settings._init['Spectra']:
            self.__setattr__(extra, settings._init['Spectra'][extra](self))

        ################
        # loading data #
        ################
        if data is not None:
            self.data = data
        return 

    ###################
    # core attributes #
    ###################
    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        ################
        # loading data #
        ################
        if value is None:
            self._data = []
        else:
            if isinstance(value, Iterable):
                for i, s in enumerate(value):
                    if isinstance(s, Spectrum) == False:
                        raise ValueError(f'All entries must be of type brixs.spectrum.\nEntry {i} is of type {type(s)}')
                # self._data = copy.deepcopy(value)
                self._data = list(value)
            else:
                raise ValueError('data must be a list.')
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    ###################################
    # computed (read-only) attributes #
    ###################################
    pass

    #########################
    # write-only attributes #
    #########################
    pass

    #######################
    # modifier attributes #
    #######################
    @property
    def shift(self):
        return [s.shift for s in self]
    @shift.setter
    def shift(self, value):
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'
        for i, s in enumerate(self):
            s.shift = -s.shift
            s.shift = value[i]
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def calib(self):
        return [s.calib for s in self]
    @calib.setter
    def calib(self, value):
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'
        for i, s in enumerate(self):
            s.calib = 1/s.calib
            s.calib = value[i]

    @calib.deleter
    def calib(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        return [s.factor for s in self]
    @factor.setter
    def factor(self, value):
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'
        for i, s in enumerate(self):
            s.factor = 1/s.factor
            s.factor = value[i]
    @factor.deleter
    def factor(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        return [s.offset for s in self]
    @offset.setter
    def offset(self, value):
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'
        for i, s in enumerate(self):
            s.offset = -s.offset
            s.offset = value[i]
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    #################
    # magic methods #
    #################
    def __str__(self):
        return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

    def __repr__(self):
        return str({i:val for i, val in enumerate(self.data)})[1:-1].replace(', ', '\n')

    def __setattr__(self, name, value):
        if name in settings._forbidden_words['Spectra']:
            raise AttributeError(f'`{name}` is a reserved word and cannot be set as an attribute')
        super().__setattr__(name, value)

    def __getitem__(self, item):
        if numanip.is_integer(item, allow_float=False, allow_str=False):
            return self._data[item]
        elif isinstance(item, slice):
            ss = Spectra(self._data[item])

            # # transfer attrs
            # for attr in self.get_attrs():
            #     value = copy.deepcopy(self.__dict__[attr])
            #     ss.__setattr__(attr, value)

            return ss 
        elif isinstance(item, Iterable):
            # print(item)
            # print([isinstance(_, bool) for _ in item])
            # print(sum([isinstance(_, bool) for _ in item]))
            # print(len(item))
            # assert sum([isinstance(_, bool) for _ in item]) == len(item), f'Error when trying indexing via list. All list elements must be either True or False'#\nbool item list:{item}'
            assert len(item) == len(self), f'Error when trying indexing via list. List must be the same size as number of spectra = {len(self)}'#All list elements must be either True or False'#\nbool item list:{item}'
            ss = Spectra()
            for i, s in enumerate(self):
                if item[i]:
                    ss.append(s)
            return ss 
        else:
            raise TypeError('Index must be int, slice, or list of bool, not {}'.format(type(item).__name__))

    def __setitem__(self, item, value):
        if isinstance(value, Spectrum) == False:
            raise ValueError(f'value must be of type brixs.spectrum, not {type(value)}')
        self._data[item] = value

    def __len__(self):
        return len(self.data)

    def __delitem__(self, item):
        self.remove(item)

    #########
    # attrs #
    #########
    def create_attr_from_spectra(self, attr, name=None):
        """create new Spectra attr based on Spectrum attrs.

        Example:
            >>> s1.a = 2
            >>> s2.a = 5
            >>> s3.a = 'g'
            >>> ss = br.Spectra(s1, s2, s3)
            >>> ss2 = ss.copy_attr_from_spectra('a')
            >>> print(ss2.a) # --> [2, 5, 'g']

        Args:
            attr (str): spectrum attribute name. All spectra inside Spectra
                object must have this attr. 
            name (str, optional): name of the new Spectra attr. If None, the same
                name is used.

        Returns:
            :py:class:`Spectra`
        """
        # assert all spectra has attr
        # check if all spectra has attr
        has = {i: hasattr(s, attr) for i, s in enumerate(self)}
        if False in [has[i] for i in has]:
            filtered = [i for i in has if has[i] == False]
            raise ValueError(f'Some spectra do not have attr: {attr}.\nSpectra without attr:{filtered}')

        # create new attr
        if name is None: name = attr
        
        # copy
        ss = self.copy()
        ss.__setattr__(name, [getattr(s, attr) for s in self])
        return ss

    def reorder_by_attr(self, attr, attrs2reorder=None, decreasing=False):
        """Reorder spectra based on a Spectra attr.

        Example:
            >>> ss.a = ['c', 'a', 'e', 'd', 'h', 'i']
            >>> ss.b = [ 3,   1,   5,   4,   8,   9]
            >>> ss2 = ss.reorder_by_attr(attr='b', attrs2reorder='a')
            >>> print(ss2.b) --> [1, 3, 4, 5, 8, 9]
            >>> print(ss2.a) --> ['a', 'c', 'd', 'e', 'h', 'i']
    
        Args:
            attr (str): name of the reference attr. The attr must be a list of 
                numbers with same lenght of number of spectra.
            attrs2reorder (list of str, optional): list of other Spectra attrs that 
                must also be sorted based on the ref attr.
            decreasing (bool, optional): if True, small attr value comes last.
        
        Returns:
            :py:class:`Spectra`
        """
        ################
        # empty object #
        ################
        if len(self) == 0:
            raise ValueError('cannot operate on empty spectra')
        
        # get ref attr
        ref = self.__getattribute__(attr)

        # check validity
        assert isinstance(ref, Iterable), 'ref attr must be an iterable type'
        assert len(ref) == len(self), f'Lenght of attr must be the same as the number of spectra.\nlenght of attr: {len(ref)}\nnumber of spectra: {len(self)}'
        
        if attrs2reorder is not None:
            for attr in attrs2reorder:
                temp = self.__getattribute__(attr)
                assert isinstance(temp, Iterable), f'{attr} must be an iterable type'
                assert len(temp) == len(self), f'Lenght of attr {attr} must be the same as the number of spectra.\nlenght of attr: {len(attr)}\nnumber of spectra: {len(self)}'

        # sort
        ss = self._copy()
        ss.data = arraymanip.sort(ref, ss.data)
        ss.__setattr__(attr, arraymanip.sort(ref, ref))
        if attrs2reorder is not None:
            for a in attrs2reorder:
                temp = self.__getattribute__(a)
                ss.__setattr__(a, arraymanip.sort(ref, temp))

        # flip order if decreasing is True
        if decreasing:
            ss.data = ss.data[::-1]
            ss.__setattr__(attr, self.__getattribute__(attr)[::-1])
            if attrs2reorder is not None:
                for a in attrs2reorder:
                    temp = self.__getattribute__(a)
                    ss.__setattr__(a, temp[::-1])

        return ss
    
    def get_by_attr(self, attr, value, closest=True, verbose=True):
        """Return spectrum with attr closest to value.

        Only works for numerical attr. If multiple spectrum have same value, it
        returns the first one.

        Args:
            attr (str): name of the attr
            value (number): value to look for
            closest (bool, optional): if False, returns spectrum only if value is exact.
            verbose (bool, optional): default is True, prints info about the spectrum.
        
        Returns:
            Spectrum
        """
        # check if all spectra has attr
        has = {i: hasattr(s, attr) for i, s in enumerate(self)}
        if False in [has[i] for i in has]:
            filtered = [i for i in has if has[i] == False]
            raise ValueError(f'Some spectra do not have attr: {attr}.\nSpectra without attr:{filtered}')

        # get attr
        temp = [getattr(s, attr) for s in self]

        if closest:
            i = arraymanip.index(temp, value)
        else:
            i = temp.index(value)

        if verbose:
            print(f'Spectrum number {i} have {attr}={temp[i]}')
                
        return self[i]

    def merge_duplicates(self, ref, limits=None, attrs2merge=None):
        """return spectra where spectrum with same attr are merged

        Args:
            ref (str or list): reference value for interpolating. If `str`, it
                will get values from attribute. If list, list must be the same 
                length as the number of Spectra.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            attrs2merge (list, optional): if not None, only spectra named here,
                will be copied to the final spectrum. The attr must be a list of
                numbers with the same length as the number of spectra. The value
                saved on the returned Spectrum for merged spectrum is a weighted avereged sum of 
                these attrs.

        Return:
            :py:class:`Spectra`
        """
        ################
        # empty object #
        ################
        if len(self) == 0:
            raise ValueError('cannot operate on empty spectra')
        
        #####################
        # check attr length #
        #####################
        if isinstance(ref, str):
            values = self.__getattribute__(ref)      
        elif isinstance(ref, Iterable):
            values = ref
        else:
            raise ValueError(f'`ref` must be type str or Iterable, not type `{type(ref)}`')

        ##########################
        # check length of values #
        ##########################
        assert len(values) == len(self), f'number of values ({len(values)}) must be the same as the number of spectra ({len(self)})'

        ##############################
        # check extra attrs to merge #
        ##############################
        if attrs2merge is not None:
            assert isinstance(attrs2merge, Iterable), 'attrs2merge must be a list'
            if isinstance(ref, str): attrs2merge += ref
            for attr in attrs2merge:
                temp = self.__getattribute__(attr)
                assert isinstance(temp, Iterable), f'{attr} must be an iterable type'
                assert len(temp) == len(self), f'Lenght of attr {attr} must be the same as the number of spectra.\nlenght of attr: {len(attr)}\nnumber of spectra: {len(self)}'
                assert sum([numanip.is_number(x) for x in temp]) == len(temp), f'{attr} must be a list of numbers'
        else:
            if isinstance(ref, str): 
                attrs2merge = [ref, ]
         
        ##########################################
        # check if spectra indeed has duplicates #
        ##########################################
        if arraymanip.has_duplicates(values) == False:
            return self.copy(limits=limits)
        
        ###############
        # new spectra #
        ###############
        ss = Spectra()
        ss._copy_attrs_from(self)
        if attrs2merge is not None:
            for attr in attrs2merge:
                ss.__setattr__(attr, [])

        ####################
        # merge duplicates #
        ####################
        for i, _x in enumerate(values):
            if list(values).count(_x) == 1:
                _s = self[i].copy(limits=limits)
                ss.append(_s)
                if attrs2merge is not None:
                    for attr in attrs2merge:
                        ss.__getattribute__(attr).append(self.__getattribute__(attr)[i])
            else:
                indexes = [j for j, x in enumerate(values) if x == _x]
                if i == indexes[0]:
                    _s = self.merge(indexes=indexes, limits=limits, attrs2merge=attrs2merge)
                    ss.append(_s)
                    if attrs2merge is not None:
                        for attr in attrs2merge:
                            ss.__getattribute__(attr).append(_s.__getattribute__(attr))

        return ss

    def interp_spectra(self, ref, start=None, stop=None, num=None, step=None, x=None, limits=None, attrs2interp=None):
        """create new averaged spectra 

        Args:
            ref (str or list): reference value for interpolating. If `str`, it
                will get values from attribute. If list, list must be the same 
                length as the number of Spectra.
            start (number, optional): The starting value of the sequence. If `None`,
                the minium attr value will be used.
            stop (number, optional): The end value of the sequence. If `None`,
                the maximum attr value will be used.
            num (int, optional): Number of samples to generate.
            step (number, optional): Spacing between values. This overwrites ``num``.
            x (list or array, optional): The values at which to
                evaluate the interpolated values. This overwrites all other arguments.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            attrs2interp (list, optional): if not None, only spectra named here,
                will be copied to the final spectra. The attr must be a list of
                numbers with the same length as the number of spectra. The value
                saved on the returned Spectra is a interpolates value of 
                these attrs.

        Return:
            :py:class:`Spectra`
        """
        #######################
        # check x is the same #
        #######################
        try:
            _ = self.check_same_x()
        except ValueError:
            raise ValueError('Cannot create new spectra. x axis are different.\nMaybe try interpolating the x axis (Spectra.interp())')
        
        #############################
        # check attrs2interp format #
        #############################
        assert isinstance(attrs2interp, Iterable) or attrs2interp is None, f'attrs2interp must be a list of strings or None'

        #####################
        # check attr length #
        #####################
        if isinstance(ref, str):
            values = self.__getattribute__(ref)
            if attrs2interp is None:
                attrs2interp = [ref, ]
            else:
                attrs2interp = list(attrs2interp) + [ref, ]
        elif isinstance(ref, Iterable):
            values = ref
        else:
            raise ValueError(f'`ref` must be type str or Iterable, not type `{type(ref)}`')

        ##########################
        # check length of values #
        ##########################
        assert len(values) == len(self), f'number of values ({len(values)}) must be the same as the number of spectra ({len(self)})'

        #############################
        # check values is monotonic #
        #############################
        assert arraymanip.check_monotonicity(values) == 1, f'values ({values}) must be increasingly monotonic. Use Spectra.merge_duplicates() and Spectra.reorder_by_attr()'

        ##############################
        # check extra attrs to merge #
        ##############################
        if attrs2interp is not None:
            for attr in attrs2interp:
                temp = self.__getattribute__(attr)
                assert isinstance(temp, Iterable), f'{attr} must be an iterable type'
                assert len(temp) == len(self), f'Lenght of attr {attr} must be the same as the number of spectra.\nlenght of attr: {len(attr)}\nnumber of spectra: {len(self)}'
                assert sum([numanip.is_number(x) for x in temp]) == len(temp), f'{attr} must be a list of numbers'

        ##########
        # sort x #
        ##########
        if x is None:
            if start is None: start = min(values)
            if stop is None:  stop  = max(values)
            if step is None and num is None:
                raise ValueError(f'step or num must be defined')
            elif step is not None:
                x = np.arange(start=start, stop=stop, step=step)
            else:
                x = np.linspace(start=start, stop=stop, num=num)
                
        ###############
        # new spectra #
        ###############
        ss = Spectra()
        if attrs2interp is not None:
            for attr in attrs2interp:
                ss.__setattr__(attr, [])

        ##################
        # interp spectra #
        ##################
        for j, _x in enumerate(x):
            i = arraymanip.index(values, _x, closest=True)
            if values[i] == _x:
                _s = self[i].copy(limits=limits)
                ss.append(_s)
                if attrs2interp is not None:
                    for attr in attrs2interp:
                        ss.__setattr__(attr, ss.__getattribute__(attr).append(self.__getattribute__(attr)[i]))
            elif values[i] < _x:
                _s1 = self[i].copy(limits=limits).set_factor(1 - (_x - values[i])/(values[i+1] - values[i]))
                _s2 = self[i+1].copy(limits=limits).set_factor(1 - (values[i+1] - _x)/(values[i+1] - values[i]))
                _s = Spectra(data=(_s1, _s2)).calculate_average()
                ss.append(_s)
                if attrs2interp is not None:
                    for attr in attrs2interp:
                        _v1 = self.__getattribute__(attr)[i]   * (1 - (_x - values[i])/(values[i+1] - values[i]))
                        _v2 = self.__getattribute__(attr)[i+1] * (1 - (values[i+1] - _x)/(values[i+1] - values[i]))
                        _v  = (_v1 + _v2)/2
                        ss.__setattr__(attr, ss.__getattribute__(attr).append(_v))
            elif values[i] > _x:
                _s1 = self[i].copy(limits=limits).set_factor(1 - (values[i] - _x)/(values[i] - values[i-1]))
                _s2 = self[i-1].copy(limits=limits).set_factor(1 - (_x - values[i-1])/(values[i] - values[i-1]))
                _s = Spectra(data=(_s1, _s2)).calculate_average()
                ss.append(_s)
                if attrs2interp is not None:
                    for attr in attrs2interp:
                        _v1 = self.__getattribute__(attr)[i]   * (1 - (values[i] - _x)/(values[i] - values[i-1]))
                        _v2 = self.__getattribute__(attr)[i-1] * (1 - (_x - values[i-1])/(values[i] - values[i-1]))
                        _v  = (_v1 + _v2)/2
                        ss.__setattr__(attr, ss.__getattribute__(attr).append(_v))
        return ss

    ###########
    # support #
    ###########
    def _check_limits(self, limits):
        """returns limits in the right format.

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            None or limits in the following format:
                ((xi_1, xf_1), (xi_2, xf_2), ...)               
        """
        #####################
        # if limits is None #
        #####################
        if limits is None:
            return None

        ##################################
        # assert that limits is Iterable #
        ##################################
        assert isinstance(limits, Iterable), f'`limits` must be an Iterable, not {type(limits)}'
        
        ################
        # empty object #
        ################
        if len(self) == 0:
            raise ValueError('cannot operate on empty spectra')

        ################################
        # get min and max range values #
        ################################
        vmin = min(min(s.x) for s in self)
        vmax = max(max(s.x) for s in self)

        ################
        # empty limits #
        ################
        if len(limits) == 0:
            return ((vmin, vmax),)

        ##############
        # fix format #
        ##############
        # one pair
        if len(limits) == 1: # ((xi, xf), )
            assert isinstance(limits[0], Iterable), f'wrong format for limits={limits}'
        # two pairs 
        elif len(limits) == 2: # (xi, xf), or ((xi1, xf1), (xi2, xf2))
            if isinstance(limits[0], Iterable) == False:
                if isinstance(limits[1], Iterable) == False:
                    limits = [limits, ]
                else:
                    raise ValueError(f'wrong format for limits={limits}')
        final = []
        # three or more pairs 
        for lim in limits:
            assert isinstance(lim, Iterable), f'wrong format for limits={limits}'
            temp = [lim[0], lim[1]]
            if temp[0] == None: temp[0] = vmin
            if temp[1] == None: temp[1] = vmax
            final.append((temp[0], temp[1]))   

        return final
    
    def _gather_ys(self, limits=None):
        """Return two lists, x and y's within a range.
        
        This structure speeds up some operations.
        
        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Raises:
            ValueError: if Spectra is empty
        
        Returns:
            x, y's
        """
        ###########
        # check x #
        ###########
        x = self.check_same_x()

        ######################
        # if object is empty #
        ######################
        if len(self) == 0:
            raise ValueError('cannot operate on empty spectra')
        
        #############
        # gather ys #
        #############
        length = len(x)
        ys = np.zeros((length, len(self)))
        for i in range(len(self)):
            ys[:, i] = self[i].y

        ##################
        # limits is None #
        ##################
        limits = self._check_limits(limits=limits)
        if limits is None:
            return x, ys
        
        ################################
        # get min and max range values #
        ################################
        vmin = min(min(s.x) for s in self)
        vmax = max(max(s.x) for s in self)
        
        ########################################
        # check if extract is really necessary #
        ########################################
        if len(limits) == 1:
            if limits[0][0] <= vmin and limits[0][1] >= vmax:
                return x, ys 

        ###########
        # extract #
        ###########
        x, ys = arraymanip.extract(x, ys, limits=limits)
        if len(ys.shape) == 1:
            ys = np.array([[_] for _ in ys])
        return x, ys
    
    ################
    # core methods #
    ################
    def _append(self, *args):
        """Return Spectra with spectrum appended to the spectra list.

        Usage:
            >>> ss = br.Spectra()
            >>> 
            >>> ss = ss.append(s)
            >>> ss = ss.append(s1, s2, s3)
            >>> ss = ss.append([s1, s2, s3])

        Args:
            *args (Spectrum or list): Spectrum object to be added or
                list of Spectrum.

        Returns:
            Spectra

        See Also:
            :py:func:`Spectra.remove`.
        """
        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input. Spectrum cannot be appended. Please, use one ' +\
                        'of the examples below:\n' +\
                        '\n' +\
                        'ss = br.Spectra()\n' +\
                        '\n' +\
                        'ss.append(s)\n' +\
                        'ss.append(s1, s2, ...)\n' +\
                        'ss.append([s1, s2, ...])\n'
        if args == ():
            raise AttributeError(error_message)
        
        ###############
        # create copy #
        ###############
        ss = self.copy()
        
        ################
        # loading data #
        ################        
        if len(args) == 1:
            if isinstance(args[0], Spectrum):
                ss._data += [args[0]]
            elif isinstance(args[0], Iterable):
                for i, s in enumerate(args[0]):
                    assert isinstance(s, Spectrum), f'Cannot append item {i}.\nAll items must be of type brixs.Spectrum.\nItem {i} is of type: {type(s)}'
                    ss._data += [s]
            else:
                raise AttributeError(error_message)
        elif len(args)>1:
            for i, s in enumerate(args):
                assert isinstance(s, Spectrum), f'Cannot append item {i}.\nAll items must be of type brixs.Spectrum.\nItem {i} is of type: {type(s)}'
                ss._data += [s]
        else:
            raise AttributeError(error_message)

        return ss
        
    def append(self, *args):
        """Append spectrum to the spectrum list.

        Usage:
            >>> ss = br.Spectra()
            >>> 
            >>> ss.append(s)
            >>> ss.append(s1, s2, s3)
            >>> ss.append([s1, s2, s3])

        Args:
            *args (Spectrum or list): Spectrum object to be added or
                list of Spectrum.

        Returns:
            None

        See Also:
            :py:func:`Spectra.remove`.
        """
        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input. Spectrum cannot be appended. Please, use one ' +\
                        'of the examples below:\n' +\
                        '\n' +\
                        'ss = br.Spectra()\n' +\
                        '\n' +\
                        'ss.append(s)\n' +\
                        'ss.append(s1, s2, ...)\n' +\
                        'ss.append([s1, s2, ...])\n'
        if args == ():
            raise AttributeError(error_message)
        
        ################
        # loading data #
        ################        
        if len(args) == 1:
            if isinstance(args[0], Spectrum):
                self._data += [args[0]]
            elif isinstance(args[0], Iterable):
                for i, s in enumerate(args[0]):
                    assert isinstance(s, Spectrum), f'Cannot append item {i}.\nAll items must be of type brixs.Spectrum.\nItem {i} is of type: {type(s)}'
                    self._data += [s]
            else:
                raise AttributeError(error_message)
        elif len(args)>1:
            for i, s in enumerate(args):
                assert isinstance(s, Spectrum), f'Cannot append item {i}.\nAll items must be of type brixs.Spectrum.\nItem {i} is of type: {type(s)}'
                self._data += [s]
        else:
            raise AttributeError(error_message)

    def remove(self, idx):
        """Remove spectrum.

        Args:
            idx (int): index of the spectrum.

        Returns:
            Spectra

        See Also:
            :py:func:`Spectra.append`.
        """
        ss = self.copy()
        del ss._data[idx]
        return ss

    def reorder(self, i1, i2):
        """reorder spectra.

        Args:
            i1, i2: Index of spectra to switch places.

        Returns:
            Spectra object with i1 and i2 switched
        """
        ss = self.copy()
        ss._data[i1] = self[i2]
        ss._data[i2] = self[i1]
        return ss
    
    def flip_order(self):
        """reorder spectra backwards:

        If ss.data = [s1, s2, s3], then ss2 = ss.flip_order() will make 
        ss2.data = [s3, s2, s1]

        Returns:
            Spectra with fliped order
        """
        ss = self.copy()
        ss._data = ss._data[::-1]
        return ss

    ########
    # copy #
    ########
    def _copy(self, limits=None):
        """Same as copy(), but attributes are not copied to the new object."""
        #############
        # copy data #
        #############
        data = copy.deepcopy(self.data)

        #####################
        # if limits is None #
        #####################
        limits = self._check_limits(limits=limits)
        if limits is None:
            ss = Spectra(data=data)
            return ss
        
        ##########################
        # check if empty spectra #
        ##########################
        if len(self) == 0:
            ss = Spectra(data=data)
            return ss
        
        ################################
        # get min and max range values #
        ################################
        vmin = min(min(s.x) for s in self)
        vmax = max(max(s.x) for s in self)

        ########################################
        # check if extract is really necessary #
        ########################################
        if len(limits) == 1:
            if limits[0][0] <= vmin and limits[0][1] >= vmax:
                ss = Spectra(data=data)
                return ss

        ###########
        # extract #
        ###########
        # if x is the same for all spectra, this operation is much faster
        ss = Spectra()
        try:
            x, ys = self._gather_ys(limits=limits)
            for i in range(len(self)):
                ss.append(Spectrum(x=x, y=ys[:, i]))
        except ValueError:
            for i, s in enumerate(self):
                ss.append(s.copy(limits=limits))
        return ss

    def copy(self, limits=None):
        """Return a copy of the object with data contained in a range.

        Usage:
            >>> # full copy
            >>> ss2 = ss1.copy()  # ss2 is now a copy of ss1
            >>>
            >>> # ss2 will have only data between x=0 and 10 and between x=90 and 100
            >>> ss2 = ss1.copy((0, 10), (90, 100))  

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:attr:`Spectra`
        """
        ss = self._copy(limits=limits)
        ss._copy_attrs_from(self)

        # extra
        for extra in settings._copy['Spectra']:
            if hasattr(ss, extra):
                ss.__setattr__(extra, self.__getattribute__(extra).copy(ss))
        for extra in settings._copy['Spectrum']:
            for s in ss:
                if hasattr(s, extra):
                    s.__setattr__(extra, self.__getattribute__(extra).copy(s))   
        return ss
    
    #################
    # save and load #
    #################
    def save(self, folderpath, prefix='spectrum_', suffix='.dat', filenames=None, zfill=None, number_of_decimal_places=None, only_data=False, verbose=False, **kwargs):
        r"""Save spectra. Wrapper for `numpy.savetxt()`_.

        Warning:
            Attrs are saved as comments if only_data is False. Saving attrs to file
            is not always reliable because requires converting variables to string. 
            Only attrs that are of type: string, number, and list of number,
             list of list of number and strings have been tested. Dictionaries are not saved.

        Args:
            folderpath (string or pathlib.Path): folderpath, folder handle. 
                If None, filepath can be inferred from a
                user defined attr 'ss.dirpath'. Default is None. Last used 
                dirpath is saved to ss.dirpath.
            prefix (string, optional): prefix used for naming the files.
            suffix (string, optional): suffix used for naming the files. If the
                filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            filenames (string, optional): overwrites prefix and suffix. Use `{<attr>}` to 
                include spectrum attr in the name. Use `{i}` to include the
                spectrum index. Example: 'spectrum_{i}_T{T}K.dat'
            zfill (int, optional): number of digits for file numbering. If `None`,
                zfill will be determined.
            number_of_decimal_places (int, optional): if not None, this argument
                defines the number of decimal places to save the data (it does
                not affect attrs, only x and y arrays). The 'fmt' argument overwrites
                number_of_decimal_places. Default is None. If None, the number
                of decimal places will be such that data will be saved with the 
                best precision necessary.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            verbose (bool, optional): turn verbose on and off. Default is `False`.
            **kwargs: kwargs are passed to ``np.savetxt()`` that saves the data.
            
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
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        ##############################
        # check for wrong parameters #
        ##############################
        if 'filepath' in kwargs:
            raise ValueError('filepath is not a valid parameter.\nPlease, use folderpath')
        
        ##############
        # folderpath #
        ##############
        folderpath = Path(folderpath)
        assert folderpath.exists(), f'folderpath does not exists.\nfolderpath: {folderpath}'
        assert folderpath.is_dir(), f'folderpath is not a directory.\nfolderpath: {folderpath}'

        #################
        # set filenames #
        #################
        if zfill is None:
            zfill = numanip.n_digits(len(self)-1)

        # saving
        # if verbose: print('saving {} files...')
        for i, s in enumerate(self.data):
            i_str = str(i).zfill(zfill)
            if filenames is not None:
                filename = filenames
                for slot in filenames.split('{'):
                    if '}' in slot:
                        name = slot.split('}')[0]
                        if '{' + name + '}' == '{i}':
                            filename = filename.replace('{i}', i_str)
                        else:
                            filename = filename.replace('{' + name + '}', str(getattr(self[i], name)))
            else:
                filename = f'{prefix}' + f'{i_str}' + f'{suffix}'

            if verbose:  print(f'{i}/{len(self)-1}: {filename}')
            s.save(filepath=folderpath/filename, number_of_decimal_places=number_of_decimal_places, only_data=only_data, check_overwrite=False, verbose=False, **kwargs)
        # if verbose: print('Done!')
    
    def load(self, folderpath, string='*', only_data=False, verbose=False, **kwargs):
        """Load data from text files. Wrapper for `numpy.genfromtxt()`_.

        Args:
            folderpath (string or path object): folderpath. If None,
                Last used folderpath is used.
            string (str, optional): file names without this string will be ignored.
                Use '*' for matching anything. Default is '*'.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is loaded. Default is False.
            verbose (bool, optional): If True, it prints the progress of loading. 
                Default is `False`.
            **kwargs: kwargs are passed to ``numpy.genfromtxt()`` that loads the data.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``#``. Attributes picked up
                from the header will be loaded too.
            usecols (tuple, optional): Default is (0, 1).

        Returns:
            Spectra with loaded spectrum

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        ##############################
        # check for wrong parameters #
        ##############################
        if 'filepath' in kwargs:
            raise ValueError('filepath is not a valid parameter.\nPlease, use folderpath')
                
        #################################################################
        # check if folder points to a folder and that folderpath exists #
        #################################################################
        folderpath = Path(folderpath)
        assert isinstance(folderpath, str) or isinstance(folderpath, Path), f'folderpath is unknown format\n{folderpath}'
        assert folderpath.exists(), f'folderpath does not exist\n{folderpath}'
        assert folderpath.is_dir(), f'folderpath must point to a file\n{folderpath}'

        ################
        # get filelist #
        ################
        fl = filemanip.filelist(dirpath=folderpath, string=string)
        
        ############
        # get data #
        ############
        ss = Spectra()
        for i, filepath in enumerate(fl):
            if verbose: print(f'    {i+1}/{len(fl)}: {filepath.name}')
            s = Spectrum().load(filepath=filepath, only_data=only_data, **kwargs)
            ss.append(s)
        return ss

    def save_all_single_file(self, filepath, number_of_decimal_places=None, only_data=False, limits=None, check_overwrite=False, verbose=False, **kwargs):
        r"""Save all Spectra in one single file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            number_of_decimal_places (int, optional): if not None, this argument
                defines the number of decimal places to save the data (it does
                not affect attrs, only x and y arrays). If list, it must have two 
                values corresponding to the number_of_decimal_places to be used for 
                to x and y separately. The 'fmt' argument overwrites
                number_of_decimal_places. Default is None. If None, the number
                of decimal places will be such that data will be saved with the 
                best precision necessary.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            limits (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                The additive factor will be calculated for spectra such the average value of
                data inside limits is the same as for the first spectrum. Spectra
                outside this range will not be saved.

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
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        ##############################
        # check for wrong parameters #
        ##############################
        if 'folderpath' in kwargs:
            raise ValueError('folderpath is not a valid parameter.\nPlease, use filepath')
                
        ######################################
        # check if filepath points to a file #
        ######################################
        filepath = Path(filepath)
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        #######################
        # check x is the same #
        #######################
        try:
            _x = self.check_same_x()
        except ValueError:
            raise ValueError('Cannot save spectra in one file. x axis are different.\nMaybe try interpolating the x axis (Spectra.interp()) or use Spectra.save() to save spectra in multiple files.')
        
        ###################
        # check overwrite #
        ###################
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if other.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')
                
        ##########
        # kwargs #
        ##########
        kwargs.setdefault('delimiter', ', ')
        kwargs.setdefault('newline', '\n')
        kwargs.setdefault('comments', '# ')

        ##############################
        # Prepare final data to save #
        ##############################
        if len(self) == 0:
            final = []
        else:
            # gather ys
            x, ys = self._gather_ys(limits=limits)

            # kwargs
            if 'fmt' in kwargs: # pick best format
                pass
            elif number_of_decimal_places is not None:
                if isinstance(number_of_decimal_places, Iterable):
                    kwargs['fmt'] = [f'%.{number_of_decimal_places[0]}f'] + [f'%.{number_of_decimal_places[1]}f']*len(temp)
                else:
                    kwargs['fmt'] = f'%.{number_of_decimal_places}f'
            else: # pick best format
                has_nan = self.check_nan()
                if has_nan:
                    raise ValueError('Cannot save spectra with NaN values. Please, use ss.remove_nan() to remove NaN values.')
                else:
                    temp = self.copy()

                number_of_decimal_places_x = max([numanip.n_decimal_places(x) for x in _x])
                number_of_decimal_places_y = 0
                for s in temp:
                    t = max([numanip.n_decimal_places(y) for y in s.y])
                    if t > number_of_decimal_places_y:
                        number_of_decimal_places_y = t
                kwargs['fmt'] = [f'%.{number_of_decimal_places_x}f'] + [f'%.{number_of_decimal_places_y}f']*len(temp)

            # final data to save
            final = np.zeros((len(_x), len(self)+1))
            final[:, 0]  = x
            final[:, 1:] = ys
            
        #####################
        # header and footer #
        #####################
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
        else:
            attrs_dict = {_:self.__getattribute__(_) for _ in settings._reserved_words['Spectra']['vars'] if _ not in ['_data', ]}
            attrs_dict.update(self.get_attrs_dict())
            header = '\n'.join(_attr2str(attrs_dict, verbose)) + '\n'

            if 'header' not in kwargs:
                kwargs['header'] = header
            else:
                if kwargs['header'] == '':
                    kwargs['header'] = header
                elif kwargs['header'][-1] != '\n':
                    kwargs['header'] += '\n'
        ########
        # save #
        ########
        np.savetxt(Path(filepath), final, **kwargs)

    def load_from_single_file(self, filepath, only_data=False, verbose=False, **kwargs):
        """load multiple spectra from file.

        Args:
            filepath (str or pathlib.Path): filepath. If None, the 
                last used filepath is used.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            verbose (book, optional): Default is False. If True, it will print
                an warning when attributes cannot be loaded from the file.
            **kwargs: kwargs are passed to ``numpy.genfromtxt()`` that loads the data.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``#``. Attributes picked up
                from the header will be loaded too.
            usecols (tuple, optional): Default is (0, 1, ...). First number in 
                usecols must point to the x axis.
            
        Returns:
            Spectra with loaded data

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        ##############################
        # check for wrong parameters #
        ##############################
        if 'folderpath' in kwargs:
            raise ValueError('folderpath is not a valid parameter.\nPlease, use filepath')
        
        ######################################
        # check if filepath points to a file #
        ######################################
        filepath = Path(filepath)
        assert filepath.is_file(), f'filepath must point to a file\n{filepath}'
        assert filepath.exists(), f'filepath does not exist\n{filepath}'

        ##########
        # kwargs #
        ##########
        kwargs.setdefault('delimiter', ', ')
        kwargs.setdefault('comments', '# ')
        kwargs.setdefault('usecols', None)
        # kwargs['usecols'] = [i for i in range(data.shape[1])]

        #############
        # read data #
        #############
        data = np.genfromtxt(Path(filepath), **kwargs)
        x    = data[:, 0]
        ys   = data[:, 1:]

        ##########
        # assign #
        ##########
        ss = Spectra()
        for i in range(ys.shape[1]):
            s = Spectrum(x=x, y=ys[:, i])
            ss.append(s)

        ###############
        # read header #
        ###############
        if only_data is False:
            # get header
            header = filemanip.load_comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            
            # remove comment flag (#)
            comment_flag_len = len(kwargs['comments'])
            if header:
                for i, line in enumerate(header):
                    header[i] = line[comment_flag_len:]

                # attrs dict
                attrs_dict = _str2attr(header[:-1], verbose=verbose)

                # set attrs
                for attr in attrs_dict:
                    ss.__setattr__(attr, attrs_dict[attr])
        return ss

    #########
    # check #
    #########
    def check_nan(self):
        """Check if at least one spectrum have non-numeric (NaN) values
        
        Returns:
            bool (True if has nan)
        """
        ################
        # empty object #
        ################
        if self.data is None:
            return False
        
        #############
        # check nan #
        #############
        for s in self:
            if s.has_nan is None: 
                s.check_nan()
            if s.has_nan:
                return True
        return False
    
    def find_nan(self):
        """Return a dict with spectra indexes and index where non-numeric (NaN) values were found
        
        Returns:
            dict
        """
        nan_places = {}
        for i, s in self:
            nan_places[i] = s.find_nan()
        return nan_places
    
    def remove_nan(self):
        """remove data points (x, y) that contain non-numeric (NaN) values
        
        Returns:
            Spectra
        """
        has_nan = self.check_nan()

        if has_nan == False:
            return self.copy()
        else:
            ss = self.copy()

            for i, s in enumerate(ss):
                ss[i] = s.remove_nan()

            return ss    
    
    def check_monotonicity(self):
        """Sets monotonicity attribute to 'increasing' or 'decreasing'.

        Raises:
            ValueError if data is not monotonic.

        Returns:
            string ('increasing' or 'decreasing')
        """
        ########################
        # check empty spectrum #
        ########################
        if len(self) == 0:
            raise ValueError('cannot check monotonicity for empty spectra')
        
        monotonicity = [None]*len(self)
        for i in range(len(self)):
            try:
                self[i].check_monotonicity()
                monotonicity[i] = self.data[i].monotonicity
            except ValueError:
                pass
        if all(x == 'increasing' for x in monotonicity):
            return 'increasing'
        elif all(x == 'decreasing' for x in monotonicity):
            return 'decreasing'
        else:
            text = ''
            for i in range(len(self)):
                text += f'spectrum: {i}, motonicity: {monotonicity[i]}\n'
            raise ValueError(f'some spectra have different monotonicity (increasing, decreasing) or no monotonicity at all (None): \n{text}')

    def fix_monotonicity(self, mode='increasing'):
        """Rearrange x, y such as x array is monotonically increasing or decreasing.

        Args:
            mode (str, optional): increasing or decreasing.

        Returns:
            Spectra
        """
        ########################
        # check empty spectrum #
        ########################
        if len(self) == 0:
            raise ValueError('cannot operate on empty spectra')
    
        ss = self.copy()
        for i, s in enumerate(ss):
            ss[i] = s.fix_monotonicity(mode=mode)
        return ss

    def check_length(self):
        """Checks if all spectra has the same length.

        If all spectra have the same length, it sets 
        :py:attr:`Spectra.length` = length.
        Otherwise, it raises an error.

        Returns:
            number (length)

        Raises:
            ValueError: spectra does not have the same length.

        See Also:
            :py:func:`Spectra.check_step`, :py:func:`Spectra.check_same_x`.
        """
        ########################
        # check empty spectrum #
        ########################
        if len(self) == 0:
            raise ValueError('cannot check length for empty spectra')
        
        # if only one spectra exists, then length is immediately defined
        if len(self) == 1:
            return len(self[0].x)
        
        # collect
        length = [None]*len(self)
        for i in range(len(self)):
            try:
                length[i] = len(self.data[i])
            except ValueError:
                pass

        # apply or raise error
        if all(x == length[0] for x in length):
            return length[0]
        else:
            text = ''
            for i in range(len(self)):
                text += f'spectrum: {i}, length: {length[i]}\n'
            raise ValueError(f'some spectra have different length: \n{text}')

    def check_step(self, max_error=0.1):
        """Check step between data points in the x-coordinates.

            If data has a well defined step size, it sets :py:attr:`Spectra.step` = step.
            Otherwise, it raises an error.

            Args:
                max_error (number, optional): percentage value 
                (of the average x-step) of the maximum allowed error. 
                Default is 0.1 %.                

            Note:
                This method checks if 1) the step between two data points is the 
                same through out the x vector for each spectrum (vector 
                uniformity). Step uniformity within each spectrum is verified by 
                the following equation (see Spectrum.check_step()):

                    (max(steps) - min(steps))/np.mean(steps) * 100 < max_error

                This method also checks 2) if the step is the same between 
                different spectra using the follwing equation

                    (step[i]-step[i+1])/avg(step) * 100 < max_error

            Returns:
                step

            Raises:
                ValueError: If condition 1, or 2 are not satisfied.

            See Also:
                :py:func:`Spectra.check_length`, :py:func:`Spectra.check_same_x`
        """
        ########################
        # check empty spectrum #
        ########################
        if len(self) == 0:
            raise ValueError('cannot check step for empty spectra')

        # if all spectra have the same x, check step becames easier
        try: 
            _x = self.check_same_x(max_error=max_error)
            temp = Spectrum(x=_x)
            try:
                temp.check_step(max_error=max_error)
            except ValueError:
                raise ValueError(f"Spectra have the same x-coordinates, but it is not uniform.")
            return temp.step
        except ValueError:
            # 1) check step uniformity
            steps = ['not uniform']*len(self)
            for idx, s in enumerate(self.data):
                if s.step is None:
                    try:
                        s.check_step(max_error=max_error)
                        steps[idx] = s.step
                    except ValueError:
                        pass
                else:
                    steps[idx] = s.step
            # raise error
            if all(x == steps[0] for x in steps):
                pass
            else:
                text = ''
                for i in range(len(self)):
                    text += f'spectrum: {i}, step: {steps[i]}\n'
                raise ValueError(f'some spectra have different step: \n{text}')

            # 2) check step between spectra
            avg_step = np.mean(steps)
            if sum([abs(steps[i]-steps[i+1]) > abs(avg_step*max_error/100) for i in range(len(self)-1)]) > 0:
                raise ValueError(f"Spectra seems to have different step size. Calculated step sizes: \n{text}")
            return avg_step

    def check_same_x(self, max_error=0.1):
        """Check if spectra have same x-coordinates.

        If data has same x-coordinates, it sets :py:attr:`Spectra.x` = x.
        Otherwise, it raises an error.

        Args:
            max_error (number, optional): percentage value (in terms of the
                average x-step) of the maximum allowed error. If None, the 
                Default value from settings will be used.

                max(s[i].x - s[i+1].x)/avg(step) * 100 < max_error

        Returns:
            x values

        Raises:
            ValueError: If any x-coodinate of any two spectrum is different.

        See Also:
            :py:func:`Spectra.check_length`, :py:func:`Spectra.check_step`.
        """
        ########################
        # check empty spectrum #
        ########################
        if len(self) == 0:
            raise ValueError('cannot check same x for empty spectra')

        # if only one spectra exists, then x is immediately defined
        if len(self) == 1:
            return self[0].x

        # if empty spectrum exist
        text  = ''
        empty = False
        for i, s in enumerate(self):
            if len(s) == 0:
                text += f'spectrum: {i}: empty\n'
                empty = True
            else:
                text += f'spectrum: {i}: ok\n'
        if empty:
            raise ValueError(f'some spectra are empty\n{text}')
        
        # check length (will raise an error if different length)
        length = self.check_length()

        # average step
        step = []
        for s in self:
            step.append(np.mean(np.diff(s.x)))
        # step = step/len(self)
        # if step == 0:

        # if self.step is None:
        #     try:
        #         self.check_step()
        #         step = self.step
        #     except ValueError:
        #         raise ValueError(f'some spectra have different x: \n{text}\n\nUse brixs.Spectra.interp() to interpolate the data and make the x axis for different spectra match.') 
        #         step = 0
        #         for s in self:
        #             step += np.mean(np.diff(s.x))
        #         step = step/len(self)
        # else:
        #     step = self.step

        # check x between spectra
        x = ['same as the next']*len(self)
        for idx in range(len(self)-1):
            if max(abs(self[idx].x - self[idx+1].x))*100/abs(step[idx]) > max_error:
                x[idx] = 'different'

        # apply or raise error
        if 'different' not in x:
            return self[0].x
        else:
            text = ''
            for i in range(len(self)):
                text += f'spectrum: {i}: {x[i]}\n'
            raise ValueError(f'some spectra have different x: \n{text}\n\nUse brixs.Spectra.interp() to interpolate the data and make the x axis for different spectra match.')

    ##################
    # BASE modifiers #
    ##################
    def set_shift(self, value):
        """Shift data recursively.

        Args:
            value (number or list): value will be added to x-coordinates. If list,
                number of values must be the same as number of Spectra.

        Returns:
            :py:class:`Spectra`

        See Also:
            :py:func:`Spectra.calculate_shift()`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_shift(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        # extra
        for extra in settings._shift['Spectra']:
            if hasattr(s, extra):
                ss.__getattribute__(extra).set_shift(value)

        return ss

    def set_calib(self, value):
        """Apply multiplicative x factor recursively.

        Args:
            value (number or list): Multiply x-coordinates by value. If 
                list, number of values must be the same as number of Spectra.

        Returns:
            :py:class:`Spectra`
        
        See Also:
            :py:func:`Spectra.calculate_calib()`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)
        
        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_calib(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        # extra
        for extra in settings._calib['Spectra']:
            if hasattr(s, extra):
                ss.__getattribute__(extra).set_shift(value)
        
        return ss

    def set_factor(self, value):
        """Apply multiplicative y factor recursively.

        Args:
            value (number or list): Multiply y-coordinates by value. If 
                list, number of values must be the same as number of Spectra.

        Returns:
            :py:class:`Spectra`

        See Also:
            :py:func:`Spectra.calculate_factor()`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_factor(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        # extra
        for extra in settings._factor['Spectra']:
            if hasattr(s, extra):
                ss.__getattribute__(extra).set_shift(value)

        return ss

    def set_offset(self, value):
        """Apply additive y factor recursively.

        Args:
            value (number or list): value will be added to x-coordinates. If list,
                number of values must be the same as number of Spectra.

        Returns:
            :py:class:`Spectra`

        See Also:
            :py:func:`Spectra.calculate_offset()`, :py:func:`Spectra.floor()`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_offset(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        # extra
        for extra in settings._factor['Spectra']:
            if hasattr(s, extra):
                ss.__getattribute__(extra).set_shift(value)

        return ss

    #############
    # modifiers #
    #############
    def set_x_via_polyval(self, p):
        """Set x to np.polyval(p, x).

        Args:
            p (array or list of arrays): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term. Or list of 1D array of polynomial coefficients with same
                length of number of Spectra.

        Returns:
            :py:class:`Spectra`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(p, Iterable) == False:
            p = [p]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(p) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(ss):
            ss[i] = s.set_x_via_polyval(p)

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        return ss

    def set_x_via_function(self, f):
        """Set x to f(x).

        Args:
            f (function): function where argument is x-coordinate elements.
                Or list of functions with same length of number of Spectra.

        Returns:
            :py:class:`Spectra`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(f, Iterable) == False:
            f = [f]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(f) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(ss):
            ss[i] = s.set_x_via_function(f)

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        return ss

    def set_roll(self, value):
        """Roll array elements of the y-coordinate [same as ss.set_y_roll()]

        Note:
            Elements that roll beyond the last position are re-introduced at the
            first.

        Args:
            value (number or list): value to roll y-coordinates. If list,
                number of values must be the same as number of Spectra.

        Returns:
            :py:class:`Spectra`

        See Also:
            :py:func:`Spectra.set_shift()`, :py:func:`Spectra.calculate_shift()`, and :py:func:`Spectra.calculate_roll()`
        """
        return self.set_y_roll(value)
    
    def set_y_roll(self, value):
        """Roll array elements of the y-coordinate

        Note:
            Elements that roll beyond the last position are re-introduced at the
            first.

        Args:
            value (number or list): value to roll y-coordinates. If list,
                number of values must be the same as number of Spectra.

        Returns:
            :py:class:`Spectra`

        See Also:
            :py:func:`Spectra.set_shift()`, :py:func:`Spectra.calculate_shift()`, and :py:func:`Spectra.calculate_roll()`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_y_roll(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        return ss

    def set_x_roll(self, value):
        """Roll array elements of the x-coordinates

        Note:
            Elements that roll beyond the last position are re-introduced at the
             first.

        Args:
            value (number or list): value to roll x-coordinates. If list,
                number of values must be the same as number of Spectra.

        Returns:
            :py:class:`Spectra`

        See Also:
            :py:func:`Spectra.set_shift()`, :py:func:`Spectra.calculate_shift()`, and :py:func:`Spectra.calculate_roll()`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_x_roll(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        return ss

    def flip_x(self):
        """Flip x axis of every spectra

        Returns:
            :py:class:`Spectra`
        """
        return self.set_calib(value=-1)
    
    def flip_y(self):
        """Flip y axis of every spectra

        Returns:
            :py:class:`Spectra`
        """
        return self.set_factor(value=-1)

    def floor(self, limits=None):
        """Sets zero value for y-coordinates (shifts data verticaly).

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectra`
        """
        ############################
        # check for empty spectrum #
        ############################
        text = ''
        flag = False
        for i, s in enumerate(self):
            if s.x is None:
                flag = True
                text += f'{i}: empty\n'
            else:
                text += f'{i}: ok\n'
        if flag: raise ValueError(f'some spectra are empty:\n{text}')

        # floor
        ss   = self.copy()
        text = ''
        flag = False
        for i, s in enumerate(self):
            try:
                ss[i] = s.floor(limits=limits)
                text += f'{i}: ok\n'
            except AssertionError:
                flag = True
                _s   = Spectrum()
                _s._copy_attrs_from(s)
                ss[i] = s
                text += f'{i}: empty\n'
        if flag:
            print(f'Warning: some spectra has no datapoints within limits={limits}\n{text}')
        return ss

    def normalize(self, value=1, limits=None):
        """Set a factor such as the average y between limits is equal to value.

        Args:
            value (number): value. Default is 1.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectra`

        See Also:
            :py:func:`Spectra.calculate_factor`
        """
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.normalize(value=value, limits=limits)
        return ss

    ############
    # advanced #
    ############
    def concatenate(self):
        """Return spectrum of concatenate spectra.
        
        Note:
            attrs are copied to the returned spectrum, but attrs from each spectrum is lost.

        Returns:
            :py:class:`Spectrum`
        """
        x = np.concatenate([s.x for s in self.data])
        y = np.concatenate([s.y for s in self.data])
        s = Spectrum(x=x, y=y)
        s._copy_attrs_from(self)

        return s

    def interp(self, start=None, stop=None, num=None, step=None, x=None):
        """Interpolate data.

        Args:
            start (number, optional): The starting value of the sequence. If `None`,
                the minium x value will be used.
            stop (number, optional): The end value of the sequence. If `None`,
                the maximum x value will be used.
            num (int, optional): Number of samples to generate.
            step (number, optional): Spacing between values. This overwrites ``num``.
            x (list or array, optional): The x-coordinates at which to
                evaluate the interpolated values. This overwrites all other arguments.

        None:
            If none arguments is given, this function will adjust all spectra to
            have the same x.

        Returns:
            :py:class:`Spectra`
        """
        # sort x
        if x is None:
            try:
                _x = self.check_same_x()
                if start is not None or stop is not None or num is not None or step is not None:
                    if start is None:
                        start = max(_x)
                    if stop is None:
                        stop = min(_x)
                    if num is None:
                        num = len(_x)
                else:
                    return self
            except ValueError:
                if start is None:
                    start = min([min(s.x) for s in self.data])
                if stop is None:
                    stop = max([max(s.x) for s in self.data])
                if num is None:
                    num = max([len(s.x) for s in self.data])              

        # new spectra object
        ss = self.copy()

        # interpolate
        for i, s in enumerate(self):
            ss[i] = s.interp(x=x, start=start, stop=stop, num=num, step=step)

        return ss
    
    def crop(self, start=None, stop=None):
        """Crop data.

        Args:
            start (number, optional): start x value. If None, the minimum value of
                x will be used.
            stop (number, optional): final x value. If None, the maximum value of
                x will be used.

        Returns:
            None
        """
        # sort start stop
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

        # new spectra object
        ss = self.copy()

        # crop
        for i, s in enumerate(self):
            ss[i] = s.crop(start, stop)

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None  

        return ss

    def switch_xy(self):
        """Switch x and y axis"""
        # copy
        ss = self.copy()

        # interpolate
        for i, s in enumerate(self):
            ss[i] = s.switch_xy()

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        return ss

    def derivative(self, order=1):
        """Returns Spectra object with the derivative of each spectrum.

        Args:
            order (number, optional): derivative order. Default is 1.

        Returns:
            Spectra object with derivate spectra
        """
        # assert spectrum is not empty
        assert len(self) > 0, 'Cannot act on a empty Spectra'

        # copy
        ss = self.copy()

        # interpolate
        for i, s in enumerate(self):
            ss[i] = s.derivative(order=order)

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        return ss

    def moving_average(self, n):
        """Returns Spectra with the moving average of each spectrum.

        Each smoothed Spectrum s inside Spectra has length given by (len(s)-n+1).

        Args:
            n (int): number of points to average.

        Returns:
            :py:class:`Spectra` with spectrum of length given by (len(s.x)-n+1)
        """
        # assert spectrum is not empty
        assert len(self) > 0, 'Cannot act on a empty Spectra'

        # copy
        ss = self.copy()

        # interpolate
        for i, s in enumerate(self):
            ss[i] = s.moving_average(n=n)

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        return ss

    def smooth(self, n, force_divisible=False):
        """Returns Spectra with the moving average of each spectrum.

        Args:
            n (int): number of points to average.
            force_divisible (bool, optional): if True, the length of the data
                must be divisible by n. Default is False.

        Warning:        
            If force_divisible=False, the last data point of each spectrum might be averaged by 
            less points then n.

        Returns:
            :py:class:`Spectra` spectrum of length given by (len(s.x)/n).
        """
        # assert spectrum is not empty
        assert len(self) > 0, 'Cannot act on a empty Spectra'

        # copy
        ss = self.copy()

        # interpolate
        for i, s in enumerate(self):
            ss[i] = s.smooth(n=n, force_divisible=force_divisible)

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        return ss


    def merge(self, indexes, weights=None, limits=None, attrs2merge=None):
        """Return averaged Spectrum 

        Args:
            indexes (list): spectra index numbers
            weights (list or None, optional): if not None, applies a weighted 
                average. Default is None.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            attrs2merge (list, optional): if not None, only spectra named here,
                will be copied to the final spectrum. The attr must be a list of
                numbers with the same length as the number of spectra. The value
                saved on the returned Spectrum is a weighted avereged sum of 
                these attrs.
        
        Return:
            :py:class:`Spectrum`
        """
        ################
        # empty object #
        ################
        if len(self) == 0:
            raise ValueError('cannot operate on empty spectra')
        
        assert isinstance(indexes, Iterable), 'indexes must be a list'

        # check weights
        if weights is None:
            weights = [1]*len(indexes)
        
        assert isinstance(weights, Iterable), 'weights must be a list'
        assert len(weights) == len(indexes), 'number of indexes must be the same as the number of weights'

        ########################
        # check attrs to merge #
        ########################
        if attrs2merge is not None:
            for attr in attrs2merge:
                temp = self.__getattribute__(attr)
                assert isinstance(temp, Iterable), f'{attr} must be an iterable type'
                assert len(temp) == len(self), f'Lenght of attr {attr} must be the same as the number of spectra.\nlenght of attr: {len(attr)}\nnumber of spectra: {len(self)}'
                assert sum([numanip.is_number(x) for x in temp]) == len(temp), f'{attr} must be a list of numbers'

        ss = Spectra()
        for j, i in enumerate(indexes):
            ss.append(self[i].copy().set_factor(weights[j]))

        if attrs2merge is not None:
            for attr in attrs2merge:
                temp = self.__getattribute__(attr)
                new = []
                for j, i in enumerate(indexes):
                    new.append(temp[i]*weights[j])
                ss.__setattr__(attr, np.mean(new))

        s = ss.calculate_average(limits=limits)
        # double check if attrs were copied (necessary for `_` attrs)
        for attr in attrs2merge:
            if hasattr(s, attr) == False:
                s.__setattr__(attr, ss.__getattribute__(attr))
        return s
    
    def merge_and_replace(self, indexes, weights=None, limits=None, attrs2merge=None):
        """return spectra with replaced spectra

        Args:
            indexes (list): spectra index numbers
            weights (list or None, optional): if not None, applies a weighted 
                average. Default is None.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            attrs2merge (list, optional): if not None, only spectra named here,
                will be copied to the final spectrum. The attr must be a list of
                numbers with the same length as the number of spectra. The value
                saved on the returned Spectrum for merged spectrum is a weighted avereged sum of 
                these attrs.

        Return:
            :py:class:`Spectra`
        """
        ##############################
        # check extra attrs to merge #
        ##############################
        if attrs2merge is not None:
            for attr in attrs2merge:
                temp = self.__getattribute__(attr)
                assert isinstance(temp, Iterable), f'{attr} must be an iterable type'
                assert len(temp) == len(self), f'Lenght of attr {attr} must be the same as the number of spectra.\nlenght of attr: {len(attr)}\nnumber of spectra: {len(self)}'
                assert sum([numanip.is_number(x) for x in temp]) == len(temp), f'{attr} must be a list of numbers'

        ss = Spectra()
        if attrs2merge is not None:
            for attr in attrs2merge:
                ss.__setattr__(attr, [])

        for i, s in enumerate(self):
            if i in indexes:
                if list(indexes).index(i) == 1:
                    _s = self.merge(indexes=indexes, weights=weights, limits=limits, attrs2merge=attrs2merge)
                    ss.append(_s)
                    if attrs2merge is not None:
                        for attr in attrs2merge:
                            ss.__setattr__(attr, ss.__getattribute__(attr).append(_s.__getattribute__(attr)))
            else:
                ss.append(s.copy(limits=limits))
                if attrs2merge is not None:
                    for attr in attrs2merge:
                        ss.__setattr__(attr, ss.__getattribute__(attr).append(self.__getattribute__(attr)[i]))
        return ss

    ########################
    # calculation and info #
    ########################
    def calculate_sum(self, limits=None):
        """Returns Spectrum object with the sum of all spectra.

        Note:
                All spectra must have the same x-coordinates. This is verified.

        Warning:
            attrs are copied to the final spectrum, but attrs from each spectrum is lost.

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectrum` object.
        """
        # gather ys
        x, ys = self._gather_ys(limits=limits)

        # calculate sum
        y = np.zeros(len(x))
        for i in range(len(self)):
            y += ys[:, i]
        s = Spectrum(x=x, y=y)
        s._copy_attrs_from(self)

        return s

    def calculate_average(self, limits=None):
        """Returns Spectrum object with the average of all spectra.

        Note:
            All spectra must have the same x-coordinates. This is verified.

        Warning:
            attrs are copied to the final spectrum, but attrs from each spectrum is lost.

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectrum`
        """
        # gather ys
        x, ys = self._gather_ys(limits=limits)

        # calculate sum
        y = np.zeros(len(x))
        for i in range(len(self)):
            y += ys[:, i]
        y = y/len(self)
        s = Spectrum(x=x, y=y)
        s._copy_attrs_from(self)

        return s


    def calculate_map(self, axis=0, limits=None):
        """Return image representation of spectra.

        Note:
            All spectra must have the same x-coordinates. This is verified.

        Warning:
            attrs are copied to the final image, but attrs from each spectrum is lost.

        Args:
            axis (int, optional): Image axis along which spectra will be laid out.
                If `axis=0`, spectra will be placed horizontally (each spectrum 
                will be a "row of pixels"). If `axis=1`, spectra will be placed
                vertically (each spectrum will be a "column"). Default is 0.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Image`.
        """
        # check axis
        if axis != 0 and axis != 1:
            raise ValueError('axis must be 0 or 1')
        
        if axis == 0:
            return self.stack_spectra_as_columns(limits=limits)
        else:
            return self.stack_spectra_as_columns(limits=limits)
        
        # # gather ys
        # y, ys = self._gather_ys(limits=limits)


        # if axis == 1:
        #     ys = ys.transpose()
        #     x = y
        #     y = centers

        # im = Image(data=ys)
        # im._copy_attrs_from(self)
        # im.x_centers = x
        # im.y_centers = y

        # return im
    
    def stack_spectra_as_columns(self, x_centers=None, limits=None):
        """Return image representation of spectra, where spectra are layed out in the vertical direction.

        Note:
            All spectra must have the same x-coordinates. This is verified.

        Warning:
            attrs are copied to the final image, but attrs from each spectrum is lost.

        Args:
            x_centers (list, optional): numerical value for each spectrum. 
                 indicating a numerical value for each spectrum. Default is None, i.e., each
                spectrum will be numbered from 1 to len(ss).
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Image`.
        """
        # check centers
        if x_centers is not None:
            assert len(x_centers) == len(self), f'centers must have the same number of items as the number of spectra.\nnumber of centers: {len(centers)}\nnumber of spectra: {len(self)}'

        # check if array is monotonic
        # This is necessary because the way images are plotted
        # The 'lower' value of the y-axis will be to the top
        if self.check_monotonicity() != 'increasing':
            raise ValueError('monotonicity of data must be strictly increasing. Use Spectra.fix_monotonicity().')

        # gather ys
        y, ys = self._gather_ys(limits=limits)
        im = Image(data=ys)
        im._copy_attrs_from(self)

        if im.data is not None:
            im.x_centers = x_centers
            im.y_centers = y

        return im
    
    def stack_spectra_as_rows(self, y_centers=None, limits=None):
        """Return image representation of spectra, where spectra are layed out in the horizontal direction.

        Note:
            All spectra must have the same x-coordinates. This is verified.

        Warning:
            attrs are copied to the final image, but attrs from each spectrum is lost.

        Args:
            y_centers (list, optional): numerical value for each spectrum. 
                 indicating a numerical value for each spectrum. Default is None, i.e., each
                spectrum will be numbered from 1 to len(ss).
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Image`.
        """
        # check centers
        if y_centers is not None:
            assert len(y_centers) == len(self), f'centers must have the same number of items as the number of spectra.\nnumber of centers: {len(centers)}\nnumber of spectra: {len(self)}'

        # check if array is monotonic
        # This is necessary because the way images are plotted
        # The 'lower' value of the x-axis will be to the left
        if self.check_monotonicity() != 'increasing':
            raise ValueError('monotonicity of data must be strictly increasing. Use Spectra.fix_monotonicity().')

        # gather ys
        y, ys = self._gather_ys(limits=limits)
        if len(y) == 0:
            im = Image()
            im._copy_attrs_from(self)
            return im

        ys = ys.transpose()
        im = Image(data=ys)
        im._copy_attrs_from(self)
        im.x_centers = y
        im.y_centers = y_centers

        return im
    

    def calculate_shift(self, mode='cc', limits=None, **kwargs):
        """Returns shift list so all spectra is aligned to the first spectrum.

        Args:
            mode (string, optional): method used. Options are: 

                 1) 'cc': align spectra via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous spectrum.

                 3) 'max': Align the max point of every spectrum. 
                 
                 4) 'peak': Fit one peak in each spectrum and align them 
                 (requires that `brixs.addons.fitting` is imported)
                 
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`
            
        Returns:
            list

        See Also:
            :py:func:`Spectra.set_shift`
        """
        ##################
        # Initialization #
        ##################
        values = np.array([0.0]*len(self))

        ####################
        # cross-corelation #
        ####################
        if mode == 'cc':
            ########
            # crop #
            ########
            ss = self._copy(limits=limits)

            ######################
            # x must be the same #
            ######################
            _x = ss.check_same_x()

            #########################################################
            # x must be uniform (same step between each data point) #
            #########################################################
            step = ss.check_step()

            values = np.array(ss.calculate_roll(mode='cc'))
            values = list(values*step)
        ###############################
        # sequential cross-corelation #
        ###############################
        elif mode == 'seq':
            ########
            # crop #
            ########
            ss = self._copy(limits=limits)

            ######################
            # x must be the same #
            ######################
            _x = ss.check_same_x()

            #########################################################
            # x must be uniform (same step between each data point) #
            #########################################################
            step = ss.check_step()

            values = np.array(ss.calculate_roll(mode='seq'))
            values = list(values*step)
        #######
        # max #
        #######
        elif mode == 'max':
            ss  = self._copy(limits=limits)
            ref = ss[0].x[np.argmax(ss[0].y)]
            for i, s in enumerate(ss):
                values[i] = -(s.x[np.argmax(s.y)] - ref)
        ########
        # peak #
        ########
        elif mode == 'peak':
            # check if fitting was imported
            if hasattr(self, 'fit_peak') == False and callable(self.fit_peak) == False:
                raise ValueError('cannot calculate shifts via `peak` because fitting functions are not imported\nPlease import fitting function via `import brixs.addons.fitting`')
            _result = self.fit_peak(limits=limits, **kwargs)
            popt = _result['popt']
            values = np.array([_[1] for _ in popt])
            values = -values + values[0]
        else:
            raise ValueError(f'mode=`{mode}` not valid. Valid modes: `cc`, `max`, `peak`')
        return values

    def calculate_roll(self, mode='cc', limits=None, **kwargs):
        """return roll values so all spectra is aligned to the first spectrum.

        Args:
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align spectra via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous spectrum.

                 3) 'max': Align the max point of every spectrum. 
                 
                 4) 'peak': Fit one peak in each spectrum and align them 
                 (requires that `brixs.addons.fitting` is imported)

            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'` 
            
        Returns:
            list

        See Also:
            :py:func:`Spectra.set_roll`
        """
        ##################
        # Initialization #
        ##################
        values = np.array([0.0]*len(self))

        ########
        # crop #
        ########
        ss = self._copy(limits=limits)

        ######################
        # x must be the same #
        ######################
        _x = ss.check_same_x()

        #########################################################
        # x must be uniform (same step between each data point) #
        #########################################################
        step = ss.check_step()

        ####################
        # cross-corelation #
        ####################
        if mode == 'cc':
            for i, s in enumerate(ss):
                cc        = np.correlate(ss[0].y, s.y, mode='full')
                values[i] = np.argmax(cc)
            values = values - values[0]
        ###############################
        # sequential cross-corelation #
        ###############################
        elif mode == 'seq':
            for i in range(1, len(ss)):
                cc        = np.correlate(ss[i-1].y, ss[i].y, mode='full')
                values[i] = np.argmax(cc) - (len(ss[i-1].y) - 1) + values[i-1]
            values = values - values[0]
        #######
        # max #
        #######
        elif mode == 'max':
            values = np.array(ss.calculate_shift(mode='max', **kwargs))
            values = list(int(round(values/step)))
        #########
        # peaks #
        #########
        elif mode == 'peak':
            values = np.array(ss.calculate_shift(mode='peak', **kwargs))
            values = list(int(round(values/step)))
        else:
            raise ValueError(f'mode={mode} not valid. Valid modes: `cc`, `max`, `peak`')
        return values

    def calculate_factor(self, mode='max', limits=None, **kwargs):
        """return mult. factor so spectra have the same height as the first spectrum

        Args:
            mode (string, optional): method used. Options are: 
                 
                 1) 'max': spectra have max same as first spectrum
                 
                 2) 'delta': spectra have y variation `[max(y)-min(y)]` is same
                 as first spectrum 

                 3) 'area': spectra have area the same as first spectrum 
                 
                 4) 'peak': Fit one peak in each spectrum and make max of that 
                  peak same a the peak in the first spectrum (requires that `brixs.addons.fitting` is imported)

            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'` 
            
        Returns:
            list

        See Also:
            :py:func:`Spectra.set_factor`
        """
        ##################
        # Initialization #
        ##################
        values = np.array([0.0]*len(self))
        
        #######
        # max #
        #######
        if mode == 'max':
            ss = self._copy(limits=limits)
            ref = max(ss[0].y)
            for i, s in enumerate(ss):
                values[i] = ref/max(s.y)
        #########
        # delta #
        #########
        elif mode == 'delta':
            ss  = self._copy(limits=limits)
            ref = max(ss[0].y) - min(ss[0].y)
            for i, s in enumerate(ss):
                values[i] = ref/(max(s.y) - min(s.y))
        ########
        # area #
        ########
        elif mode == 'area':
            ss  = self._copy(limits=limits)
            ref = ss[0].calculate_area()
            for i, s in enumerate(ss):
                values[i] = ref/s.calculate_area()
        #########                
        # peaks #
        #########                
        elif mode == 'peak':
            # check if fitting was imported
            if hasattr(self, 'fit_peak') == False and callable(self.fit_peak) == False:
                raise ValueError('cannot calculate shifts via `peaks` because fitting functions are not imported\nPlease import fitting function via `import brixs.addons.fitting`')
            _result = self.fit_peak(limits=limits, **kwargs)
            popt = _result['popt']
            values = np.array([_[0] for _ in popt])
            values = -values + values[0]
        else:
            raise ValueError(f'mode={mode} not valid. Valid modes: `max`, `delta`, `area`, `peak`')
        return values

    def calculate_offset(self, mode='average', limits=None):
        """return additive factor so spectra have the same average y value as the first spectrum

        Args:
            mode (string, optional): method used. Options are: 
                 
                 1) 'average': spectra have y average same as first spectrum
                 
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.

        Returns:
            list

        See Also:
            :py:func:`Spectra.set_offset`, :py:func:`Spectra.floor()`,
        """
        ##################
        # Initialization #
        ##################
        values = np.array([0.0]*len(self))

        ###########
        # average #
        ###########
        if mode == 'average':
            ss  = self._copy(limits=limits)
            ref = np.mean(ss[0].y)
            for i, s in enumerate(ss):
                values[i] = ref - np.mean(s.y)
        else:
            raise ValueError(f'mode={mode} not valid. Valid modes: `average`')
        return values

    def calculate_calib(self, values, mode='cc', deg=1, limits=None, **kwargs):
        """return a calibration factor via :py:func:`Spectra.calculate_shift()`.

        Note:
            The calibration factor is the shift value (x-coord) as a function of
             the `values` array. Assuming each Spectrum was collected based on a certain experimental 
            condition (i.e, with a photon energy E), than one can find a calibration
            factor to multiply the x axis so that the x-coordinates reflect this
            experimental condition.

        Args:
            values (list): value list where each element represents a for each spectrum.
            mode (string, optional): method used. Options are: 

                 1) 'cc': align spectra via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous spectrum.

                 3) 'max': Align the max point of every spectrum. 
                 
                 4) 'peak': Fit one peak in each spectrum and align them 
                 (requires that `brixs.addons.fitting` is imported)
            
            deg (int, optional): a polynomial degree order used to fit the curve
                `shift vs value`.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            :py:class:`Spectrum` -> s

            s.x     = values
            s.y     = calculated shifts
            s.fit   = speectrum type (polynomial curve)
            s.popt  = polynomial coeff. (highest power first.) that fit the calibration curve.
            s.model = function f(x) of the calibration curve.
            s.R2    = R2 factor of the fitting.
        """
        # check number of values matches the number of spectra
        if len(self) != len(values):
            raise ValueError(f'number of values ({len(values)}) do not match the number of spectra ({len(self)})')

        # CALCULATION
        x = self.calculate_shift(mode=mode, limits=limits, **kwargs)

        # return
        final  = Spectrum(x=-np.array(x), y=values)
        result = final.polyfit(deg=deg)
        final.fit   = result['fit']
        final.popt  = result['popt']
        final.model = result['model']
        final.R2    = result['R2']
        return final

    def calculate_area(self, limits=None):
        """Returns a list of the calculated area under the curve for each spectrum. Wrapper for `numpy.trapz()`_.

        Usage:
            >>> s.calculate_area()                # returns the area for the whole dataset
            >>> s.calculate_area(limits=(0, 10))  # returns the area between x=0 and 10
            
        Warning:
            Calculate_area for broken limits, i.e., `ss.calculate_area(limits=((0, 10), (90, 100)))`
             works, but will 
            typically not return a desirable value, because the ``trapz()``
            algorithm will consider the area between 10 and 90 as a rectangle.
            The correct approach in this case would be to calculated the area
            between 0 and 10 and between 90 and 100 separately.

            >>> areas = ss.calculate_area(0, 10) + ss.calculate_area(90, 100)

            this behavior might change in the future.

        Warning:
            numpy trapz does not require regular spacing between data points,
            however, a different number o points between spectra can lead to 
            area values that cannot be compared between spectra reliabily.
            To avoid this, this function will raise a error if x axis is not
            the same for all spectra.

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
        
        Returns:
            list

        .. _numpy.trapz(): https://numpy.org/doc/stable/reference/generated/numpy.trapz.html
        """
        _ = self.check_same_x(max_error=0.1)
        return [s.calculate_area(limits=limits) for s in self]

    def calculate_y_sum(self, limits=None):
        """Returns a list of the sum of y elements within a range for each spectra.
        
        Usage:
            >>> s.calculate_y_sum()  # returns the y sum for the whole dataset
            >>> s.calculate_y_sum((0, 10), (90, 100))  # returns the y sum from data between x=0 and 10 and between x=90 and 100

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
        
        Returns:
            list
        """
        return [s.calculate_y_sum(limits=limits) for s in self]

    def calculate_x_sum(self, limits=None):
        """Returns a list of the sum of x elements within a range for each spectra.
        
        Usage:
            >>> s.calculate_x_sum()  # returns the x sum for the whole dataset
            >>> s.calculate_x_sum((0, 10), (90, 100))  # returns the x sum from data between x=0 and 10 and between x=90 and 100

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
        
        Returns:
            list
        """
        return [s.calculate_x_sum(limits=limits) for s in self]

    def calculate_y_average(self, limits=None):
        """returs a list of the average x value within range for each spectrum

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
        
        Warning:
            calculating the y average does not require regular spacing between data points,
            however, a different number o points between spectra can lead to 
            wrong average values that cannot be compared between spectra reliabily.
            To avoid this, this function will raise a error if x axis is not
            the same for all spectra.
         
        Returns:
            list
        """

        _ = self.check_same_x(max_error=0.1)
        return [s.calculate_y_average(limits=limits) for s in self]

    def polyfit(self, deg, limits=None):
        """Fit data recursively with a polynomial. Wrapper for `numpy.polyfit()`_.

        Args:
            deg (int): degree of the fitting polynomial.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
         
        Returns:
            dictionary {fit, popt, R2, f(x)}

                fit (spectrum): polynomial fit spectrum with 100x more intepolated points

                popt (np.array): 1D array of polynomial coefficients 
                    (including coefficients equal to zero) from highest degree to 
                    the constant term.

                R2 (number): R2 error

                model (function): funcion f(x_centers)

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        popt  = [0]*len(self)
        model = [0]*len(self)
        fit   = [0]*len(self)
        R2    = [0]*len(self)
        for i in range(len(self)):
            _result = self[i].polyfit(deg=deg, limits=limits)
            fit[i]   = _result['fit']
            popt[i]  = _result['popt']
            model[i] = _result['model']
            R2[i]    = _result['R2']
        return {'fit': fit, 'popt': popt, 'R2': R2, 'model': model}

    ############
    # composed #
    ############
    def align(self, mode='cc', **kwargs):
        """Uses calculate_shift and set_shift.

        Args:
            same as calculate_shift() or calculate_roll()

        Returns:
            :py:class:`Spectra`

        See Also:
            :py:func:`Spectra.calculate_shift`, :py:func:`Spectra.calculate_roll`
        """
        if mode in ['cross-correlation', 'cc']:
            value = self.calculate_roll(mode=mode, **kwargs)
            return self.set_roll(value=value)
        else:
            value = self.calculate_shift(mode=mode, **kwargs)
            return self.set_shift(value=value)

    ##########################        
    # plot and visualization #
    ##########################  
    def plot(self, ax=None, smooth=1, limits=None, switch_xy=False, vi=0, hi=0, pvi=0, phi=0, verbose=True, **kwargs):
        """Plot spectra. Wrapper for `matplotlib.pyplot.plot()`_.

        Note:
            If `label` is `None` and if spectrum inside spectra have attr 
            `label`, this attr will be used as label, e.g., 
            `plt.plot(s.x, s.y, label=s.label)`.  

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            smooth (int, optional): number of points to average data. Default is 1.
            label or labels (str, number, or list, optional): if str or number, this label will be 
                applied to every spectra. If list, it must have the same length 
                as the number of spectra. If None and if 
                spectrum `s` inside spectra `ss` have attr `label`, 
                this attr will be used as label, e.g., `plt.plot(s.x, s.y, label=s.label)`.
                Default is None. 
            color, colors, or c (str, number, or list, optional): if str or number, this color will be 
                applied to every spectra. If list, it must have the same length 
                as the number of spectra and each element must be a color (a 
                color can be a str or a 3 element (RGB) list.
                Default is None. 
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            switch_xy (bool, optional): Switch x and y axis.
            hi, vi (number, optional): horizontal and vertical increments for 
                cascading plots.
            phi, pvi (number, optional): percentage wise horizontal and vertical 
                increments for cascading plots (percentage of the y-range for 
                each spectrum).
            verbose (bool, optional): if True, prints warning if ploted data has
                nun-numeric values (NaN). Default is True.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            `Line2D`_ list

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        ###################
        # figure and axes #
        ###################
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW or len(plt.get_fignums()) == 0:
                figure()

        #########
        # label #
        #########
        if 'label' in kwargs or 'labels' in kwargs:
            if 'label' in kwargs:
                _label = kwargs.pop('label')
            else:
                _label = kwargs.pop('labels')

            if isinstance(_label, Iterable) == True and isinstance(_label, str) == False:
                assert len(_label) == len(self), f'label must be a number of a list with length compatible with the number of spectra.\nnumber of labels: {len(label)}\nnumber of spectra: {len(self)}'
                label = _label
            else:
                label = [_label]*len(self)
        else:
            label = [None]*len(self)

        #########
        # color #
        #########
        if 'color' in kwargs or 'colors' in kwargs  or 'c' in kwargs:
            if 'color' in kwargs:
                _color = kwargs.pop('color')
            elif 'colors' in kwargs:
                _color = kwargs.pop('colors')
            else:
                _color = kwargs.pop('c')
            
            if isinstance(_color, Iterable) == True and isinstance(_color, str) == False:
                assert len(_color) == len(self), f'`color` must be a number of a list with length compatible with the number of spectra.\nnumber of colors: {len(_color)}\nnumber of spectra: {len(self)}'
                colors = _color
            else:
                colors = [_color]*len(self)
        else:
            colors = [None]*len(self)


        #####################################
        # vertical and horizontal increment #
        #####################################
        vi  = [vi]*len(self)
        hi  = [hi]*len(self)

        #####################################################
        # percentage wise vertical and horizontal increment #
        #####################################################
        pvi = [max(s.y)*pvi/100 for i, s in enumerate(self)]
        phi = [max(s.x)*phi/100  for i, s in enumerate(self)]

        ##########
        # offset #
        ##########
        offset = [0]*len(self)
        for i in range(len(self)):
            offset[i] = offset[i] + (vi[i] * i) + (pvi[i] * i)
        #########
        # shift #
        #########
        shift = [0]*len(self)
        for i in range(len(self)):
            shift[i] = shift[i] + (hi[i] * i) + (phi[i] * i)
        
        ########
        # plot #
        ########
        temp = [0]*len(self)
        if abs(sum(shift)) > 0:
            if abs(sum(offset)) > 0:
                for i in range(len(self)):
                    temp[i] = self[i].set_shift(shift[i]).set_offset(offset[i]).plot(ax=ax, label=label[i], color=colors[i], smooth=smooth, switch_xy=switch_xy, limits=limits, verbose=verbose, **kwargs)
            else:
                for i in range(len(self)):
                    temp[i] = self[i].set_shift(shift[i]).plot(ax=ax, label=label[i], color=colors[i], smooth=smooth, switch_xy=switch_xy, limits=limits, verbose=verbose, **kwargs)
        elif abs(sum(offset)) > 0:
            for i in range(len(self)):
                temp[i] = self[i].set_offset(offset[i]).plot(ax=ax, label=label[i], color=colors[i], smooth=smooth, switch_xy=switch_xy, limits=limits, verbose=verbose, **kwargs)
        else:
            for i in range(len(self)):
                temp[i] = self[i].plot(ax=ax, label=label[i], color=colors[i], smooth=smooth, switch_xy=switch_xy, limits=limits, verbose=verbose, **kwargs)

        return temp

# %% =============================== Image =============================== %% #
class Image(_BrixsObject, metaclass=_Meta):
    """Returns a ``spectra`` object.

    Args:
        data (list or array, optional): list of :py:class:`spectrum` objects.
        x_centers, y_centers: (list, optional): pixel center labels. `x_centers` 
            from left to right and `y_centers` from top to bottom.
        **kwargs: kwargs are passed to :py:func:`Image.load` function.

    Usage:        
            >>> im = br.Image()
            >>> im = br.Image(data)
            >>> im = br.Image(data=data)
            >>> im = br.Image(data=data, x_centers=x_centers, y_centers=y_centers)
            >>> im = br.Image().loadtxt(filepath=<filepath>)
            >>> im = br.Image().loadtxt(filepath=<filepath>, delimiter=',')
            >>> im = br.Image().loadtxt(filepath=<filepath>)
            >>> im = br.Image().loadtxt(filepath=<filepath>, delimiter=',')
            >>>
            >>> print(im.get_core_attrs()) # print list of core attrs
            >>> print(im.get_attrs())      # print list of attrs
            >>> print(im.get_methods())    # print list of methods available

    Notes:
        in numpy arrays, the start of slices are inclusive but the stop are exclusive, 
         i.e. array[start(inclusive):stop(exclusive)]. In brixs Images, both start
         and stop are inclusives, i.e. im[start(inclusive):stop(inclusive)].
    
    Notes:
        Three methods for ploting images are defined:
            im.imshow(): pixels are squares. x and y axes are given in terms of pixels
            im.plot(): Pixels are squares. x and y axes are given in terms of x and y centers 
                (x_centers and y_centers must be monotonic).
            im.pcolormesh(): Allows for irregular pixel row/columns. x and y axes are 
                set based on x and y edges (or x and y centers when edges are not available).

    Attributes:
        Every BRIXS object has 5 types of attributes: 
        `Core` , `Check`, `Modifiers`, `Labels`, `User`.

        *1. Core*
            data (array): 2D array.
        
        *2. Check*
            x_step, y_step (number): None or a number if the step between two data points 
                at x_centers or y_centers. Can only be modified by 
                im.check_x_step() and im.check_y_step() methods.
            x_monotonicity, y_monotonicity (string): None if data is not monotonic or 'increasing'
                 or 'decreasing' if x_center or y_centers is monotonic. Can only
                 be modified by im.check_x_monotonicity() and im.check_y_monotonicity()
                 methods.

        *3. Modifiers*
            factor, offset (number): absolute values of 
            modifications made to the data. 
        
        *4. Labels*
        x_centers, y_centers (array): 1D arrays representing values associated 
            with each pixel row and column. Note that, defining x_edges and 
             y_edges changes x_centers and y_centers.
        x_edges, y_edges (array): None or monotonic 1D arrays with same lenght 
            of x_centers and y_centers plus 1 representing the edges of the 
            pixel rows and columns. These are connected with x_centers and 
            y_centers (i.e., changing x_edges also changes x_centers).

        *5. User*
            anything that the user defined on the fly.        
    """
    # read only and non-removable arguments
    _read_only     = ['x_step', 'x_monotonicity', 'y_step', 'y_monotonicity', 'has_nan']
    _non_removable = []
    
    def __init__(self, data=None, x_centers=None, y_centers=None):
        """Initialize the object instance"""
        ###########################
        # Initializing attributes #
        ###########################
        # core
        self._data = None

        # check
        self._x_step = None
        self._y_step = None
        self._x_monotonicity = None
        self._y_monotonicity = None
        self._has_nan = None

        # modifiers
        self._factor = 1
        self._offset = 0 

        # labels
        self._x_centers = None
        self._y_centers = None
        self._x_edges   = None
        self._y_edges   = None

        # extra
        for extra in settings._init['Image']:
            self.__setattr__(extra, settings._init['Image'][extra](self))

        ################
        # loading data #
        ################
        if data is not None:
            self.data = data
            if x_centers is not None:
                self.x_centers = x_centers
            if y_centers is not None:
                self.y_centers = y_centers
        return

    ###################
    # core attributes #
    ###################
    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        ###########
        # if None #
        ###########
        if value is None:
            self._data = None           
        ############
        # if array #
        ############
        elif hasattr(value, 'shape'):
            if len(value.shape) != 2:
                raise ValueError(f'data must be a 2d matrix. data={value}')
            if value.shape == (0, 0):
                self._data = None
            else:
                self._data = np.array(value, dtype='float')    
        ###########
        # if list #
        ###########
        elif isinstance(value, Iterable):
            ##############
            # empty list #
            ##############
            if len(value) == 0:
                self._data = None
            else:
                ####################
                # check list is 2d #
                ####################
                for row in value:
                    if isinstance(row, Iterable) == False:
                        raise ValueError(f'data must be a 2d matrix. data={value}')
                ####################################################
                # check if all lists inside list are the same size #
                ####################################################
                for i, row in enumerate(value):
                    if len(value[0]) != len(row):
                        raise ValueError(f'data must be a rectangular or square matrix. row={i} seems to have a different number of elemets. data={value}')
                ###########################
                # check if list are empty #
                ###########################
                if len(value[0]) == 0:
                    self._data = None
                else:
                    self._data = np.array(value, dtype='float')    
        else:
            raise ValueError(f'data must be a 2d matrix. data={value}')
        self._x_step         = None
        self._y_step         = None
        self._x_monotonicity = None
        self._y_monotonicity = None
        self._has_nan = None
        self._factor         = 1
        self._offset         = 0
        self.x_centers       = None
        self.y_centers       = None  
        return
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x_centers(self):
        return self._x_centers
    @x_centers.setter
    def x_centers(self, value):
        if value is None:
            if self.data is None:
                self._x_centers = None 
                self._x_edges   = None 
                self._x_step         = None
                self._x_monotonicity = None
                return
            else:
                value = np.arange(0, self.data.shape[1])
        elif isinstance(value, Iterable):
            assert len(value) == self.shape[1], f"number of x centers ({len(value)}) must be the same as the number of pixel columns ({self.shape[1]})"
        else:
            raise ValueError(f"x centers must be None or an iterable (list, tuple, or 1D array)")
        
        # setting centers
        self._x_centers = np.array(value, dtype='float')
        self._x_edges   = None

        # reseting checks
        self._x_step         = None
        self._x_monotonicity = None
        return
    @x_centers.deleter
    def x_centers(self):
        self._x_centers = np.arange(0, self.data.shape[1])

    @property
    def y_centers(self):
        return self._y_centers
    @y_centers.setter
    def y_centers(self, value):
        if value is None:
            if self.data is None:
                self._y_centers = None 
                self._y_edges   = None 
                self._y_step         = None
                self._y_monotonicity = None
                return
            else:
                value = np.arange(0, self.data.shape[0])
        elif isinstance(value, Iterable):
            assert len(value) == self.shape[0], f"number of y centers ({len(value)}) must be the same as the number of pixel columns ({self.shape[0]})"
        else:
            raise ValueError(f"y centers must be None or an iterable (list, tuple, or 1D array)")
        
        # setting centers
        self._y_centers = np.array(value, dtype='float')
        self._y_edges  = None
        
        # reseting checks
        self._y_step         = None
        self._y_monotonicity = None
        return
    @y_centers.deleter
    def y_centers(self):
        self._y_centers  = np.arange(0, self.data.shape[0])

    @property
    def x_edges(self):
        return self._x_edges
    @x_edges.setter
    def x_edges(self, value):
        if value is None:
            if self.data is None:
                self._x_centers = None 
                self._x_edges   = None 
                return
            else:
                centers = np.arange(0, self.data.shape[1])
                temp    = list(arraymanip.moving_average(centers, 2))
                value   = [centers[0] - temp[0]] + temp + [2*centers[-1] - temp[-1]]
        elif isinstance(value, Iterable):
            assert len(value) == self.shape[1] + 1, f"number of x edges ({len(value)}) must be the same as the number of pixel columns ({self.shape[1]}) plus one ({self.shape[1] + 1})"
            monotonicity = arraymanip.check_monotonicity(value)
            assert monotonicity != 0, f"edge values must be a monotonically array (either increasing or decreasing)"
        else:
            raise ValueError(f"x edges must be None or an iterable (list, tuple, or 1D array)")
        
        self._x_edges        = np.array(value, dtype='float')
        self._x_monotonicity = 'increasing' if monotonicity == 1 else 'decreasing'
        self._x_centers      = np.array(arraymanip.moving_average(value, 2), dtype='float')
    @x_edges.deleter
    def x_edges(self):
        raise NotImplementedError('this is not implemented yet')

    @property
    def y_edges(self):
        return self._y_edges
    @y_edges.setter
    def y_edges(self, value):
        if value is None:
            if self.data is None:
                self._y_centers = None 
                self._y_edges   = None 
                return
            else: 
                centers = np.arange(0, self.data.shape[0])
                temp    = list(arraymanip.moving_average(centers, 2))
                value   = [centers[0] - temp[0]] + temp + [2*centers[-1] - temp[-1]]
        elif isinstance(value, Iterable):
            assert len(value) == self.shape[0] + 1, f"number of y edges ({len(value)}) must be the same as the number of pixel rows ({self.shape[0]}) plus one ({self.shape[0] + 1})"
            monotonicity = arraymanip.check_monotonicity(value)
            assert monotonicity != 0, f"edge values must be a monotonically array (either increasing or decreasing)"
        else:
            raise ValueError(f"y edges must be None or an iterable (list, tuple, or 1D array)")
        self._y_edges = np.array(value, dtype='float')
        self._y_monotonicity = 'increasing' if monotonicity == 1 else 'decreasing'
        self._y_centers = np.array(arraymanip.moving_average(value, 2), dtype='float')
    @y_edges.deleter
    def y_edges(self):
        raise NotImplementedError('this is not implemented yet')

    ###################################
    # computed (read-only) attributes #
    ###################################
    @property
    def shape(self):
        if self.data is None:
            return (0, 0)
        return (self.data.shape[0], self.data.shape[1])
    @shape.setter
    def shape(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @shape.deleter
    def shape(self):
        raise AttributeError('Cannot delete object.')

    @property
    def histogram(self):
        if self.data is None:
            return None
        return self.calculate_histogram()
    @histogram.setter
    def histogram(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @histogram.deleter
    def histogram(self):
        raise AttributeError('Cannot delete object.')

    @property
    def columns(self):
        if self.shape[1] > 100:
            raise ValueError('cannot return image columns with more than 100 columns. Please, use im.get_columns()')
        return self.get_columns()
    @columns.setter
    def columns(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @columns.deleter
    def columns(self):
        raise AttributeError('Cannot delete object.')

    @property
    def rows(self):
        if self.shape[0] > 100:
            raise ValueError('cannot return image rows with more than 100 rows. Please, use im.get_rows()')
        return self.get_rows()
    @rows.setter
    def rows(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @rows.deleter
    def rows(self):
        raise AttributeError('Cannot delete object.')

    #########################
    # write-only attributes #
    #########################
    pass

    #######################
    # modifier attributes #
    #######################
    @property
    def offset(self):
        return self._offset
    @offset.setter
    def offset(self, value):
        _im = self.set_offset(-self._offset).set_offset(value)
        self._data   = _im.data
        self._offset = value
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        return self._factor
    @factor.setter
    def factor(self, value):
        _im = self.set_factor(1/self._factor).set_factor(value)
        self._data   = _im.data
        self._factor = value
        self._offset = _im.offset
    @factor.deleter
    def factor(self):
        raise AttributeError('Cannot delete object.')

    #################
    # magic methods #
    #################
    def __getitem__(self, item):
        if isinstance(item, tuple):
            # assert
            assert len(item) == 2, 'indexing must be lenght 2. Tuple of numbers (x_center, y_center) or slice (x_start:x_stop, y_start:y_stop)'
            
            # from centers to index
            if isinstance(item[0], slice) or isinstance(item[1], slice):
                if isinstance(item[0], slice):
                    y_start = item[0].start
                    y_stop  = item[0].stop
                else:
                    x_start = item[0]
                    x_stop = item[0]
                if isinstance(item[1], slice):
                    x_start = item[1].start
                    x_stop  = item[1].stop
                else:
                    x_start = item[1]
                    x_stop = item[1]
            
                # get slice
                return self.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)

            else:
                return self._data[arraymanip.index(self.y_centers, item[0]), arraymanip.index(self.x_centers, item[1])]
        else:
            raise TypeError('Index must be a tuple of numbers (x_center, y_center) or slice (x_start:x_stop, y_start:y_stop), not {}'.format(type(item).__name__))

    def __setitem__(self, item, value):
        if isinstance(item, tuple):
            # assert
            assert len(item) == 2, 'indexing must be lenght 2. Tuple of numbers (x_center, y_center) or slice (x_start:x_stop, y_start:y_stop)'
            
            # from centers to index
            if isinstance(item[0], slice) or isinstance(item[1], slice):
                if isinstance(item[0], slice):
                    y_start = item[0].start
                    y_stop  = item[0].stop
                else:
                    y_start = item[0]
                    y_stop  = item[0]
                if isinstance(item[1], slice):
                    x_start = item[1].start
                    x_stop  = item[1].stop
                else:
                    x_start = item[1]
                    x_stop  = item[1]
            
                # set slice
                # assert that centers are monotonic
                if self.x_monotonicity is None:
                    try:
                        self.check_x_monotonicity()
                    except ValueError:
                        raise ValueError('x_centers must be monotonic for croping to make sense. Fix x_centers or set it to None')
                if self.y_monotonicity is None:
                    try:
                        self.check_y_monotonicity()
                    except ValueError:
                        raise ValueError('y_centers must be monotonic for croping to make sense. Fix y_centers or set it to None')
                # convert start and stop from centers to pixel
                if x_stop < x_start:
                    _x_start2 = x_stop
                    _x_stop2  = x_start
                else:
                    _x_start2 = x_start
                    _x_stop2  = x_stop
                if self.x_monotonicity.startswith('inc'):
                    x_centers = self.x_centers
                    _x_start = bisect.bisect_left(x_centers, _x_start2)
                    _x_stop  = bisect.bisect_right(x_centers, _x_stop2)
                else:
                    x_centers = self.x_centers[::-1]
                    _x_stop  = len(x_centers) - bisect.bisect_left(x_centers, _x_start2)
                    _x_start = len(x_centers) - bisect.bisect_right(x_centers, _x_stop2)

                if y_stop < y_start:
                    _y_start2 = y_stop
                    _y_stop2  = y_start
                else:
                    _y_start2 = y_start
                    _y_stop2  = y_stop
                if self.y_monotonicity.startswith('inc'):
                    y_centers = self.y_centers
                    _y_start = bisect.bisect_left(y_centers, _y_start2)
                    _y_stop  = bisect.bisect_right(y_centers, _y_stop2)
                else:
                    y_centers = self.y_centers[::-1]
                    _y_stop  = len(y_centers) - bisect.bisect_left(y_centers, _y_start2)
                    _y_start = len(y_centers) - bisect.bisect_right(y_centers, _y_stop2)

                self._data[_y_start:_y_stop, _x_start:_x_stop] = value
                return 

            else:
                assert numanip.is_number(value), f'value must be a number, not {type(value)}'
                self._data[arraymanip.index(self.y_centers, item[0]), arraymanip.index(self.x_centers, item[1])] = value
                return 
        else:
            raise TypeError('Index must be a tuple of numbers (x_center, y_center) or slice (x_start:x_stop, y_start:y_stop), not {}'.format(type(item).__name__))

    def __setattr__(self, name, value):
        if name in settings._forbidden_words['Image']:
            raise AttributeError(f'`{name}` is a reserved word and cannot be set as an attribute')
        super().__setattr__(name, value)

    def __len__(self):
        if self._data is None:
            return 0
        else:
            return len(self._data)

    def __add__(self, object):
        if isinstance(object, Image):
            if self.shape == object.shape:
                if np.issubdtype(self.data.dtype, np.integer) and np.issubdtype(object.data.dtype, np.integer):
                    dtype = 'int'
                else:
                    dtype = 'float'
                final = Image(data=np.add(self.data, object.data, dtype=dtype))
            else:
                raise ValueError(f'Shape is different.\nShape 1: {self.shape}\nShape 2: {object.shape}')
        elif numanip.is_number(object):
            if np.issubdtype(self.data.dtype, np.integer) and isinstance(object, (np.floating, float))==False:
                dtype = 'int'
            else:
                dtype = 'float'
            final = Image(data=np.add(self.data, object, dtype=dtype))
        else:
            raise ValueError(f'Cannot operate type {type(object)} with Image')

        final._copy_attrs_from(self)
        final._x_centers = self.x_centers
        final._y_centers = self.y_centers
        final._x_edges   = self.x_edges
        final._y_edges   = self.y_edges
        return final
        
    def __sub__(self, object):
        if isinstance(object, Image):
            if self.shape == object.shape:
                if np.issubdtype(self.data.dtype, np.integer) and np.issubdtype(object.data.dtype, np.integer):
                    dtype = 'int'
                else:
                    dtype = 'float'
                final = Image(data=np.subtract(self.data, object.data, dtype=dtype))
            else:
                raise ValueError(f'Shape is different.\nShape 1: {self.shape}\nShape 2: {object.shape}')
        elif isinstance(object, (np.floating, float, int)):
            if np.issubdtype(self.data.dtype, np.integer) and isinstance(object, (np.floating, float))==False:
                dtype = 'int'
            else:
                dtype = 'float'
            final = Image(data=np.subtract(self.data, object, dtype=dtype))
        else:
            raise ValueError(f'Cannot operate type {type(object)} with Image')

        final._copy_attrs_from(self)
        final._x_centers = self.x_centers
        final._y_centers = self.y_centers
        final._x_edges   = self.x_edges
        final._y_edges   = self.y_edges
        return final
    
    def __mul__(self, object):
        if isinstance(object, Image):
            if self.shape == object.shape:
                if np.issubdtype(self.data.dtype, np.integer) and np.issubdtype(object.data.dtype, np.integer):
                    dtype = 'int'
                else:
                    dtype = 'float'
                final = Image(data=np.multiply(self.data, object.data, dtype=dtype))
            else:
                raise ValueError(f'Shape is different.\nShape 1: {self.shape}\nShape 2: {object.shape}')
        elif isinstance(object, (np.floating, float, int)):
            if np.issubdtype(self.data.dtype, np.integer) and isinstance(object, (np.floating, float))==False:
                dtype = 'int'
            else:
                dtype = 'float'
            final = Image(data=np.multiply(self.data, object, dtype=dtype))
        else:
            raise ValueError(f'Cannot operate type {type(object)} with Image')
        
        final._copy_attrs_from(self)
        final._x_centers = self.x_centers
        final._y_centers = self.y_centers
        final._x_edges   = self.x_edges
        final._y_edges   = self.y_edges
        return final
    
    def __div__(self, object):
        if isinstance(object, Image):
            if self.shape == object.shape:
                if 0 in object.data:
                    raise ZeroDivisionError(f'Image contain zeros. Cannot divide by zero.')
                else:
                    dtype = 'float'
                    final = Image(data=np.divide(self.data, object.data, dtype=dtype))
                    # final = Image(data = self.data / object.data)
            else:
                raise ValueError(f'Shape is different.\nShape 1: {self.shape}\nShape 2: {object.shape}')
        elif isinstance(object, (np.floating, float, int)):
            if object == 0:
                raise ZeroDivisionError(f'Cannot divide by zero.')
            else:
                dtype = 'float'
                final = Image(data=np.divide(self.data, object, dtype=dtype))
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Image')
        
        final._copy_attrs_from(self)
        final._x_centers = self.x_centers
        final._y_centers = self.y_centers
        final._x_edges   = self.x_edges
        final._y_edges   = self.y_edges
        return final
    
    def __truediv__(self, object):
        return self.__div__(object)

    #########
    # attrs #
    #########
    pass

    ###########
    # support #
    ###########
    def _check_mask(self, mask):
        """returns mask in the right format.

        Note:
            mask must be in terms of x_center and y_center.

        Args:
            mask (None or list): list with 4 elements `(x_start, x_sto, y_start, y_stop)`, a list 
                like `((xi_1, xf_1, yi_1, yf_1), (xi_2, xf_2, yi_2, yf_2), ...)` 
                in terms of x_centers and y_centers, or None. If None, 
                this function simply returns None. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If mask = [], i.e.,
                an empty list, it assumes `limits = (None, None, None, None)`.

        Returns:
            None or limits in the following format:
                ((xi_1, xf_1, yi_1, yf_1), (xi_2, xf_2, yi_2, yf_2), ...)              
        """
        #####################
        # if limits is None #
        #####################
        if mask is None:
            return None
        
        ##################################
        # assert that limits is Iterable #
        ##################################
        assert isinstance(mask, Iterable), f'`mask` must be an Iterable, not {type(mask)}'
        
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 

        ################################
        # get min and max range values #
        ################################
        xmin = min(self.x_centers)
        xmax = max(self.x_centers)
        ymin = min(self.y_centers)
        ymax = max(self.y_centers)

        ################
        # empty limits #
        ################
        if len(mask) == 0:
            return ((xmin, xmax, ymin, ymax),)

        ##############
        # fix format #
        ##############
        # 4 elements 
        if len(mask) == 4: # (xi, xf, yi, yf), or ((xi1, xf1, yi1, yf1), ..., (xi2, xf2, yi2, yf2))
            if isinstance(mask[0], Iterable) == False:
                if isinstance(mask[1], Iterable) == False and isinstance(mask[2], Iterable) == False and isinstance(mask[3], Iterable) == False:
                    mask = [mask, ]
                else:
                    raise ValueError(f'wrong format for mask={mask}')
        # 1, 2, 3 or more than 4 elements 
        final = []
        for m in mask:
            assert isinstance(m, Iterable), f'wrong format for mask={mask}'
            temp = [m[0], m[1], m[2], m[3]]
            if temp[0] == None: temp[0] = xmin
            if temp[1] == None: temp[1] = xmax
            if temp[2] == None: temp[2] = ymin
            if temp[3] == None: temp[3] = ymax
            final.append((temp[0], temp[1], temp[2], temp[3]))   
        return final

    def _calculated_vmin_vmax(self):
        """returns optimal vmin and vmax for visualization
        
        vmin is set the the max of the intensity distribution (max of histogram)
        
        vmax is set when intensity distribution (histogram) drops below 0.01% of the max.

        Returns:
            vmin, vmax
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 

        ##############
        # vmin, vmax #
        ##############
        hist = self.calculate_histogram()
        x    = hist.x[np.argmax(hist.y)+1:]
        y    = hist.y[np.argmax(hist.y)+1:]
        filtered = np.array([[i, j] for i, j in zip(x, y) if j > abs(max(y)*0.0001)])  # clean zeros
        try:
            vmin = hist.x[np.argmax(hist.y)]
            vmax = filtered[-1, 0]
        except IndexError:  # in case the max of y is too high
            vmin = min([min(x) for x in self.data])
            vmax = max([max(x) for x in self.data])

        # in case vmax ends up being higher then vmin
        if vmax <= vmin:
            vmax = self.min() + (self.max()-self.min())/2
            vmin = self.min()
            if vmax <= vmin:
                vmax = None
                vmin = None
        return vmin, vmax

    ################
    # core methods #
    ################
    def estimate_x_edges_from_centers(self):
        """Returns copy of image with x_edges defined from averaging x_centers

        Note:
            x_centers must be at least monotonic.

        Returns:
            Image with x_edges defined
        """
        ########################
        # check if empty image #
        ########################
        if self.x_centers is None:
            return self.copy()
        
        ##################
        # check validity #
        ##################
        if self.x_monotonicity is None:
            try:
                self.check_x_monotonicity()
            except ValueError:
                raise ValueError('x_centers must be monotonic for calculating x_edges')

        #################
        # setting edges #
        #################
        im = self.copy()
        centers = im.x_centers
        if len(centers) == 1:
            im._x_edges = [centers[0] - 0.5] + [centers[0] + 0.5]
        else:
            temp = list(arraymanip.moving_average(centers, 2))
            im._x_edges = [2*centers[0] - temp[0]] + temp + [2*centers[-1] - temp[-1]]
        
        return im

    def estimate_y_edges_from_centers(self):
        """Returns copy of image with y_edges defined from averaging y_centers

        Note:
            y_centers must be at least monotonic.

        Returns:
            Image with y_edges defined
        """
        ########################
        # check if empty image #
        ########################
        if self.y_centers is None:
            return self.copy()
        
        ##################
        # check validity #
        ##################
        if self.y_monotonicity is None:
            try:
                self.check_y_monotonicity()
            except ValueError:
                raise ValueError('y_centers must be monotonic for calculating y_edges')

        #################
        # setting edges #
        #################
        im = self.copy()
        centers = im.y_centers
        if len(centers) == 1:
            im._y_edges = [centers[0] - 0.5] + [centers[0] + 0.5]
        else:
            temp = list(arraymanip.moving_average(centers, 2))
            im._y_edges = [2*centers[0] - temp[0]] + temp + [2*centers[-1] - temp[-1]]
        
        return im

    def get_rows(self, max_number_of_rows=100):
        """Returns pixel rows in a spectra object.

        Args:
            max_number_of_rows (int, optional): raises error if number of pixel 
                rows is higher than `max_number_of_rows`.

        Returns:
            Spectra
        """
        ####################################
        # raise error if object is too big #
        ####################################
        if self.shape[0] > max_number_of_rows:
            raise ValueError(f'Image has {self.shape[0]} rows. Cannot return image rows with more than {max_number_of_rows} rows. Please, increase `max_number_of_rows`')
        
        ################
        # empty object #
        ################
        if self.data is None:
            return None

        ############
        # get data #
        ############
        ss = Spectra()
        for i in range(self.shape[0]):
            ss.append(Spectrum(x=self.x_centers, y=self.data[i, :]))
        ss._copy_attrs_from(self)
        return ss
    
    def get_columns(self, max_number_of_columns=100):
        """Returns pixel columns in a spectra object.

        Args:
            max_number_of_columns (int, optional): raises error if number of pixel 
                columns is higher than `max_number_of_columns`.

        Returns:
            Spectra
        """
        ####################################
        # raise error if object is too big #
        ####################################
        if self.shape[1] > max_number_of_columns:
            raise ValueError(f'Image has {self.shape[1]} columns. Cannot return image columns with more than {max_number_of_columns} columns. Please, increase `max_number_of_columns`')

        ################
        # empty object #
        ################
        if self.data is None:
            return None
        
        ############
        # get data #
        ############
        ss = Spectra()
        for i in range(self.shape[1]):
            ss.append(Spectrum(x=self.y_centers, y=self.data[:, i]))
        ss._copy_attrs_from(self)
        return ss
    
    def index2center(self, y, x):
        """Return y and x pixel index  in terms of y_centers and x_centers
        
        Args:
            y, x (list or array): image positions in terms of pixel index
        
        Returns:
            y, x in terms of y_centers and x_centers
        """
        x = self.x_centers[x]
        y = self.y_centers[y]
        return y, x

    def center2index(self, y, x, closest=True, roundup=False):
        """Return y_centers and x_centers value in terms of y and x pixel index
        
        Args:
            y, x (list or array): image positions in terms of y_centers and x_centers
            closest (book, optional): if True, returns the index of the element in 
                array which is closest to value.
            roundup (bool, optional): if closest=True, and value is exactly midway
                between 2 items in array x, rounup=True will return the index of 
                item in x with highest value. Default is False.
        
        Returns:
            y, x in terms of pixel index
        """
        x = arraymanip.index(self.x_centers, x, closest, roundup)
        y = arraymanip.index(self.y_centers, y, closest, roundup)
        return y, x

    ########
    # copy #
    ########
    def _copy(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Same as self.copy(), but attrs are NOT copied."""
        #############
        # copy data #
        #############
        data = copy.deepcopy(self.data)

        #####################
        # if limits is None #
        #####################
        if x_start is None and x_stop is None and y_start is None and y_stop is None:
            im = Image(data=data)
            im._x_step         = self.x_step
            im._y_step         = self.y_step
            im._x_monotonicity = self.x_monotonicity
            im._y_monotonicity = self.y_monotonicity
            im._factor         = self.factor
            im._offset         = self.offset
            im._x_centers      = copy.deepcopy(self.x_centers)
            im._y_centers      = copy.deepcopy(self.y_centers)
            im._x_edges        = copy.deepcopy(self.x_edges)
            im._y_edges        = copy.deepcopy(self.y_edges)
            return im
        
        ########################
        # check if empty image #
        ########################
        if data is None:
            im = Image(data=data)
            im._x_step         = self.x_step
            im._y_step         = self.y_step
            im._x_monotonicity = self.x_monotonicity
            im._y_monotonicity = self.y_monotonicity
            im._factor         = self.factor
            im._offset         = self.offset
            im._x_centers      = copy.deepcopy(self.x_centers)
            im._y_centers      = copy.deepcopy(self.y_centers)
            im._x_edges        = copy.deepcopy(self.x_edges)
            im._y_edges        = copy.deepcopy(self.y_edges)
            return im
        
        #################
        # check if None #
        #################
        if x_start is None: x_start = min(self.x_centers)
        if x_stop  is None: x_stop  = max(self.x_centers)
        if y_start is None: y_start = min(self.y_centers)
        if y_stop  is None: y_stop  = max(self.y_centers)

        ########################################
        # check if extract is really necessary #
        ########################################
        if x_start <= min(self.x_centers) and x_stop >= max(self.x_centers):
            if y_start <= min(self.y_centers) and y_stop >= max(self.y_centers):
                im = Image(data=data)
                im._x_step         = self.x_step
                im._y_step         = self.y_step
                im._x_monotonicity = self.x_monotonicity
                im._y_monotonicity = self.y_monotonicity
                im._factor         = self.factor
                im._offset         = self.offset
                im._x_centers      = copy.deepcopy(self.x_centers)
                im._y_centers      = copy.deepcopy(self.y_centers)
                im._x_edges        = copy.deepcopy(self.x_edges)
                im._y_edges        = copy.deepcopy(self.y_edges)
                return im
            
        ##################
        # validate input #
        ##################
        # assert x_stop > x_start, f'x_start must be smaller than x_stop.'
        # assert y_stop > y_start, f'y_start must be smaller than y_stop.'        
        
        ########
        # crop #
        ########
        # assert that centers are monotonic
        if self.x_monotonicity is None:
            try:
                self.check_x_monotonicity()
            except ValueError:
                raise ValueError('x_centers must be monotonic for croping to make sense. Fix x_centers or set it to None')
        if self.y_monotonicity is None:
            try:
                self.check_y_monotonicity()
            except ValueError:
                raise ValueError('y_centers must be monotonic for croping to make sense. Fix y_centers or set it to None')
        # convert start and stop from centers to pixel
        if x_stop < x_start:
            _x_start2 = x_stop
            _x_stop2  = x_start
        else:
            _x_start2 = x_start
            _x_stop2  = x_stop
        if self.x_monotonicity.startswith('inc'):
            x_centers = self.x_centers
            _x_start = bisect.bisect_left(x_centers, _x_start2)
            _x_stop  = bisect.bisect_right(x_centers, _x_stop2)
        else:
            x_centers = self.x_centers[::-1]
            _x_stop  = len(x_centers) - bisect.bisect_left(x_centers, _x_start2)
            _x_start = len(x_centers) - bisect.bisect_right(x_centers, _x_stop2)

        if y_stop < y_start:
            _y_start2 = y_stop
            _y_stop2  = y_start
        else:
            _y_start2 = y_start
            _y_stop2  = y_stop
        if self.y_monotonicity.startswith('inc'):
            y_centers = self.y_centers
            _y_start = bisect.bisect_left(y_centers, _y_start2)
            _y_stop  = bisect.bisect_right(y_centers, _y_stop2)
        else:
            y_centers = self.y_centers[::-1]
            _y_stop  = len(y_centers) - bisect.bisect_left(y_centers, _y_start2)
            _y_start = len(y_centers) - bisect.bisect_right(y_centers, _y_stop2)

        im = Image(data=data[_y_start:_y_stop, _x_start:_x_stop])
        im.x_centers = self.x_centers[_x_start:_x_stop]
        im.y_centers = self.y_centers[_y_start:_y_stop]
 
        # edges y
        if self.y_edges is not None:
            edges = np.array(self.y_edges)
            if self.y_monotonicity is None:
                try:
                    self.check_y_monotonicity()
                except ValueError:
                    raise ValueError('y edges must be monotonic. Please, set suitable y edges (or set y_edges to None).')  # this error might never be used because if y_edges are defined, then y_edges are monotonic by definition
            if self.y_monotonicity.startswith('inc'):
                im.y_edges = [_ for _ in edges if _ >= edges[edges < im.y_centers[0]].max() and _ <= edges[edges > im.y_centers[-1]].min()]
            elif self.y_monotonicity.startswith('dec'):
                im.y_edges = [_ for _ in edges if _ >= edges[edges < im.y_centers[-1]].max() and _ <= edges[edges > im.y_centers[0]].min()]
        # edges x
        if self.x_edges is not None:
            edges = np.array(self.x_edges)
            if self.x_monotonicity is None:
                try:
                    self.check_x_monotonicity()
                except ValueError:
                    raise ValueError('x edges must be monotonic. Please, set suitable y edges (or set y_edges to None).')
            if self.x_monotonicity.startswith('inc'):
                im.x_edges = [_ for _ in edges if _ >= edges[edges < im.x_centers[0]].max() and _ <= edges[edges > im.x_centers[-1]].min()]
            elif self.x_monotonicity.startswith('dec'):
                im.x_edges = [_ for _ in edges if _ >= edges[edges < im.x_centers[-1]].max() and _ <= edges[edges > im.x_centers[0]].min()]
            
        # attrs
        im._x_step         = self.x_step
        im._y_step         = self.y_step
        im._x_monotonicity = self.x_monotonicity
        im._y_monotonicity = self.y_monotonicity
        im._factor         = self.factor
        im._offset         = self.offset
        return im

    def copy(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Return a copy of the object.

        Usage:
            >>> # full copy
            >>> im2 = im1.copy()  # im2 is now a copy of im1
            >>>
            >>> # im3 will be a croped image of im1
            >>> im3 = im1.copy(x_start, x_stop, y_start, y_stop) 

        Args:
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.

        Returns:
            :py:attr:`Image`
        """
        im = self._copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
        im._copy_attrs_from(self)

        # extra
        for extra in settings._copy['Image']:
            if hasattr(im, extra):
                im.__setattr__(extra, self.__getattribute__(extra).copy(im))

        return im
    
    #################
    # save and load #
    #################
    def savenpy(self, filepath, check_overwrite=False, verbose=False, **kwargs):
        """Save image as numpy binary (.npy). Wrapper for Wrapper for `np.save()`_.

        Warning:
            metadata is not saved to file

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user wants to overwrite file.
            verbose (bool, optional): turn verbose on and off. Default is `False`.
            **kwargs: kwargs are passed to ``np.save()`` that saves the data.

        Returns:
            None

        See Also:
            :py:func:`Image.loadnpy`

        .. _np.save(): https://numpy.org/doc/stable/reference/generated/numpy.save.html
        """
        ##############################
        # check for wrong parameters #
        ##############################
        if 'folderpath' in kwargs:
            raise ValueError('folderpath is not a valid parameter.\nPlease, use filepath')
        
        ######################################
        # check if filepath points to a file #
        ######################################
        filepath = Path(filepath)
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        ###################
        # check overwrite #
        ###################
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if other.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        ########
        # save #
        ########
        np.save(Path(filepath).with_suffix('.npy'), self._data, **kwargs)

    def loadnpy(self, filepath, verbose=False, **kwargs):
        """Load data from a numpy binary file (.npy). Wrapper for `np.load()`_.

        Warning:
            metadata is not saved to file

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first 
                decompressed. Last used filepath is saved to im.filepath.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is loaded.
            verbose (book, optional): Default is False. If True, it will print
                an warning when attributes cannot be loaded from the file.
            **kwargs: kwargs are passed to ``np.load()`` that loads the data.

        Returns:
            br.image

        See Also:
            :py:func:`Image.savenpy`

        .. _np.loadnpy(): https://numpy.org/doc/stable/reference/generated/numpy.load.html
        """
        ############
        # filepath #
        ############
        filepath = Path(filepath)
        assert filepath.is_file(), f'filepath must point to a file, {filepath}'
        assert filepath.exists(),  f'filepath does not exist, {filepath}'

        ########
        # read #
        ########
        data = np.load(Path(filepath), **kwargs)

        ##########
        # assign #
        ##########
        return Image(data=data)

    def savetiff(self, filepath, check_overwrite=False, verbose=False, **kwargs):
        """Save image as tiff. Wrapper for Wrapper for `plt.imsave()`_.

        Warning:
            metadata is not saved to file

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user wants to overwrite file.
            verbose (bool, optional): turn verbose on and off. Default is `False`.
            **kwargs: kwargs are passed to ``plt.imsave()`` that saves the data.

        Returns:
            None

        See Also:
            :py:func:`Image.loadtiff`

        .. _plt.imsave(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.imsave.html
        """
        ##############################
        # check for wrong parameters #
        ##############################
        if 'folderpath' in kwargs:
            raise ValueError('folderpath is not a valid parameter.\nPlease, use filepath')
        
        ######################################
        # check if filepath points to a file #
        ######################################
        filepath = Path(filepath)
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        ###################
        # check overwrite #
        ###################
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if other.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        ########
        # save #
        ########
        plt.imsave(Path(filepath).with_suffix('.tiff'), self._data, **kwargs)

    def savetxt(self, filepath, only_data=False,  check_overwrite=False, verbose=False, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Warning:
            Attrs are saved as comments if only_data is False. Saving attrs to file
            is not always reliable because requires converting variables to string. 
            Only attrs that are of type: string, number, and list of number,
             list of list of number and strings have been tested. Dictionaries are not saved.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user wants to overwrite file.
            verbose (bool, optional): turn verbose on and off. Default is `False`.
            **kwargs: kwargs are passed to ``np.savetxt()`` that saves the data.

        If not specified, the following parameters are passed to `numpy.savetxt()`_:

        Args:
            fmt (string, or list, optional): A single format (like ``%10.5f``), or a
                sequence of formats. See numpy's documentation for more information.
                If not specified, best fmt is calculated based on the
                number of decimal places of the data. Specifying fmt makes the code
                runs a little faster (not much tough, but it might make a difference
                if saving a lot of files).
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (", ").
            newline (str, optional): String or character separating lines.
                Default is ``\n``.
            header (bool, optional): String that will be written at the beginning of the file.
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        See Also:
            :py:func:`Image.load`

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        ##############################
        # check for wrong parameters #
        ##############################
        if 'folderpath' in kwargs:
            raise ValueError('folderpath is not a valid parameter.\nPlease, use filepath')
        
        ######################################
        # check if filepath points to a file #
        ######################################
        filepath = Path(filepath)
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        ###################
        # check overwrite #
        ###################
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if other.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        ##########
        # kwargs #
        ##########
        if 'fmt' not in kwargs: # pick best format
            if self._data != [] and self._data is not None:
                decimal = max([numanip.n_decimal_places(x) for x in arraymanip.flatten(self._data)])
                kwargs['fmt'] = f'%.{decimal}f'
        kwargs.setdefault('delimiter', ', ')
        kwargs.setdefault('newline', '\n')
        kwargs.setdefault('comments', '# ')

        #####################
        # header and footer #
        #####################
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
        else:
            attrs_dict = {_:self.__getattribute__(_) for _ in settings._reserved_words['Image']['vars'] if _ not in ['_data', ]}
            attrs_dict.update(self.get_attrs_dict())
            header = '\n'.join(_attr2str(attrs_dict, verbose)) + '\n'

            if 'header' not in kwargs:
                kwargs['header'] = header
            else:
                if kwargs['header'] == '':
                    kwargs['header'] = header
                elif kwargs['header'][-1] != '\n':
                    kwargs['header'] += '\n'

        ########
        # save #
        ########
        np.savetxt(Path(filepath), self._data, **kwargs)
        return

    def loadtxt(self, filepath, only_data=False, verbose=False, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Warning:
            This a very simple loading function that works well with matrix text
            file.

        Note:
            If file was saved by br.Image.save(), then the metadata (comments) can be 
            recovered. If not, one might get better results by setting `only_data = True`.

            
        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first 
                decompressed. Last used filepath is saved to im.filepath.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is loaded.
            verbose (book, optional): Default is False. If True, it will print
                an warning when attributes cannot be loaded from the file.
            **kwargs: kwargs are passed to ``plt.genfromtxt()`` that loads the data.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``# ``. Attributes picked up
                from the header will be loaded too.

        Returns:
            br.Image

        See Also:
            :py:func:`Image.save`

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        ############
        # filepath #
        ############
        filepath = Path(filepath)
        assert filepath.is_file(), f'filepath must point to a file, {filepath}'
        assert filepath.exists(),  f'filepath does not exist, {filepath}'

        ##########
        # kwargs #
        ##########
        kwargs.setdefault('delimiter', ', ')
        kwargs.setdefault('comments', '# ')

        ########
        # read #
        ########
        data = np.genfromtxt(Path(filepath), **kwargs)

        ##########
        # assign #
        ##########
        im = Image(data=data)

        ###############
        # read header #
        ###############
        if only_data is False:
            # get header
            header = filemanip.load_comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            
            # remove comment flag (#)
            comment_flag_len = len(kwargs['comments'])
            for i, line in enumerate(header):
                header[i] = line[comment_flag_len:]

            # attrs dict
            attrs_dict = _str2attr(header[:-1], verbose=verbose)

            # set attrs
            for attr in attrs_dict:
                im.__setattr__(attr, attrs_dict[attr])
        return im
    
    #########
    # check #
    #########
    def check_nan(self):
        """Check if data have non-numeric (NaN) values
        
        Result (True or False) is stored in im.has_nan attribute

        Returns:
            None
        """
        ################
        # empty object #
        ################
        if self.data is None:
            self._has_nan = False
            return
        
        #############
        # check nan #
        #############
        if np.isnan(self.data).any():
            self._has_nan = True
        else:
            self._has_nan = False
        return
    
    def find_nan(self):
        """Return a list with positions where non-numeric (NaN) values were found
        
        Returns:
            list
        """
        has_nan = self.check_nan()

        if has_nan:
            return np.argwhere(np.isnan(self.data))
        return []
        
    def check_x_step(self, max_error=0.1):
        """Checks vector uniformity of the x centers.

            If the step between two data points is the same through out the
            x vector, it sets :py:attr:`Image.x_step` with the value of the step size.

            Args:
                max_error (number, optional): percentage (of the x step) value 
                of the max error. Default is 0.1 %
                
            Note:
                Step uniformity is verified by the following equation:

                    (max(steps) - min(steps))/np.mean(steps) * 100 < max_error

            Returns:
                None

            Raises:
                ValueError: If x centers are not uniform.
        
            See Also:
                :py:func:`Image.check_y_step`, :py:func:`Image.check_x_monotonicity`
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 

        ##############
        # check step #
        ##############
        s = Spectrum(x=self.x_centers)
        try:
            s.check_step()
        except ValueError:
            raise ValueError(f"Step in the x centers seems not to be uniform. Set im.x_centers = None or change im.x_centers")
        self._x_step = s.step

        ###############################
        # monotonicity comes for free #
        ###############################
        self._x_monotonicity = s._monotonicity
        return

    def check_y_step(self, max_error=0.1):
        """Checks vector uniformity of the y centers.

            If the step between two data points is the same through out the
            x vector, it sets :py:attr:`Image.y_step` with the value of the step size.

            Args:
                max_error (number, optional): percentage (of the x step) value 
                of the max error. Default is 0.1 %
                
            Note:
                Step uniformity is verified by the following equation:

                    (max(steps) - min(steps))/np.mean(steps) * 100 < max_error

            Returns:
                None

            Raises:
                ValueError: If y centers are not uniform.
        
            See Also:
                :py:func:`Image.check_x_step`, :py:func:`Image.check_y_monotonicity`
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 

        s = Spectrum(x=self.y_centers)
        try:
            s.check_step()
        except ValueError:
            raise ValueError(f"Step in the y centers seems not to be uniform. Set im.y_centers = None or change im.y_centers")
        self._y_step = s.step

        ###############################
        # monotonicity comes for free #
        ###############################
        self._y_monotonicity = s._monotonicity
        return
    
    def check_x_monotonicity(self):
        """Sets x monotonicity attribute to 'increasing' or 'decreasing'.

        Raises:
            ValueError if data is not monotonic.

        Returns:
            None

        See Also:
            :py:func:`Image.check_x_step`, :py:func:`Image.fix_x_monotonicity`
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 

        s = Spectrum(x=self.x_centers)
        try:
            s.check_monotonicity()
        except ValueError:
            raise ValueError('x centers is not monotonic. Use Image.fix_x_monotonicity()')
        self._x_monotonicity = s.monotonicity
        return
    
    def check_y_monotonicity(self):
        """Sets y monotonicity attribute to 'increasing' or 'decreasing'.

        Raises:
            ValueError if data is not monotonic.

        Returns:
            None

        See Also:
            :py:func:`Image.check_y_step`, :py:func:`Image.fix_y_monotonicity`
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 

        s = Spectrum(x=self.y_centers)
        try:
            s.check_monotonicity()
        except ValueError:
            raise ValueError('y centers is not monotonic. Use Image.fix_y_monotonicity()')
        self._y_monotonicity = s.monotonicity
        return

    def fix_x_monotonicity(self, mode='increasing', attrs2reorder=None):
        """Rearrange image columns such as x centers is monotonically increasing or decreasing.

            Args:
                mode (str, optional): increasing or decreasing.
                attrs2reorder (list of str, optional): list of aditional Image attrs that 
                    must also be sorted based on the x centers.
            
            Returns:
                :py:clas:`Image`
            
            See Also:
                :py:func:`Image.check_x_monotonicity`
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 

        # check mode
        increasing = ['inc', 'i', 'up', 'increasing', 'increasingly']
        decreasing = ['dec', 'd', 'down', 'decreasing', 'decreasingly']
        if mode not in increasing and mode not in decreasing:
            raise ValueError('mode should be "decreasing" or "increasing".')
        
        if mode in increasing: mode = 'increasing'
        else: mode = 'decreasing'

        # turn array into monotonic
        if self.x_monotonicity is None:
            try:
                self.check_x_monotonicity()
                if self.x_monotonicity == mode:
                    return self.copy()
            except ValueError:
                pass
        elif self.x_monotonicity == mode:
            return self.copy()
        
        _decreasing = False
        if mode in decreasing: 
            _decreasing = True
        _im = self.copy()
        _ss = _im.columns
        _ss.__temporary__attr__123__ = self.x_centers
        if arraymanip.has_duplicates(_ss.__temporary__attr__123__):
            _ss = _ss.merge_duplicates(ref='__temporary__attr__123__', attrs2merge=attrs2reorder)
        ss = _ss.reorder_by_attr(attr='__temporary__attr__123__', attrs2reorder=attrs2reorder, decreasing=_decreasing)
        im = ss.stack_spectra_as_columns(x_centers=ss.__temporary__attr__123__)
        im.check_x_monotonicity()

        return im
    
    def fix_y_monotonicity(self, mode='increasing', attrs2reorder=None):
        """Rearrange image columns such as x centers is monotonically increasing or decreasing.

            Args:
                mode (str, optional): increasing or decreasing.
                attrs2reorder (list of str, optional): list of aditional Image attrs that 
                must also be sorted based on the y centers.
            

            Returns:
                :py:clas:`Image`
            
            See Also:
                :py:func:`Image.check_x_monotonicity`
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 

        # check mode
        increasing = ['inc', 'i', 'up', 'increasing', 'increasingly']
        decreasing = ['dec', 'd', 'down', 'decreasing', 'decreasingly']
        if mode not in increasing and mode not in decreasing:
            raise ValueError('mode should be "decreasing" or "increasing".')
        
        if mode in increasing: mode = 'increasing'
        else: mode = 'decreasing'

        # turn array into monotonic
        if self.y_monotonicity is None:
            try:
                self.check_y_monotonicity()
                if self.y_monotonicity == mode:
                    return self.copy()
            except ValueError:
                pass
        elif self.y_monotonicity == mode:
            return self.copy()
        
        _decreasing = False
        if mode in decreasing: 
            _decreasing = True
        _im = self.copy()
        _ss = _im.rows
        _ss.__temporary__attr__123__ = self.y_centers
        if arraymanip.has_duplicates(_ss.__temporary__attr__123__):
            _ss = _ss.merge_duplicates(ref='__temporary__attr__123__', attrs2merge=attrs2reorder)
        ss = _ss.reorder_by_attr(attr='__temporary__attr__123__', attrs2reorder=attrs2reorder, decreasing=_decreasing)
        im = ss.stack_spectra_as_rows(y_centers=ss.__temporary__attr__123__)
        im.check_y_monotonicity()

        return im

    #############
    # modifiers #
    #############
    def set_offset(self, value):
        """Add value to intensity (matrix elements).

        Args:
            value (value): offset value.

        Returns:
            :py:class:`Image`
        """
        im = self.copy()
        im._data   += float(value)
        im._offset += float(value)

        # extra
        for extra in settings._offset['Image']:
            if hasattr(im, extra):
                im.__getattribute__(extra).set_shift(value)

        return im

    def set_factor(self, value):
        """Multiply intensity (matrix elements) by value.

        Args:
            value (number): multiplicative factor

        Raises:
            AttributeError: if value is 0.

        Returns:
            :py:class:`Spectrum`
        """
        if value == 0:
            raise AttributeError('cannot set factor = 0.')
        im          = self.copy()
        im._data   *= float(value)
        im._factor *= float(value)
        im._offset *= float(value)

        # extra
        for extra in settings._factor['Image']:
            if hasattr(im, extra):
                im.__getattribute__(extra).set_shift(value)

        return im

    def set_horizontal_shift(self, value):
        """Roll pixels rows left and right in terms of x centers.

        Note:
            The shift value in terms of x center scale is converted to number of
             pixels to be rolled left or right. Elements that roll beyond the 
             last position are re-introduced at the other side.

        Args:
            value (number or list): shift value in terms of x center by which the data are 
                shifted. If list, then it must be of the same size as the number of
                 rows. First element will be assigned to the first row (top to 
                 bottom) and so on. 

        Returns:
            :py:class:`Image`
        
        See Also:
            :py:func:`Image.set_horizontal_roll`, :py:func:`Image.set_vertical_roll`, :py:func:`Image.set_vertical_shift`
        """
        ###################################
        # asserting validity of the input #
        ###################################
        # x axis must be uniform if mode is roll
        if self.x_step is None:
            try:
                self.check_x_step()
            except ValueError:
                raise ValueError(f'Cannot shift data, because x centers are not uniform. Use im.x_interp() to make data uniform, or shift data using im.set_roll()')
        
        #####################################
        # asserting validity of the input 2 #
        #####################################
        centers = self.y_centers
        if isinstance(value, Iterable) == False:
            value = [value]*len(centers)
        assert len(value) == len(centers), f'Number of values ({len(value)}) must be the same as the number of rows ({len(centers)})'
        
        #################
        # shift to roll #
        #################
        value = [int(round(k/self.x_step)) for k in value]
        
        ########
        # roll #
        ########
        return self.set_horizontal_roll(value=value)
    
    def set_vertical_shift(self, value):
        """Roll pixels columns up and down in terms of y centers.

        Note:
            The shift value in terms of y center scale is converted to number of
             pixels to be rolled left or right. Elements that roll beyond the 
             last position are re-introduced at the other side.

        Args:
            value (number or list): shift value in terms of y center by which the data are 
                shifted. If list, then it must be of the same size as the number of
                 rows. First element will be assigned to the first row (top to 
                 bottom) and so on. 

        Returns:
            :py:class:`Image`

        See Also:
            :py:func:`Image.set_horizontal_roll`, :py:func:`Image.set_vertical_roll`, :py:func:`Image.set_horizontal_shift`
        """
        ###################################
        # asserting validity of the input #
        ###################################
        # x axis must be uniform if mode is roll
        if self.y_step is None:
            try:
                self.check_y_step()
            except ValueError:
                raise ValueError(f'Cannot shift data, because x centers are not uniform. Use im.x_interp() to make data uniform, or shift data using im.set_roll()')
        
        #####################################
        # asserting validity of the input 2 #
        #####################################
        centers = self.x_centers
        if isinstance(value, Iterable) == False:
            value = [value]*len(centers)
        assert len(value) == len(centers), f'Number of values ({len(value)}) must be the same as the number of columns ({len(centers)})'
        
        #################
        # shift to roll #
        #################
        value = [int(round(k/self.y_step)) for k in value]
        
        ########
        # roll #
        ########
        return self.set_vertical_roll(value=value)

    def set_horizontal_roll(self, value):
        """Roll pixels rows left and right.

        Note:
            Elements that roll beyond the last position are re-introduced at the
             other side. 

        Warning:
            Roll values must be an integer. Float values will be rounded to int.

        Args:
            value (int or list): The number of pixels by which the data are 
                shifted. If list, then it must be of the same size as the number of
                 rows. First element will be assigned to the first row (top to 
                 bottom) and so on. If elements are not int, it will rounded to 
                 an integer value. 

        Returns:
            :py:class:`Image`
        
        See Also:
            :py:func:`Image.set_vertical_shift`, :py:func:`Image.set_vertical_roll`, :py:func:`Image.set_horizontal_shift`
        """
        ###################################
        # asserting validity of the input #
        ###################################
        centers = self.y_centers
        if isinstance(value, Iterable) == False:
            value = [value]*len(centers)
        assert len(value) == len(centers), f'Number of values ({len(value)}) must be the same as the number of rows ({len(centers)})'
        
        #####################
        # everything to int #
        #####################
        value = np.array([int(round(k)) for k in value])
        
        ########
        # roll #
        ########
        # I copyied this part from stackoverflow. 
        # Honestly, I don't 100% understand it, but it runs two order of magnetude
        # faster than the old implementation
        im = self.copy()   
        rows, column_indices = np.ogrid[:im.data.shape[0], :im.data.shape[1]]
        value[value < 0] += im.data.shape[1]
        column_indices = column_indices - value[:, np.newaxis]
        try:
            im._data = im.data[rows, column_indices]
        except IndexError:
            raise ValueError('Column rolling seems to be too large. Please, reduce the shift size.')

        ####################################
        # OLD IMPLEMENTATION (MUCH SLOWER) #
        ####################################
        # im = self.copy()
        # for i, v in enumerate(value):
        #     if v != 0:
        #         s = Spectrum(y=im._data[i, :]).set_roll(v)
        #         im._data[i, :] = s.y
        return im
    
    def set_vertical_roll(self, value):
        """Roll pixels columns up and down.

        Note:
            Elements that roll beyond the last position are re-introduced at the
             other side. 

        Warning:
            Roll values must be an integer. Float values will be rounded to int.

        Args:
            value (int or list): The number of pixels by which the data are 
                shifted. If list, then it must be of the same size as the number of
                 columns. First element will be assigned to the first column (left to 
                 right) and so on. If elements are not int, it will rounded to 
                 an integer value. 

        Returns:
            :py:class:`Image`
        
        See Also:
            :py:func:`Image.set_horizontal_roll`, :py:func:`Image.set_vertical_shift`, :py:func:`Image.set_horizontal_shift`
        """
        ###################################
        # asserting validity of the input #
        ###################################
        centers = self.x_centers
        if isinstance(value, Iterable) == False:
            value = [value]*len(centers)
        assert len(value) == len(centers), f'Number of values ({len(value)}) must be the same as the number of rows ({len(centers)})'
        
        #####################
        # everything to int #
        #####################
        value = np.array([int(round(k)) for k in value])
        
        ########
        # roll #
        ########
        # I copied this part from stackoverflow. 
        # Honestly, I don't 100% understand it, but it runs two order of magnitude
        # faster than the old implementation
        im = self.copy()    
        row_indices, cols = np.ogrid[:im.data.shape[0], :im.data.shape[1]]
        value[value < 0] += im.data.shape[0]
        row_indices = row_indices - value[np.newaxis, :]
        try:
            im._data = im.data[row_indices, cols]
        except IndexError:
            raise ValueError('Column rolling seems to be too large. Please, reduce the shift size.')

        ####################################
        # OLD IMPLEMENTATION (MUCH SLOWER) #
        ####################################
        # im = self.copy()    
        # for i, v in enumerate(value):
        #     if v != 0:
        #         s = Spectrum(y=im._data[:, i]).set_roll(v)
        #         im._data[:, i] = s.y
        return im

    ###############
    # modifiers 2 #
    ###############
    def set_horizontal_shift_via_polyval(self, p):
        """Set horizontal shift values to np.polyval(p, y_centers).

        Args:
            p (array): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term.

        Returns:
            :py:class:`Image`
        """
        f = lambda y: np.polyval(p, y)
        return self.set_horizontal_shift_via_function(f)
    
    def set_vertical_shift_via_polyval(self, p):
        """Set vertical shift values to np.polyval(p, x_centers).

        Args:
            p (array): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term.

        Returns:
            :py:class:`Image`
        """
        f = lambda x: np.polyval(p, x)
        return self.set_vertical_shift_via_function(f)

    def set_vertical_shift_via_function(self, f):
        """Set vertical shift values to f(x_center).

        Args:
            f (function): function where argument is x centers elements

        Returns:
            :py:class:`Image`
        """
        value = np.array([f(x) for x in self.x_centers])
        return self.set_vertical_shift(value=value)

    def set_horizontal_shift_via_function(self, f):
        """Set horizontal shift values to f(y_center).

        Args:
            f (function): function where argument is y centers elements

        Returns:
            :py:class:`Image`
        """
        value = np.array([f(x) for x in self.y_centers])
        return self.set_horizontal_shift(value=value)


    def floor(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Set intensity to zero inside limits to zero.

        Args:
            x_start, x_stop, y_start, y_stop (int): data range where average will
                be set to zero. Pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.

        Returns:
            Copy of image
        """
        temp  = self.copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
        value = temp.calculate_average()
        return self.set_offset(-value)

    ############
    # advanced #
    ############
    def crop(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Crop Image.

        Args:
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.

        Returns:
            :py:class:`Image`
        """
        return self.copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)

    def transpose(self):
        """Transpose image

        Returns:
            :py:class:`Image`
        """
        im = Image(data=self.data.transpose())

        ##################
        # transfer attrs #
        ##################
        im._copy_attrs_from(self)
        im._x_step         = self.y_step
        im._y_step         = self.x_step
        im._x_monotonicity = self.y_monotonicity
        im._y_monotonicity = self.x_monotonicity
        im._factor          = self.factor
        im._offset          = self.offset
        im._x_centers       = copy.deepcopy(self.y_centers)
        im._y_centers       = copy.deepcopy(self.x_centers)       
        im._x_edges         = copy.deepcopy(self.y_edges)
        im._y_edges         = copy.deepcopy(self.x_edges) 
        return im

    def x_interp(self, start=None, stop=None, num=None, step=None, x=None):
        """return image with interpolated x centers

        Args:
            start (number, optional): The starting value for x centers. If `None`,
                the minium attr value will be used.
            stop (number, optional): The end value for x centers. If `None`,
                the maximum attr value will be used.
            num (int, optional): Number of x centers values.
            step (number, optional): Spacing between x centers. This overwrites ``num``.
            x (list or array, optional): The values at which to
                evaluate the interpolated values for x centers. This overwrites all other arguments.

        Return:
            :py:class:`Image`
        """
        ###############
        # get columns #
        ###############
        cols = self.columns

        ##########
        # interp #
        ##########
        cols.ref = self.x_centers
        ss = cols.interp_spectra(ref='ref', start=start, stop=stop, num=num, step=step, x=x)
        im = ss.stack_spectra_as_columns()

        ##################
        # transfer attrs #
        ##################
        im._copy_attrs_from(self)
        im._x_step         = None
        im._y_step         = self.y_step
        im._x_monotonicity = None
        im._y_monotonicity = self.y_monotonicity
        im._factor          = self.factor
        im._offset          = self.offset
        im._x_centers       = ss.ref
        im._y_centers       = copy.deepcopy(self.y_centers)     
        im._x_edges         = None
        im._y_edges         = copy.deepcopy(self.y_edges)
        return im

    def y_interp(self, start=None, stop=None, num=None, step=None, y=None):
        """return image with interpolated y centers

        Args:
            start (number, optional): The starting value for y centers. If `None`,
                the minium attr value will be used.
            stop (number, optional): The end value for y centers. If `None`,
                the maximum attr value will be used.
            num (int, optional): Number of x centers values.
            step (number, optional): Spacing between y centers. This overwrites ``num``.
            y (list or array, optional): The values at which to
                evaluate the interpolated values for y centers. This overwrites all other arguments.

        Return:
            :py:class:`Image`
        """
        ############
        # get rows #
        ############
        rows = self.rows

        ##########
        # interp #
        ##########      
        rows.ref = self.y_centers
        ss = rows.interp_spectra(ref='ref', start=start, stop=stop, num=num, step=step, x=y)
        im = ss.stack_spectra_as_rows()
    
        ##################
        # transfer attrs #
        ##################
        im._copy_attrs_from(self)
        im._x_step         = self.x_step
        im._y_step         = None
        im._x_monotonicity = self.x_monotonicity
        im._y_monotonicity = None
        im._factor          = self.factor
        im._offset          = self.offset
        im._x_centers       = copy.deepcopy(self.x_centers)
        im._y_centers       = ss.ref
        im._x_edges         = copy.deepcopy(self.x_edges)
        im._y_edges         = None
        return im

    def binning(self, ncols=None, nrows=None):
        """Compute the 2D histogram of the data (binning of the data).

        Args:
            ncols, nrows (int or None, optional): number of columns and rows of the returned Image.
                If None, no binning is applied. Default is None. If the number of pixels in the image
                cannot be divided by the selected number of bins, it will raise 
                an error. 

        Returns:
            :py:class:`Image` binned image
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if ncols is None:
            ncols = self.shape[1]
        if nrows is None:
            nrows = self.shape[0]
        if numanip.is_integer(ncols) == False or numanip.is_integer(nrows)  == False or ncols < 0 or nrows < 0:
            raise ValueError("Number of bins must be a positive integer.")

        # is divisible
        assert self.shape[1] % ncols == 0, f"The {self.shape[1]} pixels in a row is not evenly divisible by {ncols}\nPlease, pick one of the following numbers: {np.sort(list(numanip.factors(self.shape[1])))}"
        assert self.shape[0] % nrows == 0, f"The {self.shape[0]} pixels in a column is not evenly divisible by {nrows}\nPlease, pick one of the following numbers: {np.sort(list(numanip.factors(self.shape[0])))}"

        # is uniform
        if self.x_step is None:
            try:
                self.check_x_step()
            except ValueError:
                raise ValueError('x_centers not uniform. Binning only makes sense for uniform images. Fix x_centers or set it to None')
        if self.y_step is None:
            try:
                self.check_y_step()
            except ValueError:
                raise ValueError('y_centers not uniform. Binning only makes sense for uniform images. Fix y_centers or set it to None')

        ###############
        # Calculation #
        ###############
        _bins_size = np.array((self.shape[0]/nrows, self.shape[1]/ncols))
        reduced    = Image(np.add.reduceat(np.add.reduceat(self._data, list(map(float, np.arange(0, self.shape[0], _bins_size[0]))), axis=0), list(map(float, np.arange(0, self.shape[1], _bins_size[1]))), axis=1))
        
        # x and y centers
        reduced._x_centers = Spectrum(x=self.x_centers).smooth(int(self.shape[1]/ncols), force_divisible=True).x
        reduced._y_centers = Spectrum(x=self.y_centers).smooth(int(self.shape[0]/nrows), force_divisible=True).x
        # _x_edges = np.arange(0, self.shape[1]+_bins_size[1]/2, _bins_size[1])
        # _y_edges = np.arange(0, self.shape[0]+_bins_size[0]/2, _bins_size[0])
        # reduced._x_centers = arraymanip.moving_average(_x_edges, n=2)
        # reduced._y_centers = arraymanip.moving_average(_y_edges, n=2)
        # reduced._x_edges = _x_edges
        # reduced._y_edges = _y_edges
        reduced = reduced.estimate_x_edges_from_centers()
        reduced = reduced.estimate_y_edges_from_centers()

        ##################
        # transfer attrs #
        ##################
        reduced._copy_attrs_from(self)
        reduced._x_step         = None
        reduced._y_step         = None
        reduced._x_monotonicity = self.x_monotonicity
        reduced._y_monotonicity = self.y_monotonicity
        reduced._factor          = self.factor
        reduced._offset          = self.offset
        # reduced._x_centers       = None
        # reduced._y_centers       = None
        # reduced._x_edges         = None
        # reduced._y_edges         = None
        return reduced
    
    def moving_average_x(self, n):
        """Returns an Image object with moving average on rows (x direction).

        Note:
            moving average is also applied to x_centers.

        Args:
            n (int): number of points to average.

        Returns:
            :py:class:`Image` with number of columns given by (number_of_columns - n + 1)
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 
        
        ###########
        # has nan #
        ###########
        if self.has_nan is None:
            self.check_nan()
        if self.has_nan:
            raise ValueError('Image has non-numeric values (NaN). Please, use im.find_nan() and remove nan values')
        
        ##################
        # moving average #
        ##################
        _data = sliding_window_view(self.data, window_shape=(1, n)).mean(axis=2).mean(axis=2)
        final = Image(data=_data)

        ##################
        # transfer attrs #
        ##################
        final._copy_attrs_from(self)
        final._x_step         = None
        final._y_step         = self.y_step
        final._x_monotonicity = self.x_monotonicity
        final._y_monotonicity = self.y_monotonicity
        final._factor         = self.factor
        final._offset         = self.offset
        final._x_centers      = Spectrum(y=self.x_centers).moving_average(n).y
        final._y_centers      = copy.deepcopy(self.y_centers)
        final._x_edges        = None
        final._y_edges        = copy.deepcopy(self.y_edges)
        return final

    def moving_average_y(self, n):
        """Returns an Image object with moving average on columns (y direction)

        Note:
            moving average is also applied to y_centers.

        Args:
            n (int): number of points to average.

        Returns:
            :py:class:`Image` with number of rows given by (number_of_rows - n + 1)
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 
        
        ###########
        # has nan #
        ###########
        if self.has_nan is None:
            self.check_nan()
        if self.has_nan:
            raise ValueError('Image has non-numeric values (NaN). Please, use im.find_nan() and remove nan values')
        
        ##################
        # moving average #
        ##################
        _data = sliding_window_view(self.data, window_shape=(n, 1)).mean(axis=2).mean(axis=2)
        final = Image(data=_data)
                   
        ##################
        # transfer attrs #
        ##################
        final._copy_attrs_from(self)
        final._x_step         = self.x_step
        final._y_step         = None
        final._x_monotonicity = self.x_monotonicity
        final._y_monotonicity = self.y_monotonicity
        final._factor         = self.factor
        final._offset         = self.offset
        final._x_centers      = copy.deepcopy(self.x_centers)
        final._y_centers      = Spectrum(y=self.y_centers).moving_average(n).y
        final._x_edges        = copy.deepcopy(self.x_edges)
        final._y_edges        = None
        return final

    def moving_average(self, n):
        """Returns an Image object with 2d moving average.

        Note:
            moving average is also applied to x_centers and y_centers.

        Args:
            n (int): number of points to average.

        Returns:
            :py:class:`Image` with number of rows given by (number_of_rows - n + 1) 
            and columns given by (number_of_columns - n + 1)
        """
        ################
        # empty object #
        ################
        if self.data is None:
            raise ValueError('cannot operate on empty image') 
        
        ###########
        # has nan #
        ###########
        if self.has_nan is None:
            self.check_nan()
        if self.has_nan:
            raise ValueError('Image has non-numeric values (NaN). Please, use im.find_nan() and remove nan values')
        
        ##################
        # moving average #
        ##################
        _data = sliding_window_view(self.data, window_shape=(n, n)).mean(axis=2).mean(axis=2)
        final = Image(data=_data)

        ##################
        # transfer attrs #
        ##################
        final._copy_attrs_from(self)
        final._x_step         = None
        final._y_step         = None
        final._x_monotonicity = self.x_monotonicity
        final._y_monotonicity = self.y_monotonicity
        final._factor         = self.factor
        final._offset         = self.offset
        final._x_centers      = Spectrum(y=self.x_centers).moving_average(n).y
        final._y_centers      = Spectrum(y=self.y_centers).moving_average(n).y
        final._x_edges        = None
        final._y_edges        = None
        return final
    
    def rotate(self, direction='clockwise'):
        """returns 90 degrees rotated image

        Args:
            direction (str, optional): rotation direction. Options are 'clockwise'
            and 'counterclockwise'. Default is clockwise.

        Returns:
            rotated image
        """
        assert direction in ['clockwise', 'counterclockwise'], f"invalid direction `{direction}`. Valid options are 'clockwise', 'counterclockwise'"
        
        im = self.copy()
        if direction == 'clockwise':
            im._data = np.rot90(self.data, k=3)
            im.x_centers = self.y_centers[::-1]
            im.y_centers = self.x_centers
        else:
            im._data = np.rot90(self.data, k=1)
            im.x_centers = self.y_centers
            im.y_centers = self.x_centers[::-1]
        return im

    def flipx(self):
        """flip data (in relation to x_centers)

        Returns:
            flipped image
        """    
        im = self.copy()
        im._data = np.flip(self.data, axis=1)
        return im

    def flipy(self):
        """flip data (in relation to y_centers)

        Returns:
            flipped image
        """    
        im = self.copy()
        im._data = np.flip(self.data, axis=0)
        return im

    ########################
    # calculation and info #
    ########################
    def calculate_average(self):
        """returns the average intensity value"""
        return self.data.mean()
        # return np.mean([s.calculate_y_average() for s in self.columns])

    def calculate_sigma(self):
        """returns the standard deviation of intensity values"""
        return self.data.std()

    def max(self):
        """return max intensity value"""
        return np.max(self.data)
    
    def min(self):
        """return min intensity value"""
        return np.min(self.data)

    def mean(self):
        """returns the average intensity value"""
        return self.data.mean()

    def std(self):
        """returns the standard deviation of intensity values"""
        return self.data.std()

    def argmax(self, coordinates='centers'):
        """returns the y, x position of the pixel with max intensity
        
        Args:
            coordinates (str, optional): If `pixels`, (y, x) is given in pixel 
                coordinates. If `centers`, (y, x) is given in terms of x_centers
                and y_centers. Default is `centers`.

        Returns:
            y, x
        """
        assert coordinates in ['pixels', 'centers'], f'coordinates must be `pixels` or `centers`, not `{coordinates}`'
        candidates = [(y, x) for y, x in enumerate(self.data.argmax(axis=1))]
        y, x = candidates[np.argmax([self.data[y, x] for y, x in candidates])]
        if coordinates == 'centers':
            y, x = self.index2center(y=y, x=x)
        return y, x

    def multiply(self, im):
        """Element-wise multiplication of two images. Wrapper for np.multiply()
        
        Args:
            im (Image): image to be multiplied.

        Returns:
            multiplied image
        """
        im2 = self.copy()
        im2.data = np.multiply(self.data, im.data)
        return im2
    

    def calculate_histogram(self, nbins=None):
        """Compute the intensity histogram of the data. Wrapper for `numpy.histogram()`_.

        Args:
            nbins (int, optional): number of bins. If not specified, it will
                be set to 1000. If 1000 is too high (maximum value of the histogram
                is less than 5 % of the total integrated intensity) bins will be
                reduced 10 by 10 until the criteria is satisfied or until nbins 
                = 50.

        Returns
            brixs.Spectrum

        .. _numpy.histogram(): https://numpy.org/doc/stable/reference/generated/numpy.histogram.html
        """
        if nbins is None:
            nbins = 1000
            _flat  = arraymanip.flatten(self._data)
            flat = _flat[~np.isnan(_flat)]

            hist, bin_edges = np.histogram(flat, bins=nbins)
            while max(hist) < self.shape[0]*self.shape[1]*0.05:
                nbins -= 10
                hist, bin_edges = np.histogram(flat, bins=nbins)
                if nbins < 50:
                    break
        elif numanip.is_integer(nbins):
            _flat  = arraymanip.flatten(self._data)
            flat = _flat[~np.isnan(_flat)]
            hist, bin_edges = np.histogram(flat, bins=int(nbins))
        else:
            raise TypeError('nbins must be a integer')

        x = arraymanip.moving_average(bin_edges, 2)
        s = Spectrum(x=x, y=hist)

        ##################
        # transfer attrs #
        ##################
        s._copy_attrs_from(self)
 
        return s


    def integrated_rows_vs_y_centers(self):
        """return spectra with integrated row values vs y centers

        Returns:
            :py:class:`Spectrum`
        """
        s = Spectrum(x=self.y_centers, y=np.sum(self._data, axis=1))
        s._copy_attrs_from(self)
        return s
    
    def integrated_columns_vs_x_centers(self):
        """return spectra with integrated column values vs x centers

        Returns:
            :py:class:`Spectrum`
        """
        s = Spectrum(x=self.x_centers, y=np.sum(self._data, axis=0))
        s._copy_attrs_from(self)
        return s
    
    def calculate_spectrum(self, axis=1):
        """Integrate data in one direction (sum columns or rows).

        Args:
            axis (int or string, optional): Axis along which elements are integrated.
                If axis = 0, spectra will be integrated vertically (each spectrum 
                datapoint is the sum of a "column of pixels"). If axis = 1, 
                spectra will be integrated
                horizontally (each spectrum datapoint is a "row"). 
                Default is 1.

        Returns:
            :py:class:`Spectrum`.
        """
        # check axis
        if axis != 0 and axis != 1:
            raise ValueError('axis must be 0 or 1')

        # calculation
        if axis == 0:
            return self.integrated_columns_vs_x_centers()
        elif axis == 1:
            return self.integrated_rows_vs_y_centers()


    def possible_nbins(self):
        """return possible values for nbins in the y (nrows) and x (ncols) directions."""
        return np.sort(list(numanip.factors(self.shape[0]))), np.sort(list(numanip.factors(self.shape[1])))


    def calculate_horizontal_shift(self, mode='cc', xlimits=None, limit_size=100, **kwargs):
        """Calculate intensity misalignments in terms of x centers.

        Args:
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align rows via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous row.

                 3) 'max': Align the max point of every row. 
                 
                 4) 'peak': Fit one peak in each row and align them 
                 (requires that `brixs.addons.fitting` is imported)

            xlimits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            limit_size (int, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. It 
                ensures that the number of rows is not bigger than 
                limit_size. Default is 100.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'` 
            
        Returns:
            list
        """
        ss = self.get_rows(max_number_of_rows=limit_size)
        values = ss.calculate_shift(mode=mode, limits=xlimits, **kwargs)
        return values
    
    def calculate_vertical_shift(self, mode='cc', ylimits=None, limit_size=100, **kwargs):
        """Calculate intensity misalignments in terms of y centers.

        Args:
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align columns via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous column.

                 3) 'max': Align the max point of every column. 
                 
                 4) 'peak': Fit one peak in each column and align them 
                 (requires that `brixs.addons.fitting` is imported)

            ylimits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            limit_size (int, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. It 
                ensures that the number of colums is not bigger than 
                limit_size. Default is 100.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'` 
            
        Returns:
            list
        """
        ss = self.get_columns(max_number_of_columns=limit_size)
        values = ss.calculate_shift(mode=mode, limits=ylimits, **kwargs)
        return values
    
    ############
    # composed #
    ############
    def calculate_vertical_shift_curvature(self, deg=2, mode='cc', ylimits=None, limit_size=100, **kwargs):
        """Calculate vertical shift values to fix curvature.

        Note:
            This function calculates the SHIFT values, and not the polynomial 
            curve that fits the curvature. The image will be corrected uppon
            shifting pixel columns or rows by these shift values. 

        Args:
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align columns via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous column.

                 3) 'max': Align the max point of every column. 
                 
                 4) 'peak': Fit one peak in each column and align them 
                 (requires that `brixs.addons.fitting` is imported)

            ylimits (None or list): y center pair of values `(y_start, y_stop)`, a list 
                of pairs `((yi_1, yf_1), (yi_2, yf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `y_start = None` or `y_stop = None` to indicate the minimum or 
                maximum y value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            limit_size (int, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. It 
                ensures that the number of colums is not bigger than 
                limit_size. Default is 100.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            Spectrum: spectrum with shift values vs x centers with the following
                attributes:
            
                    s.fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

                    s.popt (np.array): 1D array of polynomial coefficients 
                        (including coefficients equal to zero) from highest degree to 
                        the constant term.

                    s.R2 (number): R2 error

                    s.model (function): funcion f(x_centers)
        """
        im = self.copy()
        values = im.calculate_vertical_shift(mode=mode, ylimits=ylimits, limit_size=limit_size, **kwargs)

        # calculate poly
        final = Spectrum(x=self.x_centers, y=values)
        result = final.polyfit(deg=deg)
        final.fit   = result['fit']
        final.popt  = result['popt']
        final.model = result['model']
        final.R2    = result['R2']
        return final
    
    def calculate_horizontal_shift_curvature(self, deg=2, mode='cc', xlimits=None, limit_size=100, **kwargs):
        """Calculate horizontal shift values to fix curvature.

        Note:
            This function calculates the SHIFT values, and not the polynomial 
            curve that fits the curvature. The image will be corrected uppon
            shifting pixel columns or rows by these shift values. 

        Args:
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align rows via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous row.

                 3) 'max': Align the max point of every row. 
                 
                 4) 'peak': Fit one peak in each row and align them 
                 (requires that `brixs.addons.fitting` is imported)

            xlimits (None or list): x center pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            limit_size (int, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. It 
                ensures that the number of rows is not bigger than 
                limit_size. Default is 100.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            Spectrum: spectrum with shift values vs y centers with the following
                attributes:
            
                    s.fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

                    s.popt (np.array): 1D array of polynomial coefficients 
                        (including coefficients equal to zero) from highest degree to 
                        the constant term.

                    s.R2 (number): R2 error

                    s.model (function): funcion f(y_centers)
        """
        im = self.copy()
        values = im.calculate_horizontal_shift(mode=mode, xlimits=xlimits, limit_size=limit_size, **kwargs)

        # calculate poly
        final = Spectrum(x=self.y_centers, y=values)
        result = final.polyfit(deg=deg)
        final.fit   = result['fit']
        final.popt  = result['popt']
        final.model = result['model']
        final.R2    = result['R2']
        return final
    

    def fix_vertical_shift_curvature(self, deg=2, mode='cc', ylimits=None, limit_size=100, **kwargs):
        """Roll column of pixels to fix curvature.

        Args:
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align columns via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous column.

                 3) 'max': Align the max point of every column. 
                 
                 4) 'peak': Fit one peak in each column and align them 
                 (requires that `brixs.addons.fitting` is imported)

            ylimits (None or list): y center pair of values `(y_start, y_stop)`, a list 
                of pairs `((yi_1, yf_1), (yi_2, yf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `y_start = None` or `y_stop = None` to indicate the minimum or 
                maximum y value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            limit_size (int, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. It 
                ensures that the number of colums is not bigger than 
                limit_size. Default is 100.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            Image
        """
        _result = self.calculate_vertical_shift_curvature(deg=deg, mode=mode, ylimits=ylimits, limit_size=limit_size, **kwargs)

        return self.set_vertical_shift_via_polyval(_result['popt']) 
    
    def fix_horizontal_shift_curvature(self, deg=2, mode='cc', xlimits=None, limit_size=100, **kwargs):
        """Roll row of pixels to fix curvature.

        Args:
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align rows via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous row.

                 3) 'max': Align the max point of every row. 
                 
                 4) 'peak': Fit one peak in each row and align them 
                 (requires that `brixs.addons.fitting` is imported)

            xlimits (None or list): x center pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            limit_size (int, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. It 
                ensures that the number of rows is not bigger than 
                limit_size. Default is 100.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            Image
        """
        _result = self.calculate_horizontal_shift_curvature(deg=deg, mode=mode, xlimits=xlimits, limit_size=limit_size, **kwargs)

        return self.set_horizontal_shift_via_polyval(_result['popt']) 

    
    def get_spot(self, y, x, ny, nx, coordinates='centers'):
        """return a square image with side (ny*2+1, nx*2+1) with pixels surrounding x, y

        Note:
            if x, y is too close to the edges, spot is considered up the edge and 
            final spot may not be a square array and may have sides shorter than  
            n*2+1

        Args:
            y, x (tuple or list): (y, x) position in terms of pixel value (if 
                coordinates='pixels') or in terms of y_centers and x_centers (if
                coordinates='centers').
            ny, nx (int): number of neighbors to include, e.g., if ny=nx=1, only first neighbors 
                will be considered and the spot will be a 3x3 array.
            coordinates (str, optional): `pixels` if pixel positions are given in 
                pixel coordinates or `centers` if positions are given in terms of 
                y_centers and x_centers.

        Returns:
            spot, yslice, xslice, is_near_edge = im.get_spot()

            spot is a image

            yslice, xslice are the slices that created spot, i.e., spot=image[yslice, xslice] 
            always in terms of pixels.

            is_near_edge is a dict indicates if the spot touched the edges of the image.
        """
        assert coordinates in ['pixels', 'centers'], f'coordinates must be `pixels` or `centers`, not `{coordinates}`'
        assert numanip.is_integer(nx), f'n must be of type integer, not {type(nx)}'
        assert numanip.is_integer(ny), f'n must be of type integer, not {type(ny)}'
        assert nx >= 0, f'nx must be an integer number higher or equal to 0'
        assert ny >= 0, f'ny must be an integer number higher or equal to 0'

        # returns the full image if n is larger that the figure
        if ny > self.shape[0] and nx > self.shape[1]:
            is_near_edge = {'left':True, 'right':True, 'bottom':True, 'top':True}
            return self.copy(), slice(0, self.shape[0]), slice(0, self.shape[1]), is_near_edge

        # get square
        if coordinates == 'centers':
            y, x = self.center2index(y=y, x=x)
        x_start = int(x-int(nx))
        x_stop  = int(x+int(nx))+1
        y_start = int(y-int(ny))
        y_stop  = int(y+int(ny))+1

        # get edges
        is_near_edge = {'left':True, 'right':True, 'bottom':True, 'top':True}
        if x_start <= 0: 
            x_start = 0
            is_near_edge['left'] = False
        if y_start <= 0: 
            y_start = 0
            is_near_edge['bottom'] = False
        if x_stop >= self.shape[1]: 
            x_stop = self.shape[1]
            is_near_edge['right'] = False
        if y_stop >= self.shape[0]: 
            y_stop = self.shape[0]
            is_near_edge['top'] = False

        # get spot
        spot = Image(data=self.data[y_start:y_stop, x_start:x_stop])
        # if coordinates == 'centers':
        spot.x_centers = self.x_centers[x_start:x_stop]
        spot.y_centers = self.y_centers[y_start:y_stop]

        return spot, slice(y_start, y_stop), slice(x_start, x_stop), is_near_edge

    def patch(self, pos, ny, nx, value=None, coordinates='centers'):
        """Return an image with all y, x pixel positions patched by a ny x nx square
        
        The pixels inside a (ny x nx) square around all y, x pixel positions in pos
        will be replaced by a value. 

        Args:
            pos (tuple or list): (y, x) pixel position or a list of positions 
                [(y1, x1), (y2, x2), ...]. 
            ny, nx (int): number of neighbors around (y, x) to be patched, e.g., if ny=nx=1,
                only first neighbors will be considered and the patch will be a 3x3
                array.
            value (number, optional): patch value. If None, this value will be 
                the average of all 1st neighboring pixels around the n x n patch.
            coordinates (str, optional): `pixels` if pixel positions are given in pixel 
                coordinates or `centers` if positions are given in terms of x_centers
                and y_centers. Default is `centers`.

        Return:
            Image
        """
        # check if arguments make sense
        assert coordinates in ['pixels', 'centers'], f'coordinates must be `pixels` or `centers`, not `{coordinates}`'
        assert numanip.is_integer(nx), f'n must be of type integer, not {type(nx)}'
        assert numanip.is_integer(ny), f'n must be of type integer, not {type(ny)}'
        assert nx >= 0, f'nx must be an integer number higher or equal to 0'
        assert ny >= 0, f'ny must be an integer number higher or equal to 0'

        assert isinstance(pos, Iterable), f'pos must be an Iterable (list or tuple), not type {type(pos)}'
        if isinstance(pos[0], Iterable) == False:
            pos = [pos, ]
        else:
            assert any([isinstance(_, Iterable) for _ in pos]), f'pos must be a pixel position (y, x) or a list of positions [(y1, x1), (y2, x2), ...]'
            assert any([len(_) != 2 for _ in pos])==False, f'pos must be a pixel position (y, x) or a list of positions [(y1, x1), (y2, x2), ...]'
            if coordinates == 'pixels':
                xhigh = [x > self.shape[1] for _, x in pos]
                xlow  = [x < 0 for _, x in pos]
                yhigh = [y > self.shape[0] for y, _ in pos]
                ylow  = [y < 0 for y, _ in pos]
            else:
                xhigh = [x > max(self.x_centers) for _, x in pos]
                xlow  = [x < min(self.x_centers) for _, x in pos]
                yhigh = [y > max(self.y_centers) for y, _ in pos]
                ylow  = [y < min(self.y_centers) for y, _ in pos]
            assert any(xhigh) == False, f'Some x positions are outide of the image. Indexes of bad positions: {np.where(xhigh)[0]}'
            assert any(xlow)  == False, f'Some x positions are outide of the image. Indexes of bad positions: {np.where(xlow)[0]}'
            assert any(yhigh) == False, f'Some y positions are outide of the image. Indexes of bad positions: {np.where(yhigh)[0]}'
            assert any(ylow)  == False, f'Some y positions are outide of the image. Indexes of bad positions: {np.where(ylow)[0]}'

        # patch image
        im = self._copy()
        for y, x in pos:

            # get spot around x, y
            spot, yslice, xslice, is_near_edge = self.get_spot(y=y, x=x, ny=ny, nx=nx, coordinates=coordinates)

            # if value is None, get average value around patch
            x_start, x_stop, y_start, y_stop = xslice.start, xslice.stop, yslice.start, yslice.stop
            if value is None:
                _positions_around = []
                if is_near_edge['left']:   _positions_around += [(x_start-1, _) for _ in range(y_start, y_stop)]
                if is_near_edge['right']:  _positions_around += [(x_stop, _)    for _ in range(y_start, y_stop)]
                if is_near_edge['bottom']: _positions_around += [(_, y_start-1) for _ in range(x_start, x_stop)]
                if is_near_edge['top']:    _positions_around += [(_, y_stop)    for _ in range(x_start, x_stop)]
                if is_near_edge['left']  and is_near_edge['bottom']: _positions_around += [(x_start-1, y_start-1)]
                if is_near_edge['left']  and is_near_edge['top']:    _positions_around += [(x_start-1, y_stop)]
                if is_near_edge['right'] and is_near_edge['bottom']: _positions_around += [(x_stop, y_start-1)]
                if is_near_edge['right'] and is_near_edge['top']:    _positions_around += [(x_stop, y_stop)]
                _values = []
                for _x, _y in _positions_around:
                    _values.append(self.data[_y, _x])
                valuefinal = np.mean(_values)
            else:
                valuefinal = value
                
            # replace value
            im.data[y_start:y_stop, x_start:x_stop] = np.ones((abs(y_start-y_stop), abs(x_start-x_stop)))*valuefinal
        return im


    ##########################        
    # plot and visualization #
    ########################## 
    def pcolorfast(self, ax=None, x_start=None, x_stop=None, y_start=None, y_stop=None, colorbar=False, **kwargs):
        """Fast alternative to im.pcolormesh. Wrapper for `ax.pcolorfast()`_.

        If x_edges and y_edges are not defined and x_centers and y_centers have 
            irregular pixel separation, pcolormesh does its best to defined pixel 
            edges so centers labels correspond to the real centers (nearest possible).

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
            **kwargs: kwargs are passed to `matplotlib.pyplot.pcolormesh()`_.

        If not specified, the following parameters are passed to `matplotlib.pyplot.pcolormesh()`_:

        Args:
            cmap: The Colormap instance. Default is 'jet'.
            vmin: Minimum intensity that the colormap covers. The intensity histogram is
                calculated and vmin is set on the position of the maximum.
            vmax: Maximmum intensity that the colormap covers.  The intensity histogram is
                calculated and vmax is set to the value where the 
                intensity drops below 0.01 % of the maximum.

        Returns:
            `matplotlib.collections.QuadMesh`_

        .. _matplotlib.pyplot.pcolormesh(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.pcolormesh.html
        .. _matplotlib.collections.QuadMesh: https://matplotlib.org/3.5.0/api/collections_api.html#matplotlib.collections.QuadMesh
        """
        ##########
        # limits #
        ##########
        im = self.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)

        ##########
        # kwargs #
        ##########
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'
        if 'vmin' not in kwargs or 'vmax' not in kwargs:
            vmin, vmax = self._calculated_vmin_vmax()
            if 'vmin' not in kwargs:
                kwargs['vmin'] = vmin
            if 'vmax' not in kwargs:
                kwargs['vmax'] = vmax

        ####################
        # colorbar divider #
        ####################
        divider = False
        if ax is not None and colorbar is True:
            divider = True
        
        ###################
        # figure and axes #
        ###################
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW or len(plt.get_fignums()) == 0:
                figure()

        #############
        # get edges #
        #############
        if im.x_edges is None: 
            if im.x_monotonicity is None:
                im.check_x_monotonicity()
            im = im.estimate_x_edges_from_centers()
        if im.y_edges is None: 
            if im.y_monotonicity is None:
                im.check_y_monotonicity()
            im = im.estimate_y_edges_from_centers()

        ########
        # plot #
        ########
        if ax == plt:
            ax = plt.gca()
        X, Y = np.meshgrid(im.x_edges, im.y_edges)
        pos  = ax.pcolorfast(X, Y, im.data, **kwargs)

        ###########################################
        # show x, y, z values upon mouse hovering #
        ###########################################
        def format_coord(x, y):
            xarr = X[0,:]
            yarr = Y[:,0]
            if ((x > xarr.min()) & (x <= xarr.max()) & 
                (y > yarr.min()) & (y <= yarr.max())):
                col = np.searchsorted(xarr, x)-1
                row = np.searchsorted(yarr, y)-1
                z = im.data[row, col]
                return f'x={x:1.4f}, y={y:1.4f}, z={z:1.4f}   [{row},{col}]'
            else:
                return f'x={x:1.4f}, y={y:1.4f}'
        if ax == plt:
            ax.gca().format_coord = format_coord
        else:
            ax.format_coord = format_coord

        ############
        # colorbar #
        ############
        if colorbar:
            if divider:
                divider = make_axes_locatable(ax)
                ax_cb = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(pos, cax=ax_cb)
            else:
                plt.colorbar(pos, aspect=50)

        # add edges and centers to quadmesh
        pos.x_edges   = pos.get_coordinates()[0][:, 0]
        pos.y_edges   = pos.get_coordinates()[:, 1][:, 1]
        pos.x_centers = im.x_centers
        pos.y_centers = im.y_centers

        return pos
        
    def pcolormesh(self, ax=None, x_start=None, x_stop=None, y_start=None, y_stop=None, colorbar=False, **kwargs):
        """Display data as a mesh. Wrapper for `matplotlib.pyplot.pcolormesh()`_.

        If x_edges and y_edges are not defined and x_centers and y_centers have 
            irregular pixel separation, pcolormesh does its best to defined pixel 
            edges so centers labels correspond to the real centers (nearest possible).

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
            **kwargs: kwargs are passed to `matplotlib.pyplot.pcolormesh()`_.

        If not specified, the following parameters are passed to `matplotlib.pyplot.pcolormesh()`_:

        Args:
            cmap: The Colormap instance. Default is 'jet'.
            vmin: Minimum intensity that the colormap covers. The intensity histogram is
                calculated and vmin is set on the position of the maximum.
            vmax: Maximmum intensity that the colormap covers.  The intensity histogram is
                calculated and vmax is set to the value where the 
                intensity drops below 0.01 % of the maximum.

        Returns:
            `matplotlib.collections.QuadMesh`_

        .. _matplotlib.pyplot.pcolormesh(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.pcolormesh.html
        .. _matplotlib.collections.QuadMesh: https://matplotlib.org/3.5.0/api/collections_api.html#matplotlib.collections.QuadMesh
        """
        ##########
        # limits #
        ##########
        im = self.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)

        ##########
        # kwargs #
        ##########
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'
        if 'vmin' not in kwargs or 'vmax' not in kwargs:
            vmin, vmax = self._calculated_vmin_vmax()
            if 'vmin' not in kwargs:
                kwargs['vmin'] = vmin
            if 'vmax' not in kwargs:
                kwargs['vmax'] = vmax

        ####################
        # colorbar divider #
        ####################
        divider = False
        if ax is not None and colorbar is True:
            divider = True
        
        ###################
        # figure and axes #
        ###################
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW or len(plt.get_fignums()) == 0:
                figure()

        #############
        # get edges #
        #############
        if im.x_edges is None: 
            if im.x_monotonicity is None:
                im.check_x_monotonicity()
            assert im.x_monotonicity == 'increasing', 'x edges are not increasingly monotonic. Use im.fix_x_monotonicity() of adjust im.x_centers or im.x_edges directly'
            im = im.estimate_x_edges_from_centers()
        else:
            s = Spectrum(x=self.x_edges)
            try:
                s.check_monotonicity()
                assert s.monotonicity == 'increasing', 'x edges are not increasingly monotonic. Use im.fix_x_monotonicity() of adjust im.x_centers or im.x_edges directly'
            except ValueError:
                raise ValueError('x edges are not increasingly monotonic. Use im.fix_x_monotonicity() of adjust im.x_centers or im.x_edges directly')
        
        if im.y_edges is None: 
            if im.y_monotonicity is None:
                im.check_y_monotonicity()
            assert im.y_monotonicity == 'increasing', 'y edges are not increasingly monotonic. Use im.fix_y_monotonicity() of adjust im.y_centers or im.y_edges directly'
            im = im.estimate_y_edges_from_centers()
        else:
            s = Spectrum(x=self.y_edges)
            try:
                s.check_monotonicity()
                assert s.monotonicity == 'increasing', 'y edges are not increasingly monotonic. Use im.fix_y_monotonicity() of adjust im.y_centers or im.y_edges directly'
            except ValueError:
                raise ValueError('y edges are not increasingly monotonic. Use im.fix_y_monotonicity() of adjust im.y_centers or im.y_edges directly')

        ########
        # plot #
        ########
        X, Y = np.meshgrid(im.x_edges, im.y_edges)
        pos  = ax.pcolormesh(X, Y, im.data, **kwargs)

        ###########################################
        # show x, y, z values upon mouse hovering #
        ###########################################
        def format_coord(x, y):
            xarr = X[0, :]
            yarr = Y[:, 0]
            if ((x > xarr.min()) & (x <= xarr.max()) & 
                (y > yarr.min()) & (y <= yarr.max())):
                col = np.searchsorted(xarr, x)-1
                row = np.searchsorted(yarr, y)-1
                z = im.data[row, col]
                return f'x={x:1.4f}, y={y:1.4f}, z={z:1.4f}   [{row},{col}]'
            else:
                return f'x={x:1.4f}, y={y:1.4f}'
        if ax == plt:
            ax.gca().format_coord = format_coord
        else:
            ax.format_coord = format_coord

        ############
        # colorbar #
        ############
        if colorbar:
            if divider:
                divider = make_axes_locatable(ax)
                ax_cb = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(pos, cax=ax_cb)
            else:
                plt.colorbar(pos, aspect=50)

        # add edges and centers to quadmesh
        pos.x_edges   = pos.get_coordinates()[0][:, 0]
        pos.y_edges   = pos.get_coordinates()[:, 1][:, 1]
        pos.x_centers = im.x_centers
        pos.y_centers = im.y_centers

        return pos
        
    def imshow(self, ax=None, x_start=None, x_stop=None, y_start=None, y_stop=None, colorbar=False, origin='upper', verbose=True, **kwargs):
        """Display data as an image in terms of pixels. Wrapper for `matplotlib.pyplot.imshow()`_.

        Warning:
            Pixels are squares. For irregular pixel row/columns, see Image.pcolormesh().
            For image with axis in terms of x and y centers, use Image.plot().        

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
                (str, optional): Location of the [0, 0] index. Options are
                `upper` and `lower`. Default is 'lower'.
            origin (str, optional): Place the [0, 0] index of the array in the 
                upper left (`origin='upper'`) or lower left corner (`origin='lower'`) of the Axes. 
                Default is 'lower'.
            verbose (bool, optional): if True, a warning will show up if data
                has iregular pixel sizes. Default is true.
            **kwargs: kwargs are passed to `matplotlib.pyplot.imshow()`_.

        If not specified, the following parameters are passed to `matplotlib.pyplot.imshow()`_:

        Args:
            cmap: The Colormap instance. Default is 'jet'.
            aspect: The aspect ratio of the Axes. Default is 'auto'. If 'equal',
                an aspect ratio of 1 will be used (pixels will be square).
            interpolation: The interpolation method used. Default is 'none'.
                Supported values are 'none', 'antialiased', 'nearest', 'bilinear',
                'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite',
                'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell',
                'sinc', 'lanczos', 'blackman'.
            vmin: Minimum intensity that the colormap covers. The intensity histogram is
                calculated and vmin is set on the position of the maximum.
            vmax: Maximmum intensity that the colormap covers. The intensity histogram is
                calculated and vmax is set to the value where the 
                intensity drops below 0.01 % of the maximum.

        Returns:
            `matplotlib.image.AxesImage`_

        .. _matplotlib.pyplot.imshow(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.imshow.html
        .. _matplotlib.image.AxesImage: https://matplotlib.org/3.5.0/api/image_api.html#matplotlib.image.AxesImage
        """
        ##########
        # limits #
        ##########
        im = self.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)

        ##########
        # kwargs #
        ##########
        kwargs['origin'] = origin
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'
        if 'aspect' not in kwargs:
            kwargs['aspect'] = 'auto'
        if 'interpolation' not in kwargs:
            kwargs['interpolation'] = 'none'
        if 'vmin' not in kwargs or 'vmax' not in kwargs:
            vmin, vmax = self._calculated_vmin_vmax()
            if 'vmin' not in kwargs:
                kwargs['vmin'] = vmin
            if 'vmax' not in kwargs:
                kwargs['vmax'] = vmax

        ####################
        # colorbar divider #
        ####################
        divider = False
        if ax is not None and colorbar is True:
            divider = True 

        ###################
        # figure and axes #
        ###################
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW or len(plt.get_fignums()) == 0:
                figure()

        # plot
        pos = ax.imshow(im.data, **kwargs)

        # colorbar
        if colorbar:
            if divider:
                divider = make_axes_locatable(ax)
                ax_cb = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(pos, cax=ax_cb)
            else:
                plt.colorbar(pos, aspect=50)

        # add edges and centers to axesimage object
        pos.x_edges = np.linspace(-0.5, len(im.x_centers)+0.5, len(im.x_centers)+1)
        pos.y_edges = np.linspace(-0.5, len(im.y_centers)+0.5, len(im.y_centers)+1)
        pos.x_centers = arraymanip.moving_average(pos.x_edges, 2)
        pos.y_centers = arraymanip.moving_average(pos.y_edges, 2)
        
        return pos
    
    def plot(self, ax=None, x_start=None, x_stop=None, y_start=None, y_stop=None, colorbar=False, origin='upper', verbose=True, **kwargs):
        """Display data as an image with axis based on x and y centers. Wrapper for `matplotlib.pyplot.imshow()`_.

        Warning:
            Pixels are squares. For irregular pixel row/columns, see Image.pcolormesh()

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
                (str, optional): Location of the [0, 0] index. Options are
                `upper` and `lower`. Default is 'lower'.
            origin (str, optional): Place the [0, 0] index of the array in the 
                upper left (`origin='upper'`) or lower left corner (`origin='lower'`) of the Axes. 
                Default is 'lower'.
            verbose (bool, optional): if True, a warning will show up if data
                has iregular pixel sizes. Default is true.
            **kwargs: kwargs are passed to `matplotlib.pyplot.imshow()`_.

        If not specified, the following parameters are passed to `matplotlib.pyplot.imshow()`_:

        Args:
            cmap: The Colormap instance. Default is 'jet'.
            aspect: The aspect ratio of the Axes. Default is 'auto'. If 'equal',
                an aspect ratio of 1 will be used (pixels will be square).
            interpolation: The interpolation method used. Default is 'none'.
                Supported values are 'none', 'antialiased', 'nearest', 'bilinear',
                'bicubic', 'spline16', 'spline36', 'hanning', 'hamming', 'hermite',
                'kaiser', 'quadric', 'catrom', 'gaussian', 'bessel', 'mitchell',
                'sinc', 'lanczos', 'blackman'.
            extent: minimun and maximum x and y values. Default will be given by
                the Image.x and Image.y attributes.
            vmin: Minimum intensity that the colormap covers. The intensity histogram is
                calculated and vmin is set on the position of the maximum.
            vmax: Maximmum intensity that the colormap covers. The intensity histogram is
                calculated and vmax is set to the value where the 
                intensity drops below 0.01 % of the maximum.

        Returns:
            `matplotlib.image.AxesImage`_

        .. _matplotlib.pyplot.imshow(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.imshow.html
        .. _matplotlib.image.AxesImage: https://matplotlib.org/3.5.0/api/image_api.html#matplotlib.image.AxesImage
        """
        ##########
        # limits #
        ##########
        im = self.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)

        ##########
        # kwargs #
        ##########
        assert origin == 'lower' or origin == 'upper', f'origin can only be `lower` or `upper`, not `{origin}`'
        kwargs['origin'] = origin
        kwargs.setdefault('cmap', 'jet')
        kwargs.setdefault('aspect', 'auto')
        kwargs.setdefault('interpolation', 'none')
        if 'vmin' not in kwargs or 'vmax' not in kwargs:
            vmin, vmax = self._calculated_vmin_vmax()
            if 'vmin' not in kwargs:
                kwargs['vmin'] = vmin
            if 'vmax' not in kwargs:
                kwargs['vmax'] = vmax

        ######################
        # check monotonicity #
        ######################
        if self.x_monotonicity is None:
            try:
                self.check_x_monotonicity()
            except ValueError:
                raise ValueError(f'x centers is not monotonic. Plot image using im.imshow() or maybe use im.fix_x_monotonicity()')
        if self.y_monotonicity is None:
            try:
                self.check_y_monotonicity()
            except ValueError:
                raise ValueError(f'y centers is not monotonic. Plot image using im.imshow() or maybe use im.fix_y_monotonicity()')
        
        ##############
        # check step #
        ##############
        if self.x_step is None or self.y_monotonicity is None:
            try:
                self.check_x_step()
                self.check_y_step()
            except ValueError:
                if verbose:
                    print('Warning: Data seems to have irregular pixel size. Maybe plot it using Image.pcolormesh().\nTo turn off this warning set verbose to False.')

        # extent
        if 'extent' not in kwargs:
            # x  = np.linspace(im.x_centers[0], im.x_centers[-1], len(im.x_centers))
            # dx = np.mean(np.diff(x))
            # extent_x = [x[0]-dx/2, x[-1]+dx/2]

            if im.x_edges is None:
                im = im.estimate_x_edges_from_centers()
            # extent_x = [min(im.x_edges), max(im.x_edges)]
            extent_x = [im.x_edges[0], im.x_edges[-1]]

            if im.y_edges is None:
                im = im.estimate_y_edges_from_centers()
            if origin == 'upper':
                # extent_y = [max(im.y_edges), min(im.y_edges)]
                extent_y = [im.y_edges[-1], im.y_edges[0]]
            else:
                # extent_y = [min(im.y_edges), max(im.y_edges)]
                extent_y = [im.y_edges[0], im.y_edges[-1]]
            
            # y  = np.linspace(im.y_centers[0], im.y_centers[-1], len(im.y_centers))
            # dy = np.mean(np.diff(y))
            # if origin == 'upper':
            #     extent_y = [y[-1]+dy/2, y[0]-dy/2]
            # else:
            #     extent_y = [y[0]-dy/2, y[-1]+dy/2]

            kwargs['extent'] = np.append(extent_x, extent_y)

        ####################
        # colorbar divider #
        ####################
        divider = False
        if ax is not None and colorbar is True:
            divider = True 

        ###################
        # figure and axes #
        ###################
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW or len(plt.get_fignums()) == 0:
                figure()

        # plot
        pos = ax.imshow(im.data, **kwargs)

        # colorbar
        if colorbar:
            if divider:
                divider = make_axes_locatable(ax)
                ax_cb = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(pos, cax=ax_cb)
            else:
                plt.colorbar(pos, aspect=50)

        # add edges and centers to axesimage object
        pos.x_edges = np.linspace(pos._extent[0], pos._extent[1], len(im.x_centers)+1)
        pos.y_edges = np.linspace(pos._extent[2], pos._extent[3], len(im.y_centers)+1)
        pos.x_centers = arraymanip.moving_average(pos.x_edges, 2)
        pos.y_centers = arraymanip.moving_average(pos.y_edges, 2)
        
        return pos

# %% ============================ PhotonEvents =========================== %% #
class PhotonEvents(_BrixsObject, metaclass=_Meta):
    """Returns a ``Photon events`` object.

    Args:
        x, y (list or array, optional): array with x, y photon events coordinates
        xlim, ylim (list, optional): two element tuple with min and max possible 
            x and y coordinates. Used for defining binning limits
             and for plotting.

    Usage:        
            >>> pe = br.PhotonEvents()
            >>> pe = br.PhotonEvents(x, y)
            >>> pe = br.PhotonEvents(x=x, y=y)
            >>> pe = br.PhotonEvents(x, y, xlim=(0, 10), ylim=(0, 10))
            >>> pe = br.PhotonEvents().load(filepath=<filepath>)
            >>> pe = br.PhotonEvents().load(filepath=<filepath>, delimiter=',')
            >>>
            >>> print(pe.get_core_attrs()) # print list of core attrs
            >>> print(pe.get_attrs())      # print list of attrs
            >>> print(pe.get_methods())    # print list of methods available

    Attributes:
        Every BRIXS object has 5 types of attributes: 
        `Core` , `Check`, `Modifiers`, `Labels`, `User`.

        *1. Core*
            x, y (array): 1D arrays with uncorrelated x, y data
        
        *2. Check*
            None

        *3. Modifiers*
            None
        
        *4. Labels*
            xlim, ylim (number): two element tuple with min and max possible 
                x and y coordinates. Used for defining binning limits
                and for plotting.

        *5. User*
            anything that the user defined on the fly.        
    """
    _read_only = ['has_nan']
    _non_removable = []
    
    def __init__(self, x=[], y=[], xlim=None, ylim=None, filepath=None, **kwargs): 
        """Initialize the object instance"""
        ###########################
        # Initializing attributes #
        ###########################
        # core
        self._x = None
        self._y = None

        # check
        self._has_nan = None

        # modifiers
        pass

        # labels
        self._xlim = None
        self._ylim = None

        # extra
        for extra in settings._init['PhotonEvents']:
            self.__setattr__(extra, settings._init['PhotonEvents'][extra](self))

        ################
        # loading data #
        ################
        if filepath is not None:
            self.load(filepath=filepath, **kwargs)
        else:
            if not isinstance(y, Iterable):
                raise TypeError(f'The y-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a y-array which is not an Iterable.\nThe type of the variable you passed is: {type(y)}\nAccepted types are: list, array, ...')
            if len(y) > 0:
                try:
                    _ = sum(y)
                except:
                    raise TypeError(f'The y-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a y-array which is Iterable, but is NOT ENTIRELY made of numbers because we fail to sum all the elements of the array you are passing.')
            self._y = y
            self.x = x

            self.xlim = xlim
            self.ylim = ylim       
        return

    ###################
    # core attributes #
    ###################
    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, value):
        if value is None:
            self._x = []
            self._y = []

            self._has_nan = None
            return
        ###################################
        # asserting validity of the input #
        ###################################
        # check type
        if not isinstance(value, Iterable):
            raise TypeError(f'The x-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a x-array which is not an Iterable.\nThe type of the variable you passed is: {type(value)}\nAccepted types are: list, array, ...')
    
        # check if iterable is made of numbers
        if len(value) > 0:
            try:
                _ = sum(value)
            except:
                raise TypeError(f'The x-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a x-array which is Iterable, but is NOT ENTIRELY made of numbers because we fail to sum all the elements of the array you are passing.')
        
        # check length
        if self.y is not None:
            assert len(value) == len(self.y), f'Length of x array (len={len(value)}) you are trying to set is not compatible with current length of the y array (len={len(self.y)}).'

        # empty array
        if len(value) == 0:
            self._x = []
            self._has_nan = None
            return
        
        #################
        # set attribute #
        #################
        self._x = np.array(value, dtype='float')

        self._has_nan = None
    @x.deleter
    def x(self):
        raise AttributeError('Cannot delete object.')

    @property
    def y(self):
        return self._y
    @y.setter
    def y(self, value):
        # check None
        if value is None:
            self._x = []
            self._y = []
            self._has_nan = None

            return
        ###################################
        # asserting validity of the input #
        ###################################
        # check type
        if not isinstance(value, Iterable):
            raise TypeError(f'The y-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a y-array which is not an Iterable.\nThe type of the variable you passed is: {type(value)}\nAccepted types are: list, array, ...')
    
        # check if iterable is made of numbers
        if len(value) > 0:
            try:
                _ = sum(value)
            except:
                raise TypeError(f'The y-array must be an Iterable (list or array) of numbers.\nYou are trying to set up a y-array which is Iterable, but is NOT ENTIRELY made of numbers because we fail to sum all the elements of the array you are passing.')
        
        # check length
        if self.x is not None:
            assert len(value) == len(self.x), f'Length of y-array (len={len(value)}) you are trying to set is not compatible with current length of the x-array (len={len(self.x)}).'

        # empty array
        if len(value) == 0:
            self._y = []
            self._has_nan = None
            return
        
        #################
        # set attribute #
        #################
        self._y = np.array(value, dtype='float')

        self._has_nan = None
    @y.deleter
    def y(self):
        raise AttributeError('Cannot delete object.')

    @property
    def xlim(self):
        return self._xlim
    @xlim.setter
    def xlim(self, value):
        ###################################
        # asserting validity of the input #
        ###################################
        if value is None:
            self._xlim = None#np.array((min(self.x), max(self.x)), dtype='float')
        else:
            # check type
            if isinstance(value, Iterable) == False:
                raise TypeError(f'xlim must be an Iterable (list or array) of numbers.\nYou are trying to set up a xlim which is not an Iterable.\nThe type of the variable you passed is: {type(value)}\nAccepted types are: list, array, ...')
        
            # check if iterable is made of numbers
            try:
                _ = sum(value)
            except:
                raise TypeError(f'xlim must be an Iterable (list or array) of numbers.\nYou are trying to set up a xlim which is Iterable, but is NOT ENTIRELY made of numbers because we fail to sum all the elements of the array you are passing.')
            
            # check length
            assert len(value) == 2, f'Length of xlim must be 2 and not: {len(value)}'

            # check min and max
            assert value[1] > value[0], f'max value must be bigger than min value'

            # check if data is within limits
            # if self.x is not None:
            #     assert max(self.x) <= value[-1] and min(self.x) >= value[0], f'xlim={value} not allowed because x coordinates (from {min(self.x)} to {max(self.x)}) are outside xlim. xlim set to None'

            #################
            # set attribute #
            #################
            self._xlim = np.array(value, dtype='float')
    @xlim.deleter
    def xlim(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def ylim(self):
        return self._ylim
    @ylim.setter
    def ylim(self, value):
        ###################################
        # asserting validity of the input #
        ###################################
        if value is None:
            self._ylim = None#np.array((min(self.y), max(self.y)), dtype='float')
        else:
            # check type
            if not isinstance(value, Iterable):
                raise TypeError(f'ylim must be an Iterable (list or array) of numbers.\nYou are trying to set up a ylim which is not an Iterable.\nThe type of the variable you passed is: {type(value)}\nAccepted types are: list, array, ...')
        
            # check if iterable is made of numbers
            try:
                _ = sum(value)
            except:
                raise TypeError(f'ylim must be an Iterable (list or array) of numbers.\nYou are trying to set up a ylim which is Iterable, but is NOT ENTIRELY made of numbers because we fail to sum all the elements of the array you are passing.')
            
            # check length
            assert len(value) == 2, f'Length of ylim must be 2 and not: {len(value)}'

            # check min and max
            assert value[1] > value[0], f'max value must be bigger than min value'

            # check if data is within limits
            # if self.y is not None:
            #     assert max(self.y) <= value[-1] and min(self.y) >= value[0], f'ylim={value} not allowed bacause y coordinates (from {min(self.y)} to {max(self.y)}) are outside ylim. ylim set to None'

            #################
            # set attribute #
            #################
            self._ylim = np.array(value, dtype='float')
    @ylim.deleter
    def ylim(self):
        raise AttributeError('Cannot delete object.')

    ###################################
    # computed (read-only) attributes #
    ###################################
    @property
    def data(self):
        if self.x is not None:
            return np.vstack((self.x, self.y)).transpose()
    @data.setter
    def data(self, value):
        raise AttributeError('Attribute is "read only".')
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    #########################
    # write-only attributes #
    #########################
    pass

    #######################
    # modifier attributes #
    #######################
    pass

    #################
    # magic methods #
    #################
    def __setattr__(self, name, value):
        if name in settings._forbidden_words['PhotonEvents']:
            raise AttributeError(f'`{name}` is a reserved word and cannot be set as an attribute')
        super().__setattr__(name, value)

    def __len__(self):
        if self.x is None:
            return 0
        else:
            return len(self.x)

    def __add__(self, object):
        if isinstance(object, PhotonEvents):
            if self.x is None:
                final = PhotonEvents(x=object.x, y=object.y)
            elif object.x is None:
                final = PhotonEvents(x=self.x, y=self.y)
            else:
                final = PhotonEvents(x=list(self.x) + list(object.x), y=list(self.y) + list(object.y))
            final._copy_attrs_from(self)

            # fix limits
            try:
                final.xlim = (min([self.xlim[0], object.xlim[0]]), max([self.xlim[1], object.xlim[1]]))
            except TypeError:
                pass
            try:
                final.ylim = (min([self.ylim[0], object.ylim[0]]), max([self.ylim[1], object.ylim[1]]))
            except TypeError:
                pass
        else:
            raise ValueError(f'Cannot operate type {type(object)} with Photon Events')

        return final
        
    def __sub__(self, object):
        raise ValueError('Photon Events cannot be subtracted')
    
    def __mul__(self, object):
        raise ValueError('Photon Events cannot be multiplied')
    
    def __div__(self, object):
        raise ValueError('Photon Events cannot be divided')
    
    def __truediv__(self, object):
        return self.__div__(object)

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.x[item], self.y[item]
        elif isinstance(item, slice):
            x = self.x[item]
            y = self.y[item]

            pe = PhotonEvents(x=x, y=y)

            # transfer attrs
            for attr in self.get_attrs():
                value = copy.deepcopy(self.__dict__[attr])
                pe.__setattr__(attr, value)

            return pe
        else:
            raise TypeError('Index must be int or a slice, not {}'.format(type(item).__name__))

    def __delitem__(self, item):
        if isinstance(item, int) or isinstance(item, slice):
            self._x = np.delete(self.x, item)
            self._y = np.delete(self.y, item)
        else:
            raise TypeError('Index must be int or a slice, not {}'.format(type(item).__name__))
              
    #########
    # attrs #
    #########
    def append(self, x, y):
        """append new photon event
        
        Args:
            x, y (number): x, y positions of the new photon event 
            
        Returns:
            None
        """
        self._x = np.append(self._x, x)
        self._y = np.append(self._y, y)
        self._has_nan = None
        return

    ###########
    # support #
    ###########
    def _check_mask(self, mask):
        """returns mask in the right format.

        Args:
            mask (None or list): list with 4 elements `(x_start, x_sto, y_start, y_stop)`, a list 
                like `((xi_1, xf_1, yi_1, yf_1), (xi_2, xf_2, yi_2, yf_2), ...)` 
                in terms of x_centers and y_centers, or None. If None, 
                this function simply returns None. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If mask = [], i.e.,
                an empty list, it assumes `limits = (None, None, None, None)`.

        Returns:
            None or limits in the following format:
                ((xi_1, xf_1, yi_1, yf_1), (xi_2, xf_2, yi_2, yf_2), ...)              
        """
        #####################
        # if limits is None #
        #####################
        if mask is None:
            return None
        
        ################################
        # assert that mask is Iterable #
        ################################
        assert isinstance(mask, Iterable), f'`limits` must be an Iterable, not {type(limits)}'
        
        ################
        # empty object #
        ################
        if self.x is None:
            raise ValueError('cannot operate on empty photon events')

        ################################
        # get min and max range values #
        ################################
        if self.xlim is None:
            xmin = min(self.x)
            xmax = max(self.x)
        else:
            xmin = self.xlim[0]
            xmax = self.xlim[1]
        if self.ylim is None:
            ymin = min(self.y)
            ymax = max(self.y)
        else:
            ymin = self.ylim[0]
            ymax = self.ylim[1]

        ################
        # empty limits #
        ################
        if len(mask) == 0:
            return ((xmin, xmax, ymin, ymax),)

        ##############
        # fix format #
        ##############
        # 4 elements 
        if len(mask) == 4: # (xi, xf, yi, yf), or ((xi1, xf1, yi1, yf1), ..., (xi2, xf2, yi2, yf2))
            if isinstance(mask[0], Iterable) == False:
                if isinstance(mask[1], Iterable) == False and isinstance(mask[2], Iterable) == False and isinstance(mask[3], Iterable) == False:
                    mask = [mask, ]
                else:
                    raise ValueError(f'wrong format for mask={mask}')
        # 1, 2, 3 or more than 4 elements 
        final = []
        for m in mask:
            assert isinstance(m, Iterable), f'wrong format for mask={mask}'
            temp = [m[0], m[1], m[2], m[3]]
            if temp[0] == None: temp[0] = xmin
            if temp[1] == None: temp[1] = xmax
            if temp[2] == None: temp[2] = ymin
            if temp[3] == None: temp[3] = ymax
            final.append((temp[0], temp[1], temp[2], temp[3]))   
        return final

    ################
    # core methods #
    ################
    pass

    ########
    # copy #
    ########
    def _copy(self, x_start=None, x_stop=None, y_start=None, y_stop=None, attrs2crop=None):
        """Same as copy(), but attributes are not copied to the new object.
        
        xlim and ylim only change if boundaries are smaller then the original xlim and ylim"""
        #############
        # copy data #
        #############
        x = copy.deepcopy(self.x)
        y = copy.deepcopy(self.y)

        #####################
        # if limits is None #
        #####################
        if x_start is None and x_stop is None and y_start is None and y_stop is None:
            return PhotonEvents(x=x, y=y, xlim=copy.deepcopy(self.xlim), ylim=copy.deepcopy(self.ylim))

        ###############################
        # check if empty PhotonEvents #
        ###############################
        if self.x is None:
            return PhotonEvents(x=x, y=y, xlim=copy.deepcopy(self.xlim), ylim=copy.deepcopy(self.ylim))
            
        #################
        # check if None #
        #################
        if x_start is None: 
            if self.xlim is None:
                x_start = min(self.x)
            else:
                x_start = self.xlim[0]
        if x_stop  is None: 
            if self.xlim is None:
                x_stop  = max(self.x)
            else:
                x_stop = self.xlim[1]
        if y_start is None: 
            if self.ylim is None:
                y_start = min(self.y)
            else:
                y_start = self.ylim[0]
        if y_stop  is None: 
            if self.ylim is None:
                y_stop  = max(self.y)
            else:
                y_stop = self.ylim[1]

        ########################################
        # check if extract is really necessary #
        ########################################
        if x_start <= min(self.x) and x_stop >= max(self.x):
            if y_start <= min(self.y) and y_stop >= max(self.y):
                return PhotonEvents(x=x, y=y, xlim=copy.deepcopy(self.xlim), ylim=copy.deepcopy(self.ylim))

        #################################
        # assert validity of attrs2clip #
        #################################
        if attrs2crop is not None:
            _attrs2crop = {}
            assert isinstance(attrs2crop, Iterable), f'attrs2crop must be a list (Iterable), not type `{type(attrs2crop)}`'
            for attr in attrs2crop:
                assert hasattr(self, attr), f'attrs2crop cannot find attr `{attr}`'
                assert isinstance(self.__getattribute__(attr), Iterable), f'cannot clip attr `{attr}`. It must be an Iterable'
                assert len(self.__getattribute__(attr)), f'cannot clip attr `{attr}` (length={len(self.__getattribute__(attr))}) because it does not have the same lenght as the number of photon events ({len(self)})'
                _attrs2crop[attr] = []

        ##################
        # validate input #
        ##################
        assert x_stop > x_start, f'x_start ({x_start}) must be smaller than x_stop ({x_stop})'
        assert y_stop > y_start, f'y_start ({y_start}) must be smaller than y_stop ({y_stop})'     
        
        ########
        # crop #
        ########
        if attrs2crop is not None:  # exection is much slower with attrs
            x = []
            y = []
            for i in range(len(self)):
                _x = self.x[i]
                _y = self.y[i]
                if (_x > x_start and _x < x_stop) and (_y > y_start and _y < y_stop):
                    x.append(_x)
                    y.append(_y)
                    for attr in _attrs2crop:
                        _attrs2crop[attr].append(self.__getattribute__(attr)[i])
            # _pe = PhotonEvents(x=x, y=y, xlim=(x_start, x_stop), ylim=(y_start, y_stop))
            # _pe = PhotonEvents(x=x, y=y, xlim=copy.deepcopy(self.xlim), ylim=copy.deepcopy(self.ylim))
            
            # fix limits
            xlim = [self.xlim[0], self.xlim[1]]
            ylim = [self.ylim[0], self.ylim[1]]
            if x_start > self.xlim[0]:
                xlim[0] = x_start
            if x_stop < self.xlim[1]:
                xlim[1] = x_stop
            if y_start > self.ylim[0]:
                ylim[1] = y_start
            if y_stop < self.ylim[1]:
                ylim[1] = y_stop
            _pe = PhotonEvents(x=x, y=y, xlim=xlim, ylim=ylim)
            for attr in _attrs2crop:
                _pe.__setattr__(attr, _attrs2crop[attr])
            return _pe
        else:
            # if len(self.x) == 0:
            #     return PhotonEvents(xlim=copy.deepcopy(self.xlim), ylim=copy.deepcopy(self.ylim))
                
            temp = np.array([(x, y) for x, y in zip(self.x, self.y) if ((x > x_start and x < x_stop) and (y > y_start and y < y_stop))])
            
            # fix limits
            if self.xlim is not None:
                xlim = [self.xlim[0], self.xlim[1]]
                if x_start > self.xlim[0]:
                    xlim[0] = x_start
                if x_stop < self.xlim[1]:
                    xlim[1] = x_stop
            else:
                xlim = None#[x_start, x_stop]
            if self.ylim is not None:
                ylim = [self.ylim[0], self.ylim[1]]
                if y_start > self.ylim[0]:
                    ylim[1] = y_start
                if y_stop < self.ylim[1]:
                    ylim[1] = y_stop
            else:
                ylim = None#[y_start, y_stop]

            if len(temp) == 0: # no photons inside limits
                # return PhotonEvents(xlim=copy.deepcopy(self.xlim), ylim=copy.deepcopy(self.ylim))
                return PhotonEvents(xlim=xlim, ylim=ylim)
            else:
                # return PhotonEvents(x=list(temp[:, 0]), y=list(temp[:, 1]), xlim=copy.deepcopy(self.xlim), ylim=copy.deepcopy(self.ylim))
                return PhotonEvents(x=list(temp[:, 0]), y=list(temp[:, 1]), xlim=xlim, ylim=ylim)
        
    def copy(self, x_start=None, x_stop=None, y_start=None, y_stop=None, attrs2crop=None):
        """Return a copy of the object.

        Note:
            pe.xlim[0] only changes if x_start is bigger than pe.xlim[0] and
            pe.xlim[1] only changes if x_stop is smaller than pe.xlim[1]. Otherwise
            pe.xlim is kept. Same for pe.ylim.


        Args:
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.
            attrs2crop (None or list, optional): list of attrs to be crop 
                toghether with the data (x and y). attrs must be a list of same
                length as the number of photon events.

        Returns:
            :py:attr:`PhotonEvents`
        """
        pe = self._copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, attrs2crop=attrs2crop)
        if attrs2crop is not None:
            _attr2crop = {attr: pe.__getattribute__(attr) for attr in attrs2crop}
        pe._copy_attrs_from(self)
        if attrs2crop is not None:
            for attr in attrs2crop:
                pe.__setattr__(attr, _attr2crop[attr])

        # extra
        for extra in settings._copy['PhotonEvents']:
            if hasattr(pe, extra):
                pe.__setattr__(extra, self.__getattribute__(extra).copy(pe))

        return pe

    def clip(self, mask, attrs2clip=None):
        """Return a masked copy of the object.

        Args:
            mask (list): list with rectangular coordinates `(x_start, x_stop, y_start, y_stop)`
                or a list with multiple rectangular coordinates, i.e., `[(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])`
            attrs2clip (None or list, optional): list of attrs to be cliped 
                toghether with the data (x and y). attrs must be a list of same
                length as the number of photon events.

        Returns:
            :py:attr:`PhotonEvents`
        """
        ###############################
        # check if empty PhotonEvents #
        ###############################
        if self.x is None:
            return PhotonEvents(x=None, y=None, xlim=copy.deepcopy(self.xlim), ylim=copy.deepcopy(self.ylim))

        ######################
        # assert mask format #
        ######################
        assert isinstance(mask, Iterable), 'mask must be iterable'
        if len(mask) == 4:
            if isinstance(mask[0], Iterable) == False:
                mask = [mask, ] 
        for m in mask:
            assert len(m) == 4, 'mask must have the format: [(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])'
        
        #################################
        # assert validity of attrs2clip #
        #################################
        if attrs2clip is not None:
            _attrs2clip = {}
            assert isinstance(attrs2clip, Iterable), f'attrs2clip must be a list (Iterable), not type `{type(attrs2clip)}`'
            for attr in attrs2clip:
                assert hasattr(self, attr), f'attrs2clip cannot find attr `{attr}`'
                assert isinstance(self.__getattribute__(attr), Iterable), f'cannot clip attr `{attr}`. It must be an Iterable'
                assert len(self.__getattribute__(attr)), f'cannot clip attr `{attr}` (length={len(self.__getattribute__(attr))}) because it does not have the same lenght as the number of photon events ({len(self)})'
                _attrs2clip[attr] = []

        ########
        # clip #
        ########
        x = []
        y = []
        xlim = None
        ylim = None
        for x_start, x_stop, y_start, y_stop in mask: 
            _pe = self._copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, attrs2crop=attrs2clip)
            if _pe.x is not None:
                x += list(_pe.x)
                y += list(_pe.y)
                if attrs2clip is not None:
                    for attr in _attrs2clip:
                        _attrs2clip[attr] += _pe.__getattribute__(attr)
            
            # # fix xlim and ylim
            # if _pe.xlim is not None:
            #     if xlim is None: xlim = list(_pe.xlim)
            #     if _pe.xlim[0] < xlim[0]: xlim[0] = _pe.xlim[0]
            #     if _pe.xlim[1] > xlim[1]: xlim[1] = _pe.xlim[1]
            # if _pe.ylim is not None:
            #     if ylim is None: ylim = list(_pe.ylim)
            #     if _pe.ylim[0] < ylim[0]: ylim[0] = _pe.ylim[0]
            #     if _pe.ylim[1] > ylim[1]: ylim[1] = _pe.ylim[1]
        
        # check xlim and ylim
        if xlim is None:
            xlim = copy.deepcopy(self.xlim)
        if ylim is None:
            ylim = copy.deepcopy(self.ylim)

        # final
        pe = PhotonEvents(x=x, y=y, xlim=xlim, ylim=ylim)
        pe._copy_attrs_from(self)
        if attrs2clip is not None:
            for attr in _attrs2clip:
                pe.__setattr__(attr, _attrs2clip[attr])

        return pe

    #################
    # save and load #
    #################
    def save(self, filepath=None, only_data=False, check_overwrite=False, verbose=False, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Warning:
            Attrs are saved as comments if only_data is False. Saving attrs to file
            is not always reliable because requires converting variables to string. 
            Only attrs that are of type: string, number, and list of number,
             list of list of number and strings have been tested. Dictionaries are not saved.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format. If None, filepath can be inferred from a
                user defined attr 'im.filepath'. Default is None.  Last used 
                filepath is saved to pe.filepath.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            check_overwrite (bool, optional): if True, it will check if file exists
                and ask if user wants to overwrite file.
            **kwargs: kwargs are passed to ``np.savetxt()`` that saves the data.

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
                Note that, object attributes will be saved at the beginning of the file.
                If only_data=True, header and footer is ignored.
            comments (str, optional): String that will be prepended to the
                header and footer strings, to mark them as comments. Default is "# ".

        Returns:
            None

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        ##############################
        # check for wrong parameters #
        ##############################
        if 'folderpath' in kwargs:
            raise ValueError('folderpath is not a valid parameter.\nPlease, use filepath')
        
        ######################################
        # check if filepath points to a file #
        ######################################
        filepath = Path(filepath)
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        ###################
        # check overwrite #
        ###################
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if other.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        ##########
        # kwargs #
        ##########
        if 'fmt' not in kwargs: # pick best format
            if len(self._x) > 0:
                decimal = max([numanip.n_decimal_places(x) for x in arraymanip.flatten(self.data)])
                kwargs['fmt'] = f'%.{decimal}f'
        kwargs.setdefault('delimiter', ', ')
        kwargs.setdefault('newline', '\n')
        kwargs.setdefault('comments', '# ')

        #####################
        # header and footer #
        #####################
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
        else:
            attrs_dict = {_:self.__getattribute__(_) for _ in settings._reserved_words['PhotonEvents']['vars'] if _ not in ['_x', '_y']}
            attrs_dict.update(self.get_attrs_dict())
            header = '\n'.join(_attr2str(attrs_dict, verbose)) + '\n'

            if 'header' not in kwargs:
                kwargs['header'] = header
            else:
                if kwargs['header'] == '':
                    kwargs['header'] = header
                elif kwargs['header'][-1] != '\n':
                    kwargs['header'] += '\n'
        
        ########
        # save #
        ########
        if self.x is None:
            np.savetxt(Path(filepath), [], **kwargs)
        else:
            np.savetxt(Path(filepath), self.data, **kwargs)
        return 
    
    def load(self, filepath, only_data=False, verbose=False, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first 
                decompressed. Last used filepath is saved to im.filepath.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is loaded.
            verbose (book, optional): Default is False. If True, it will print
                an warning when attributes cannot be loaded from the file.
            **kwargs: kwargs are passed to ``plt.genfromtxt()`` that loads the data.

        If not specified, the following parameters are passed to `numpy.genfromtxt()`_:

        Args:
            delimiter (str, optional): String or character separating columns.
                Use ``\\t`` for tab. Default is comma (', ').
            comments (str, optional): The character used to indicate the start
                of a comment. Default is ``# ``. Attributes picked up
                from the header will be loaded too.

        Returns:
            br.PhotonEvents

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        ############
        # filepath #
        ############
        filepath = Path(filepath)
        assert filepath.is_file(), f'filepath must point to a file, {filepath}'
        assert filepath.exists(),  f'filepath does not exist, {filepath}'

        ##########
        # kwargs #
        ##########
        kwargs.setdefault('delimiter', ', ')
        kwargs.setdefault('comments', '# ')

        ########
        # read #
        ########
        data = np.genfromtxt(Path(filepath), **kwargs)

        ##########
        # assign #
        ##########
        pe = PhotonEvents()
        if len(data) != 0:
            pe._x = data[:, 0]
            pe._y = data[:, 1]
        else:
            pe._x = None
            pe._x = None
        
        ##########################
        # reset check attributes #
        ##########################
        # pe._has_nan      = None

        ###############
        # read header #
        ###############
        if only_data is False:
            # get header
            header = filemanip.load_comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            
            # remove comment flag (#)
            comment_flag_len = len(kwargs['comments'])
            for i, line in enumerate(header):
                header[i] = line[comment_flag_len:]

            # attrs dict
            attrs_dict = _str2attr(header[:-1], verbose=verbose)

            # set attrs
            for attr in attrs_dict:
                pe.__setattr__(attr, attrs_dict[attr])
        return pe
    
    #########
    # check #
    #########
    def check_nan(self):
        """Check if x or y have non-numeric (NaN) values

        Result (True or False) is stored in im.has_nan attribute

        Returns:
            None
        """
        ################
        # empty object #
        ################
        if self.x is None:
            self._has_nan = False
            return
        
        #############
        # check nan #
        #############
        if np.isnan(self.x).any():
            self._has_nan = True
        elif np.isnan(self.y).any():
            self._has_nan = True
        else:
            self._has_nan = False
        return None
        
    def remove_nan(self):
        """remove data points (x, y) that contain non-numeric (NaN) values
        
        Returns:
            Spectrum
        """
        if self.has_nan is None:
            self.check_nan()

        if self.has_nan:
            return self.copy()
        else:
            index2remove = list(np.argwhere(np.isnan(self.x))) + list(np.argwhere(np.isnan(self.y)))
            x = np.delete(self.x, index2remove)
            y = np.delete(self.y, index2remove)

            # copy
            pe = self.copy()
            pe._x = x
            pe._y = y
            pe._has_nan = False
            return pe
    
    #############
    # modifiers #
    #############
    def set_horizontal_shift(self, value):
        """Shift photon events horizontally (x).

        Args:
            value (number or list): shift value by which the data are 
                shifted. If list, then it must be of the same size as the number of
                 photon events. First element will be assigned to the first photon event and so on.

        Warning:
            pe.xlim is set to None.

        Returns:
            :py:class:`PhotonEvents`
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if isinstance(value, Iterable):
            assert len(value) == len(self), f'value length ({len(value)}) must have the same length as the data\nLength of the data = {len(self)}.'
        else:
            value = [value]*len(self)

        #########
        # shift #
        #########
        pe = self.copy()
        pe.xlim = None
        for i, v in enumerate(value):
            if v != 0:
                pe.x[i] += v                   
        return pe
    
    def set_vertical_shift(self, value):
        """Shift photon events vertically (y).

        Args:
            value (number or list): shift value by which the data are 
                shifted. If list, then it must be of the same size as the number of
                 photon events. First element will be assigned to the first photon event and so on.

        Warning:
            pe.ylim is set to None.

        Returns:
            :py:class:`PhotonEvents`
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if isinstance(value, Iterable):
            assert len(value) == len(self), f'value length ({len(value)}) must have the same length as the data\nLength of the data = {len(self)}.'
        else:
            value = [value]*len(self)

        #########
        # shift #
        #########
        pe = self.copy()
        pe.ylim = None
        for i, v in enumerate(value):
            if v != 0:
                pe.y[i] += v
        return pe

    ###############
    # modifiers 2 #
    ###############
    def set_horizontal_shift_via_polyval(self, p):
        """Set horizontal shift values to np.polyval(p, y).

        Args:
            p (array): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term.

        Warning:
            pe.xlim is set to None.

        Returns:
            :py:class:`Image`
        """
        f = lambda y: np.polyval(p, y)
        return self.set_horizontal_shift_via_function(f)
    
    def set_vertical_shift_via_polyval(self, p):
        """Set vertical shift values to np.polyval(p, x).

        Args:
            p (array): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term.
            
        Warning:
            pe.ylim is set to None.

        Returns:
            :py:class:`Image`
        """
        f = lambda x: np.polyval(p, x)
        return self.set_vertical_shift_via_function(f)

    def set_vertical_shift_via_function(self, f):
        """Set vertical shift values to f(x).

        Args:
            f (function): function where argument is x coordinates

        Warning:
            pe.ylim is set to None.

        Returns:
            :py:class:`Image`
        """
        value = np.array([f(x) for x in self.x])
        return self.set_vertical_shift(value=value)

    def set_horizontal_shift_via_function(self, f):
        """Set horizontal shift values to f(y).

        Args:
            f (function): function where argument is y coordinates

        Warning:
            pe.xlim is set to None.

        Returns:
            :py:class:`Image`
        """
        value = np.array([f(y) for y in self.y])
        return self.set_horizontal_shift(value=value)

    ############
    # advanced #
    ############
    def crop(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Crop photon events out.

        Args:
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.

        Returns:
            PhotonEvents
        """
        return self.copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)

    def transpose(self):
        """Switch x and y positions.

        Returns:
            PhotonEvents
        """
        # create new PhotonEvents with switched x and y
        pe = PhotonEvents(x=self.y, y=self.x)
        pe._copy_attrs_from(self)
        return pe

    def binning(self, ncols, nrows):
        """Compute the 2D histogram of the data (binning of the data).

        Args:
            ncols, nrows (int or None): number of columns and rows of the returned Image.

        Returns:
            :py:class:`Image` binned image
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if numanip.is_integer(ncols) == False or numanip.is_integer(nrows) == False or ncols < 0 or nrows < 0:
            raise ValueError("Number of bins must be a positive integer.")
        ncols = int(ncols)
        nrows = int(nrows)

        ###########
        # binning #
        ###########
        temp, _x_edges, _y_edges = np.histogram2d(self.x, self.y, bins=(ncols, nrows), range=(self.xlim, self.ylim))

        #########
        # Image #
        #########
        im = Image(temp.transpose())
        im.x_edges = _x_edges
        im.y_edges = _y_edges
        im._copy_attrs_from(self)

        return im

    ########################
    # calculation and info #
    ########################
    def integrated_rows_vs_y_centers(self, nrows):
        """Integrate data in the horizontal direction.

        Args:
            nrows (int, optional): number of y bins.

        Returns:
            :py:class:`Spectrum`.
        """
        im = self.binning(ncols=1, nrows=nrows)
        s  = im.integrated_rows_vs_y_centers()
        s._copy_attrs_from(self)
        return s
    
    def integrated_columns_vs_x_centers(self, ncols):
        """Integrate data in the vertical direction.

        Args:
            ncols (int, optional): number of x bins.

        Returns:
            :py:class:`Spectrum`.
        """
        im = self.binning(ncols=ncols, nrows=1)
        s  = im.integrated_columns_vs_x_centers()
        s._copy_attrs_from(self)
        return s

    def calculate_spectrum(self, nbins, axis=1):
        """Integrate data in one direction (sum columns or rows).

        Args:
            nbins (int, optional): number of bins.
            axis (int or string, optional): Axis along which elements are integrated.
                If axis = 0, spectra will be integrated vertically (each spectrum 
                datapoint is the sum of a "column of pixels"). If axis = 1, 
                spectra will be integrated
                horizontally (each spectrum datapoint is a "row"). 
                Default is 1.

        Returns:
            :py:class:`Spectrum`.
        """
        # check axis
        if axis != 0 and axis != 1:
            raise ValueError('axis must be 0 or 1')
        
        if axis == 0:
            return self.integrated_columns_vs_x_centers(ncols=nbins)
        else: 
            return self.integrated_rows_vs_y_centers(nrows=nbins)

    ############
    # composed #
    ############
    def calculate_vertical_shift(self, ncols, nrows, mode='cc', ylimits=None, limit_size=1000, **kwargs):
        """Calculate intensity misalignments in terms of y centers.

        Args:
            ncols, nrows (int or None): number of columns and rows for binning.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align columns via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous column.

                 3) 'max': Align the max point of every column. 
                 
                 4) 'peak': Fit one peak in each column and align them 
                 (requires that `brixs.addons.fitting` is imported)

            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'` 

        Returns:
            list
        """
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.calculate_vertical_shift(mode=mode, ylimits=ylimits, limit_size=limit_size, **kwargs)

    def calculate_horizontal_shift(self, ncols, nrows, mode='cc', xlimits=None, limit_size=1000, **kwargs):
        """Calculate intensity misalignments in terms of x centers.

        Args:
            ncols, nrows (int or None): number of columns and rows for binning.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align rows via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous row.

                 3) 'max': Align the max point of every row. 
                 
                 4) 'peak': Fit one peak in each row and align them 
                 (requires that `brixs.addons.fitting` is imported)

            xlimits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'` 
            
        Returns:
            list
        """
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.calculate_horizontal_shift(mode=mode, xlimits=xlimits, limit_size=limit_size, **kwargs)


    def calculate_vertical_shift_curvature(self, ncols, nrows, deg=2, mode='cc', ylimits=None, limit_size=1000, **kwargs):
        """Calculate vertical shift values to fix curvature.
        
        Note:
            This function calculates the SHIFT values, and not the polynomial 
            curve that fits the curvature. The events will be corrected uppon
            shifting events by these shift values. 

        Args:
            ncols, nrows (int or None): number of columns and rows for binning.
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align columns via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous column.

                 3) 'max': Align the max point of every column. 
                 
                 4) 'peak': Fit one peak in each column and align them 
                 (requires that `brixs.addons.fitting` is imported)

            ylimits (None or list): y center pair of values `(y_start, y_stop)`, a list 
                of pairs `((yi_1, yf_1), (yi_2, yf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `y_start = None` or `y_stop = None` to indicate the minimum or 
                maximum y value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            Spectrum: spectrum with shift values vs x centers with the following
                attributes:
            
                    s.fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

                    s.popt (np.array): 1D array of polynomial coefficients 
                        (including coefficients equal to zero) from highest degree to 
                        the constant term.

                    s.R2 (number): R2 error

                    s.model (function): funcion f(x_centers)
        """
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.calculate_vertical_shift_curvature(deg=deg, mode=mode, ylimits=ylimits, limit_size=limit_size, **kwargs)
     
    def calculate_horizontal_shift_curvature(self, ncols, nrows, deg=2, mode='cc', xlimits=None, limit_size=1000, **kwargs):
        """Calculate horizontal shift values to fix curvature.
        
        Note:
            This function calculates the SHIFT values, and not the polynomial 
            curve that fits the curvature. The events will be corrected uppon
            shifting events by these shift values. 

        Args:
            ncols, nrows (int or None): number of columns and rows for binning.
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align rows via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous row.

                 3) 'max': Align the max point of every row. 
                 
                 4) 'peak': Fit one peak in each row and align them 
                 (requires that `brixs.addons.fitting` is imported)

            xlimits (None or list): x center pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            Spectrum: spectrum with shift values vs y centers with the following 
                attributes:
            
                    s.fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

                    s.popt (np.array): 1D array of polynomial coefficients 
                        (including coefficients equal to zero) from highest degree to 
                        the constant term.

                    s.R2 (number): R2 error

                    s.model (function): funcion f(y_centers)
        """
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.calculate_horizontal_shift_curvature(deg=deg, mode=mode, xlimits=xlimits, limit_size=limit_size, **kwargs)
     

    def fix_vertical_shift_curvature(self, ncols, nrows, deg=2, mode='cc', ylimits=None, limit_size=1000, **kwargs):
        """shift photon events vertically to fix curvature.

        Args:
            ncols, nrows (int or None): number of columns and rows for binning.
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align columns via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous column.

                 3) 'max': Align the max point of every column. 
                 
                 4) 'peak': Fit one peak in each column and align them 
                 (requires that `brixs.addons.fitting` is imported)

            limits (None or list): y center pair of values `(y_start, y_stop)`, a list 
                of pairs `((yi_1, yf_1), (yi_2, yf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `y_start = None` or `y_stop = None` to indicate the minimum or 
                maximum y value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            PhotonEvents
        """   
        _result = self.calculate_vertical_shift_curvature(ncols=ncols, nrows=nrows, deg=deg, mode=mode, ylimits=ylimits, limit_size=limit_size)
        return self.set_vertical_shift_via_polyval(p=_result['popt'])

    def fix_horizontal_shift_curvature(self, ncols, nrows, deg=2, mode='cc', xlimits=None, limit_size=1000, **kwargs):
        """shift photon events horizontally to fix curvature.

        Args:
            ncols, nrows (int or None): number of columns and rows for binning.
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align rows via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous row.

                 3) 'max': Align the max point of every row. 
                 
                 4) 'peak': Fit one peak in each row and align them 
                 (requires that `brixs.addons.fitting` is imported)

            limits (None or list): x center pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.
            **kwargs (dict)
                kwargs to be passed to ss.fit_peak() function when `mode='peak'`

        Returns:
            PhotonEvents
        """   
        _result = self.calculate_horizontal_shift_curvature(ncols=ncols, nrows=nrows, deg=deg, mode=mode, xlimits=xlimits, limit_size=limit_size)
        return self.set_horizontal_shift_via_polyval(p=_result['popt'])

    ##########################        
    # plot and visualization #
    ##########################     
    def plot(self, ax=None, s=0.1, show_limits=False, x_start=None, x_stop=None, y_start=None, y_stop=None, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.scatter()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            s (number, optional): The marker size in points**2. Default is 0.1.
            show_limits (bool, optional): if True and pe.xlim and pe.ylim are
                defined, it will draw a grey rectangle around xlim and ylim.
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.
            **kwargs: kwargs are passed to ``plt.scatter()`` that plots the data.            

        Returns:
            `matplotlib.image.AxesImage`_

        .. _matplotlib.pyplot.scatter(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.scatter.html
        .. _matplotlib.image.AxesImage: https://matplotlib.org/3.5.0/api/image_api.html#matplotlib.image.AxesImage
        """
        ###################
        # figure and axes #
        ###################
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW or len(plt.get_fignums()) == 0:
                figure()

        ########
        # plot #
        ########
        pe  = self._copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
        pos = ax.scatter(pe.x, pe.y, s=s, **kwargs)

        #############
        # check lim #
        #############
        if show_limits:
            if self.xlim is not None and self.ylim is not None:
                if hasattr(ax, 'gca'):
                    figmanip.rectangle(xlim=self.xlim, ylim=self.ylim, edgecolor='grey', lw=1)
                    ax.gca().autoscale_view()
                else:
                    figmanip.rectangle(ax=ax, xlim=self.xlim, ylim=self.ylim, edgecolor='grey', lw=1)
                    ax.autoscale_view()
        return pos

# %% =============================== Dummy =============================== %% #
class Dummy(_BrixsObject, metaclass=_Meta):
    
    def __init__(self, data=None):
        self._data  = []
        if data is not None:
            self.data = data

        # extra
        for extra in settings._init['Dummy']:
            self.__setattr__(extra, settings._init['Dummy'][extra](self))


    ###################
    # core attributes #
    ###################
    @property
    def data(self):
        return self._data
    @data.setter
    def data(self, value):
        if value is None:
            self._data = []
        elif isinstance(value, Iterable):
            self._data = value
        else:
            raise ValueError('value must be an iterable')
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    ###################################
    # computed (read-only) attributes #
    ###################################
    pass

    #########################
    # write-only attributes #
    #########################
    pass

    #######################
    # modifier attributes #
    #######################
    pass

    #################
    # magic methods #
    #################
    # def __setattr__(self, name, value):
    #     super().__setattr__(name, value)

    # def __getattr__(self, name):
    #     super().__getattr__(name)

    def __getitem__(self, item):
        return self._data[item]

    def __setitem__(self, item, value):
        self._data[item] = value

    def __len__(self):
        return len(self.data)

    def __delitem__(self, item):
        del self._data[item]

    #########
    # attrs #
    #########
    pass

    ###########
    # support #
    ###########
    pass

    ################
    # core methods #
    ################
    def append(self, value):
        """Append something to dummy.

        Args:
            value: something to append to Dummy

        Returns:
            None

        See Also:
            :py:func:`Dummy.remove`
        """        
        self.data.append(value)

    ########
    # copy #
    ########
    def _copy(self):
        """Same as self.copy(), but attrs are NOT copied."""
        #############
        # copy data #
        #############
        data = copy.deepcopy(self.data)
        return Dummy(data=data)

    def copy(self):
        """Return a copy of the object.

        Usage:
            >>> # full copy
            >>> dummy2 = dummy1.copy()  # dummy2 is now a copy of dummy1

        Returns:
            :py:attr:`Dummy`
        """
        dummy = self._copy()
        dummy._copy_attrs_from(self)

        # extra
        for extra in settings._copy['Dummy']:
            if hasattr(dummy, extra):
                dummy.__setattr__(extra, self.__getattribute__(extra).copy(dummy))

        return dummy
    

    #################
    # save and load #
    ################# 
    pass

    #########
    # check #
    #########
    pass

    #############
    # modifiers #
    #############
    def set_shift(self, value):
        """Calls set_shift for each object inside br.Dummy

        Args:
            value (number or list): value will be passed to set_shift. If list,
                number of values must be the same as the length of br.Dummy

        Returns:
            :py:class:`Dummy`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the lenght of Dummy.\nnumber of values: {len(value)}\nnumber of objects: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_shift(value=value[i])

        return ss

    def set_calib(self, value):
        """Calls set_calib for each object inside br.Dummy

        Args:
            value (number or list): value will be passed to set_calib. If list,
                number of values must be the same as the length of br.Dummy

        Returns:
            :py:class:`Dummy`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)
        
        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the lenght of Dummy.\nnumber of values: {len(value)}\nnumber of objects: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_calib(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        # extra
        for extra in settings._calib['Spectra']:
            if hasattr(s, extra):
                ss.__getattribute__(extra).set_shift(value)
        
        return ss

    def set_factor(self, value):
        """Calls set_factor for each object inside br.Dummy

        Args:
            value (number or list): value will be passed to set_calib. If list,
                number of values must be the same as the length of br.Dummy

        Returns:
            :py:class:`Dummy`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the lenght of Dummy.\nnumber of values: {len(value)}\nnumber of objects: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_factor(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        # ss._x            = None
        # ss._monotonicity = None

        # extra
        for extra in settings._factor['Spectra']:
            if hasattr(s, extra):
                ss.__getattribute__(extra).set_shift(value)

        return ss

    def set_offset(self, value):
        """Calls set_offset for each object inside br.Dummy

        Args:
            value (number or list): value will be passed to set_calib. If list,
                number of values must be the same as the length of br.Dummy

        Returns:
            :py:class:`Dummy`
        """
        ##############################
        # check if value is a number #
        ##############################
        if isinstance(value, Iterable) == False:
            value = [value]*len(self)

        ##################################
        # value must be the right length #
        ##################################
        assert len(value) == len(self), f'value must have the same number of items as the lenght of Dummy.\nnumber of values: {len(value)}\nnumber of objects: {len(self)}'

        ##############
        # set values #
        ##############
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.set_offset(value=value[i])

        return ss


    ###############
    # modifiers 2 #
    ###############
    pass

    ############
    # advanced #
    ############
    pass

    ########################
    # calculation and info #
    ########################
    def calculate_sum(self):
        """Tries to return object which is the sum of all objects in the list"""
        if isinstance(self[0], Image):
            im = Image(data=np.sum([_.data for _ in self], axis=0))
            im._copy_attrs_from(self[0])
            return im
        else:
            obj = self[0]
            for _obj in self[1:]:
                obj += _obj
            return obj

    def calculate_average(self, limits=None):
        """Tries to return object which is the average of all objects in the list"""
        return self.calculate_sum().set_factor(1/len(self))



    ############
    # composed #
    ############
    pass

    ##########################        
    # plot and visualization #
    ##########################
    def plot(self, ax=None, **kwargs):
        """Plot items one by one. Wrapper for `matplotlib.pyplot.plot()`_.

        Note:
            If `label` is `None` and if spectrum inside spectra have attr 
            `label`, this attr will be used as label, e.g., 
            `plt.plot(s.x, s.y, label=s.label)`.  

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            smooth (int, optional): number of points to average data. Default is 1.
            label or labels (str, number, or list, optional): if str or number, this label will be 
                applied to every spectra. If list, it must have the same length 
                as the number of spectra. If None and if 
                spectrum `s` inside spectra `ss` have attr `label`, 
                this attr will be used as label, e.g., `plt.plot(s.x, s.y, label=s.label)`.
                Default is None. 
            color, colors, or c (str, number, or list, optional): if str or number, this color will be 
                applied to every spectra. If list, it must have the same length 
                as the number of spectra and each element must be a color (a 
                color can be a str or a 3 element (RGB) list.
                Default is None. 
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Limits are
                inclusive. Use `x_start = None` or `x_stop = None` to indicate 
                the minimum or maximum x value of the data, respectively. If 
                limits = [], i.e., an empty list, it assumes `limits = (None, None)`.
            switch_xy (bool, optional): Switch x and y axis.
            hi, vi (number, optional): horizontal and vertical increments for 
                cascading plots.
            phi, pvi (number, optional): percentage wise horizontal and vertical 
                increments for cascading plots (percentage of the y-range for 
                each spectrum).
            verbose (bool, optional): if True, prints warning if ploted data has
                nun-numeric values (NaN). Default is True.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            `Line2D`_ list

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        ###################
        # figure and axes #
        ###################
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW or len(plt.get_fignums()) == 0:
                figure()
        
        ########
        # plot #
        ########
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].plot(ax=ax, **kwargs)

        return temp
