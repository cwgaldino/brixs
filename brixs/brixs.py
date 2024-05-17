#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core brixs module. Defines the main objects: Spectrum, Spectra, PhotonEvents, Image"""

# %% ------------------------- Standard Imports --------------------------- %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import warnings
import copy

# %% -------------------------- Special Imports --------------------------- %% #
from collections.abc import Iterable, MutableMapping
from mpl_toolkits.axes_grid1 import make_axes_locatable

# %% ----------------------------- backpack ------------------------------- %% #
from .backpack import filemanip, arraymanip, figmanip, numanip, query

# %% ------------------------------ settings ------------------------------ %% #
from .config import settings
# %%

# %% ============================= support class ========================== %% #
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
# %%

# %% ========================= common support functions =================== %% #
def _attr2str(s, attrs, verbose):
    """returns a list with strings for each attr and attr value
    
        Warning:
            eval() must be able to run the attr value string in order for the attr value
            to be later readable by load() functions
        
        Args:
            s (brixs object): Spectrum, Spectra, Image, PhotonEvents
            attrs (list): list of attr names
            verbose (bool): if True, message is printed when attr cannot be 
                converted to string.
        
        Returns:
            list
    """
    final = []

    #######################
    # collect attr values #
    #######################
    for name in attrs:
        try:
            ##############
            # type: None #
            ##############
            if s.__dict__[name] is None:
                final.append(f'{name}: None')
            #############
            # type: str #
            #############
            elif isinstance(s.__dict__[name], str):
                temp2 = str(s.__dict__[name]).replace('\n','\\n')
                final.append(f'{name}: \"{temp2}\"')
            ##############
            # type: dict #
            ##############
            elif isinstance(s.__dict__[name], dict) or isinstance(s.__dict__[name], MutableMapping):
                if verbose:
                    type_ = str(type(s.__dict__[name]))
                    print('\nWarning: Cannot save attr of type: ' + type_ + '\nattr name: '+ name + '\nTo turn off this warning, set verbose to False.')
            ########################
            # type: list and tuple #
            ########################
            elif isinstance(s.__dict__[name], Iterable):
                final.append(f'{name}: {list(s.__dict__[name])}')
            ################
            # type: number #
            ################
            elif numanip.is_number(s.__dict__[name]):
                tosave = str(s.__dict__[name])
                if tosave[-1] == '\n':
                    tosave = tosave[:-1]
                final.append(f'{name}: {tosave}')
            else:
                temp2 = str(s.__dict__[name]).replace('\n','\\n')
                final.append(f'{name}: \"{temp2}\"')
        except:
            if verbose:
                type_ = str(type(s.__dict__[name]))
                print('\nWarning: Cannot save attr of type: ' + type_ + '\nattr name: '+ name + '\nTo turn off this warning, set verbose to False.')
    return final
# %%

# %% ====================== modified figure function ======================= %% #
def figure(*args, **kwargs):
    """Create figure object. Wrapper for `plt.figure()`_.

    The following br.settings affect figure:

        br.settings.FIGURE_POSITION
        br.settings.FIGURE_FORCE_ON_TOP
        br.settings.FIGURE_DPI
        br.settings.FIGURE_SIZE
        br.settings.FIGURE_GRID


    Mouse click behavior:

        Right click:
            x value is copied to the clipboard.
        Left click OR (y + Right click):
            y value is copied to the clipboard.
        Middle click:
            copies cursor position in terms of figure coordinates.
    Args:
        *args, **kwargs: args and kwargs are passed to `plt.figure()`.

    Note:
        This function overwrites the behavior of `figsize` parameters. In
        plt.figure(figsize=(w, h)), w and h must be given in inches. However,
        this function gets `w` and `h` in cm. 
    
    Returns:
        figure object
    
    .. _plt.figure(): https://matplotlib.org/stable/api/figure_api.html
    """
    fig = figmanip.figure(*args, **kwargs)

    ############
    # position #
    ############
    if settings.FIGURE_POSITION is not None:
        figmanip.set_window_position(settings.FIGURE_POSITION)

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
    """Create a figure and a set of subplots. Wrapper for `plt.subplots()`_.

    The following br.settings affect figure:

        br.settings.FIGURE_POSITION
        br.settings.FIGURE_FORCE_ON_TOP
        br.settings.FIGURE_DPI
        br.settings.FIGURE_SIZE
        br.settings.FIGURE_GRID


    Mouse click behavior:

        Right click:
            x value is copied to the clipboard.
        Left click OR (y + Right click):
            y value is copied to the clipboard.
        Middle click:
            copies cursor position in terms of figure coordinates.
    
    Args:
        nrows, ncols (int): Number of rows/columns of the subplot grid.
        sharex, sharey (bool or str, optional): Share the x or y axis between 
            axes. Options are: True, False, 'all', 'row', 'col'. Default is False
        hspace, wspace (float, optional): The amount of height/width reserved for space 
            between subplots, expressed as a fraction of the average axis 
            height/width. Default is 0.3.
        width_ratios, height_ratios (list, optional): Defines the relative heights/widths of the 
            ows/columns. Each row/column gets a relative height/width of 
            ratios[i] / sum(ratios). If not given, all rows/columns will have 
            the same width. 
        **fig_kw: All additional keyword arguments are passed to the plt.figure call


    Note:
        This function overwrites the behavior of `figsize` parameters. In
        plt.figure(figsize=(w, h)), w and h must be given in inches. However,
        this function gets `w` and `h` in cm. 
    
    Returns:
        fig, axes
    
    .. _plt.subplots(): https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.subplots.html
    """
    fig, axes = figmanip.subplots(*args, **kwargs)

    ############
    # position #
    ############
    if settings.FIGURE_POSITION is not None:
        figmanip.set_window_position(settings.FIGURE_POSITION)

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

# %% ============================== Spectrum ============================== %% #
class Spectrum(metaclass=_Meta):
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
            >>> s = br.Spectrum(filepath=<filepath>)
            >>> s = br.Spectrum(filepath=<filepath>, delimiter=',')
            >>>
            >>> print(s.get_core_attrs()) # print list of core attrs
            >>> print(s.get_attrs())      # print list of attrs
            >>> print(s.get_methods())    # print list of methods available
    """
    _read_only     = ['step', 'monotonicity']
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

        # modifiers
        self._calib  = 1
        self._factor = 1
        self._shift  = 0
        self._offset = 0

        # labels
        pass

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
       
        #################
        # set attribute #
        #################
        self._x = np.array(value, dtype='float')

        ##########################
        # reset check attributes #
        ##########################
        self._step         = None
        self._monotonicity = None
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
            self._x = np.arange(0, len(value))

        #################
        # set attribute #
        #################
        self._y = np.array(value, dtype='float')
    @y.deleter
    def y(self):
        raise AttributeError('Cannot delete object.')

    ###################################
    # computed (read-only) attributes #
    ###################################
    @property
    def data(self):
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
    @property
    def calib(self):
        return self._calib
    @calib.setter
    def calib(self, value):
        _s = self.set_calib(1/self._calib).set_calib(value)
        self.x      = _s.x
        self._calib = value
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

    def __iter__(self):
        # this is what makes possible max(s) and min(s)
        for y in self.y:
            yield y

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
        else:
            raise TypeError('Index must be int or a slice, not {}'.format(type(item).__name__))
        
    #########
    # attrs #
    #########
    def get_core_attrs(self): 
        """return a list of core attrs"""
        return settings._reserved_words['Spectrum']['pseudovars']
    
    def get_attrs(self):
        """return a list of user defined attrs"""
        return [key for key in self.__dict__.keys() if key.startswith('_') == False and key not in settings._reserved_words]

    def get_methods(self):
        """return a list of methods available"""
        return [key for key in self.__dir__() if key.startswith('_') == False and key not in self.get_attrs() + self.get_core_attrs()]
    
    def remove_attrs(self):
        """Delete all user defined attrs."""
        for attr in [key for key in self.__dict__.keys() if key.startswith('_') == False and key not in settings._reserved_words]:
            self.__delattr__(attr)

    def copy_attrs_from(self, s):
        """Copy user defined attributes from another brixs object.

        Args:
            s (brixs object): Either a Spectrum, Spectra, Image, or PhotonEvents
                to copy user defined attributes from.
        
        Returns:
            None
        """
        # check type
        if isinstance(s, Spectrum) or isinstance(s, Spectra) or isinstance(s, Image) or isinstance(s, PhotonEvents):
            pass
        else:
            raise TypeError(f'type {type(s)} not valid\nCan only copy user attrs from type br.Spectrum, br.Spectra, br.Image, or br.PhotonEvents')

        # transfer attrs
        for attr in s.get_attrs():
            value = copy.deepcopy(s.__dict__[attr])
            self.__setattr__(attr, value)

    ###########
    # attrs 2 #
    ###########
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            None or limits in the following format:
                ((xi_1, xf_1), (xi_2, xf_2), ...)               
        """
        #####################
        # if limits is None #
        #####################
        if limits is None:
            return None
        
        ######################
        # if object is empty #
        ######################
        if len(self) == 0:
            return None
        
        ##################################
        # assert that limits is Iterable #
        ##################################
        assert isinstance(limits, Iterable), f'`limits` must be an Iterable, not {type(limits)}'
        
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
        final = ()
        # one pair
        if len(limits) == 1: # ((xi, xf), )
            assert isinstance(limits[0], Iterable), f'wrong format for limits={limits}'
        # two pairs 
        elif len(limits) == 2: # (xi, xf), or ((xi1, xf1), (xi2, xf2))
            if isinstance(limits[0], Iterable) == False:
                if isinstance(limits[1], Iterable) == False:
                    final.append(limits)
                else:
                    raise ValueError(f'wrong format for limits={limits}')
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
    pass

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
            return Spectrum(x=x, y=y)
        
        ########################################
        # check if extract is really necessary #
        ########################################
        if len(limits) == 1:
            if limits[0][0] <= min(self.x) and limits[0][1] >= max(self.x):
                return Spectrum(x=x, y=y)
        
        ###########################
        # check if empty spectrum #
        ###########################
        if len(self) == 0:
            print(f'Warning: No datapoints within limits={limits}')
            return Spectrum()
        
        ###########
        # extract #
        ###########
        try:
            x, y = arraymanip.extract(x, y, limits)
            return Spectrum(x=x, y=y)
        except RuntimeError:
            print(f'Warning: No datapoints within limits={limits}')
            return Spectrum()

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectrum`
        """
        s = self._copy(limits=limits)

        ##################
        # transfer attrs #
        ##################
        s.copy_attrs_from(self)
        if limits == None: 
            s._step         = self.step
            s._monotonicity = self.monotonicity
        s._shift        = self.shift
        s._calib        = self.calib
        s._offset       = self.offset
        s._factor       = self.factor

        return s

    #################
    # save and load #
    #################     
    def save(self, filepath, only_data=False, check_overwrite=False, verbose=False, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Warning:
            Attrs are saved as comments if only_data is False. Saving attrs to file
            is not always reliable because requires converting variables to string. 
            Only attrs that are of type: string, number, and list of number and strings are 
            saved somewhat correctly. Dictionaries are not saved. 

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
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
                    if query.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        ##########
        # kwargs #
        ##########
        if 'fmt' not in kwargs: # pick best format
            decimal = max([numanip.n_decimal_places(x) for x in arraymanip.flatten(self.data)])
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        #####################
        # header and footer #
        #####################
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
        else:
            attrs  = ['_step', '_monotonicity']
            attrs += ['_shift', '_factor', '_calib', '_offset']
            attrs += self.get_attrs()
            header = '\n'.join(_attr2str(self, attrs, verbose)) + '\n'

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
        np.savetxt(Path(filepath), self.data, **kwargs)

    def load(self, filepath, only_data=False, verbose=False, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

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
            None

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
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '#'
        if 'usecols' not in kwargs:
            kwargs['usecols'] = (0, 1)

        ########
        # read #
        ########
        data = np.genfromtxt(Path(filepath), **kwargs)
        
        ##############
        # check data #
        ##############
        x = data[:, 0]
        y = data[:, 1]
        assert len(x) == len(y), f'Length of x array (len={len(x)}) is not compatible with y array (len={len(y)}).'
        
        ##########
        # assign #
        ##########
        self._x = x
        self._y = y

        ##########################
        # reset check attributes #
        ##########################
        self._step         = None
        self._monotonicity = None

        ###################
        # reset modifiers #
        ###################
        self._calib  = 1
        self._factor = 1
        self._shift  = 0
        self._offset = 0

        ###############
        # read header #
        ###############
        if only_data is False:
            header = filemanip.load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            if header:
                for line in header:
                    if ':' not in line:
                        pass
                    else:
                        # extract name and value
                        name = line[1:-1].split(':')[0].strip()
                        try:
                            value = eval(str(':'.join(line[1:-1].split(':')[1:])).strip())
                        except:
                            value = str(':'.join(line[1:-1].split(':')[1:])).strip()
                        try:
                            setattr(self, name, value)
                        except Exception as e:
                            if verbose:
                                print(f'cannot read attr ({name}: {value})\nAttribute not set\n{e}\n')

    #########
    # check #
    #########
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

            Returns:
                None

            Raises:
                ValueError: If x-coordinates are not uniform.
        
            See Also:
                :py:func:`Spectrum.check_monotonicity`
        """
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

    def check_monotonicity(self):
        """Sets monotonicity attribute to 'increasing' or 'decreasing'.

        Raises:
            ValueError if data is not monotonic.

        Returns:
            None

        See Also:
                :py:func:`Spectrum.check_step`, :py:func:`Spectrum.fix_monotonicity`
        """
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

            Returns:
                None
            
            See Also:
                :py:func:`Spectrum.check_monotonicity`
        """
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

    #############
    # modifiers #
    #############
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
        return s

    ###############
    # modifiers 2 #
    ###############
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

    def set_roll(self, value):
        """Roll array elements for the x-coordinates

        Note:
            Elements that roll beyond the last position are re-introduced at the
             first. y-coordinates are fully preserved.

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectrum`
        """
        if limits is None: limits = (None, None)
        temp  = self._copy(limits=limits)
        assert len(temp.y) > 0, 'no data points within limits = {limits}'
        value = -np.mean(temp.y)
        return self.set_offset(value)
    
    def normalize(self, value=1, limits=None):
        """Set a factor such as the average y between limits is equal to value.

        Args:
            value (number): value. Default is 1.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectrum`
        """
        if limits is None: limits = (None, None)
        temp = self._copy(limits=limits)
        assert len(temp.y) > 0, 'no data points within limits = {limits}'
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
                temp, _ = self._copy(self.x, self.y, (start, stop))
                x = np.linspace(start, stop, num=len(temp))

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
        
        return s

    def crop(self, start=None, stop=None):
        """Crop edges of the dataset.

        Args:
            start (number, optional): start x value. If None, the minimum value of
                x will be used.
            stop (number, optional): final x value. If None, the maximum value of
                x will be used.

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            number

        .. _numpy.trapz(): https://numpy.org/doc/stable/reference/generated/numpy.trapz.html
        """
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
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
            >>> fit, popt, R2, f = polyfit(deg)

        Args:
            deg (int): degree of the fitting polynomial.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
        
        Returns:
            fit (spectrum), popt, R2, f(x)

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
        _x = np.arange(start, stop, len(x)*100)
        arr100 = Spectrum(x, y=model(_x))

        return arr100, popt, R2, model
   
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

    ############
    # composed #
    ############
    pass

    ##########################        
    # plot and visualization #
    ##########################        
    def plot(self, ax=None, offset=0, shift=0, roll=0, factor=1, calib=1, smooth=1, limits=None, switch_xy=False, **kwargs):
        """Plot spectrum. Wrapper for `matplotlib.pyplot.plot()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            calib, shift (number, optional): multiplicative and additive factor
                 on the x-coordinates. calib is applied first.
            factor, offset (number, optional): multiplicative and additive factor
                 on the y-coordinates. Factor is applied first.
            roll (int, optional): Roll value of array elements of the x-coordinates
            smooth (int, optional): number of points to average data. Default is 1.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            switch_xy (bool, optional): Switch x and y axis.
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
            self = self._copy(limits=limits)
        x = self.x
        y = self.y

        #############
        # switch xy #
        #############
        if switch_xy:
            _x = x
            x = y
            y = _x

        #############
        # modifiers #
        #############
        x = (x*calib) + shift
        y = (y*factor) + offset

        ########
        # roll #
        ########
        if roll != 0:
            assert numanip.is_integer(roll), 'roll must be an integer'
            x    = np.roll(x, roll)
            x, y = arraymanip.sort(x, x, y)

        ##########
        # smooth #
        ##########
        if smooth > 1:
            ids = np.arange(len(x))//int(smooth)
            x = np.bincount(ids, x)/np.bincount(ids)

            ids = np.arange(len(y))//int(smooth)
            y = np.bincount(ids, y)/np.bincount(ids)
        
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
        line[0].offset = offset
        line[0].shift  = shift
        line[0].roll   = roll
        line[0].calib  = calib
        line[0].factor = factor
        line[0].smooth = smooth

        return line[0]

# %% =============================== Spectra ============================== %% #
class Spectra(metaclass=_Meta):
    """Returns a ``spectra`` object.

    Args:
        data (list or array, optional): list of :py:class:`spectrum` objects.
        folderpath (string or path, optional): folderpath with spectrum files.
        filepath (str or Path, optional): filepath with multiple spectra to read.
        **kwargs: kwargs are passed to :py:func:`Spectra.load` function or 
            :py:func:`Spectra.load_from_single_file`

    Usage:        
            >>> ss = br.Spectra()
            >>> ss = br.Spectra([s1, s2, ...])
            >>> ss = br.Spectra(data=[s1, s2, ...])
            >>> ss = br.Spectrum(folderpath=<folderpath>)
            >>> ss = br.Spectrum(folderpath=<folderpath>, delimiter=',')
            >>> ss = br.Spectrum(filepath=<filepath>)
            >>> ss = br.Spectrum(filepath=<filepath>, delimiter=',')
            >>>
            >>> print(ss.get_core_attrs()) # print list of core attrs
            >>> print(ss.get_attrs())      # print list of attrs
            >>> print(ss.get_methods())    # print list of methods available
    """
    _read_only     = ['step', 'length', 'x', 'monotonicity']
    _non_removable = []

    def __init__(self, data=None, folderpath=None, filepath=None, **kwargs):
        """Initialize the object instance"""
        ###########################
        # Initializing attributes #
        ###########################
        # core
        self._data  = []

        # check
        self._length       = None
        self._step         = None
        self._x            = None
        self._monotonicity = None

        # modifiers
        pass

        # labels
        pass

        ################
        # loading data #
        ################
        if data is not None:
            self.data = data
        elif folderpath is not None:
            self.load(folderpath=folderpath, **kwargs)
        elif filepath is not None:
            self.load_from_single_file(filepath=filepath, **kwargs)
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
        
        ##########################
        # reset check attributes #
        ##########################
        self._length       = None
        self._step         = None
        self._x            = None
        self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        self._calculated_calib   = None
        self._calculated_factor = None
        self._calculated_offset = None
        self._calculated_shift  = None
        self._calculated_roll   = None
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
        if isinstance(item, int):
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
        
        ##########################
        # reset check attributes #
        ##########################
        self._length       = None
        self._step         = None
        self._x            = None
        self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        self._calculated_calib   = None
        self._calculated_factor = None
        self._calculated_offset = None
        self._calculated_shift  = None
        self._calculated_roll   = None

    def __len__(self):
        return len(self.data)

    def __delitem__(self, item):
        self.remove(item)

    #########
    # attrs #
    #########
    def get_core_attrs(self): 
        """return a list of core attrs"""
        return settings._reserved_words['Spectra']['pseudovars']
    
    def get_attrs(self):
        """return attrs that are user defined.""" 
        return [key for key in self.__dict__.keys() if key.startswith('_') == False and key not in settings._reserved_words]

    def get_methods(self):
        """return a list of methods available"""
        return [key for key in self.__dir__() if key.startswith('_') == False and key not in self.get_attrs() + self.get_core_attrs()]
    
    def remove_attrs(self):
        """Delete all user defined attrs."""
        for attr in self.get_attrs():
            self.__delattr__(attr)

    def copy_attrs_from(self, s):
        """Copy user defined attributes from another brixs object.

        Args:
            s (brixs object): Either a Spectrum, Spectra, Image, or PhotonEvents
                to copy user defined attributes from.
        
        Returns:
            None
        """
        # check type
        if isinstance(s, Spectrum) or isinstance(s, Spectra) or isinstance(s, Image) or isinstance(s, PhotonEvents):
            pass
        else:
            raise TypeError(f'type {type(s)} not valid\nCan only copy user attrs from type br.Spectrum, br.Spectra, br.Image, or br.PhotonEvents')

        # transfer attrs
        for attr in s.get_attrs():
            value = copy.deepcopy(s.__dict__[attr])
            self.__setattr__(attr, value)

    ###########
    # attrs 2 #
    ###########
    def create_attr_from_spectra(self, attr, name=None):
        """create new Spectra attr based on Spectrum attrs.

        Example:
            >>> s1.a = 2
            >>> s2.a = 5
            >>> s3.a = 'g'
            >>> ss = br.Spectra(s1, s2, s3)
            >>> ss.copy_attr_from_spectra('a')
            >>> print(ss.a) # --> [2, 5, 'g']

        Args:
            attr (str): spectrum attribute name. All spectra inside Spectra
                object must have this attr. 
            name (str, optional): name of the new Spectra attr. If None, the same
                name is used.

        Returns:
            None
        """
        # assert all spectra has attr
        # check if all spectra has attr
        has = {i: hasattr(s, attr) for i, s in enumerate(self)}
        if False in [has[i] for i in has]:
            filtered = [i for i in has if has[i] == False]
            raise ValueError(f'Some spectra do not have attr: {attr}.\nSpectra without attr:{filtered}')

        # create new attr
        if name is None: name = attr
        self.__setattr__(name, [getattr(s, attr) for s in self])

    def reorder_by_attr(self, attr, attrs2reorder=None, decreasing=False):
        """Reorder spectra based on a Spectra attr.

        Example:
            >>> ss.a = ['c', 'a', 'e', 'd', 'h', 'i']
            >>> ss.b = [ 3,   1,   5,   4,   8,   9]
            >>> ss.reorder_by_attr(attr='b', attrs2reaorder='a')
            >>> print(ss.b) --> [1, 3, 4, 5, 8, 9]
            >>> print(ss.a) --> ['a', 'c', 'd', 'e', 'h', 'i']
    
        Args:
            attr (str): name of the reference attr. The attr must be a list of 
                numbers with same lenght of number of spectra.
            attrs2reorder (list of str, optional): list of other Spectra attrs that 
                must also be sorted based on the ref attr.
            decreasing (bool, optional): if True, small attr value comes last.
        
        Returns:
            None
        """
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
        self.data = arraymanip.sort(ref, self.data)
        self.__setattr__(attr, arraymanip.sort(ref, ref))
        if attrs2reorder is not None:
            for a in attrs2reorder:
                temp = self.__getattribute__(a)
                self.__setattr__(a, arraymanip.sort(ref, temp))

        # flip order if decreasing is True
        if decreasing:
            self.data = self.data[::-1]
            self.__setattr__(attr, self.__getattribute__(attr)[::-1])
            if attrs2reorder is not None:
                for a in attrs2reorder:
                    temp = self.__getattribute__(a)
                    self.__setattr__(a, self.__getattribute__(a)[::-1])

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

    ###########
    # support #
    ###########
    def _check_limits(self, limits):
        """returns limits in the right format.

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            None or limits in the following format:
                ((xi_1, xf_1), (xi_2, xf_2), ...)               
        """
        #####################
        # if limits is None #
        #####################
        if limits is None:
            return None
        
        ######################
        # if object is empty #
        ######################
        if len(self) == 0:
            return None

        ##################################
        # assert that limits is Iterable #
        ##################################
        assert isinstance(limits, Iterable), f'`limits` must be an Iterable, not {type(limits)}'
        
        ################################
        # get min and max range values #
        ################################
        if self.x is not None:
            vmin = min(self.x)
            vmax = max(self.x)
        else:
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
        final = ()
        # one pair
        if len(limits) == 1: # ((xi, xf), )
            assert isinstance(limits[0], Iterable), f'wrong format for limits={limits}'
        # two pairs 
        elif len(limits) == 2: # (xi, xf), or ((xi1, xf1), (xi2, xf2))
            if isinstance(limits[0], Iterable) == False:
                if isinstance(limits[1], Iterable) == False:
                    final.append(limits)
                else:
                    raise ValueError(f'wrong format for limits={limits}')
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Raises:
            ValueError: if Spectra is empty
        
        Returns:
            x, y's
        """
        ###########
        # check x #
        ###########
        if self.x is None:
            self.check_same_x()

        ######################
        # if object is empty #
        ######################
        if len(self) == 0:
            raise ValueError('cannot operate on empty spectra')
        
        ###################
        # validate limits #
        ###################
        limits = self._check_limits(limits=limits)

        ys = np.zeros((self.length, len(self)))
        for i in range(len(self)):
            ys[:, i] = self[i].y

        try:
            x, ys = arraymanip.extract(self.x, ys, limits=limits)
        except RuntimeError:
            x  = []
            ys = []
            warnings.warn(f'It seems like all spectra has no data points within range: {limits}.\nPlease, fix limits so all spectra have at least one data point within range.')
        
        return x, ys
    
    ################
    # core methods #
    ################
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
        
        ##########################
        # reset check attributes #
        ##########################
        self._length       = None
        self._step         = None
        self._x            = None
        self._monotonicity = None

    def remove(self, idx):
        """Remove spectrum.

        Args:
            idx (int): index of the spectrum.

        Returns:
            None

        See Also:
            :py:func:`Spectra.append`.
        """
        del self._data[idx]

    def reorder(self, i1, i2):
        """reorder spectra.

        Args:
            i1, i2: Index of spectra to switch places.

        Returns:
            None
        """
        temp           = self[i1]
        self._data[i1] = self[i2]
        self._data[i2] = temp
    
    def flip_order(self):
        """reorder spectra backwards:

        If ss.data = [s1, s2, s3], then flip_order() would make it 
        ss.data = [s3, s2, s1]

        Returns:
            None
        """
        self._data = self._data[::-1]

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
            return Spectra(data=data)
        
        ################################
        # get min and max range values #
        ################################
        if self.x is not None:
            vmin = min(self.x)
            vmax = max(self.x)
        else:
            vmin = min(min(s.x) for s in self)
            vmax = max(max(s.x) for s in self)

        ########################################
        # check if extract is really necessary #
        ########################################
        if len(limits) == 1:
            if limits[0][0] <= vmin and limits[0][1] >= vmax:
                return Spectra(data=data)
            
        ##########################
        # check if empty spectra #
        ##########################
        if len(self) == 0:
            print(f'Warning: No spectra to copy')
            return Spectra()

        ###########
        # extract #
        ###########
        # if x is the same for all spectra, this operation is much faster
        if self.x is None:
            try:
                self.check_same_x()
                if self.x is not None:
                    ss = Spectra()
                    x, ys = self._gather_ys(limits=limits)
                    for i in range(len(self)):
                        ss.append(Spectrum(x=x, y=ys[i]))
                    return ss
            except:
                pass
        # if x is not the same, extract data recursively
        ss = Spectra()
        for i, s in enumerate(self):
            # try:
            ss.append(s._copy(limits=limits))
            # except RuntimeError:
            #     raise RuntimeError(f'It seems like spectrum number {i} has no data points within range: {limits}.\nPlease, fix limits (or delete spectrum) so all spectra have at least one data point within range.')
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:attr:`Spectra`
        """
        ss = self._copy(limits=limits)

        ##################
        # transfer attrs #
        ##################
        ss.copy_attrs_from(self)
        if limits == None: 
            ss._x            = self.x
            ss._step         = self.step
            ss._monotonicity = self.monotonicity
            ss._length       = self.length

        return ss
    
    #################
    # save and load #
    #################
    def save(self, folderpath, prefix='spectrum_', suffix='.dat', filenames=None, zfill=None, only_data=False, verbose=False, **kwargs):
        r"""Save spectra. Wrapper for `numpy.savetxt()`_.

        Warning:
            Attrs are saved as comments if only_data is False. Saving attrs to file
            is not always reliable because requires converting variables to string. 
            Only attrs that are of type: string, number, and list of number and strings are 
            saved somewhat correctly. Dictionaries are not saved. 

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
            zfill = figmanip.n_digits(len(self)-1)

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
            s.save(filepath=folderpath/filename, only_data=only_data, check_overwrite=False, verbose=False, **kwargs)
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
            None

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
        self._data = []
        for i, filepath in enumerate(fl):
            # if verbose: print(f'Loading: {folderpath}')
            if verbose: print(f'    {i+1}/{len(fl)}: {filepath.name}')
            s = Spectrum()
            s.load(filepath=filepath, only_data=only_data, **kwargs)
            self.append(s)
        # if verbose: print('Done!')

        ##########################
        # reset check attributes #
        ##########################
        self._length       = None
        self._step         = None
        self._x            = None
        self._monotonicity = None

    def save_all_single_file(self, filepath, only_data=False, limits=None, check_overwrite=False, verbose=False, **kwargs):
        r"""Save all Spectra in one single file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
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
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        #######################
        # check x is the same #
        #######################
        if self.x is None:
            try:
                self.check_same_x()
            except ValueError:
                raise ValueError('Cannot save spectra in one file. x axis are different.\nMaybe try interpolating the x axis (Spectra.interp()) or use Spectra.save() to save spectra in multiple files.')
        
        ###################
        # check overwrite #
        ###################
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')
                
        ##########
        # kwargs #
        ##########
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        ##############################
        # Prepare final data to save #
        ##############################
        if len(self) == 0:
            final = []
        else:
            # gather ys
            x, ys = self._gather_ys(limits=limits)

            # kwargs
            if 'fmt' not in kwargs: # pick best format
                decimal = max([numanip.n_decimal_places(x) for x in self.x])
                temp_decimal = max([numanip.n_decimal_places(y) for y in arraymanip.flatten(ys)])
                if temp_decimal > decimal:
                    decimal = copy.deepcopy(temp_decimal)
                kwargs['fmt'] = f'%.{decimal}f'

            # final data to save
            final = np.zeros((self.length, len(self)+1))
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
            np.savetxt(Path(filepath), final, **kwargs)
        else:
            attrs  = ['_step', '_monotonicity', '_x', '_length']
            attrs += self.get_attrs()
            header = '\n'.join(_attr2str(self, attrs, verbose)) + '\n'

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
            None

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
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '#'
        if 'usecols' not in kwargs:
            kwargs['usecols'] = None
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
        self._data = []
        for i in range(ys.shape[1]):
            s = Spectrum(x=x, y=ys[:, i])
            self.append(s)

        ##########################
        # reset check attributes #
        ##########################
        self._step         = None
        self._monotonicity = None
        self._x            = None
        self._length       = None

        ###############
        # read header #
        ###############
        if only_data is False:
            header = filemanip.load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            if header:
                for line in header:
                    if ':' not in line:
                        pass
                    else:
                        # extract name and value
                        name = line[1:-1].split(':')[0].strip()
                        try:
                            value = eval(str(':'.join(line[1:-1].split(':')[1:])).strip())
                        except:
                            value = str(':'.join(line[1:-1].split(':')[1:])).strip()
                        
                        try:
                            setattr(self, name, value)
                        except Exception as e:
                            if verbose:
                                print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')

    #########
    # check #
    #########
    def check_monotonicity(self):
        """Sets monotonicity attribute to 'increasing' or 'decreasing'.

        Raises:
            ValueError if data is not monotonic.

        Returns:
            None
        """
        monotonicity = [None]*len(self)
        for i in range(len(self)):
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
            for i in range(len(self)):
                text += f'spectrum: {i}, motonicity: {monotonicity[i]}\n'
            raise ValueError(f'some spectra have different monotonicity (increasing, decreasing) or no monotonicity at all (None): \n{text}')

    def fix_monotonicity(self, mode='increasing'):
        """Rearrange x, y such as x array is monotonically increasing or decreasing.

        Args:
            mode (str, optional): increasing or decreasing.

        Returns:
            None
        """
        for s in self.data:
            s.fix_monotonicity(mode=mode)
        self.check_monotonicity()

    def check_length(self):
        """Checks if all spectra has the same length.

        If all spectra have the same length, it sets 
        :py:attr:`Spectra.length` = length.
        Otherwise, it raises an error.

        Returns:
            None

        Raises:
            ValueError: spectra does not have the same length.

        See Also:
            :py:func:`Spectra.check_step`, :py:func:`Spectra.check_same_x`.
        """
        # if zero spectra exists, then length is immediately defined
        # and a warn is raised
        if len(self) == 0:
            self._length = 0
            warnings.warn('no spectra found')
            return
        
        # if only one spectra exists, then length is immediately defined
        if len(self) == 1:
            self._length = len(self.x)
            return
        
        # collect
        length = [None]*len(self)
        for i in range(len(self)):
            try:
                length[i] = len(self.data[i])
            except ValueError:
                pass

        # apply or raise error
        if all(x == length[0] for x in length):
            self._length = length[0]
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
                max_error (number, optional): percentage value (of the average x 
                step) of the maximum allowed error. Default is 0.1 %.                

            Note:

                This method checks if the step between two data points is the 
                same through out the x vector for each spectrum (vector 
                uniformity). Step uniformity within each spectrum is verified by 
                the following equation (see Spectrum.check_step()):

                    (max(steps) - min(steps))/np.mean(steps) * 100 < max_error

                This method also checks if the step is the same between 
                different spectra using the follwing equation

                    (step[i]-step[i+1])/avg(step) * 100 < max_error

            Returns:
                None

            Raises:
                ValueError: If condition 1, or 2 are not satisfied.

            See Also:
                :py:func:`Spectra.check_length`, :py:func:`Spectra.check_same_x`
        """
        # check spectra exists
        if len(self) == 0:
            raise ValueError('no spectra found')

        if self.x is None:
            # try and see if spectra have the same x
            try:
                self.check_same_x(max_error=max_error)

                # check step uniformity
                temp = Spectrum(x=self.x, y=self.x)
                try:
                    temp.check_step(max_error=max_error)
                except ValueError:
                    raise ValueError(f"Spectra have the same x-coordinates, but it is not uniform.")
                self._step = temp.step
                return

            # if spectra have different x-coordinates
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
                    self._step = None
                    raise ValueError(f"Spectra seems to have different step size. Calculated step sizes = {steps}")
                self._step = avg_step
                return
        
        # if all spectra have the same x, check step becames easier
        else:
            # check step uniformity
            temp = Spectrum(x=self.x, y=self.x)
            try:
                temp.check_step(max_error=max_error)
            except ValueError:
                raise ValueError(f"Spectra have the same x-coordinates, but it is not uniform.")
            self._step = temp.step
            return

    def check_same_x(self, max_error=0.1):
        """Check if spectra have same x-coordinates.

        If data has same x-coordinates, it sets :py:attr:`Spectra.x` = x.
        Otherwise, it raises an error.

        Args:
            max_error (number, optional): percentage value (in terms of the
                average x step) of the maximum allowed error. If None, the 
                Default value from settings will be used.

                max(s[i].x - s[i+1].x)/avg(step) * 100 < max_error

        Returns:
            None

        Raises:
            ValueError: If any x-coodinate of any two spectrum is different.

        See Also:
            :py:func:`Spectra.check_length`, :py:func:`Spectra.check_step`.
        """
        # if zero spectra exists, then x is immediately defined
        # and a warn is raised
        if len(self) == 0:
            warnings.warn('no spectra found')
            self._x      = []
            self._length = 0
            return
        
        # if only one spectra exists, then x is immediately defined
        if len(self) == 1:
            self._x      = self[0].x
            self._length = len(self.x)
            return
        
        # check length
        self.check_length()

        # average step
        if self.step is None:
            step = 0
            for s in self:
                step += np.mean(np.diff(s.x))
            step = step/len(self)
        else:
            step = self.step

        # check x between spectra
        x = ['same as the previous']*len(self)
        for idx in range(len(self)-1):
            if max(abs(self[idx].x - self[idx+1].x))*100/abs(step) > max_error:
                x[idx] = 'different'

        # apply or raise error
        if 'different' not in x:
            self._x = self[0].x
        else:
            text = ''
            for i in range(len(self)):
                text += f'spectrum: {i}: {x[i]}\n'
            raise ValueError(f'some spectra have different x: \n{text}\n\nUse brixs.Spectra.interp() to interpolate the data and make the x axis for different spectra match.')

        # update length
        self._length = len(self.x)

    #############
    # modifiers #
    #############
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
        ss._x            = None
        # ss._monotonicity = None

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
        ss._step         = None
        ss._x            = None
        ss._monotonicity = None
        
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

        return ss

    ###############
    # modifiers 2 #
    ###############
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
        ss._step         = None
        ss._x            = None
        ss._monotonicity = None

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
        ss._step         = None
        ss._x            = None
        ss._monotonicity = None

        return ss

    def set_roll(self, value):
        """Roll y data recursively.

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
            ss[i] = s.set_roll(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # ss._length       = None
        # ss._step         = None
        ss._x            = None
        ss._monotonicity = None

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Spectra`
        """
        ss = self.copy()
        for i, s in enumerate(self):
            ss[i] = s.floor(limits=limits)
        return ss

    def normalize(self, value=1, limits=None):
        """Set a factor such as the average y between limits is equal to value.

        Args:
            value (number): value. Default is 1.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
        s.copy_attrs_from(self)

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
            if self.x is None:
                if start is None:
                    start = min([min(s.x) for s in self.data])
                if stop is None:
                    stop = max([max(s.x) for s in self.data])
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
                    return self

        # new spectra object
        ss = self.copy()

        # interpolate
        for i, s in enumerate(self):
            ss[i] = s.interp(x=x, start=start, stop=stop, num=num, step=step)

        ##########################
        # reset check attributes #
        ##########################
        ss._length       = None
        ss._step         = None
        ss._x            = None
        ss._monotonicity = None

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
        if self.x is None:
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
        ss._length       = None
        # ss._step         = None
        ss._x            = None
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
        ss._step         = None
        ss._x            = None
        ss._monotonicity = None

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
        ss._length       = None
        ss._step         = None
        ss._x            = None
        ss._monotonicity = None

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
        ss._length       = None
        ss._step         = None
        ss._x            = None
        ss._monotonicity = None

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            attrs2merge (list, optional): if not None, only spectra named here,
                will be copied to the final spectrum. The attr must be a list of
                numbers with the same length as the number of spectra. The value
                saved on the returned Spectrum is a weighted avereged sum of 
                these attrs.
        
        Return:
            :py:class:`Spectrum`
        """
        assert isinstance(indexes, Iterable), 'indexes must be a list'

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
        for i in indexes:
            ss.append(self[i].copy().set_factor(weights[i]))

        if attrs2merge is not None:
            for attr in attrs2merge:
                temp = self.__getattribute__(attr)
                new = []
                for i in indexes:
                    temp.append(temp[i]*weights[i])
                ss.__setattr__(attr, np.mean(new))

        return ss.calculate_average(limits=limits)
    
    def merge_and_replace(self, indexes, weights=None, limits=None, attrs2merge=None):
        """return spectra with replaced spectra

        Args:
            indexes (list): spectra index numbers
            weights (list or None, optional): if not None, applies a weighted 
                average. Default is None.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
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

    def merge_duplicates(self, ref, limits=None, attrs2merge=None):
        """return spectra where spectrum with same attr are merged

        Args:
            ref (str or list): reference value for interpolating. If `str`, it
                will get values from attribute. If list, list must be the same 
                length as the number of Spectra.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            attrs2merge (list, optional): if not None, only spectra named here,
                will be copied to the final spectrum. The attr must be a list of
                numbers with the same length as the number of spectra. The value
                saved on the returned Spectrum for merged spectrum is a weighted avereged sum of 
                these attrs.

        Return:
            :py:class:`Spectra`
        """
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
            for attr in attrs2merge:
                temp = self.__getattribute__(attr)
                assert isinstance(temp, Iterable), f'{attr} must be an iterable type'
                assert len(temp) == len(self), f'Lenght of attr {attr} must be the same as the number of spectra.\nlenght of attr: {len(attr)}\nnumber of spectra: {len(self)}'
                assert sum([numanip.is_number(x) for x in temp]) == len(temp), f'{attr} must be a list of numbers'

        ##########################################
        # check if spectra indeed has duplicates #
        ##########################################
        if arraymanip.has_duplicates(values) == False:
            return self.copy(limits=limits)
        
        ###############
        # new spectra #
        ###############
        ss = Spectra()
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
                        ss.__setattr__(attr, ss.__getattribute__(attr).append(self.__getattribute__(attr)[i]))
            else:
                indexes = [i for i, x in enumerate(values) if x == _x]
                if i == indexes[0]:
                    _s = self.merge(indexes=indexes, limits=limits, attrs2merge=attrs2merge)
                    ss.append(_s)
                    if attrs2merge is not None:
                        for attr in attrs2merge:
                            ss.__setattr__(attr, ss.__getattribute__(attr).append(_s.__getattribute__(attr)))
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
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
        if self.x is None:
            try:
                self.check_same_x()
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
        s.copy_attrs_from(self)

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
        s.copy_attrs_from(self)

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
        # im.copy_attrs_from(self)
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Image`.
        """
        # check centers
        if x_centers is not None:
            assert len(x_centers) == len(self), f'centers must have the same number of items as the number of spectra.\nnumber of centers: {len(centers)}\nnumber of spectra: {len(self)}'

        # gather ys
        y, ys = self._gather_ys(limits=limits)

        im = Image(data=ys)
        im.copy_attrs_from(self)
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

        Returns:
            :py:class:`Image`.
        """
        # check centers
        if y_centers is not None:
            assert len(y_centers) == len(self), f'centers must have the same number of items as the number of spectra.\nnumber of centers: {len(centers)}\nnumber of spectra: {len(self)}'

        # gather ys
        y, ys = self._gather_ys(limits=limits)
        ys = ys.transpose()

        im = Image(data=ys)
        im.copy_attrs_from(self)
        im.x_centers = y
        im.y_centers = y_centers

        return im
    

    def calculate_shift(self, mode='max', limits=None, **kwargs):
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
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
        if mode in 'cc':
            values = list(np.array(self.calculate_roll(mode='cc', limits=limits))*self.step)
        ###############################
        # sequential cross-corelation #
        ###############################
        if mode in 'seq':
            values = list(np.array(self.calculate_roll(mode='seq', limits=limits))*self.step)
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
            fit, popt, err, fs = self.fit_peak(limits=limits, **kwargs)
            values = np.array([_[1] for _ in popt])
            values = -values + values[0]
        else:
            raise ValueError(f'mode={mode} not valid. Valid modes: `cc`, `max`, `peak`')
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
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

        ######################
        # x must be the same #
        ######################
        if self.x is None:
            self.check_same_x()

        #########################################################
        # x must be uniform (same step between each data point) #
        #########################################################
        self.check_step()

        ####################
        # cross-corelation #
        ####################
        if mode == 'cc':
            ss = self._copy(limits=limits)
            for i, s in enumerate(self):
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
            values = np.array(self.calculate_shift(mode='max', limits=limits, **kwargs))
            values = list(int(round(values/self.step)))
        #########
        # peaks #
        #########
        elif mode == 'peak':
            values = np.array(self.calculate_shift(mode='peak', limits=limits, **kwargs))
            values = list(int(round(self.calculated_shift/self.step)))
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
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
            fit, popt, err, fs = self.fit_peak(limits=limits, **kwargs)
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.

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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
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
        # check number of values matches the number of spectra
        if len(self) != len(values):
            raise ValueError(f'number of values ({len(values)}) do not match the number of spectra ({len(self)})')

        # CALCULATION
        x = self.calculate_shift(mode=mode, limits=limits, **kwargs)

        # return
        final = Spectrum(x=-np.array(x), y=values)
        popt, model, R2 = final.polyfit(deg=deg)
        final.popt      = popt
        final.model     = model
        final.R2        = R2
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

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
        
        Returns:
            list

        .. _numpy.trapz(): https://numpy.org/doc/stable/reference/generated/numpy.trapz.html
        """
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
        
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
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
        
        Returns:
            list
        """
        return [s.calculate_x_sum(limits=limits) for s in self]

    def calculate_y_average_per_spectrum(self, limits=None):
        """returs a list of the average x value within range for each spectrum

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
         
        Returns:
            list
        """
        return [s.calculate_x_average(limits=limits) for s in self]

    def calculate_y_average_per_spectrum(self, limits=None):
        """returs a list of the average y value within range for each spectrum

        Args:
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
         
        Returns:
            list
        """
        return [s.calculate_y_average(limits=limits) for s in self]

    def polyfit(self, deg, limits=None):
        """Fit data recursively with a polynomial. Wrapper for `numpy.polyfit()`_.

        Args:
            deg (int): degree of the fitting polynomial.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
         
        Returns:
            popt, f(x), R2
            list with polynomial coefficients, highest power first.
            list with Model function f(x).
            list with R2.

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        popt  = [0]*len(self)
        model = [0]*len(self)
        R2    = [0]*len(self)
        for i in range(len(self)):
            popt[i], model[i], R2[i] = self[i].polyfit(deg=deg, limits=limits)

        return popt, model, R2

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
    def plot(self, ax=None, offset=0, shift=0, roll=0, factor=1, calib=1, smooth=1, limits=None, switch_xy=False, vi=0, hi=0, pvi=0, phi=0, **kwargs):
        """Plot spectra. Wrapper for `matplotlib.pyplot.plot()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            calib, shift (number or list, optional): multiplicative and additive factor
                 on the x-coordinates. calib is applied first. If list, it must 
                 have the same length as the number of spectra.
            factor, offset (number or list, optional): multiplicative and additive factor
                 on the y-coordinates. Factor is applied first. If list, it must 
                 have the same length as the number of spectra.
            roll (int or list, optional): Roll value of array elements of the x-coordinates.
                If list, it must have the same length as the number of spectra.
            smooth (int, optional): number of points to average data. Default is 1.
            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
                of pairs `((xi_1, xf_1), (xi_2, xf_2), ...)`, or None. If None, 
                this function simply returns None. If pairs, each pair 
                represents the start and stop of a data range from x. Use 
                `x_start = None` or `x_stop = None` to indicate the minimum or 
                maximum x value of the data, respectively. If limits = [], i.e.,
                an empty list, it assumes `limits = (None, None)`.
            switch_xy (bool, optional): Switch x and y axis.
            hi, vi (number, optional): horizontal and vertical increments for 
                cascading plots.
            phi, pvi (number, optional): percentage wise horizontal and vertical 
                increments for cascading plots (percentage of the y-range for 
                each spectrum).
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

        #####################################
        # vertical and horizontal increment #
        #####################################
        vi  = [hi]*len(self)
        hi  = [hi]*len(self)

        #########
        # calib #
        #########
        if isinstance(calib, Iterable):
            assert len(calib) == len(self), f'calib must be a number of a list with length compatible with the number of spectra.\nnumber of calib: {len(calib)}\nnumber of spectra: {len(self)}'
        else:
            calib = [calib]*len(self)
        ##########
        # factor #
        ##########
        if isinstance(factor, Iterable):
            assert len(factor) == len(self), f'factor must be a number of a list with length compatible with the number of spectra.\nnumber of factor: {len(factor)}\nnumber of spectra: {len(self)}'
        else:
            factor = [factor]*len(self)

        #####################################################
        # percentage wise vertical and horizontal increment #
        #####################################################
        pvi = [max(s.y)*factor[i]*pvi/100 for i, s in enumerate(self)]
        phi = [max(s.x)*calib[i]*phi/100  for i, s in enumerate(self)]

        ##########
        # offset #
        ##########
        if isinstance(offset, Iterable):
            assert len(offset) == len(self), f'offset must be a number of a list with length compatible with the number of spectra.\nnumber of offsets: {len(offset)}\nnumber of spectra: {len(self)}'
        else:
            offset = [offset]*len(self)
            for i in range(len(self)):
                offset[i] = offset[i] + (vi[i] * i) + (pvi[i] * i)
        #########
        # shift #
        #########
        if isinstance(shift, Iterable):
            assert len(shift) == len(self), f'shift must be a number of a list with length compatible with the number of spectra.\nnumber of shift: {len(shift)}\nnumber of spectra: {len(self)}'
        else:
            shift = [shift]*len(self)
            for i in range(len(self)):
                shift[i] = shift[i] + (hi[i] * i) + (phi[i] * i)
        ########
        # roll #
        ########
        if isinstance(roll, Iterable):
            assert len(roll) == len(self), f'roll must be a number of a list with length compatible with the number of spectra.\nnumber of shift: {len(shift)}\nnumber of spectra: {len(self)}'
        else:
            roll = [roll]*len(self)
        
        ########
        # plot #
        ########
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self.data[i].plot(ax=ax, offset=offset[i], shift=shift[i], factor=factor[i], calib=calib[i], smooth=smooth, switch_xy=switch_xy, limits=limits, **kwargs)

        return temp

# %% ================================ Image =============================== %% #
class Image(metaclass=_Meta):
    """Returns a ``spectra`` object.

    Args:
        data (list or array, optional): list of :py:class:`spectrum` objects.
        filepath (str or Path, optional): filepath.
        x_centers, y_centers: (list, optional): pixel center labels.
        **kwargs: kwargs are passed to :py:func:`Image.load` function.

    Usage:        
            >>> im = br.Image()
            >>> im = br.Image(data)
            >>> im = br.Image(data=data)
            >>> im = br.Image(data=data, x_centers=x_centers, y_centers=y_centers)
            >>> im = br.Image(filepath=<filepath>)
            >>> im = br.Image(filepath=<filepath>, delimiter=',')
            >>>
            >>> print(im.get_core_attrs()) # print list of core attrs
            >>> print(im.get_attrs())      # print list of attrs
            >>> print(im.get_methods())    # print list of methods available
    """
    # read only and non-removable arguments
    _read_only     = ['x_step', 'x_monotonicity', 'y_step', 'y_monotonicity']
    _non_removable = []
    
    def __init__(self, data=None, filepath=None, x_centers=None, y_centers=None, **kwargs):
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

        # modifiers
        self._factor = 1
        self._offset = 0 

        # labels
        self._x_centers = None
        self._y_centers = None
        self._x_edges   = None
        self._y_edges   = None

        ################
        # loading data #
        ################
        if data is not None:
            self.data = data
            if x_centers is not None:
                self.x_centers = x_centers
            if y_centers is not None:
                self.y_centers = y_centers
        elif filepath is not None:
            self.load(filepath, **kwargs)
        return

    ###################
    # core attributes #
    ###################
    @property
    def data(self):
        return copy.deepcopy(self._data)
    @data.setter
    def data(self, value):
        self._data  = np.array(value, dtype='float')
        self.x_centers = None
        self.y_centers = None    
        self._x_edges  = None
        self._y_edges  = None    
    @data.deleter
    def data(self):
        raise AttributeError('Cannot delete object.')

    @property
    def x_centers(self):
        return self._x_centers
    @x_centers.setter
    def x_centers(self, value):
        if value is None:
            value = np.arange(0, self.data.shape[1])
        elif isinstance(value, Iterable):
            assert len(value) == self.shape[1], f"number of x centers ({len(value)}) must be the same as the number of pixel columns ({self.data.shape[1]})"
            # assert arraymanip.check_monotonicity(value) == 1, f"x centers must be a monotonically increasing array"
        else:
            raise ValueError(f"x centers must be None or an iterable (list, tuple, or 1D array)")
        self._x_centers = np.array(value, dtype='float')
        temp            = list(arraymanip.moving_average(value, 2))
        self.x_edges    = [value[0] - temp[0]] + temp + [2*value[-1] - temp[-1]]
    @x_centers.deleter
    def x_centers(self):
        self._x_centers = np.arange(0, self.data.shape[1])

    @property
    def y_centers(self):
        return self._y_centers
    @y_centers.setter
    def y_centers(self, value):
        if value is None:
            value = np.arange(0, self.data.shape[0])
        elif isinstance(value, Iterable):
            assert len(value) == self.shape[0], f"number of y centers ({len(value)}) must be the same as the number of pixel columns ({self.data.shape[0]})"
            # assert arraymanip.check_monotonicity(value) == 1, f"y centers must be a monotonically increasing array"
        else:
            raise ValueError(f"y centers must be None or an iterable (list, tuple, or 1D array)")
        self._y_centers = np.array(value, dtype='float')
        temp            = list(arraymanip.moving_average(value, 2))
        self.y_edges    = [value[0] - temp[0]] + temp + [2*value[-1] - temp[-1]]
    @y_centers.deleter
    def y_centers(self):
        self._y_centers  = np.arange(0, self.data.shape[0])

    @property
    def x_edges(self):
        return self._x_edges
    @x_edges.setter
    def x_edges(self, value):
        if value is None:
            centers = np.arange(0, self.data.shape[1])
            temp    = list(arraymanip.moving_average(centers, 2))
            value   = [centers[0] - temp[0]] + temp + [2*centers[-1] - temp[-1]]
        elif isinstance(value, Iterable):
            assert len(value) == self.shape[1] + 1, f"number of x edges ({len(value)}) must be the same as the number of pixel columns plus one ({self.data.shape[1] + 1})"
            # assert arraymanip.check_monotonicity(value) == 1, f"x edges must be a monotonically increasing array"
        else:
            raise ValueError(f"x edges must be None or an iterable (list, tuple, or 1D array)")
        self._x_edges   = np.array(value, dtype='float')
        self._x_centers = np.array(arraymanip.moving_average(value, 2), dtype='float')
    @x_edges.deleter
    def x_edges(self):
        raise NotImplementedError('this is not implemented yet')

    @property
    def y_edges(self):
        return self._y_edges
    @y_edges.setter
    def y_edges(self, value):
        if value is None:
            centers = np.arange(0, self.data.shape[0])
            temp    = list(arraymanip.moving_average(centers, 2))
            value   = [centers[0] - temp[0]] + temp + [2*centers[-1] - temp[-1]]
        elif isinstance(value, Iterable):
            assert len(value) == self.shape[0] + 1, f"number of y edges ({len(value)}) must be the same as the number of pixel rows plus one ({self.data.shape[0] + 1})"
            assert arraymanip.check_monotonicity(value) == 1, f"y edges must be a monotonically increasing array"
            self._y_centers = np.array(arraymanip.moving_average(value, 2), dtype='float')
        else:
            raise ValueError(f"y edges must be None or an iterable (list, tuple, or 1D array)")
        self._y_edges = np.array(value, dtype='float')
    @y_edges.deleter
    def y_edges(self):
        raise NotImplementedError('this is not implemented yet')

    ###################################
    # computed (read-only) attributes #
    ###################################
    @property
    def shape(self):
        if self._data is None:
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
        if self._data is None:
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
        if self._data is None:
            return None
        ss = Spectra(n=self.shape[1])
        for i in range(self.shape[1]):
            ss[i] = Spectrum(x=self.y_centers, y=self.data[:, i])
        ss.copy_attrs_from(self)
        return ss
    @columns.setter
    def columns(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @columns.deleter
    def columns(self):
        raise AttributeError('Cannot delete object.')

    @property
    def rows(self):
        if self._data is None:
            return None
        ss = Spectra(n=self.shape[0])
        for i in range(self.shape[0]):
            ss[i] = Spectrum(x=self.x_centers, y=self.data[i, :])
        ss.copy_attrs_from(self)
        return ss
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
    @factor.deleter
    def factor(self):
        raise AttributeError('Cannot delete object.')

    #################
    # magic methods #
    #################
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

        final.copy_attrs_from(self)
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

        final.copy_attrs_from(self)
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
        
        final.copy_attrs_from(self)
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
        
        final.copy_attrs_from(self)
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
    def get_core_attrs(self): 
        """return a list of core attrs"""
        return settings._reserved_words['Image']['pseudovars']
    
    def get_attrs(self):
        """return attrs that are user defined.""" 
        return [key for key in self.__dict__.keys() if key.startswith('_') == False and key not in settings._reserved_words]
    
    def get_methods(self):
        """return a list of methods available"""
        return [key for key in self.__dir__() if key.startswith('_') == False and key not in self.get_attrs() + self.get_core_attrs()]
    
    def remove_attrs(self):
        """Delete all user defined attrs."""
        for attr in self.get_attrs():
            self.__delattr__(attr)

    def copy_attrs_from(self, s):
        """Copy user defined attributes from another brixs object.

        Args:
            s (brixs object): Either a Spectrum, Spectra, Image, or PhotonEvents
                to copy user defined attributes from.
        
        Returns:
            None
        """
        # check type
        if isinstance(s, Spectrum) or isinstance(s, Spectra) or isinstance(s, Image) or isinstance(s, PhotonEvents):
            pass
        else:
            raise TypeError(f'type {type(s)} not valid\nCan only copy user attrs from type br.Spectrum, br.Spectra, br.Image, or br.PhotonEvents')

        # transfer attrs
        for attr in s.get_attrs():
            value = copy.deepcopy(s.__dict__[attr])
            self.__setattr__(attr, value)

    ###########
    # attrs 2 #
    ###########
    pass

    ###########
    # support #
    ###########
    def _calculated_vmin_vmax(self):
        """returns optimal vmin and vmax for visualization
        
        vmin is set the the max of the intensity distribution (max of histogram)
        
        vmax is set when intensity distribution (histogram) drops below 0.01% of the max.

        Returns:
            vmin, vmax
        """
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
        return vmin, vmax
    
    ################
    # core methods #
    ################
    pass

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
            im.x_centers = copy.deepcopy(self.x_centers)
            im.y_centers = copy.deepcopy(self.y_centers)
            return im
        
        ########################################
        # check if extract is really necessary #
        ########################################
        if x_start <= min(self.x_centers) and x_stop >= max(self.x_centers):
            if y_start <= min(self.y_centers) and y_stop >= max(self.y_centers):
                im = Image(data=data)
                im.x_centers = copy.deepcopy(self.x_centers)
                im.y_centers = copy.deepcopy(self.y_centers)
                return im
            
        #################
        # check if None #
        #################
        if x_start is None: x_start = min(self.x_centers)
        if x_stop  is None: x_stop  = max(self.x_centers)
        if y_start is None: y_start = min(self.y_centers)
        if y_stop  is None: y_stop  = max(self.y_centers)

        ##################
        # validate input #
        ##################
        assert x_stop > x_start, f'x_start must be smaller than x_stop.'
        assert y_stop > y_start, f'y_start must be smaller than y_stop.'        
        
        ########################
        # check if empty image #
        ########################
        if self.data is None:
            print(f'Warning: copying empty Image')
            return Image()
        if len(self.data) == 0:
            print(f'Warning: copying empty Image')
            return Image()
        
        ########
        # crop #
        ########
        y_start = int(arraymanip.index(self.y_centers, y_start))
        y_stop  = int(arraymanip.index(self.y_centers, y_stop))
        x_start = int(arraymanip.index(self.x_centers, x_start))
        x_stop  = int(arraymanip.index(self.x_centers, x_stop))
        im      = Image(data=self.data[y_start:y_stop+1, x_start:x_stop+1])
        im.x_centers = self.x_centers[x_start:x_stop+1]
        im.y_centers = self.y_centers[y_start:y_stop+1]
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

        ##################
        # transfer attrs #
        ##################
        im.copy_attrs_from(self)
        return im
    
    #################
    # save and load #
    #################
    def save(self, filepath=None, only_data=False,  check_overwrite=False, verbose=False, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Warning:
            Attrs are saved as comments if only_data is False. Saving attrs to file
            is not always reliable because requires converting variables to string. 
            Only attrs that are of type: string, number, and list of number and strings are 
            saved somewhat correctly. Dictionaries are not saved. 

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format. If None, filepath can be inferred from a
                user defined attr 'im.filepath'. Default is None.  Last used 
                filepath is saved to im.filepath.
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
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        ###################
        # check overwrite #
        ###################
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        ##########
        # kwargs #
        ##########
        if 'fmt' not in kwargs: # pick best format
            decimal = max([numanip.n_decimal_places(x) for x in arraymanip.flatten(self._data)])
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        #####################
        # header and footer #
        #####################
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
        else:
            attrs  = ['_x_step', '_y_step', '_x_monotonicity', '_y_monotonicity']
            attrs += ['_factor', '_offset']
            attrs += ['_x_centers', '_y_centers', '_x_edges', '_y_edges']
            attrs += self.get_attrs()
            header = '\n'.join(_attr2str(self, attrs, verbose)) + '\n'

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

    def load(self, filepath, only_data=False, verbose=False, **kwargs):
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
            None

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
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        ########
        # read #
        ########
        data = np.genfromtxt(Path(filepath), **kwargs)

        ##########
        # assign #
        ##########
        self.data = data

        ###############
        # read header #
        ###############
        if only_data is False:
            header = filemanip.load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            if header:
                for line in header:
                    if ':' not in line:
                        pass
                    else:
                        # extract name and value
                        name = line[1:-1].split(':')[0].strip()
                        try:
                            value = eval(str(':'.join(line[1:-1].split(':')[1:])).strip())
                        except:
                            value = str(':'.join(line[1:-1].split(':')[1:])).strip()
                        try:
                            setattr(self, name, value)
                        except Exception as e:
                            if verbose:
                                print(f'cannot read attr ({name}: {value})\nAttribute not set\n{e}\n')

    #########
    # check #
    #########
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
        s = Spectrum(x=self.x_centers, y=self.x_centers)
        try:
            s.check_step()
        except ValueError:
            raise ValueError(f"Step in the x centers seems not to be uniform. Set im.x_centers = None or change im.x_centers")
        self._x_step = s.step
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
        s = Spectrum(x=self.y_centers, y=self.y_centers)
        try:
            s.check_step()
        except ValueError:
            raise ValueError(f"Step in the y centers seems not to be uniform. Set im.y_centers = None or change im.y_centers")
        self._y_step = s.step
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
        s = Spectrum(x=self.x_centers, y=self.x_centers)
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
        s = Spectrum(x=self.y_centers, y=self.y_centers)
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
        # check mode
        increasing = ['inc', 'i', 'up', 'increasing', 'increasingly']
        decreasing = ['dec', 'd', 'down', 'decreasing', 'decreasingly']
        if mode not in increasing and mode not in decreasing:
            raise ValueError('mode should be "decreasing" or "increasing".')
        
        if mode in increasing: mode = 'increasing'
        else: mode = 'decreasing'

        # turn array into monotonic
        if self.monotonicity is None:
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
        im._data   += float(value).astype(self.data.dtype)
        im._offset += float(value).astype(self.data.dtype)
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
        im._data   *= float(value).astype(self.data.dtype)
        im._factor *= float(value).astype(self.data.dtype)
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
        """Roll pixels rows left and right in terms of x centers.

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
        value = [int(round(k)) for k in value]
        
        ########
        # roll #
        ########
        im = self.copy()
        for i, v in enumerate(value):
            if v != 0:
                s = Spectrum(im._data[:, i]).set_roll(v)
                im._data[:, i] = s.y
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
        value = [int(round(k)) for k in value]
        
        ########
        # roll #
        ########
        im = self.copy()
        for i, v in enumerate(value):
            if v != 0:
                s = Spectrum(im._data[i, :]).set_roll(v)
                im._data[i, :] = s.y
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


    def set_horizontal_roll_via_polyval(self, p):
        """Set horizontal roll to np.polyval(p, y pixel).

        Args:
            p (array): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term.

        Returns:
            :py:class:`Image`
        """
        f = lambda y: np.polyval(p, y)
        return self.set_horizontal_roll_via_function(f)
    
    def set_vertical_roll_via_polyval(self, p):
        """Set vertical roll to np.polyval(p, x pixel).

        Args:
            p (array): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term.

        Returns:
            :py:class:`Image`
        """
        f = lambda x: np.polyval(p, x)
        return self.set_vertical_roll_via_function(f)

    def set_vertical_roll_via_function(self, f):
        """Set vertical roll to f(x pixel).

        Args:
            f (function): function where argument is x centers elements

        Returns:
            :py:class:`Image`
        """
        value = np.array([f(x) for x in range(len(self.x_centers))])
        return self.set_vertical_roll(value=value)

    def set_horizontal_roll_via_function(self, f):
        """Set horizontal roll to f(y pixel).

        Args:
            f (function): function where argument is y centers elements

        Returns:
            :py:class:`Image`
        """
        value = np.array([f(x) for x in range(len(self.y_centers))])
        return self.set_horizontal_roll(value=value)


    def floor(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Set intensity to zero inside limits to zero.

        Args:
            x, y (int, optional): x and y position to sample background intensity.
            n, nx, ny (int, optional): size of the pixel window around x, y.

        Returns:
            None
        """
        temp = self.copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
        value = temp.calculate_average()
        im = self.copy()
        return im.set_factor(value)

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
        im.copy_attrs_from(self)
        im.x_centers = copy.deepcopy(self.y_centers)
        im.y_centers = copy.deepcopy(self.x_centers)       

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
        cols = self.columns

        cols.ref = self.x_centers
        ss = cols.interp_spectra(ref='ref', start=start, stop=stop, num=num, step=step, x=x)
        
        im = ss.stack_spectra_as_columns()
        im.copy_attrs_from(self)
        im.x_centers = ss.ref
        im.y_centers = self.y_centers
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
        rows = self.rows
        rows.ref = self.y_centers
        ss = rows.interp_spectra(ref='ref', start=start, stop=stop, num=num, step=step, x=y)
        
        im = ss.stack_spectra_as_rows()
        im.copy_attrs_from(self)
        im.x_centers = self.x_centers
        im.y_centers = ss.ref
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
            ncols = self.shape[0]
        if nrows is None:
            nrows = self.shape[1]
        if numanip.is_integer(ncols) == False or numanip.is_integer(nrows)  == False or ncols < 0 or nrows < 0:
            raise ValueError("Number of bins must be a positive integer.")

        # is divisible
        assert self.shape[1] % nrows == 0, f"The {self.shape[1]} pixels in a row is not evenly divisible by {nrows}\nPlease, pick one of the following numbers: {np.sort(list(numanip.factors(self.shape[1])))}"
        assert self.shape[0] % ncols == 0, f"The {self.shape[0]} pixels in a column is not evenly divisible by {ncols}\nPlease, pick one of the following numbers: {np.sort(list(numanip.factors(self.shape[0])))}"

        ###############
        # Calculation #
        ###############
        _bins_size = np.array((self.shape[0]/ncols, self.shape[1]/nrows))
        reduced    = Image(np.add.reduceat(np.add.reduceat(self._data, list(map(float, np.arange(0, self.shape[0], _bins_size[0]))), axis=0), list(map(float, np.arange(0, self.shape[1], _bins_size[1]))), axis=1))
        
        # x and y centers
        _x_edges = np.arange(0, self.shape[1]+_bins_size[1]/2, _bins_size[1])
        _y_edges = np.arange(0, self.shape[0]+_bins_size[0]/2, _bins_size[0])
        
        reduced._x_centers = arraymanip.moving_average(_x_edges, n=2)
        reduced._y_centers = arraymanip.moving_average(_y_edges, n=2)
        reduced._x_edges = _x_edges
        reduced._y_edges = _y_edges

        ##################
        # transfer attrs #
        ##################
        reduced.copy_attrs_from(self)

        return reduced
    
    ########################
    # calculation and info #
    ########################
    def calculate_average(self):
        """returns the average intensity (matrix element) value

        Returns:
            number
        """
        return np.mean([s.calculate_average() for s in self.columns])

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
            hist, bin_edges = np.histogram(arraymanip.flatten(self._data), bins=nbins)
            while max(hist) < self.shape[0]*self.shape[1]*0.05:
                nbins -= 10
                hist, bin_edges = np.histogram(arraymanip.flatten(self._data), bins=nbins)
                if nbins < 50:
                    break
        elif numanip.is_integer(nbins):
            hist, bin_edges = np.histogram(arraymanip.flatten(self._data), bins=int(nbins))
        else:
            raise TypeError('nbins must be a integer')

        x = arraymanip.moving_average(bin_edges, 2)
        s = Spectrum(x=x, y=hist)

        ##################
        # transfer attrs #
        ##################
        s.copy_attrs_from(self)
 
        return s


    def integrated_rows_vs_y_centers(self):
        """return spectra with integrated row values vs y centers

        Returns:
            :py:class:`Spectrum`
        """
        s = Spectrum(x=self.y_centers, y=np.sum(self._data, axis=1))
        s.copy_attrs_from(self)
        return s
    
    def integrated_columns_vs_x_centers(self):
        """return spectra with integrated column values vs x centers

        Returns:
            :py:class:`Spectrum`
        """
        s = Spectrum(x=self.y_centers, y=np.sum(self._data, axis=1))
        s.copy_attrs_from(self)
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


    def calculate_horizontal_shift(self, mode='cc', limits=None, limit_size=1000, **kwargs):
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

            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
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
            list
        """
        ss = self.rows
        if limit_size:
            if len(ss) > limit_size:
                raise ValueError(f'Number of rows is bigger than limit_size.\nImage is seems to be too big.\nAre you sure you want to calculate shifts for such a big image.\nIf so, either set limit_size to False or a higher value.\nNumber of columns: {len(self.y_centers)}\nlimit size: {limit_size}')

        # calculate
        values = ss.calculate_shift(mode=mode, limits=limits, **kwargs)
        return values
    
    def calculate_vertical_shift(self, mode='cc', limits=None, limit_size=1000, **kwargs):
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

            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
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
            list
        """
        ss = self.columns
        if limit_size:
            if len(ss) > limit_size:
                raise ValueError(f'Number of columns is bigger than limit_size.\nImage is seems to be too big.\nAre you sure you want to calculate shifts for such a big image.\nIf so, either set limit_size to False or a higher value.\nNumber of columns: {len(self.x_centers)}\nlimit size: {limit_size}')

        # calculate
        values = ss.calculate_shift(mode=mode, limits=limits, **kwargs)
        return values
    
    def calculate_horizontal_roll(self, mode='cc', limits=None, limit_size=1000, **kwargs):
        """Calculate intensity misalignments in terms of pixel roll.

        Args:
            mode (string, optional): method used. Options are: 
                 
                 1) 'cc': align rows via cross-correlation (cc), where cc for
                 all spectra is calculated against the frist spectrum. 
                 
                 2) 'seq': align via 'cros-correlation' (cc), where cc is 
                 calculated against previous row.

                 3) 'max': Align the max point of every row. 
                 
                 4) 'peak': Fit one peak in each row and align them 
                 (requires that `brixs.addons.fitting` is imported)

            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
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
            list
        """
        ss = self.rows
        if limit_size:
            if len(ss) > limit_size:
                raise ValueError(f'Number of rows is bigger than limit_size.\nImage is seems to be too big.\nAre you sure you want to calculate shifts for such a big image.\nIf so, either set limit_size to False or a higher value.\nNumber of columns: {len(self.y_centers)}\nlimit size: {limit_size}')

        # change from x centers to x pixels
        for s in ss:
            s._x = np.arange(len(self.x_centers))

        # calculate
        values = ss.calculate_roll(mode=mode, limits=limits, **kwargs)
        return values
    
    def calculate_vertical_roll(self, mode='cc', limits=None, limit_size=1000, **kwargs):
        """Calculate intensity misalignments in terms of pixel roll.

        Args:
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
            list
        """
        ss = self.columns
        if limit_size:
            if len(ss) > limit_size:
                raise ValueError(f'Number of columns is bigger than limit_size.\nImage is seems to be too big.\nAre you sure you want to calculate shifts for such a big image.\nIf so, either set limit_size to False or a higher value.\nNumber of columns: {len(self.x_centers)}\nlimit size: {limit_size}')

        # change from y centers to y pixels
        for s in ss:
            s._x = np.arange(len(self.y_centers))

        # calculate
        values = ss.calculate_roll(mode=mode, limits=limits, **kwargs)
        return values
    
    ############
    # composed #
    ############
    def calculate_vertical_shift_curvature(self, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
        """Calculate vertical shift values to fix curvature.

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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs x centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(x_centers)
        """
        im = self.copy()
        im.floor()
        values = im.calculate_vertical_shift(mode=mode, limits=limits, limit_size=limit_size, **kwargs)

        # calculate poly
        s = Spectrum(x=self.x_centers, y=values)
        fit, popt, R2, model = s.polyfit(deg=deg)
        return s, fit, popt, R2, model
    
    def calculate_horizontal_shift_curvature(self, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
        """Calculate horizontal shift values to fix curvature.

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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs y centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(y_centers)
        """
        im = self.copy()
        im.floor()
        values = im.calculate_horizontal_roll(mode=mode, limits=limits, limit_size=limit_size, **kwargs)

        # calculate poly
        s = Spectrum(x=self.y_centers, y=values)
        fit, popt, R2, model = s.polyfit(deg=deg)
        return s, fit, popt, R2, model
    
    def calculate_vertical_roll_curvature(self, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
        """Calculate vertical roll values to fix curvature.

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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs x centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(x_centers)
        """
        im = self.copy()
        im.floor()
        values = im.calculate_vertical_roll(mode=mode, limits=limits, limit_size=limit_size, **kwargs)

        # calculate poly
        s = Spectrum(x=np.arange(len(values)), y=values)
        fit, popt, R2, model = s.polyfit(deg=deg)
        return s, fit, popt, R2, model
    
    def calculate_horizontal_roll_curvature(self, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
        """Calculate horizontal roll values to fix curvature.

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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs y centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(y_centers)
        """
        im = self.copy()
        im.floor()
        values = im.calculate_horizontal_roll(mode=mode, limits=limits, limit_size=limit_size, **kwargs)

        # calculate poly
        s = Spectrum(x=np.arange(len(values)), y=values)
        fit, popt, R2, model = s.polyfit(deg=deg)
        return s, fit, popt, R2, model
    
    
    def fix_vertical_shift_curvature(self, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs x centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(x_centers)
        """
        s, fit, popt, R2, model = self.calculate_vertical_shift_curvature(deg=deg, mode=mode, limits=limits, limit_size=limit_size, **kwargs)

        return self.set_vertical_shift_via_polyval(popt) 
    
    def fix_horizontal_shift_curvature(self, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs y centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(y_centers)
        """
        s, fit, popt, R2, model = self.calculate_horizontal_shift_curvature(deg=deg, mode=mode, limits=limits, limit_size=limit_size, **kwargs)

        return self.set_horizontal_shift_via_polyval(popt) 
    
    def fix_vertical_roll_curvature(self, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs x centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(x_centers)
        """
        s, fit, popt, R2, model = self.calculate_vertical_roll_curvature(deg=deg, mode=mode, limits=limits, limit_size=limit_size, **kwargs)

        return self.set_vertical_roll_via_polyval(popt) 
    
    def fix_horizontal_roll_curvature(self, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs y centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(y_centers)
        """
        s, fit, popt, R2, model = self.calculate_horizontal_roll_curvature(deg=deg, mode=mode, limits=limits, limit_size=limit_size, **kwargs)

        return self.set_horizontal_roll_via_polyval(popt) 
    
    ##########################        
    # plot and visualization #
    ########################## 
    def pcolormesh(self, ax=None, x_start=None, x_stop=None, y_start=None, y_stop=None, colorbar=False, **kwargs):
        """Display data as a mesh. Wrapper for `matplotlib.pyplot.pcolormesh()`_.

        If x_centers and y_centers have irregular pixel separation, pcolormesh
            does its best to defined pixel edges so centers labels correspond 
            to the real centers (nearest possible).

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

        ######################
        # check monotonicity #
        ######################
        if self.x_monotonicity is None:
            self.check_x_monotonicity()
        if self.y_monotonicity is None:
            self.check_y_monotonicity()

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
        X, Y = np.meshgrid(im.x_centers, im.y_centers)
        pos  = ax.pcolormesh(X, Y, im.data, **kwargs)

        # colorbar
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
    
    def imshow(self, ax=None, x_start=None, x_stop=None, y_start=None, y_stop=None, colorbar=False, origin='lower', verbose=True, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.imshow()`_.

        Warning:
            Pixels are always square. For irregular pixel row/columns, see Image.pcolormesh().
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
    
    def plot(self, ax=None, x_start=None, x_stop=None, y_start=None, y_stop=None, colorbar=False, origin='lower', verbose=True, **kwargs):
        """Display data as an image with axis based on x and y centers. Wrapper for `matplotlib.pyplot.imshow()`_.

        Warning:
            Pixels are always square. For irregular pixel row/columns, see Image.pcolormesh()

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

        ######################
        # check monotonicity #
        ######################
        if self.x_monotonicity is None:
            self.check_x_monotonicity()
        if self.y_monotonicity is None:
            self.check_y_monotonicity()

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
            x  = np.linspace(im.x_centers[0], im.x_centers[-1], len(im.x_centers))
            dx = np.mean(np.diff(x))
            extent_x = [x[0]-dx/2, x[-1]+dx/2]

            y  = np.linspace(im.y_centers[0], im.y_centers[-1], len(im.y_centers))
            dy = np.mean(np.diff(y))
            extent_y = [y[0]-dy/2, y[-1]+dy/2]

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

# %% ============================= PhotonEvents =========================== %% #
class PhotonEvents(metaclass=_Meta):
    """Returns a ``Photon events`` object.

    Args:
        x, y (list or array, optional): array with x, y photon events coordinates
        xlim, ylim (list, optional): two element tuple with min and max possible 
            x and y coordinates.
        filepath (str or Path, optional): filepath.
        **kwargs: kwargs are passed to :py:func:`PhotonEvents.load` function.

    Usage:        
            >>> pe = br.PhotonEvents()
            >>> pe = br.PhotonEvents(x, y)
            >>> pe = br.PhotonEvents(x=x, y=y)
            >>> pe = br.PhotonEvents(x, y, xlim=(0, 10), ylim=(0, 10))
            >>> pe = br.PhotonEvents(filepath=<filepath>)
            >>> pe = br.PhotonEvents(filepath=<filepath>, delimiter=',')
            >>>
            >>> print(pe.get_core_attrs()) # print list of core attrs
            >>> print(pe.get_attrs())      # print list of attrs
            >>> print(pe.get_methods())    # print list of methods available
    """
    _read_only = []
    _non_removable = []
    
    def __init__(self, x=None, y=None, xlim=None, ylim=None, filepath=None, **kwargs): 
        """Initialize the object instance"""
        ###########################
        # Initializing attributes #
        ###########################
        # core
        self._x    = None
        self._y    = None

        # modifiers
        pass

        # labels
        self._xlim = None
        self._ylim = None

        ################
        # loading data #
        ################
        # keyword arguments
        if y is not None:
            if x is None:
                raise ValueError('missing x coordinates')
            self.x = x
            self.y = y

            self.xlim = xlim
            self.ylim = ylim       
        elif filepath is not None:
            self.load(filepath=filepath, **kwargs)
        return

    ###################
    # core attributes #
    ###################
    @property
    def x(self):
        return self._x
    @x.setter
    def x(self, value):
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
       
        #################
        # set attribute #
        #################
        self._x = np.array(value, dtype='float')
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

        #################
        # set attribute #
        #################
        self._y = np.array(value, dtype='float')
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
            self._xlim = np.array((min(self.x), max(self.x)), dtype='float')
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
            assert max(self.x) <= value[-1] or min(self.x) >= value[0], f'x coordinates (from {min(self.x)} to {max(self.x)}) are outside xlim ({value})'

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
            self._tylim = np.array((min(self.y), max(self.y)), dtype='float')
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
            assert max(self.y) <= value[-1] or min(self.y) >= value[0], f'y coordinates (from {min(self.y)} to {max(self.y)}) are outside ylim ({value})'


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
            final = PhotonEvents(x=list(self.x) + list(object.x), y=list(self.y) + list(object.y))
            final.copy_attrs_from(self)
            final.xlim = (min([self.xlim[0], object.xlim[0]]), max([self.xlim[1], object.xlim[1]]))
            final.ylim = (min([self.ylim[0], object.ylim[0]]), max([self.ylim[1], object.ylim[1]]))
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
    def get_core_attrs(self): 
        """return a list of core attrs"""
        return settings._reserved_words['PhotonEvents']['pseudovars']
       
    def get_attrs(self):
        """return attrs that are user defined.""" 
        return [key for key in self.__dict__.keys() if key.startswith('_') == False and key not in settings._reserved_words]
    
    def get_methods(self):
        """return a list of methods available"""
        return [key for key in self.__dir__() if key.startswith('_') == False and key not in self.get_attrs() + self.get_core_attrs()]
    
    def remove_attrs(self):
        """Delete all user defined attrs."""
        for attr in self.get_attrs():
            self.__delattr__(attr)

    def copy_attrs_from(self, s):
        """Copy user defined attributes from another brixs object.

        Args:
            s (brixs object): Either a Spectrum, Spectra, Image, or PhotonEvents
                to copy user defined attributes from.
        
        Returns:
            None
        """
        # check type
        if isinstance(s, Spectrum) or isinstance(s, Spectra) or isinstance(s, Image) or isinstance(s, PhotonEvents):
            pass
        else:
            raise TypeError(f'type {type(s)} not valid\nCan only copy user attrs from type br.Spectrum, br.Spectra, br.Image, or br.PhotonEvents')

        # transfer attrs
        for attr in s.get_attrs():
            value = copy.deepcopy(s.__dict__[attr])
            self.__setattr__(attr, value)

    ###########
    # attrs 2 #
    ###########
    pass

    ###########
    # support #
    ###########
    pass

    ################
    # core methods #
    ################
    pass

    ########
    # copy #
    ########
    def _copy(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Same as copy(), but attributes are not copied to the new object."""
        #############
        # copy data #
        #############
        x = copy.deepcopy(self.x)
        y = copy.deepcopy(self.y)

        #####################
        # if limits is None #
        #####################
        if x_start is None and x_stop is None and y_start is None and y_stop is None:
            return PhotonEvents(x=x, y=y, xlim=self.xlim, ylim=self.ylim)

        ########################################
        # check if extract is really necessary #
        ########################################
        if x_start <= min(self.x) and x_stop >= max(self.x):
            if y_start <= min(self.y) and y_stop >= max(self.y):
                return PhotonEvents(x=x, y=y, xlim=self.xlim, ylim=self.ylim)
            
        #################
        # check if None #
        #################
        if x_start is None: x_start = min(self.x)
        if x_stop  is None: x_stop  = max(self.x)
        if y_start is None: y_start = min(self.y)
        if y_stop  is None: y_stop  = max(self.y)

        ##################
        # validate input #
        ##################
        assert x_stop > x_start, f'x_start must be smaller than x_stop.'
        assert y_stop > y_start, f'y_start must be smaller than y_stop.' 

        ########################
        # check if empty image #
        ########################
        if self.x is None:
            print(f'Warning: copying empty PhotonEvents')
            return PhotonEvents()
        if len(self.x) == 0:
            print(f'Warning: copying empty PhotonEvents')
            return PhotonEvents()       
        
        ########
        # crop #
        ########
        temp = np.array([(x, y) for x, y in zip(self.x, self.y) if ((x > x_start and x < x_stop) and (y > y_start and y < y_stop))])
        return PhotonEvents(x=list(temp[:, 0]), y=list(temp[:, 1]), xlim=(x_start, x_stop), ylim=(y_start, y_stop))
    
    def copy(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Return a copy of the object.

        Args:
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.

        Returns:
            :py:attr:`PhotonEvents`
        """
        pe = self._copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
        pe.copy_attrs_from(self)
        return pe

    def clip(self, mask):
        """Return a masked copy of the object.

        Args:
            mask (list): list with rectangular coordinates `(x_start, x_stop, y_start, y_stop)`
                or a list with multiple rectangular coordinates, i.e., `[(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])`

        Returns:
            :py:attr:`PhotonEvents`
        """
        ######################
        # assert mask format #
        ######################
        assert isinstance(mask, Iterable), 'mask must be iterable'
        if len(mask) == 4:
            if isinstance(mask[0], Iterable) == False:
                mask = [mask, ] 
        for m in mask:
            assert len(m) == 4, 'mask must have the format: [(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])'
        
        ########
        # clip #
        ########
        pe = PhotonEvents(x=[], y=[])
        for x_start, x_stop, y_start, y_stop in mask: 
            pe += self._copy(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
        pe.copy_attrs_from(self)

        return pe

    #################
    # save and load #
    #################
    def save(self, filepath=None, only_data=False, check_overwrite=False, verbose=False, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Warning:
            Attrs are saved as comments if only_data is False. Saving attrs to file
            is not always reliable because requires converting variables to string. 
            Only attrs that are of type: string, number, and list of number and strings are 
            saved somewhat correctly. Dictionaries are not saved. 

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
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        ###################
        # check overwrite #
        ###################
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        ##########
        # kwargs #
        ##########
        if 'fmt' not in kwargs: # pick best format
            decimal = max([numanip.n_decimal_places(x) for x in arraymanip.flatten(self.data)])
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        #####################
        # header and footer #
        #####################
        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
        else:
            attrs = ['_xlim', '_ylim']
            attrs += self.get_attrs()
            header = '\n'.join(_attr2str(self, attrs, verbose)) + '\n'
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
        np.savetxt(Path(filepath), self.data, **kwargs)

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
            None

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
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        ########
        # read #
        ########
        data = np.genfromtxt(Path(filepath), **kwargs)

        ##########
        # assign #
        ##########
        self._x = data[:, 0]
        self._y = data[:, 1]

        ###############
        # read header #
        ###############
        if only_data is False:
            header = filemanip.load_Comments(Path(filepath), comment_flag=kwargs['comments'], stop_flag=kwargs['comments'])
            if header:
                for line in header:
                    if ':' not in line:
                        pass
                    else:
                        # extract name and value
                        name = line[1:-1].split(':')[0].strip()
                        try:
                            value = eval(str(':'.join(line[1:-1].split(':')[1:])).strip())
                        except:
                            value = str(':'.join(line[1:-1].split(':')[1:])).strip()
                        try:
                            setattr(self, name, value)
                        except Exception as e:
                            if verbose:
                                print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')

    #########
    # check #
    #########
    pass

    #############
    # modifiers #
    #############
    def set_horizontal_shift(self, value):
        """Shift photon events horizontally (x).

        Args:
            value (number or list): shift value by which the data are 
                shifted. If list, then it must be of the same size as the number of
                 photon events. First element will be assigned to the first photon event and so on.

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

        Returns:
            :py:class:`Image`
        """
        f = lambda x: np.polyval(p, x)
        return self.set_vertical_shift_via_function(f)

    def set_vertical_shift_via_function(self, f):
        """Set vertical shift values to f(x).

        Args:
            f (function): function where argument is x coordinates

        Returns:
            :py:class:`Image`
        """
        value = np.array([f(x) for x in self.x])
        return self.set_vertical_shift(value=value)

    def set_horizontal_shift_via_function(self, f):
        """Set horizontal shift values to f(y).

        Args:
            f (function): function where argument is y coordinates

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
        pe.copy_attrs_from(self)
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

        ###########
        # binning #
        ###########
        temp, _x_edges, _y_edges = np.histogram2d(self.x, self.y, bins=(nrows, ncols), range=(self.xlim, self.ylim))

        #########
        # Image #
        #########
        im = Image(temp.transpose())
        im.x_edges = _x_edges
        im.y_edges = _y_edges
        im.copy_attrs_from(self)

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
        s.copy_attrs_from(self)
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
        s.copy_attrs_from(self)
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
    def calculate_vertical_shift(self, ncols, nrows, mode='cc', limits=None, limit_size=1000, **kwargs):
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
            list
        """
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.calculate_vertical_shift(mode=mode, limits=limits, limit_size=limit_size, **kwargs)

    def calculate_horizontal_shift(self, ncols, nrows, mode='cc', limits=None, limit_size=1000, **kwargs):
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

            limits (None or list): a pair of values `(x_start, x_stop)`, a list 
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
            list
        """
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.calculate_horizontal_shift(mode=mode, limits=limits, limit_size=limit_size, **kwargs)


    def calculate_vertical_shift_curvature(self, ncols, nrows, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
        """Calculate vertical shift values to fix curvature.
        
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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs x centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(x_centers)
        """
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.calculate_vertical_shift_curvature(deg=deg, mode=mode, limits=limits, limit_size=limit_size, **kwargs)
     
    def calculate_horizontal_shift_curvature(self, ncols, nrows, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
        """Calculate horizontal shift values to fix curvature.
        
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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs y centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(y_centers)
        """
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.calculate_horizontal_shift_curvature(deg=deg, mode=mode, limits=limits, limit_size=limit_size, **kwargs)
     

    def fix_vertical_shift_curvature(self, ncols, nrows, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs x centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(x_centers)
        """   
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.fix_vertical_shift_curvature(deg=deg, mode=mode, limits=limits, limit_size=limit_size, **kwargs)

    def fix_horizontal_shift_curvature(self, ncols, nrows, deg=2, mode='cc', limits=None, limit_size=1000, **kwargs):
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
            s, fit, popt, R2, model

            s (spectrum): spectrum with shift values vs y centers
            
            fit (spectrum): polynomial fit spectrum of s with 100x more intepolated points

            popt (np.array): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term.

            R2 (number): R2 error

            model (function): funcion f(y_centers)
        """   
        im = self.binning(ncols=ncols, nrows=nrows)
        return im.fix_horizontal_shift_curvature(deg=deg, mode=mode, limits=limits, limit_size=limit_size, **kwargs)


    ##########################        
    # plot and visualization #
    ########################## 
    def plot(self, ax=None, s=0.1, x_start=None, x_stop=None, y_start=None, y_stop=None, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.scatter()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            s (number, optional): The marker size in points**2. Default is 0.1.
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

        return pos

# %% ============================= Dummy ================================== %% #
class Dummy():
    
    def __init__(self, *args, **kwargs):
        self._data  = []
        if len(args) > 0:
            self._data = list(args)

    def __setattr__(self, name, value):
        super().__setattr__(name, value)

    def __getattr__(self, name):
        super().__getattr__(name)

    def __getitem__(self, item):
        return self._data[item]

    def __setitem__(self, item, value):
        self._data[item] = value

    def __len__(self):
        return len(self.data)

    def __delitem__(self, item):
        del self._data[item]

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
