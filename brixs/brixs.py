#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core brixs module

Defines the main objects: Spectrum, Spectra, PhotonEvents, Image

TODO:
- Spectrum.__init__ without args and kwargs
- Spectra.__init__ without args and kwargs
- Image.__init__ without args and kwargs
- PhotonEvents.__init__ without args and kwargs

- maybe remove runtimeerror fom ss._gather_ys()
- implement s.calculate_calib()
- ss.plot(), should pvi scale for every spectrum? for now yes
- how methods deal with empty objects

- import all brixs addons at once
"""
# %% ------------------------- Standard Imports --------------------------- %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import matplotlib
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
        filepath (str or pathlib.Path, optional): filepath to read.

    How to initialize a Spectrum object:
        **Empty**
        
            >>> s = br.Spectrum()

        **From array**

            >>> s = br.Spectrum(x, y)
            >>> s = br.Spectrum(x=x, y=y)
            >>> s = br.Spectrum(y)
            >>> s = br.Spectrum(y=y)

        where ``x`` and ``y`` are 1D arrays (or list). If only y data is passed, a dummy 
        x-array will be set.

        **From file**

            >>> s = br.Spectrum(filepath=<filepath>)
            >>> s = br.Spectrum(<filepath>)

        where ``filepath`` must point to a xy-type file, where comments 
        must be marked with `#` and columns must be separated by `,` (comma).

    Attributes:
        x (array): vector with x-coordinate values
        y (array): vector with y-coordinate values
        filepath (str or pathlib.Path): filepath associated with data.

    **Computed attributes (also read-only):**
        data (array)
            2 column data (x, y). Computed when calling s.data

        step (number or None)
            step size of the x-array. Computed via s.check_step()

        monotonicity (str or None)
            monotonicity of the x-array. Computed via s.check_monotonicity()

    **Write-only attributes:**
        calib
            calls s.set_calib()

        shift
            calls s.set_shift()

        roll
            calls s.set_roll()

        factor
            calls s.set_factor()

        offset
            calls s.set_offset()
    """
    _read_only     = ['step', 'monotonicity']
    _non_removable = []

    def __init__(self, x=None, y=None, filepath=None, **kwargs):
        """Initialize the object instance.

        Raises:
            AttributeError: if kwargs and args cannot be read.
        """
        ###########################
        # Initializing attributes #
        ###########################
        # core
        self._x    = None
        self._y    = None

        # check
        self._step         = None
        self._monotonicity = None

        # modifiers
        self._calib  = 1
        self._factor = 1
        self._shift  = 0
        self._offset = 0

        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input for creating Spectrum object\n' +\
                        'Examples:\n' +\
                        '\n' +\
                        's = br.Spectrum()\n' +\
                        '\n' +\
                        's = br.Spectrum(x, y)\n' +\
                        's = br.Spectrum(x=x, y=y)\n' +\
                        '\n' +\
                        's = br.Spectrum(y)\n' +\
                        's = br.Spectrum(y=y)\n' +\
                        '\n' +\
                        's = br.Spectrum(filepath=<filepath>, **kwargs)\n' +\
                        's = br.Spectrum(<filepath>, **kwargs)\n' +\
                        '\n' +\
                        'where x and y must be a 1D array (or list), filepath ' +\
                        'must be a string or pathlib.Path object, and kwargs are ' +\
                        'passed to the load function (see br.Spectrum.load())'
        if kwargs != {} and args != ():
            raise AttributeError(error_message)
        if any([item not in ['x', 'y', 'filepath'] for item in kwargs.keys()]):
            raise AttributeError(error_message)
        if len(args) > 2 or len(kwargs) > 2:
            raise AttributeError(error_message)
        if 'x' in kwargs and 'y' not in kwargs:
            raise AttributeError(error_message)
        if ('x' in kwargs or 'y' in kwargs) and 'filepath' in kwargs:
            raise AttributeError(error_message)

        ################
        # loading data #
        ################
        # if x is available, x must be set first.
        # Otherwise, when y is set, it will create a fake x axis using np.arrange()

        # keyword arguments
        if 'y' in kwargs:
            if 'x' in kwargs:
                self.x = kwargs['x']
            self.y = kwargs['y']
            return
        if 'filepath' in kwargs:
            self.load(kwargs['filepath'])
            return

        # positional arguments
        if len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                self.load(args[0])
                return
            else:
                self.y = args[0]
                return
        elif len(args) == 2:
            self.x = args[0]
            self.y = args[1]
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

        ############################
        # reset special attributes #
        ############################
        # self.peaks.clear()
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

        ##########################
        # reset check attributes #
        ##########################
        # self._step         = None
        # self._monotonicity = None

        ############################
        # reset special attributes #
        ############################
        # self.peaks.clear()
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
        
        # turn array into monotonic
        if self.monotonicity is None:
            try:
                self.check_monotonicity()
            except ValueError:
                unqa, ID, counts = np.unique(self.x, return_inverse=True, return_counts=True)
                data        = np.column_stack((unqa, np.bincount(ID,self.y)/counts))
                
                self._x = data[:, 0]
                self._y = data[:, 1]
                # check
                self.check_monotonicity()
            
        # make it decreasing or increasing
        if self.monotonicity != mode:
            self._x = self.x[::-1]
            self._y = self.y[::-1]
        self.check_monotonicity()

        ##########################
        # reset check attributes #
        ##########################
        self._step         = None
        # self._monotonicity = None

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
                raise ValueError(f'Cannot roll data, because x-coordinates are not uniform. Use s.interp() to make data uniform.\nOr shift data using s.set_shift()')
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

    ###########
    # methods #
    ###########
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
            popt, model, R2 = polyfit(deg)

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

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        s = self._copy(limits=limits)
        x = s.x
        y = s.y

        # fit
        popt  = np.polyfit(x, y, deg=deg)
        model = lambda x: np.polyval(popt, x)
        R2 =  1 - (sum((self.y-model(self.x))**2)/sum((self.y-np.mean(self.y))**2))
    
        return popt, model, R2
   
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
        n (int, optional): array preallocation.
        data (list or array, optional): list of :py:class:`spectrum` objects.
        folderpath (string, path object): folderpath.
        filepath (str or pathlib.Path, optional): filepath to read.

    How to initialize a Spectra object:
        **Empty**

            >>> ss = br.Spectra()

        **Pre-allocation**

            >>> ss = br.Spectra(n=10)
            >>> ss = br.Spectra(10)

        **From Spectrum**

            >>> ss = br.Spectra(s1, s2, ...)
            >>> ss = br.Spectra([s1, s2, ...])
            >>> ss = br.Spectra(data=[s1, s2, ...])

        where s1, s2, ... are :py:attr:`Spectrum` type.
        
        **From file**

            >>> ss = br.Spectra(folderpath=<folderpath>)
            >>> ss = br.Spectra(<folderpath>)
            >>> ss = br.Spectra(filepath=<filepath>)
            >>> ss = br.Spectra(<filepath>)

        where filepath must point to 
        a xy-type file, where comments must be marked with `#` and columns 
        must be separated by `,` (comma). Folderpath must point to 
        a folderpath populated only with files similar xy-type files.

    Attributes:
        data (list of :py:attr:`Spectrum`): list with :py:attr:`Spectrum` objects.
        folderpath (str or pathlib.Path): folderpath associated with data.
        filepath (str or pathlib.Path): filepath associated with data.        

    **Computed (read-only) attributes:**
        step (number)
            step size for x-coordinates. Computed via ss.check_step()

        length (number)
            length of x-coordinates vector. Comuted via ss.check_length()

        x (number) 
            x-coordinates. Computed via ss.check_same_x()

        monotonicity (str or None)
            monotonicity of the x-array. Computed via ss.check_monotonicity()

        calculated_shift (list)
            calculated x additive factor. Computed via ss.calculate_shift()

        calculated_offset (list)
            calculated y additive factor. Computed via ss.calculate_offset()

        calculated_roll (list)
            calculated x roll. Computed via ss.calculate_roll()

        calculated_calib (list)
            calculated x multiplicative factor. Computed via ss.calculate_calib()

        calculated_factor (list)
            calculated y multiplicative factor. Computed via ss.calculate_factor()

    **Write-only attributes:**
        calib
            calls s.set_calib()

        shift
            calls s.set_shift()

        roll
            calls s.set_roll()

        factor
            calls s.set_factor()

        offset
            calls s.set_offset()
    """
    _read_only     = ['step', 'length', 'x', 'monotonicity',
                      'calculated_calib', 'calculated_factor',
                      'calculated_offset', 'calculated_shift',
                      'calculated_roll']
    _non_removable = []

    def __init__(self, *args, **kwargs):
        """Initialize the object instance.
        
        Raises:
            AttributeError: if kwargs and args cannot be read.
        """
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
        self._calib  = 1
        self._factor = 1
        self._shift  = 0
        self._offset = 0

        # # calculated
        # self._calculated_calib  = None
        # self._calculated_factor = None
        # self._calculated_offset = None
        # self._calculated_shift  = None
        # self._calculated_roll   = None

        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input. Spectra object cannot be created. Please, use one ' +\
                        'of the examples below to create a spectra object:\n' +\
                        '\n' +\
                        'ss = br.Spectra()\n' +\
                        'ss = br.Spectra(10)    \# pre-allocation\n' +\
                        'ss = br.Spectra(n=10)  \# pre-allocation\n' +\
                        '\n' +\
                        'ss = br.Spectra(s1, s2, ...)\n' +\
                        'ss = br.Spectra([s1, s2, ...])\n' +\
                        'ss = br.Spectra(data=[s1, s2, ...])\n' +\
                        '\n' +\
                        'ss = br.Spectra(folderpath=<folderpath>)\n' +\
                        'ss = br.Spectra(<folderpath>)\n' +\
                        '\n' +\
                        'ss = br.Spectra(filepath=<filepath>)\n' +\
                        'ss = br.Spectra(<filepath>)\n' +\
                        '\n' +\
                        's1, s2, ... must be Spectrum, filepath ' +\
                        'must be a string or pathlib.Path object that points'  +\
                        'to a file and folderpath to a folder'
        if kwargs != {} and args != ():
            raise AttributeError(error_message)
        if any([item not in ['folderpath', 'n', 'data', 'filepath'] for item in kwargs.keys()]):
            raise AttributeError(error_message)
        if len(kwargs) >= 2:
            raise AttributeError(error_message)


        ################
        # loading data #
        ################
        # keyword arguments
        if 'folderpath' in kwargs:
            self.load(kwargs['folderpath'])
            return
        elif 'filepath' in kwargs:
            self.load_from_single_file(kwargs['filepath'])
            return
        elif 'n' in kwargs:
            n = int(kwargs['n'])
            self._data = [Spectrum()]*n
            return
        elif 'data' in kwargs:
            self.data = kwargs['data']
            return

        # positional arguments
        if len(args) == 0:
            return
        elif len(args) == 1:
            if isinstance(args[0], int):
                n = int(args[0])
                self._data = [-1]*n
                return
            elif isinstance(args[0], str) or isinstance(args[0], Path):
                temp = Path(args[0])
                assert temp.exists(), f'filepath (or folderpath) do not exist\n{temp}'
                if temp.is_file():
                    self.load_from_single_file(temp)
                elif temp.is_dir():
                    self.load(temp)
                return
            elif isinstance(args[0], Iterable):
                self.data = args[0]
                return
            else:
                raise AttributeError(error_message)
        elif len(args) > 1:
            self.data = args
            return
        else:
            raise AttributeError(error_message)

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
            for a in attrs2reorder:
                temp = self.__getattribute__(attr)
                assert isinstance(temp, Iterable), f'{a} must be an iterable type'
                assert len(temp) == len(self), f'Lenght of attr {a} must be the same as the number of spectra.\nlenght of attr: {len(a)}\nnumber of spectra: {len(self)}'

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

    ###########
    # methods #
    ###########
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

    def calculate_map(self, axis=0, centers=None, limits=None):
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
            centers (list, optional): numerical value for each spectrum. 
                If `axis=0`, this is equivalent as using im.x_centers=centers. 
                If `axis=1`, im.y_centers=centers. Default is None, i.e., each
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
        # check axis
        if axis != 0 and axis != 1:
            raise ValueError('axis must be 0 or 1')
        
        # check centers
        if centers is not None:
            assert len(centers) == len(self), f'centers must have the same number of items as the number of spectra.\nnumber of centers: {len(centers)}\nnumber of spectra: {len(self)}'

        # gather ys
        y, ys = self._gather_ys(limits=limits)
        x = centers

        if axis == 1:
            ys = ys.transpose()
            x = y
            y = centers

        im = Image(data=ys)
        im.copy_attrs_from(self)
        im.x_centers = x
        im.y_centers = y

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
    """Returns a ``Image`` object.

    Args:
        data (2D array, optional): Image.
        filepath (string or path object, optional): filepath.

    How to initialize a Image object:
        **Empty**

            >>> im = br.Image()

        **From 2D array**

            >>> im = br.Image(data)
            >>> im = br.Image(data=data)

        **From file**
            
            >>> im = br.Image(<filepath>)
            >>> im = br.Image(filepath=<filepath>)

        where ``filepath`` must point to a 2D xy-type file, where comments 
        must be marked with `#` and columns must be separated by `,` (comma).


    Attributes:
        data (2D np.array): This is where we store the Image.
        x_centers, y_centers (np.array): pixel center labels.

    Computed (read-only) attributes:
        shape (tuple, read only): Shape of data (vertical size, horizontal size),
            i.e., number of pixels rows and number of columns, respectively.
        histogram (brixs.Spectrum)
            Data intensity histogram.
        calculated_roll (list)
            Calculated rolls.
        calculated_shift (list)
            same as calculated_shift, but in terms of x or y centers.
        columns, rows (brixs.Spectra)
            Spectra obtained from each pixel columns or rows, respectively.

    Write-only attributes:
        None
    """
    # read only and non-removable arguments
    _read_only     = ['calculated_roll', 'calculated_shift']
    _non_removable = []
    
    def __init__(self, *args, **kwargs):
        """Initialize the object instance.
        
        Raises:
            AttributeError: if kwargs and args cannot be read.
        """
        ###########################
        # Initializing attributes #
        ###########################
        # basic
        self._data      = None
        self._x_centers = None
        self._y_centers = None
        self._x_edges   = None
        self._y_edges   = None
        self._filepath  = ''

        self._calculated_roll  = None
        self._calculated_shift = None

        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input. Image object cannot be created. Please, use one ' +\
                        'of the examples below to create a Image object:\n' +\
                        '\n' +\
                        'im = br.Image()\n' +\
                        '\n' +\
                        'im = br.Image(data)\n' +\
                        'im = br.Image(data=data)\n' +\
                        '\n' +\
                        'im = br.Image(filepath=<filepath>)\n' +\
                        'im = br.Image(<filepath>)\n' +\
                        '\n' +\
                        'data must be a 2D array and filepath ' +\
                        'must be a string or pathlib.Path object'
        if kwargs != {} and args != ():
            raise AttributeError(error_message)
        if any([item not in ['data', 'filepath'] for item in kwargs.keys()]):
            raise AttributeError(error_message)
        if len(args) > 1 or len(kwargs) > 1:
            raise AttributeError(error_message)

        ################
        # loading data #
        ################
        # keyword arguments
        if 'data' in kwargs:
            self.data = kwargs['data']
        elif 'filepath' in kwargs:
            self.load(kwargs['filepath'])
        # positional arguments
        elif len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                self.load(args[0])
            else:
                self.data = args[0]

    ##############
    # attributes #
    ##############
    @property
    def data(self):
        return copy.deepcopy(self._data)
    @data.setter
    def data(self, value):
        # basic attr
        self._data  = np.array(value, dtype='float')
        self.x_centers = None
        self.y_centers = None        
        self._calculated_roll  = None
        self._calculated_shift = None
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
            if len(value) != self.shape[1]:
                raise ValueError(f"Length of x must be the same as the number of pixels in the horizontal direction.\nNumber of pixels = {self.data.shape[1]}\nLength of the array = {len(value)}")
        else:
            raise ValueError(f"x must be None or an iterable (list, tuple, or 1D array)")
        self._x_centers = np.array(value, dtype='float')
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
            if len(value) != self.shape[0]:
                raise ValueError(f"Length of y must be the same as the number of pixels in the vertical direction.\nNumber of pixels = {self.data.shape[0]}\nLength of the array = {len(value)}")
        else:
            raise ValueError(f"y must be None or an iterable (list, tuple, or 1D array)")
        self._y_centers = np.array(value, dtype='float')
    @y_centers.deleter
    def y_centers(self):
        self._y_centers  = np.arange(0, self.data.shape[0])

    @property
    def x_edges(self):
        return self._x_edges
    @x_edges.setter
    def x_edges(self, value):
        raise NotImplementedError('this is not implemented yet')
    @x_edges.deleter
    def x_edges(self):
        raise NotImplementedError('this is not implemented yet')

    @property
    def y_edges(self):
        return self._y_edges
    @y_edges.setter
    def y_edges(self, value):
        raise NotImplementedError('this is not implemented yet')
    @y_edges.deleter
    def y_edges(self):
        raise NotImplementedError('this is not implemented yet')

    @property
    def filepath(self):
        return self._filepath
    @filepath.setter
    def filepath(self, value):
        if value is None:
            value = ''
        elif isinstance(value, str) or isinstance(value, Path):
            self._filepath = value
        else:
            raise TypeError(r'Invalid type ' + str(type(value)) + 'for filepath\nfilepath can only be str or pathlib.Path type.')
    @filepath.deleter
    def filepath(self):
        raise AttributeError('Cannot delete object.')  
    
    ###################################
    # computed (read-only) attributes #
    ###################################
    @property
    def shape(self):
        if self._data is None:
            return None
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
        return ss
    @rows.setter
    def rows(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @rows.deleter
    def rows(self):
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

        # ##################
        # # transfer attrs #
        # ##################
        # for attr in self.get_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     final.__setattr__(attr, value)

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

        # ##################
        # # transfer attrs #
        # ##################
        # for attr in self.get_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     final.__setattr__(attr, value)

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
        
        # ##################
        # # transfer attrs #
        # ##################
        # for attr in self.get_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     final.__setattr__(attr, value)

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
        
        # ##################
        # # transfer attrs #
        # ##################
        # for attr in self.get_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     final.__setattr__(attr, value)

        return final
    
    def __truediv__(self, object):
        return self.__div__(object)

    
    #########
    # attrs #
    #########
    def get_attrs(self):
        """return attrs that are user defined.""" 
        return [key for key in self.__dict__.keys() if key.startswith('_') == False and key not in settings._reserved_words]
    
    def _get_center_attrs(self):
        """return center attrs"""
        return ['_x_centers', '_y_centers', '_x_edges', '_y_edges']

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

    def copy_centers_from(self, s):
        """Copy center attributes from another brixs.Image object.

        Args:
            s (brixs object): Image to copy center attributes from.
        
        Returns:
            None
        """
        # check type
        if isinstance(s, Image):
            pass
        else:
            raise TypeError(f'type {type(s)} not valid\nCan only copy user attrs from type br.Image to type br.Image')

        # transfer attrs
        for attr in s._get_center_attrs():
            value = copy.deepcopy(s.__dict__[attr])
            self.__setattr__(attr, value)

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
    
    #################
    # save and load #
    #################
    def _create_header(self, verbose=False):
        """Gather attrs to be saved to a file."""
        header = ''
        attrs  = self.get_attrs()
        attrs += self._get_center_attrs()
        for name in attrs:
            try:
                if self.__dict__[name] is None:
                    header += f'{name}: None'  + '\n'
                elif isinstance(self.__dict__[name], str):
                    temp2 = str(self.__dict__[name]).replace('\n','\\n')
                    header += f'{name}: \"{temp2}\"'  + '\n'
                elif isinstance(self.__dict__[name], dict) or isinstance(self.__dict__[name], MutableMapping):
                    if verbose:
                        type_ = str(type(self.__dict__[name]))
                        print('\nWarning: Cannot save attr of type: ' + type_ + '\nattr name: '+ name + '\nTo turn off this warning, set verbose to False.')
                elif isinstance(self.__dict__[name], Iterable):
                    header += f'{name}: {list(self.__dict__[name])}'  + '\n'
                elif numanip.is_number(self.__dict__[name]):
                    tosave = str(self.__dict__[name])
                    if tosave[-1] == '\n':
                        tosave = tosave[:-1]
                    header += f'{name}: {tosave}'  + '\n'
                else:
                    temp2 = str(self.__dict__[name]).replace('\n','\\n')
                    header += f'{name}: \"{temp2}\"'  + '\n'
            except:
                if verbose:
                    type_ = str(type(self.__dict__[name]))
                    print('\nWarning: Cannot save attr of type: ' + type_ + '\nattr name: '+ name + '\nTo turn off this warning, set verbose to False.')
        return header[:-1]
    
    def save(self, filepath=None, only_data=False,  check_overwrite=False, verbose=False, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Attrs are saved as comments if only_data is False. Saving attrs to file
        is tricky because requires converting variables to string. Only attrs 
        that are of type: string, number, and list of number and strings are 
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

        .. _numpy.savetxt(): https://numpy.org/doc/stable/reference/generated/numpy.savetxt.html
        """
        # filepath
        if filepath is None:
            if filepath == '':
                raise TypeError("Missing 1 required argument: 'filepath'")
            else:
                filepath = self.filepath               
        filepath = Path(filepath)

        # check if filepath points to a file
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'
        
        # check overwrite
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([numanip.n_decimal_places(x) for x in arraymanip.flatten(self._data)])
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
        else:
            if 'header' not in kwargs:
                kwargs['header'] = self._create_header(verbose=verbose)
            else:
                if kwargs['header'] == '':
                    kwargs['header'] = self._create_header(verbose=verbose)
                elif kwargs['header'][-1] != '\n':
                    kwargs['header'] += '\n'
        np.savetxt(Path(filepath), self._data, **kwargs)

        # save filepath
        self.filepath = filepath

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
        ###########
        # default #
        ###########
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # check if filepath points to a file
        filepath = Path(filepath)
        assert filepath.is_file(), f'filepath must point to a file\n{filepath}'
        assert filepath.exists(), f'filepath does not exist\n{filepath}'

        #############
        # read data #
        #############
        data = np.genfromtxt(Path(filepath), **kwargs)
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
                                print(f'Error loading attribute: {name}\nvalue: {value}\nAttribute not set.\n{e}\n')
        
        # save filepath
        self.filepath = str(filepath)

    ##################
    # base modifiers #
    ##################
    def set_roll(self, value=None, p=None, f=None, axis=0):
        """Roll array of pixels along a given axis.

        Elements that roll beyond the edge are re-introduced at the first position.

        If value, p, and f are None, values are read from im.calculated_roll.

        Roll values must be an integer. Float values will be rounded to int.

        Args:
            value (int or tuple): The number of pixels by which the data are shifted.
                If a tuple, then it must be of the same size as the number of
                columns (``for axis=0``) or rows (``for axis=1``). If elements are
                not int, it will rounded to an integer value. Default is None.
            p (np.array, read only): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term. Polynomial as a function of x_centers for 
                axis=0 or y_centers for axis=1. Default is None.
            f (function, read only): funcion f(x_centers) for axis=0 or 
                f(y_centers) for axis=1. Default is None.
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.

        Returns:
            None
        """
        # check axis and define centers
        if axis == 0:
            centers = self.x_centers
        elif axis == 1:
            centers = self.y_centers
        else:
            raise ValueError('axis must be 0 or 1')

        # calculate value
        if value is not None:
            if isinstance(value, Iterable):
                assert len(value) == len(centers), f'Number of values ({len(value)}) must be the same as the number of columns ({self.shape[1]}) (for axis=1) or rows ({self.shape[0]}) (for axis=0)'
            else:
                value = [value]*len(centers)
        elif p is not None:
            value = np.polyval(p, centers)
        elif f is not None:
            value = f(centers)
        else:
            value = self.calculated_roll
        value = [int(round(k)) for k in value]

        # roll
        if axis == 0:
            for i, v in enumerate(value):
                if v != 0:
                    s = Spectrum(self._data[:, i])
                    s.roll = v
                    self._data[:, i] = s.y
        elif axis == 1:
            for i, h in enumerate(value):
                if h != 0:
                    s = Spectrum(self._data[i, :])
                    s.roll = h
                    self._data[i, :] = s.y

    #############
    # modifiers #
    #############
    def floor(self, x=0, y=0, n=30, nx=None, ny=None):
        """Set background intensity to zero.

        Args:
            x, y (int, optional): x and y position to sample background intensity.
            n, nx, ny (int, optional): size of the pixel window around x, y.

        Returns:
            None
        """
        assert x >= 0 and numanip.is_integer(x) and x<self.shape[1], f'x must be a positive integer smaller than {self.shape[1]}.'
        assert y >= 0 and numanip.is_integer(y) and y<self.shape[0], f'y must be a positive integer smaller than {self.shape[0]}.'

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

        self._data -= np.mean(self._data[int(y_start):int(y_stop), int(x_start):int(x_stop)]).astype(self.data.dtype)

    def fix_curvature(self, deg=2, axis=0):
        """Roll array of pixels along a given axis to fix curvature via cc.

        Elements that roll beyond the edge are re-introduced at the first position.

        Args:
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.

        Returns:
            p (np.array, read only): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term. Polynomial as a function of x_centers for 
                axis=0 or y_centers for axis=1. Default is None.
            f (function, read only): funcion f(x_centers) for axis=0 or 
                f(y_centers) for axis=1.
        """
        # check axis and define centers
        if axis == 0:
            centers = self.x_centers
        elif axis == 1:
            centers = self.y_centers
        else:
            raise ValueError('axis must be 0 or 1')

        # calculate shifts
        im = self.copy()
        im.floor()
        im.calculate_roll(axis=axis)

        # calculate poly
        popt = np.polyfit(x=centers, y=im.calculated_roll, deg=deg)
        model = lambda x: np.polyval(popt, x)

        # set roll
        self.set_roll(p=popt, axis=axis)

        return popt, model
    
    ##############
    # extractors #
    ##############
    def copy(self, *args, **kwargs):
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
        error_message = 'Wrong input. Please check input.\n' +\
                        f'args passed: {args}\n'  +\
                        f'kwargs passed: {kwargs}\n\n'  +\
                        ' === Usage ===\n'  +\
                        'im2 = im1.copy()  # im2 is now a copy of im1\n'  +\
                        'im3 = im1.copy(x_start, x_stop, y_start, y_stop)\n'
        ###################################
        # check if crop limits are passed #
        ###################################
        x_start = False
        if 'x_start' in kwargs or 'x_stop' in kwargs or 'y_start' in kwargs or 'y_stop' in kwargs:
            assert 'x_start' in kwargs and 'x_stop' in kwargs and 'y_start' in kwargs and 'y_stop' in kwargs, error_message
            x_start = kwargs['x_start']
            x_stop  = kwargs['x_stop']
            y_start = kwargs['y_start']
            y_stop  = kwargs['y_stop']
        elif len(args) == 4:
            x_start = args[0]
            x_stop  = args[1]
            y_start = args[2]
            y_stop  = args[3]
        elif len(args) != 4 and len(args) != 0:
            raise ValueError(error_message)
        elif len(kwargs) != 0:
            for x in list(kwargs.keys()):
                if x not in ['x_start', 'x_stop', 'y_start', 'y_stop']:
                    raise ValueError(f'{x} is not a recognized input for copy function\n'+error_message)
        
        ##################
        # identical copy #
        ##################
        if x_start == False:
            im = Image(data=self.data)
            im.copy_centers_from(self)
        #############################
        # if crop limits are passed #
        #############################
        else:
            im = self.copy()
            im.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)

        ##################
        # transfer attrs #
        ##################
        # im.copy_centers_from(self)
        im.copy_attrs_from(self)

        return im
    
    def binning(self, *args, **kwargs):
        """Compute the 2D histogram of the data (binning of the data).

        Usage:

            >>> # 10 rows, 10 columns
            >>> im.binning(nbins=10)       
            >>> im.binning(nbins=(10))
            >>> im.binning(10)           
            >>> im.binning((10))  
            >>>
            >>> # 10 rows, 5 columns
            >>> im.binning(nbins=(10, 5)) 
            >>> im.binning((10, 5))       
            >>> im.binning(10, 5)   
            >>>
            >>> # No binning, 5 columns
            >>> im.binning(None, 5)    

        Args:
            nbins (int or tuple): number of bins. Must be non-zero positive number.
                If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the number of rows and number of
                columns, respectively. If the number of pixels in the image
                cannot be divided by the selected number of bins, it will raise 
                an error. Use None to idicate no binning. 

        Returns:
            binned image
        """
        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input. Please, see examples below:\n' +\
                        '\n' +\
                        'im.binning(nbins=10)\n' +\
                        'im.binning(nbins=(10))\n' +\
                        'im.binning(10)\n' +\
                        'im.binning((10))\n' +\
                        '\n' +\
                        'im.binning(nbins=(10, 5))\n' +\
                        'im.binning((10, 5))\n' +\
                        'im.binning(10, 5)\n' +\
                        'im.binning(None, 5)\n'
        if kwargs != {} and args != ():
            raise AttributeError(error_message)
        if any([item not in ['nbins', ] for item in kwargs.keys()]):
            raise AttributeError(error_message)
        if len(args) > 2 or len(kwargs) > 1:
            raise AttributeError(error_message)
       
        #################
        # sorting input #
        #################
        # keyword arguments
        if 'nbins' in kwargs:
            args = kwargs['nbins']
        # positional arguments
        if isinstance(args, Iterable):
            if len(args) == 1:
                if isinstance(args[0], Iterable):
                    nbins = list(args[0])
                else:
                    nbins = [args[0], args[0]]
            else:
                nbins = list(args)
        else:
            nbins = [args, args]

        # fix format
        if nbins[0] is None:
            nbins[0] = self.shape[0]
        if nbins[1] is None:
            nbins[1] = self.shape[1]
        if numanip.is_integer(nbins[0]) == False or numanip.is_integer(nbins[1])  == False or nbins[0] < 0 or nbins[1] < 0:
            raise ValueError("Number of bins must be a positive integer.")

        # is divisible
        assert self.shape[1] % nbins[1] == 0, f"The {self.shape[1]} pixels in a row is not evenly divisible by {nbins[1]}\nPlease, pick one of the following numbers: {np.sort(list(numanip.factors(self.shape[1])))}"
        assert self.shape[0] % nbins[0] == 0, f"The {self.shape[0]} pixels in a column is not evenly divisible by {nbins[0]}\nPlease, pick one of the following numbers: {np.sort(list(numanip.factors(self.shape[0])))}"

        ###############
        # Calculation #
        ###############
        _bins_size = np.array((self.shape[0]/nbins[0], self.shape[1]/nbins[1]))
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
        # reduced.copy_centers_from(self)
        reduced.copy_attrs_from(self)

        return reduced
    
    def transpose(self):
        """Transpose image

        Returns:
            None
        """
        im = Image(data=self.data.transpose())

        ##################
        # transfer attrs #
        ##################
        # im.copy_centers_from(self)
        im.copy_attrs_from(self)
        im.x_centers = copy.deepcopy(self.y_centers)
        im.y_centers = copy.deepcopy(self.x_centers)       

        return im

    def crop(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Crop Image.

        Args:
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.

        Returns:
            None
        """
        #################
        # check if None #
        #################
        if x_start is None: x_start = self.x_centers[0]
        if x_stop  is None: x_stop  = self.x_centers[-1]
        if y_start is None: y_start = self.y_centers[0]
        if y_stop  is None: y_stop  = self.y_centers[-1]

        ####################
        # centers to index #
        ####################
        x_start = arraymanip.index(self.x_centers, x_start)
        x_stop  = arraymanip.index(self.x_centers, x_stop) + 1
        y_start = arraymanip.index(self.y_centers, y_start)
        y_stop  = arraymanip.index(self.y_centers, y_stop) + 1

        ##################
        # validate input #
        ##################
        assert x_stop > x_start, f'x_start must be smaller than x_stop.'
        assert y_stop > y_start, f'y_start must be smaller than y_stop.'

        ########
        # crop #
        ########
        im = Image(data=self.data[int(y_start):int(y_stop), int(x_start):int(x_stop)])

        ##################
        # transfer attrs #
        ##################
        # im.copy_centers_from(self)
        im.copy_attrs_from(self)
        im._x_centers = self.x_centers[x_start:x_stop]
        im._y_centers = self.y_centers[y_start:y_stop]

        return im

    def calculate_histogram(self, nbins=None):
        """Compute the histogram of data. Wrapper for `numpy.histogram()`_.

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
            s = Spectrum(x=self.x_centers, y=np.sum(self._data, axis=0))
        elif axis == 1:
            s = Spectrum(x=self.y_centers, y=np.sum(self._data, axis=1))
        
        ##################
        # transfer attrs #
        ##################
        s.copy_attrs_from(self)

        return s

    ########################
    # calculation and info #
    ########################
    def possible_nbins(self):
        """return possible values for nbins in the y (nrows) and x (ncols) directions."""
        return np.sort(list(numanip.factors(self.shape[0]))), np.sort(list(numanip.factors(self.shape[1])))

    def calculate_roll(self, axis=0, limit_size=1000):
        """Calculate intensity misalignments via cross-correlation.

        Args:
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.

        Returns:
            None. Values are save at: im.calculated_roll
        """
        # check axis
        if axis != 0 and axis != 1:
            raise ValueError('axis must be 0 or 1')

        # select axis
        if axis == 0:
            ss = self.columns
            if limit_size:
                if len(ss) > limit_size:
                    raise ValueError(f'Number of columns is bigger than limit_size.\nImage is seems to be too big.\nAre you sure you want to calculate shifts for such a big image.\nIf so, either set limit_size to False or a higher value.\nNumber of columns: {len(self.x_centers)}\nlimit size: {limit_size}')
        elif axis == 1:
            ss = self.rows
            if limit_size:
                if len(ss) > limit_size:
                    raise ValueError(f'Number of rows is bigger than limit_size.\nImage is seems to be too big.\nAre you sure you want to calculate shifts for such a big image.\nIf so, either set limit_size to False or a higher value.\nNumber of columns: {len(self.y_centers)}\nlimit size: {limit_size}')

        # calculate
        ss.calculate_roll(mode='cc')
        self._calculated_roll  = ss.calculated_roll
        self._calculated_shift = ss.calculated_roll*ss.step

    def calculate_curvature(self, deg=2, axis=0):
        """Calculate roll values along a given axis to fix curvature via cc.

        Args:
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.

        Returns:
            p (np.array, read only): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term. Polynomial as a function of x_centers for 
                axis=0 or y_centers for axis=1. Default is None.
            f (function, read only): funcion f(x_centers) for axis=0 or 
                f(y_centers) for axis=1.
        """
        # check axis and define centers
        if axis == 0:
            centers = self.x_centers
        elif axis == 1:
            centers = self.y_centers
        else:
            raise ValueError('axis must be 0 or 1')

        # calculate shifts
        im = self.copy()
        im.floor()
        im.calculate_roll(axis=axis)

        # calculate poly
        popt = np.polyfit(x=centers, y=im.calculated_roll, deg=deg)
        model = lambda x: np.polyval(popt, x)

        return popt, model
    
    ##########################        
    # plot and visualization #
    ########################## 
    def pcolormesh(self, ax=None, colorbar=False, **kwargs):
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
        # colorbar divider
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

        # kwargs
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'
        if 'vmin' not in kwargs or 'vmax' not in kwargs:
            vmin, vmax = self._calculated_vmin_vmax()
            if 'vmin' not in kwargs:
                kwargs['vmin'] = vmin
            if 'vmax' not in kwargs:
                kwargs['vmax'] = vmax

        # fix monotonicity of labels x
        if arraymanip.check_monotonicity(self.x_centers) != 1:
            x, ordering = arraymanip.fix_monotonicity(self.x_centers, np.arange(len(self.x_centers)), mode='increasing')
            assert len(x)==len(self.x_centers), f'Cannot plot when Image.x have repeated elements.\nEither fix Image.x_centers or set it to None.\nx: {self.x_centers}'
            ordering = [int(i) for i in ordering]
            data = copy.deepcopy(self.data)
            for i in ordering:
                if i != ordering[i]:
                    data[:, i] = self.data[:, ordering[i]]
        else:
            x    = self.x_centers
            data = self.data

        # fix monotonicity of labels y
        if arraymanip.check_monotonicity(self.y_centers) != 1:
            y, ordering = arraymanip.fix_monotonicity(self.y_centers, np.arange(len(self.y_centers)), mode='increasing')
            assert len(y)==len(self.y_centers), f'Cannot plot when Image.y have repeated elements.\nEither fix Image.y_centers or set it to None.\ny: {self.y_centers}'
            ordering = [int(i) for i in ordering]
            data2 = copy.deepcopy(data)
            for i in ordering:
                if i != ordering[i]:
                    data2[i, :] = data[ordering[i], :]
        else:
            y     = self.y_centers
            data2 = data

        # plot
        X, Y = np.meshgrid(x, y)
        pos  = ax.pcolormesh(X, Y, data2, **kwargs)

        # colorbar
        if colorbar:
            if divider:
                divider = make_axes_locatable(ax)
                ax_cb = divider.append_axes("right", size="5%", pad=0.05)
                plt.colorbar(pos, cax=ax_cb)
            else:
                plt.colorbar(pos, aspect=50)

        # add edges and centers to quadmesh
        pos.x_edges = pos.get_coordinates()[0][:, 0]
        pos.y_edges = pos.get_coordinates()[:, 1][:, 1]
        pos.x_centers = self.x_centers
        pos.y_centers = self.y_centers

        return pos

    def imshow(self, ax=None, xlim=None, ylim=None, colorbar=False, origin='lower', verbose=True, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.imshow()`_.

        Warning:
            Pixels are always square. For irregular pixel row/columns, see Image.pcolormesh()

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            colorbar (bool, optional): if True, colorbar is shown on the right side.
                 (str, optional): Location of the [0, 0] index. Options are
                `upper` and `lower`. Default is 'lower'.
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
        # colorbar divider
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
        
        # ylim, xlim
        if xlim is not None and ylim is None:
            im = self.crop(x_start=xlim[0], x_stop=xlim[1])
        elif ylim is not None and xlim is None:
            im = self.crop(y_start=ylim[0], y_stop=ylim[1])
        elif xlim is not None and ylim is not None:
            im = self.crop(x_start=xlim[0], x_stop=xlim[1], y_start=ylim[0], y_stop=ylim[1])
        else:
            im = self.copy()

        # default arguments
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

        # check monotonicity x
        if arraymanip.check_monotonicity(im.x_centers) != 1:
            x, ordering = arraymanip.fix_monotonicity(im.x_centers, np.arange(len(im.x_centers)), mode='increasing')
            assert len(x)==len(im.x_centers), f'Cannot plot when Image.x have repeated elements.\nEither fix Image.x or set it to None.\nx: {self.x_centers}'
            ordering = [int(i) for i in ordering]
            data     = copy.deepcopy(im.data)
            for i in ordering:
                if i != ordering[i]:
                    data[:, i] = im.data[:, ordering[i]]
                extent_x = [min(x), max(x)]
        else:
            x = im.x_centers
            x = np.linspace(x[0], x[-1], len(x))
            dx  = np.mean(np.diff(x))
            extent_x = [x[0]-dx/2, x[-1]+dx/2]
            data = im.data

        # check monotonicity y
        if arraymanip.check_monotonicity(im.y_centers) != 1:
            y, ordering = arraymanip.fix_monotonicity(im.y_centers, np.arange(len(im.y_centers)), mode='increasing')
            assert len(y)==len(im.y_centers), f'Cannot plot when Image.y have repeated elements.\nEither fix Image.x or set it to None.\ny: {self.y_centers}'
            ordering = [int(i) for i in ordering]
            data2 = copy.deepcopy(data)
            for i in ordering:
                if i != ordering[i]:
                    data2[i, :] = data[ordering[i], :]
            extent_y = [min(y), max(y)]
        else:
            y = im.y_centers
            y = np.linspace(y[0], y[-1], len(y))
            dy  = np.mean(np.diff(y))
            extent_y = [y[0]-dy/2, y[-1]+dy/2]
            data2 = data

        if 'extent' not in kwargs:
            kwargs['extent'] = np.append(extent_x, extent_y)

        # check irregular spacing (issues a warning)
        if verbose:
            sx = Spectrum(x=im.x_centers, y=im.x_centers)
            sy = Spectrum(x=im.y_centers, y=im.y_centers)
            try:
                sx.check_step()
                sy.check_step()
            except ValueError:
                print('Data seems to have irregular pixel size. Maybe plot it using Image.pcolormesh().\nTo turn off this warning set verbose to False.')

        # plot
        pos = ax.imshow(data2, **kwargs)

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

    def plot(self, *args, **kwargs):
        """Same as Image.imshow.

        See also:
            :py:func:`Image.imshow`
        """
        return self.imshow(*args, **kwargs)

    # experimental 
    def linecuts(self, axis=0, xlim=None, ylim=None, ilim=None, keep=None, **kwargs):
        """[EXPERIMENTAL] Plot image and flip trhough linecuts with keyboard arrows.

        Args:
            axis (int or string, optional): Axis along linecuts.
                By default, linecuts are in the vertical (0) direction.
            xlim, ylim, ilim (tuple, optional): Image ploting limits.
            keep (int or str, optional): index of the spectrum to keep on screen. 
                Use 'previous' or 'next' to show previous or next spectrum.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            None
        """
        # axis
        if axis == 0:
            ss = self.columns
        elif axis == 1:
            raise NotImplementedError('Not implemented yet for axis = 1')
            # ss = self.rows
        else:
            raise ValueError('axis must be 0 or 1')
    
        # change keybindings
        try:
            matplotlib.rcParams['keymap.back'].remove('left')
            matplotlib.rcParams['keymap.forward'].remove('right')
        except ValueError:
            pass
        
        # vars
        self.__i    = 0
        self.__axes = None
        self.__xlim = xlim
        self.__ylim = ylim
        self.__ilim = ilim

        # lims
        if self.__xlim is None:
            self.__xlim = [min(self._x_centers), max(self._x_centers)]
        if self.__ylim is None:
            self.__ylim = [min(self._y_centers), max(self._y_centers)]
        if self.__ilim is None:
            self.__ilim     = [min([min(l) for l in self.data]), max([max(l) for l in self.data])]

        # keep
        if keep is not None:
            if isinstance(keep, str):
                assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
            else:
                assert abs(keep) < len(ss), 'keep must be a valid spectrum index, or "previous/next"'

        # kwargs
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        # update function
        def update(ss):
            # plot spectrum
            ax = self.__axes[0]
            if axis == 0:
                ax.set_title(f'{self.__i}: {self.x_centers[self.__i]} eV')
            if keep is not None:
                if keep == 'next':
                    if self.__i+1 < len(self):
                        ss[self.__i+1].plot(color='red', alpha=0.5)
                elif keep == 'previous':
                    if self.__i-1 >= 0:
                        ss[self.__i-1].plot(color='red', alpha=0.5)
                else:
                    ss[keep].plot(color='red', alpha=0.5)
            ss[self.__i].plot(ax=ax, **kwargs)

            # lim
            if self.__ylim is not None:
                ax.set_xlim(self.__ylim)
            if self.__ilim is not None:
                ax.set_ylim(self.__ilim)



            # plot map (shouldn't we put this outside of the update function?)
            ax = self.__axes[1]
            
            # lim
            if self.__ilim is not None:
                self.pcolormesh(ax=ax, vmin=self.__ilim[0], vmax=self.__ilim[1])
            else:
                self.pcolormesh(ax=ax)
            if self.__xlim is not None:
                if axis == 0:
                    ax.set_xlim(self.__xlim)
                else:
                    pass
                    # ax.set_xlim(self.__ylim)
            if self.__ylim is not None:
                if axis == 0:
                    ax.set_ylim(self.__ylim)
                else:
                    pass
                    # ax.set_ylim(self.__ylim)

            # vline
            if axis == 0:
                if keep is not None:
                    if keep == 'next':
                        if self.__i+1 < len(self):
                            ax.plot([self.x_centers[self.__i+1], self.x_centers[self.__i+1]], [self.__ylim[0], self.__ylim[1]], color='red', alpha=0.5)
                    elif keep == 'previous':
                        if self.__i-1 >= 0:
                            ax.plot([self.x_centers[self.__i-1], self.x_centers[self.__i-1]], [self.__ylim[0], self.__ylim[1]], color='red', alpha=0.5)
                    else:
                        ax.plot([self.x_centers[keep], self.x_centers[keep]], [self.__ylim[0], self.__ylim[1]], color='red', alpha=0.5)
                ax.plot([self.x_centers[self.__i], self.x_centers[self.__i]], [self.__ylim[0], self.__ylim[1]], color='white')
                # E = self.x_centers[self.__i]
                # figmanip.vlines(E, ax=ax, color='red', lw=1)
            else:
                E = self.y_centers[self.__i]
                figmanip.hlines(E, ax=ax, color='red', lw=1)

            

        def _update(ss):
            if self.__i >= len(ss):
                self.__i = len(ss) - 1
            elif self.__i < 0:
                self.__i = 0

            update(ss)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'right' or event.key == 'up':
                self.__i = self.__i + 1

                for ax in self.__axes:
                    # ax.cla()
                    ax.lines.clear()
                _update(ss)
                plt.draw()
            elif event.key == 'left' or event.key == 'down':
                self.__i = self.__i - 1
                
                for ax in self.__axes:
                    # ax.cla()
                    ax.lines.clear()
                _update(ss)
                plt.draw()

        # axis zoom changes
        def on_1_xlims_change(event_ax):
            # print("updated xlims: ", event_ax.get_xlim())
            self.__xlim = event_ax.get_xlim()

        def on_1_ylims_change(event_ax):
            # print("updated ylims: ", event_ax.get_ylim())
            self.__ylim = event_ax.get_ylim()

        def on_0_xlims_change(event_ax):
            # print("updated xlims: ", event_ax.get_xlim())
            self.__ylim = event_ax.get_xlim()

        def on_0_ylims_change(event_ax):
            # print("updated ylims: ", event_ax.get_ylim())
            self.__ilim = event_ax.get_ylim()

        # figure
        fig, self.__axes = figmanip.subplots(nrows=2, ncols=1)
        fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=ss))
        _update(ss)

        self.__axes[0].callbacks.connect('xlim_changed', on_0_xlims_change)
        self.__axes[0].callbacks.connect('ylim_changed', on_0_ylims_change)
        self.__axes[1].callbacks.connect('xlim_changed', on_1_xlims_change)
        self.__axes[1].callbacks.connect('ylim_changed', on_1_ylims_change)

    # experimental 
    def linecuts2(self, axis=0, xlim=None, ylim=None):
        """[EXPERIMENTAL] Plot image and flip trhough linecuts with keyboard arrows.

        Warning: old version, with auto rescaling.

        Args:
            axis (int or string, optional): Axis along linecuts.
                By default, linecuts are in the vertical (0) direction.
            xlim, ylim (tuple, optional): spectra ploting limits (not image 
                ploting limits).

        Returns:
            None
        """
        # axis
        if axis == 0:
            ss = self.columns
        elif axis == 1:
            ss = self.rows
        else:
            raise ValueError('axis must be 0 or 1')
        
        # update function
        def update(ss):
            # plot
            ax = axes[0]
            ss[ss.__i].plot(ax=ax, color='black', marker='o')

            # lim
            if xlim is not None:
                ax.set_xlim(xlim)
            if ylim is not None:
                ax.set_ylim(ylim)
            # br.label_rixs(ax=ax)

            # plot map (shouldn't we put this outside of the update function?)
            ax = axes[1]
            
            if ylim is not None:
                self.pcolormesh(ax=ax, vmin=ylim[0], vmax=ylim[1])
            else:
                self.pcolormesh(ax=ax)
            if xlim is not None:
                if axis == 0:
                    ax.set_ylim(xlim)
                else:
                    ax.set_xlim(xlim)

            # vline
            if axis == 0:
                E = self.x_centers[ss.__i]
                figmanip.vlines(E, ax=ax, color='red', lw=1)
            else:
                E = self.y_centers[ss.__i]
                figmanip.hlines(E, ax=ax, color='red', lw=1)

            # title
            axes[0].set_title(f'{ss.__i}: {E} eV')

        def _update(ss):
            if ss.__i >= len(ss):
                ss.__i = len(ss) - 1
            elif ss.__i < 0:
                ss.__i = 0

            update(ss)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'right' or event.key == 'up':
                ss.__i = ss.__i + 1

                for ax in axes:
                    ax.cla()
                _update(ss)
                # for ax in axes:
                plt.draw()
            elif event.key == 'left' or event.key == 'down':
                ss.__i = ss.__i - 1
                
                for ax in axes:
                    ax.cla()
                _update(ss)
                # for ax in axes:
                plt.draw()

        # figure
        fig, axes = figmanip.subplots(nrows=2, ncols=1)
        fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=ss))
        _update(ss)

    # experimental 
    def roll_plot(self, axis=0, vlines=None, hlines=None, **kwargs):
        """[EXPERIMENTAL] Display data as an image. Wrapper for `matplotlib.pyplot.imshow()`_.

        Warning:
            Pixels are always square. For irregular pixel row/columns, see Image.pcolormesh()

        Args:
            axis (int or string, optional): Axis along roll.
                By default, columns are rolled down or up (0).
            vlines, hlines (list or number, optional): vertical or horizontal
                dashed lines for reference, default is None.
            **kwargs: kwargs are passed to `matplotlib.pyplot.imshow()`_.

        If not specified, the following parameters are passed to `matplotlib.pyplot.imshow()`_:

        Args:
            cmap: The Colormap instance. Default is 'jet'.
            aspect: The aspect ratio of the Axes. Default is 'auto'. If 'equal',
                an aspect ratio of 1 will be used (pixels will be square).
            origin: Location of the [0, 0] index. default is 'lower'.
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
        """
        # axis
        assert axis == 0 or axis == 1, 'axis must be 0 or 1'

        # default arguments
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'
        if 'aspect' not in kwargs:
            kwargs['aspect'] = 'auto'
        if 'origin' not in kwargs:
            kwargs['origin'] = 'lower'
        if 'interpolation' not in kwargs:
            kwargs['interpolation'] = 'none'
        if 'vmin' not in kwargs or 'vmax' not in kwargs:
            vmin, vmax = self._calculated_vmin_vmax()
            if 'vmin' not in kwargs:
                kwargs['vmin'] = vmin
            if 'vmax' not in kwargs:
                kwargs['vmax'] = vmax

        # vars
        self.__i = 0
        self.__xlim = None
        self.__ylim = None
        self.__ax   = None
        self.__temp = self.copy()
        if axis == 0:
            self._calculated_roll = np.array([0]*self.shape[1])
        else:
            raise NotImplementedError('axis = 1 nor implemented yet')

        # change keybindings
        try:
            matplotlib.rcParams['keymap.back'].remove('left')
            matplotlib.rcParams['keymap.forward'].remove('right')
        except ValueError:
            pass

        # lims
        if self.__xlim is None:
            self.__xlim = [min(self._x_centers), max(self._x_centers)]
        if self.__ylim is None:
            self.__ylim = [min(self._y_centers), max(self._y_centers)]

        # vlines hlines warning
        if vlines is not None or hlines is not None:
            raise NotImplementedError('vlines is not implemented yet')
        
        # core update function
        update_function = None
        if update_function is None:
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                plt.title(f'{self.__i}: {self._calculated_roll[self.__i]}')
                im = self.__temp.imshow(ax=self.__ax, verbose=False)

                if axis == 0:
                    self.__ax.plot([self.x_centers[self.__i], self.x_centers[self.__i]], [self.__ylim[0], self.__ylim[1]], color='white')
                else:
                    raise NotImplementedError('not implemented yet')
                
                # if vlines is not None:
                #     if isinstance(vlines, Iterable) == False:
                #         vlines = [vlines, ]
                #     for vline in vlines:
                #         self.__ax.plot([vline, vline], [self.__ylim[0], self.__ylim[1]], color='red', ls='--', lw=1)
                # if hlines is not None:
                #     if isinstance(hlines, Iterable) == False:
                #         hlines = [hlines, ]
                #     for hline in hlines:
                #         self.__ax.plot([self.__xlim[0], self.__xlim[1]], [hline, hline], color='red', ls='--', lw=1)

                if self.__xlim is not None:
                    plt.xlim(self.__xlim)
                if self.__ylim is not None:
                    plt.ylim(self.__ylim)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if axis == 0:
                if event.key == 'up':
                    self._calculated_roll[self.__i] += 1

                    if axis == 0:
                        rolls = np.zeros(self.shape[1])
                        rolls[self.__i] = 1
                        self.__temp.set_roll(rolls)
                    _update(ss)
                    plt.draw()
                elif event.key == 'down':
                    self._calculated_roll[self.__i] -= 1

                    if axis == 0:
                        rolls = np.zeros(self.shape[1])
                        rolls[self.__i] = -1
                        self.__temp.set_roll(rolls)
                    _update(ss)
                    plt.draw()
                elif event.key == 'right':
                    self.__i = self.__i + 1

                    # self.__ax.cla()
                    self.__ax.lines.clear()
                    _update(ss)
                    plt.draw()
                elif event.key == 'left':
                    self.__i = self.__i - 1

                    # self.__ax.cla()
                    self.__ax.lines.clear()
                    _update(ss)
                    plt.draw()

        # axis zoom changes
        def on_xlims_change(event_ax):
            # print("updated xlims: ", event_ax.get_xlim())
            self.__xlim = event_ax.get_xlim()

        def on_ylims_change(event_ax):
            # print("updated ylims: ", event_ax.get_ylim())
            self.__ylim = event_ax.get_ylim()

        # plotting
        fig, self.__ax = plt.subplots(1, 1)
        _update(self)

        # register callbacks
        fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
        # fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))
        self.__ax.callbacks.connect('xlim_changed', on_xlims_change)
        self.__ax.callbacks.connect('ylim_changed', on_ylims_change)

        # _update(self)
        return

# %% ============================= PhotonEvents =========================== %% #
class PhotonEvents(metaclass=_Meta):
    """Returns a ``Photon events`` object.

    Args:
        x (array): vector with x-coordinate of photon events.
        y (array): vector with y-coordinate of photon events.
        filepath (string or path object, optional): filepath.

    How to initialize a Photon events object:
        **Empty**
            
            >>> pe = br.PhotonEvents()

        **From array**
            >>> pe = br.PhotonEvents(x, y)
            >>> pe = br.PhotonEvents(x=x, y=y)
        
        where ``x`` and ``y`` are 1D array or list.

        **From file**
            >>> pe = br.PhotonEvents(<filepath>)
            >>> pe = br.PhotonEvents(filepath=<filepath>)
        
        where ``filepath`` must point to a xy-type file, where comments 
        must be marked with `#` and columns must be separated by `,` (comma).

    Attributes:
        x (array): vector with x-coordinate of photon events.
        y (array): vector with y-coordinate of photon events.
        filepath (str or pathlib.Path): filepath associated with data.
        xlim, ylim (tuple): two element tuple with min and max possible 
            x and y coordinates. This is only used for plotting. Usually, its
            okay if this is not set. 

    Computed (read-only) attributes:
        data (array)
            2 column data (x, y)
        calculated_shift (list)
            Calculated shifts.

    Write-only attributes:
        None
    """

    _read_only = ['calculated_shift', 'x_edges', 'y_edges', 'x_centers', 'y_centers']
    _non_removable = []
    
    def __init__(self, *args, **kwargs): 
        """Initialize the object instance.
         
        Raises:
            AttributeError: if kwargs and args cannot be read.
        """
        ###########################
        # Initializing attributes #
        ###########################
        # basic
        self._x        = None
        self._y        = None
        self._filepath = ''

        self._xlim = None
        self._ylim = None

        self._x_centers = None
        self._y_centers = None
        self._x_edges   = None
        self._y_edges   = None

        self._calculated_shift = None

        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input. Photon events object cannot be created.' +\
                        'Please, use one ' +\
                        'of the examples below to create an object:\n' +\
                        '\n' +\
                        'pe = br.PhotonEvents()\n' +\
                        '\n' +\
                        'pe = br.PhotonEvents(x, y)\n' +\
                        'pe = br.PhotonEvents(x=x, y=y)\n' +\
                        '\n' +\
                        'pe = br.PhotonEvents(filepath=<filepath>)\n' +\
                        'pe = br.PhotonEvents(<filepath>)\n' +\
                        '\n' +\
                        'x and y must be a 1D array (or list) and filepath ' +\
                        'must be a string or pathlib.Path object'
        if kwargs != {} and args != ():
            raise AttributeError(error_message)
        if any([item not in ['x', 'y', 'filepath'] for item in kwargs.keys()]):
            raise AttributeError(error_message)
        if len(args) > 2 or len(kwargs) > 2:
            raise AttributeError(error_message)
        if 'x' in kwargs and 'y' not in kwargs:
            raise AttributeError(error_message)
        if ('x' in kwargs or 'y' in kwargs) and 'filepath' in kwargs:
            raise AttributeError(error_message)

        ################
        # loading data #
        ################
        # keyword arguments
        if 'y' in kwargs:
            self.x = kwargs['x']
            self.y = kwargs['y']
            return
        if 'filepath' in kwargs:
            self.load(kwargs['filepath'])
            return

        # positional arguments
        if len(args) == 1:
            if isinstance(args[0], str) or isinstance(args[0], Path):
                self.load(args[0])
                return
            else:
                raise AttributeError(error_message)
        elif len(args) == 2:
            self.x = args[0]
            self.y = args[1]
            return

    ##############
    # attributes #
    ##############
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
        # check type
        if not isinstance(value, Iterable):
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

        #################
        # set attribute #
        #################
        self._ylim = np.array(value, dtype='float')
    @ylim.deleter
    def ylim(self):
        raise AttributeError('Cannot delete object.')

    @property
    def filepath(self):
        return self._filepath
    @filepath.setter
    def filepath(self, value):
        if value is None:
            value = ''
        elif isinstance(value, str) or isinstance(value, Path):
            self._filepath = value
        else:
            raise TypeError(r'Invalid type ' + str(type(value)) + 'for filepath\nFilepath can only be str or pathlib.Path type.')
    @filepath.deleter
    def filepath(self):
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

    #################
    # magic methods #
    #################
    def __len__(self):
        if self.x is None:
            return 0
        else:
            return len(self.x)

    def __setattr__(self, name, value):
        if name in settings._forbidden_words['PhotonEvents']:
            raise AttributeError(f'`{name}` is a reserved word and cannot be set as an attribute')
        super().__setattr__(name, value)

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
    def get_attrs(self):
        """return attrs that are user defined.""" 
        return [key for key in self.__dict__.keys() if key.startswith('_') == False and key not in settings._reserved_words]
    
    def _get_center_attrs(self):
        """return center attrs"""
        return ['_x_centers', '_y_centers', '_x_edges', '_y_edges']

    def _get_lim_attrs(self):
        """return lim attrs"""
        return ['_xlim', '_ylim']

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

    def copy_lims_from(self, s):
        """Copy lim attributes from another brixs.PhotonEvents object.

        Args:
            s (brixs object): PhotonEvents to copy lims attributes from.
        
        Returns:
            None
        """
        # check type
        if isinstance(s, PhotonEvents):
            pass
        else:
            raise TypeError(f'type {type(s)} not valid\nCan only copy user attrs from type br.PhotonEvents to br.PhotonEvents')

        # transfer attrs
        for attr in s._get_lim_attrs():
            value = copy.deepcopy(s.__dict__[attr])
            self.__setattr__(attr, value)

    def copy_center_from(self, s):
        """Copy center attributes from another brixs.PhotonEvents object.

        Args:
            s (brixs object): PhotonEvents to copy lims attributes from.
        
        Returns:
            None
        """
        # check type
        if isinstance(s, PhotonEvents):
            pass
        else:
            raise TypeError(f'type {type(s)} not valid\nCan only copy user attrs from type br.PhotonEvents to br.PhotonEvents')

        # transfer attrs
        for attr in s._get_center_attrs():
            value = copy.deepcopy(s.__dict__[attr])
            self.__setattr__(attr, value)

    ###########
    # support #
    ###########
    pass

    #################
    # save and load #
    #################
    def _create_header(self, verbose=False):
        """Gather attrs to be saved to a file."""
        header = ''
        attrs  = self.get_attrs()
        attrs += self._get_lim_attrs()
        attrs += self._get_center_attrs()
        for name in attrs:
            if self.__dict__[name] is None:
                header += f'{name}: None'  + '\n'
            elif isinstance(self.__dict__[name], str):
                temp2 = str(self.__dict__[name]).replace('\n','\\n')
                header += f'{name}: \"{temp2}\"'  + '\n'
            elif isinstance(self.__dict__[name], dict) or isinstance(self.__dict__[name], MutableMapping):
                if verbose:
                    type_ = str(type(self.__dict__[name]))
                    print('\nWarning: Cannot save attr of type: ' + type_ + '\nattr name: '+ name + '\nTo turn off this warning, set verbose to False.')
            elif isinstance(self.__dict__[name], Iterable):
                header += f'{name}: {list(self.__dict__[name])}'  + '\n'
            elif numanip.is_number(self.__dict__[name]):
                tosave = str(self.__dict__[name])
                if tosave[-1] == '\n':
                    tosave = tosave[:-1]
                header += f'{name}: {tosave}'  + '\n'
            else:
                temp2 = str(self.__dict__[name]).replace('\n','\\n')
                header += f'{name}: \"{temp2}\"'  + '\n'
        return header[:-1]
    
    def save(self, filepath=None, only_data=False,  check_overwrite=False, verbose=False, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Attrs are saved as comments if only_data is False. Saving attrs to file
        is tricky because requires converting variables to string. Only attrs 
        that are of type: string, number, and list of number and strings are 
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
        # filepath
        if filepath is None:
            try: 
                filepath = Path(self.filepath)
            except AttributeError:
                raise TypeError("Missing 1 required argument: 'filepath'")
        else:
            filepath = Path(filepath)

        # check if filepath points to a file
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        # check overwrite
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if query.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')
                
        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([numanip.n_decimal_places(x) for x in arraymanip.flatten(self.data)])
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
        else:
            if 'header' not in kwargs:
                kwargs['header'] = self._create_header(verbose=verbose)
            else:
                if kwargs['header'] == '':
                    kwargs['header'] = self._create_header(verbose=verbose)
                elif kwargs['header'][-1] != '\n':
                    kwargs['header'] += '\n'
        np.savetxt(Path(filepath), self.data, **kwargs)

        # save filepath
        self.filepath = filepath

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
        ###########
        # default #
        ###########
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        # check if filepath points to a file
        filepath = Path(filepath)
        assert filepath.is_file(), f'filepath must point to a file\n{filepath}'
        assert filepath.exists(), f'filepath does not exist\n{filepath}'

        #############
        # read data #
        #############
        data = np.genfromtxt(Path(filepath), **kwargs)
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
        
        # save filepath
        self.filepath = str(filepath)

    ##################
    # base modifiers #
    ##################
    def set_shift(self, value=None, p=None, f=None, axis=0):
        """Shift photon events along a given axis.

        If value, p, and f are None, values are read from pe.calculated_shift.

        Args:
            value (int or tuple): The number by which the data are shifted.
                If a tuple, then it must be of the same size as the number of
                columns (``for axis=0``) or rows (``for axis=1``). Default is None.
            p (np.array, read only): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term. Polynomial as a function of x for 
                axis=0 or y for axis=1. Default is None.
            f (function, read only): funcion f(x) for axis=0 or 
                f(y) for axis=1. Default is None.
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.

        Returns:
            None
        """
        # check axis
        if axis == 0:
            arr = self.x
        elif axis == 1:
            arr = self.y
        else:
            raise ValueError('axis must be 0 or 1')

        # calculate value
        if value is not None:
            if isinstance(value, Iterable):
                assert len(value) == len(self), f'value must have the same length as the data\nLength of the data = {len(self)}.'
            else:
                value = [value]*len(self)
        elif p is not None:
            value = np.polyval(p, arr)
        elif f is not None:
            value = f(arr)
        else:
            value = self.calculated_shift

        # shift
        if axis == 0:
            for i, v in enumerate(value):
                if v != 0:
                    self.y[i] += v
        elif axis == 1:
            for i, v in enumerate(value):
                if v != 0:
                    self.x[i] += v

    #############
    # modifiers #
    #############
    def fix_curvature(self, nbins, deg=2, axis=0, limit_size=1000):
        """Shift datapoints along a given axis to fix curvature via cc.

        Elements that roll beyond the edge are re-introduced at the first position.

        Args:
            nbins (int or tuple): number of bins. Must be non-zero positive number.
                If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the number of rows and number of
                columns, respectively.
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of bin columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.

        Returns:
            p (np.array, read only): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term. Polynomial as a function of x_centers for 
                axis=0 or y_centers for axis=1. Default is None.
            f (function, read only): funcion f(x_centers) for axis=0 or 
                f(y_centers) for axis=1.
        """       
        popt, model = self.calculate_shift(nbins=nbins, deg=deg, axis=axis, limit_size=limit_size)
        self.set_shift(value=self.calculated_shift)

        return popt, model

    ##############
    # extractors #
    ##############
    def _copy(self, *args, **kwargs):
        """Same as copy(), but attributes are not copied to the new object."""
        
        # check if extract is really necessary
        if kwargs == {} and args == ():
            return copy.deepcopy(self)
        else:
            limits = self._check_limits(*args, **kwargs)
        if (limits[0][0] <= min([min(s.x) for s in self]) and limits[0][1] >= max([max(s.x) for s in self]) ) and len(limits) == 1:
                return copy.deepcopy(self)

        # if x is the same for all spectra, this operation is much faster
        if self.x is None:
            try:
                self.check_same_x()
                if self.x is not None:
                    ss = Spectra(n=len(self))
                    x, ys = self._gather_ys(*args, **kwargs)
                    for i in range(len(self)):
                        ss[i].copy(self[i])
                        ss[i]._x = x
                        ss[i]._x = ys[i]
                    return ss
            except:
                pass
        
        # if x is not the same, extract data recursively
        limits = self._check_limits(*args, **kwargs)
        ss = Spectra(n=len(self))
        for i, s in enumerate(self.data):
            try:
                ss[i] = s.copy(limits=limits)
            except RuntimeError:
                raise RuntimeError(f'It seems like spectrum number {i} has no data points within range: {limits}.\nPlease, fix limits (or delete spectrum) so all spectra have at least one data point within range.')
        return ss
    
    def copy(self, *args, **kwargs):
        """Return a copy of the object.

        Usage:
            >>> # full copy
            >>> pe2 = pe1.copy()  # pe2 is now a copy of ss1
            >>>
            >>> # partial copy
            >>> pe3 = pe1.copy(x_start, x_stop, y_start, y_stop)
            >>> pe3 = pe1.copy((x1_start, x1_stop, y1_start, y1_stop))
            >>> pe3 = pe1.copy([(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])

        Args:
            pe (Image, optional): PhotonEvents to be is copied. See usage.
            mask (list): list with rectangles coordinates (x_start, x_stop, y_start, y_stop).

        Returns:
            :py:attr:`PhotonEvents`
        """
        error_message = 'Wrong input. Please check input.\n' +\
                        f'args passed: {args}\n'  +\
                        f'kwargs passed: {kwargs}\n\n'  +\
                        ' === Usage ===\n'  +\
                        'pe2 = pe1.copy()  # pe2 is now a copy of pe1\n'  +\
                        'pe3 = pe1.copy(x_start, x_stop, y_start, y_stop)\n' +\
                        'pe3 = pe1.copy([(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])'
        ############################
        # check if mask are passed #
        ############################
        mask = None
        if 'mask' in kwargs:
            mask = kwargs['mask']
        elif len(args) == 4:
            mask = [[args[0], args[1], args[2], args[3]], ]
        elif len(args) != 4 and len(args) != 0:
            raise ValueError(error_message)
        elif len(kwargs) != 0:
            for x in list(kwargs.keys()):
                if x not in ['mask']:
                    raise ValueError(f'{x} is not a recognized input for copy function\n'+error_message)

        #########################################
        # if mask is passed, clip photon events #
        #########################################
        if mask is not None:
            pe = self.copy()
            pe.clip(mask=mask)

        ####################################
        # Otherwise, return a copy of self #
        ####################################
        else:
            pe = PhotonEvents(x=self.x, y=self.y)
            pe.copy_lims_from(self)

        ##################
        # transfer attrs #
        ##################
        # pe.copy_lims_from(self)
        pe.copy_attrs_from(self)

        return pe

    def binning(self, *args, **kwargs):
        """Compute the 2D histogram of the data (binning of the data).

        Usage:
            >>> # 10 rows, 10 columns
            >>> pe.binning(nbins=10)       
            >>> pe.binning(nbins=(10))   
            >>> pe.binning(10)            
            >>> pe.binning((10))          
            >>> 
            >>> # 10 rows, 5 columns
            >>> pe.binning(nbins=(10, 5)) 
            >>> pe.binning((10, 5))      
            >>> pe.binning(10, 5)      

        Args:
            nbins (int or tuple): number of bins. Must be non-zero positive number.
                If one value is given,
                this is used for both x and y directions. If two values are given,
                they are used separetely for the number of rows and number of
                columns, respectively.

        Returns:
            Image
        """
        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input. Please, see examples below:\n' +\
                        '\n' +\
                        'im.binning(nbins=10)\n' +\
                        'im.binning(nbins=(10))\n' +\
                        'im.binning(10)\n' +\
                        'im.binning((10))\n' +\
                        '\n' +\
                        'im.binning(nbins=(10, 5))\n' +\
                        'im.binning((10, 5))\n' +\
                        'im.binning(10, 5)\n' +\
                        '\n' +\
                        'numbers have to be positive interger'
        if kwargs != {} and args != ():
            raise AttributeError(error_message)
        if any([item not in ['nbins', ] for item in kwargs.keys()]):
            raise AttributeError(error_message)
        if len(args) > 2 or len(kwargs) > 1:
            raise AttributeError(error_message)
        
        #################
        # sorting input #
        #################
        # keyword arguments
        if 'nbins' in kwargs:
            args = kwargs['nbins']
        # positional arguments
        if isinstance(args, Iterable):
            if len(args) == 1:
                if isinstance(args[0], Iterable):
                    nbins = list(args[0])
                else:
                    nbins = [args[0], args[0]]
            else:
                nbins = list(args)
        else:
            nbins = [args, args]

        # assert format
        if nbins[0] is None or nbins[1] is None:
            raise AttributeError(error_message)
        if numanip.is_integer(nbins[0]) == False or numanip.is_integer(nbins[1])  == False or nbins[0] < 0 or nbins[1] < 0:
            raise ValueError("Number of bins must be a positive integer.")

        ###############
        # Calculation #
        ###############
        if self.xlim is None:
            xlim = (min(self.x), max(self.x))
        else:
            xlim = self.xlim
        if self.ylim is None:
            ylim = (min(self.y), max(self.y))
        else:
            ylim = self.ylim
        temp, _x_edges, _y_edges = np.histogram2d(self.x, self.y, bins=nbins[::-1], range=(xlim, ylim))

        # Image
        im = Image(temp.transpose())
        _x_centers = arraymanip.moving_average(_x_edges, n=2)
        _y_centers = arraymanip.moving_average(_y_edges, n=2)
        
        for obj in (self, im):
            obj._x_centers = _x_centers
            obj._y_centers = _y_centers
            obj._x_edges   = _x_edges
            obj._y_edges   = _y_edges

        ##################
        # transfer attrs #
        ##################
        im.copy_attrs_from(self)

        return im

    def transpose(self):
        """Switch x and y positions.

        Returns:
            PhotonEvents
        """
        # create new PhotonEvents with switched x and y
        pe = PhotonEvents(x=self.y, y=self.x)
  
        ##################
        # transfer attrs #
        ##################
        # pe.copy_lims_from(self)
        pe.copy_attrs_from(self)

        return pe

    def crop(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Crop photon events out.

        Args:
            x_start, x_stop, y_start, y_stop (int): pixel range in terms of
                x_centers and y_centers. Interval is inclusive. Use None to 
                indicate the edge of the image.

        Returns:
            PhotonEvents
        """
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

        ########
        # crop #
        ########
        temp = np.array([(x, y) for x, y in zip(self.x, self.y) if ((x > x_start and x < x_stop) and (y > y_start and y < y_stop))])
        pe = PhotonEvents(x=list(temp[:, 0]), y=list(temp[:, 1]))

        ##################
        # transfer attrs #
        ##################
        pe.copy_lims_from(self)
        pe.copy_attrs_from(self)

        return pe

    def clip(self, mask):
        """Clip photon events.

        Usage:
            >>> pe.clip((xmin, xmax, ymin, ymax))
            >>> pe.clip((0, 12, 3, 12))
            >>> pe.clip([(0, 12, 3, 12), (1, 4, 15, 18)])

        Args:
            mask (list): list with rectangles coordinates (x_start, x_stop, y_start, y_stop).

        Returns:
            PhotonEvents
        """
        # assert mask is the right format
        assert isinstance(mask, Iterable), 'mask must be iterable'
        if len(mask) == 4:
            if isinstance(mask[0], Iterable) == False:
                mask = [mask, ]

        ########
        # clip #
        ########
        x = []
        y = []
        for r in mask:
            temp = np.array([(x, y) for x, y in zip(self.x, self.y) if ((x > r[0] and x < r[1]) and (y > r[2] and y < r[3]))])
            x += list(temp[:, 0])
            y += list(temp[:, 1])

        #########
        # final #
        #########
        pe = PhotonEvents(x=x, y=y)

        ##################
        # transfer attrs #
        ##################
        pe.copy_lims_from(self)
        pe.copy_attrs_from(self)

        return pe

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
        
        pe = self.copy()
        im = pe.binning(nbins=nbins)
        s  = im.calculate_spectrum(axis=axis)

        ##################
        # transfer attrs #
        ##################
        s.copy_attrs_from(self)

        return s

    ########################
    # calculation and info #
    ########################
    def calculate_shift(self, nbins, deg=2, axis=0, limit_size=1000):
        """Calculate intensity misalignments via cross-correlation.

        Args:
            nbins (int, optional): number of bins.
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of bin columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.

        Returns:
            p (np.array, read only): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term. Polynomial as a function of x_centers for 
                axis=0 or y_centers for axis=1. Default is None.
            f (function, read only): funcion f(x_centers) for axis=0 or 
                f(y_centers) for axis=1.
        """
        # check axis
        if axis != 0 and axis != 1:
            raise ValueError('axis must be 0 or 1')
        
        # calculate shifts
        pe = self.copy()
        im = pe.binning(nbins=nbins)
        im.calculate_roll(axis=axis, limit_size=limit_size)
        if axis == 0:
            centers = im.x_centers
            y = im.calculated_roll*np.mean(np.diff(im.y_centers))
        elif axis == 1:
            centers = im.y_centers
            y = im.calculated_roll*np.mean(np.diff(im.x_centers))

        # calculate poly
        popt  = np.polyfit(x=centers, y=y, deg=deg)
        model = lambda x: np.polyval(popt, x)

        # store
        if axis == 0:
            centers = self.x
        elif axis == 1:
            centers = self.y
        self._calculated_shift = model(centers)

        return popt, model

    def calculate_curvature(self, nbins, deg=2, axis=0, limit_size=1000):
        """Same as br.pe.calculate_shift()
        
        Calculate intensity misalignments via cross-correlation.

        Args:
            nbins (int, optional): number of bins.
            deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
            axis (int or string, optional): Axis along which elements are shifted.
                By default, data is shifted in the vertical (0) direction.
            limit_size (int or False, optional): prevents from mistakenly calculating
                cross-corelation for unusualy big images. If axis = 0 (1), it 
                ensures that the number of bin columns (rows) is not bigger than 
                limit_size. Default is 1000. Set to False to bypass this limit.

        Returns:
            p (np.array, read only): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term. Polynomial as a function of x_centers for 
                axis=0 or y_centers for axis=1. Default is None.
            f (function, read only): funcion f(x_centers) for axis=0 or 
                f(y_centers) for axis=1.
        """
        return self.calculate_shift(nbins=nbins, deg=deg, axis=axis, limit_size=limit_size)
     
    ##########################        
    # plot and visualization #
    ########################## 
    def plot(self, ax=None, set_limits=True, **kwargs):
        """Display data as an image. Wrapper for `matplotlib.pyplot.scatter()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            set_limits (bool, optional): If True, plt.xlim and plt.ylim are
                called and the limits are set to pe.xlim and pe.ylim.
            **kwargs: kwargs are passed to ``plt.scatter()`` that plots the data.

        If not specified, the following parameters are passed to `matplotlib.pyplot.scatter()`_:

        Args:
            s: The marker size in points**2. Default is 0.1.

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

        # kwargs
        if 's' not in kwargs:
            kwargs['s'] = 0.1

        # plot
        pos = ax.scatter(self.x, self.y, **kwargs)

        # set limits
        if set_limits:
            try:
                if self.xlim is not None:
                    ax.xlim(self.xlim)
                if self.ylim is not None:
                    ax.ylim(self.ylim)
            except AttributeError:
                if self.xlim is not None:
                    ax.set_xlim(self.xlim)
                if self.ylim is not None:
                    ax.set_ylim(self.ylim)

        return pos

    def plot_edges(self, ax=None, **kwargs):
        """[Not implemented]Plot edges after binning

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            **kwargs: kwargs are passed to ``plt.scatter()`` that plots the data.

        Returns:
            None
        """
        raise NotImplementedError('this is not implemented yet')

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
