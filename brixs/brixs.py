#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core brixs module"""

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# special libraries
from collections.abc import Iterable, MutableMapping
from mpl_toolkits.axes_grid1 import make_axes_locatable

# backpack
import brixs.backpack.filemanip  as filemanip
import brixs.backpack.arraymanip as arraymanip
import brixs.backpack.figmanip   as figmanip
import brixs.backpack.numanip    as numanip
import brixs.backpack.interact   as interact
from .peaks import Peaks

# common definitions ===========================================================
from .config import settings

# support class ================================================================
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
            return getter, setter, deleter#, 'read only attribute'

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

# BRIXS ========================================================================
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
        step (number or None, read only): step size of the x-array.
        monotonicity (str or None, read only): monotonicity of the x-array.
        peaks (:py:attr:`Peaks`): each entry represents a parameter of a peak.

    Computed (read-only) attributes:
        data (array)
            2 column data (x, y)

        area (number) 
            area under the curve.

    Write-only attributes:
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
    _read_only     = ['step', 'monotonicity', 'peaks']
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
        self._x    = None
        self._y    = None
        self._filepath = ''

        # check
        self._step         = None
        self._monotonicity = None

        # special
        self._peaks = Peaks()

        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input. Spectrum object cannot be created. Please, use one ' +\
                        'of the examples below to create a spectrum object:\n' +\
                        '\n' +\
                        's = br.Spectrum()\n' +\
                        '\n' +\
                        's = br.Spectrum(x, y)\n' +\
                        's = br.Spectrum(x=x, y=y)\n' +\
                        '\n' +\
                        's = br.Spectrum(y)\n' +\
                        's = br.Spectrum(y=y)\n' +\
                        '\n' +\
                        's = br.Spectrum(filepath=<filepath>)\n' +\
                        's = br.Spectrum(<filepath>)\n' +\
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
        # if x is available, x must be set first.
        # Otherwise, when y is set, it will create a fake x axis using np.arange()

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

        ##########################
        # reset check attributes #
        ##########################
        self._step         = None
        self._monotonicity = None

        ############################
        # reset special attributes #
        ############################
        self.peaks.clear()
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
        self.peaks.clear()
    @y.deleter
    def y(self):
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

    # @property
    # def area(self):
    #     return self.calculate_area()
    # @area.setter
    # def area(self, value):
    #     raise AttributeError('Attribute is "read only". Cannot set attribute.')
    # @area.deleter
    # def area(self):
    #     raise AttributeError('Cannot delete object.')

    #########################
    # write-only attributes #
    #########################
    @property
    def calib(self):
        pass
    @calib.setter
    def calib(self, value):
        self.set_calib(value)
    @calib.deleter
    def calib(self):
        raise AttributeError('Cannot delete object.')

    @property
    def shift(self):
        pass
    @shift.setter
    def shift(self, value):
        self.set_shift(value)
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def roll(self):
        pass
    @roll.setter
    def roll(self, value):
        self.set_roll(value)
    @roll.deleter
    def roll(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        pass
    @offset.setter
    def offset(self, value):
        self.set_offset(value)
    @offset.deleter
    def offset(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        pass
    @factor.setter
    def factor(self, value):
        self.set_factor(value)
    @factor.deleter
    def factor(self):
        raise AttributeError('Cannot delete object.')

    #################
    # magic methods #
    #################
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
            final = Spectrum(x=self.x, y=self.y + object.y)
            
        elif isinstance(object, (np.floating, float, int)):
            final = Spectrum(x=self.x, y=self.y + object)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')
        
        ##################
        # transfer attrs #
        ##################
        final._step         = self.step
        final._monotonicity = self.monotonicity

        # for attr in self.get_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     final.__setattr__(attr, value)

        return final

    def __sub__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot operate on spectra. x axis is different.\nMaybe try interpolating the x axis.')
            final = Spectrum(x=self.x, y=self.y - object.y)
        elif isinstance(object, (np.floating, float, int)):
            final = Spectrum(x=self.x, y=self.y - object)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')
        
        ##################
        # transfer attrs #
        ##################
        final._step = self.step
        final._monotonicity = self.monotonicity
        
        # for attr in self.get_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     final.__setattr__(attr, value)

        return final

    def __mul__(self, object):
        if isinstance(object, Spectrum):
            ss = Spectra([self, object])
            try:
                ss.check_same_x()
            except ValueError:
                raise ValueError('Cannot operate on spectra. x axis is different.\nMaybe try interpolating the x axis.')
            final = Spectrum(x=self.x, y=self.y * object.y)
        elif isinstance(object, (np.floating, float, int)):
            final = Spectrum(x=self.x, y=self.y * object)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')
        
        ##################
        # transfer attrs #
        ##################
        final._step = self.step
        final._monotonicity = self.monotonicity
        
        # for attr in self.get_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     final.__setattr__(attr, value)

        return final

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
                final = Spectrum(x=self.x, y=self.y/object.y)
        elif isinstance(object, (np.floating, float, int)):
            if object == 0:
                raise ZeroDivisionError(f'Cannot divide by zero.')
            else:
                final = Spectrum(x=self.x, y=self.y / object)
        else:
            raise ValueError(f'Cannot operate type {type(object)} with type Spectrum')
        
        ##################
        # transfer attrs #
        ##################
        final._step = self.step
        final._monotonicity = self.monotonicity
        
        # for attr in self.get_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     final.__setattr__(attr, value)

        return final

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
        
    ###########
    # support #
    ###########
    def _default_attrs(self):
        """return list with default attrs."""
        return ['_x', '_y', '_filepath', '_step', '_monotonicity', '_peaks']
    
    def get_attrs(self):
        """return attrs that are user defined.""" 
        return [key for key in self.__dict__.keys() if key not in self._default_attrs() and key.startswith('_') == False]
    
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

    def _validate_ranges(self, *args, **kwargs):
        """Check if ranges is the right format.

        Possible input formats are:

        xi, xf
        (xi, xf)
        (xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3)
        ((xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3))

        Input can be passed as positional arguments or with the kwarg `ranges`.

        If xi=None, this item is replaced by the min value of the dataset.
        If xf=None, this item is replaced by the max value of the dataset.

        Returns:
            ranges in the following format:
                ((xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3))
        """
        error_message = 'Wrong input. Please check ranges input.\n' +\
                        f'args passed: {args}\n' +\
                        f'kwargs passed: {kwargs}\n\n' +\
                        'ranges should be a pair (x_init, x_final) or a list of pairs.\n' +\
                        'See examples below:\n' +\
                        '\n' +\
                        'ranges = (x_init, x_final)\n' +\
                        '\n' +\
                        'ranges = ((x_init_1, x_final_1), (x_init_2, x_final_2), ...)\n' +\
                        '\n' +\
                        'Use None to indicate the minimum or maximum x value of the data.\n'
        if kwargs != {} and args != ():
            raise AttributeError(error_message)
        if any([item not in ['ranges', ] for item in kwargs.keys()]):
            raise AttributeError(error_message)
        if 'ranges' in kwargs:
            args = kwargs['ranges']

        # get min and max range values
        vmin = min(self.x)
        vmax = max(self.x)

        ##############
        # fix format #
        ##############
        # empty input
        if len(args) == 0:
            return ((vmin, vmax),)
        
        # (xi, xf)
        # ((xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3))
        elif len(args) == 1:
            if isinstance(args[0][0], Iterable):
                if len(args[0][0]) == 2:
                    args = args[0]
                else:
                    raise AttributeError(error_message)

        # xi, xf
        # (xi_, xf_1), (xi_2, xf_2)
        elif len(args) == 2:
            if isinstance(args[0], Iterable):
                args = list(args)
            else:
                args = [args, ]
        
        # (xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3)
        args = list(args)
        for i, pair in enumerate(args):
            pair    = list(pair)
            args[i] = pair
            if len(pair) == 2:
                if pair[0] == None: pair[0] = vmin
                if pair[1] == None: pair[1] = vmax
            else:
                raise AttributeError(error_message)        
        return args

    #################
    # save and load #
    #################
    def _create_header(self, verbose=False):
        """Gather attrs to be saved to a file."""
        header = ''
        attrs = self.get_attrs()
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

    def save(self, filepath=None, only_data=False, check_overwrite=False, verbose=False, **kwargs):
        r"""Save data to a text file. Wrapper for `numpy.savetxt()`_.

        Attrs are saved as comments if only_data is False. Saving attrs to file
        is not always reliable because requires converting variables to string. 
        Only attrs that are of type: string, number, and list of number and strings are 
        saved somewhat correctly. Dictionaries are not saved. 


        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format. If None, filepath can be inferred from a
                user defined attr 's.filepath'. Default is None. Last used 
                dirpath is saved to an attr ss.filepath.
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
                    if interact.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')

        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([figmanip.n_decimal_places(x) for x in arraymanip.flatten(self.data)])
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
        self.filepath = str(filepath)

    def load(self, filepath, only_data=False, verbose=False, **kwargs):
        """Load data from a text file. Wrapper for `numpy.genfromtxt()`_.

        This a very simple loading function that works well with two column text
        files. If file has more columns than two columns, the first two columns
        will be loaded. Use `usecols` to select columns, for example:
        
        usecols = (1, 4) 

        If file was saved by br.Spectrum.save(), then the metadata (comments) can be 
        recovered. If not, only_data must be set to True.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename extension is .gz or .bz2, the file is first 
                decompressed. Last used filepath is saved to an attr s.filepath.
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
                of a comment. Default is ``#``. Attributes picked up
                from the header will be loaded too.
            usecols (tuple, optional): Default is (0, 1).

        Returns:
            None

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '#'
        if 'usecols' not in kwargs:
            kwargs['usecols'] = (0, 1)

        # check if filepath points to a file
        filepath = Path(filepath)
        assert filepath.is_file(), f'filepath must point to a file\n{filepath}'
        assert filepath.exists(), f'filepath does not exist\n{filepath}'

        # read data
        data = np.genfromtxt(Path(filepath), **kwargs)
        
        # check data
        x = data[:, 0]
        y = data[:, 1]
        assert len(x) == len(y), f'Length of x array (len={len(x)}) is not compatible with y array (len={len(y)}).'
        
        # assign
        self._x = x
        self._y = y

        ##########################
        # reset check attributes #
        ##########################
        self._step         = None
        self._monotonicity = None

        ############################
        # reset special attributes #
        ############################
        self.peaks.clear()

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
        
        #################
        # save filepath #
        #################
        self.filepath = str(filepath)

    #################
    # check methods #
    #################
    def check_step(self, max_error=None):
        """Checks vector uniformity of the x-coordinates.

            If the step between two data points is the same through out the
            x vector, it sets :py:attr:`step` with the value of the step size.

            Args:
                max_error (number, optional): percentage (of the x step) value 
                of the max error.
                
                (max(steps) - min(steps))/np.mean(steps) * 100 < max_error

            Returns:
                None

            Raises:
                ValueError: If x-coordinates are not uniform.
        """
        if max_error is None:
            max_error = settings.MAX_ERROR_STEP_X

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

        # special
        self.peaks._step = self._step

    def check_monotonicity(self):
        """Sets monotonicity attribute to 'increasing' or 'decreasing'.

        Raises:
            ValueError if data is not monotonic.

        Returns:
            None
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

    ##################
    # base modifiers #
    ##################
    def set_calib(self, value):
        """Set calibration value.

        Args:
            value (number, list, function): calibration value. 
                If number, x-coordinates will be multiplied by this value.
                If list, polynomial coeff. (highest power first) is expected and
                a function new_x = f(x) will be set using np.polyval().
                If function, a function of type new_x = f(x) is expected.

        Returns:
            None
        """
        #######
        # run #
        #######
        if isinstance(value, Iterable):
            f = lambda x: np.polyval(value, x)
            self._x = np.array([f(x) for x in self.x])
        elif callable(value):
            self._x = np.array([f(x) for x in self.x])
        else:
            if value == 0:  # calib cannot be zero
                raise ValueError('cannot set calib = 0.0')
            elif value == 1:   # if calib is 1, do nothing
                return
            self._x = self.x*value

        ##########################
        # reset check attributes #
        ##########################
        self._step         = None
        self._monotonicity = None

        ###############
        # fix special #
        ###############
        self.peaks.set_calib(value=value)

    def set_shift(self, value):
        """Set shift value.

        Adds a value to the x-coordinates. y-coordinates are fully preserved.

        Args:
            value (float or int): shift value (value will be added to x-coordinates).

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if value == 0:
            return
        
        self._x += value

        ###############
        # fix special #
        ###############
        self.peaks.set_shift(value=value)

    def set_roll(self, value):
        """Set roll value.

        Roll array elements for the x-coordinates. Elements that roll beyond the
        last position are re-introduced at the first. y-coordinates are fully 
        preserved.

        Value must be lower than the length of the array.

        Args:
            value (float or int): roll value.

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if value == 0:
            return
        assert value < len(self), 'roll value must be smaller than the array length\nUse s.set_shift() for shifting the x array' 
        
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
        
        #######
        # run #
        #######
        self._y = np.roll(self.y, int(value))

        ###############
        # fix special #
        ###############
        self.peaks.set_shift(value=value*self.step)

    def set_offset(self, value):
        """Set offset value.

        Args:
            value (value): offset value (value will be added to y-coordinates).

        Returns:
            None
        """
        ###################################
        # asserting validity of the input #
        ###################################
        if value == 0:
            return
        
        self._y = self.y + value

        ###############
        # fix special #
        ###############
        self.peaks.set_offset(value=value)

    def set_factor(self, value):
        """Set y multiplicative factor.

        Args:
            value (number): multiplicative factor (y-coordinates will be
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
    
        self._y = self.y*value

        ###############
        # fix special #
        ###############
        self.peaks.set_factor(value=value)

    #############
    # modifiers #
    #############
    def floor(self, *args, **kwargs):
        """Sets zero value for y-coordinates.

        Usage:
            >>> # to bring the avg of all data points to zero.
            >>> s.floor(None, None)  
            >>>
            >>> # Brings the avg between x=0 and 10 and between x=90 and 100 to zero
            >>> s.floor((0, 10), (90, 100))  

        Args:
            ranges (list, optional): Pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                Use None to indicate the minimum or maximum x value of the data.
                The mean value of the y-coordinates inside the defined range
                will be set to zero.

        Returns:
            None
        """
        temp  = self._extract(*args, **kwargs)
        value = -np.mean(temp.y)
        self.set_offset(value)

    def flip(self):
        """Flip sign of x axis.

        Returns:
            None
        """
        self.set_calib(value=-1)

    def normalize(self, value=1, *args, **kwargs):
        """Set a factor such as the average y between ranges is equal to value.

        Args:
            value (number): value. Default is 1.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data. If None, the maximum
                y value of the data is set to 1.

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_factor`
        """
        s = self._extract(*args, **kwargs)
        self.set_factor(value/np.mean(s.y))

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
            None
        """
        self.check_monotonicity()
        if self.monotonicity != 'increasing':
            raise ValueError('monotonicity of data must be strictly increasing. Use Spectrum.fix_monotonicity().')

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
                temp, _ = self._extract(self.x, self.y, (start, stop))
                x = np.linspace(start, stop, num=len(temp))

        self._y = np.interp(x, self.x, self.y)
        self._x = x

        # check
        self._step         = None
        self._monotonicity = None
        # self._reset_modifiers()
        #special
        # self.peaks.clear()

    # def remove(self, *args, **kwargs):
    #     """Remove datapoints within a range.

    #     Args:
    #         ranges (list, optional): a pair of values or a list of pairs. Each pair represents
    #             the start and stop of a data range from x. Use None to indicate
    #             the minimum or maximum x value of the data.

    #     Returns:
    #         None
    #     """
    #     ranges = self._validate_ranges(*args, **kwargs)
    #     x, y   = arraymanip.extract(self.x, self.y, ranges, invert=True)
    #     self._x = x
    #     self._y = y

    def crop(self, start=None, stop=None):
        """Crop edges of the dataset.

        Args:
            start (number, optional): start x value. If None, the minimum value of
                x will be used.
            stop (number, optional): final x value. If None, the maximum value of
                x will be used.

        Returns:
            None
        """
        if start is None and stop is None:
            return
        
        if start is None:
            start = min(self.x)
        if stop is None:
            stop = max(self.x)
        
        if start <= min(self.x) and stop >= max(self.x):
            return
        
        s = self._extract(start, stop)
        self._x = s.x
        self._y = s.y

    def switch(self):
        """Switch x and y axis.
        
        Returns:
            None"""
        x = copy.deepcopy(self.x)
        y = copy.deepcopy(self.y)

        self.x = y
        self.y = x

    ##############
    # extractors #
    ##############
    def _extract(self, *args, **kwargs):
        """Same as self.copy(), but attrs are NOT copied."""

        # check if extract is really necessary
        if kwargs == {} and args == ():
            return copy.deepcopy(self)
        else:
            ranges = self._validate_ranges(*args, **kwargs)
            if (ranges[0][0] <= min(self.x) and ranges[0][1] >= max(self.x)) and len(ranges) == 1:
                return copy.deepcopy(self)

        # extract
        x, y   = arraymanip.extract(self.x, self.y, ranges)
        s      = Spectrum(x=x, y=y)
        return s
    
    def copy(self, *args, **kwargs):
        """Return a copy of the data within a range.

        Usage:
            >>> s2.copy(s1)     # s2 is now a copy of s1
            >>> s2 = s1.copy()  # s2 is now a copy of s2
            >>>
            >>> # s2 will have only data between x=0 and 10 and between x=90 and 100
            >>> s2 = s1.copy((0, 10), (90, 100))  

        Args:
            s (Spectrum, optional): Spectrum is copied. See usage.
            ranges (list, optional): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            :py:attr:`Spectrum`
        """
        ##############################
        # check if input is Spectrum #
        ##############################
        s = None
        if 's' in kwargs:
            s = kwargs['s']
        elif len(args) == 1:
            if isinstance(args[0], Spectrum):
                s = args[0]

        ########################################
        # if Spectrum is passed, copy Spectrum #
        ########################################
        if s is not None:
            if isinstance(s, Spectrum):
                self._x    = copy.deepcopy(s.x)
                self._y    = copy.deepcopy(s.y)
                self._filepath = s.filepath
                self._step         = s.step
                self._monotonicity = s.monotonicity

                # special
                self._peaks = s.peaks

                # user defined attrs
                for attr in self.get_attrs():
                    self.__delattr__(attr)
                for attr in s.get_attrs():
                    self.__setattr__(attr, s.__dict__[attr])
            else:
                raise TypeError('Only type br.Spectrum can be copied to type br.Spectrum')
            return
        #######################################
        # Otherwise, a type ranges is assumed #
        #######################################
        else:
            s = self._extract(*args, **kwargs)

            # transfer attrs
            for attr in self.get_attrs():
                value = copy.deepcopy(self.__dict__[attr])
                s.__setattr__(attr, value)

        return s

    def derivative(self, order=1):
        """Returns the derivative of y-coordinates as a function of x-coordinates.

        Args:
            order (number, optional): derivative order. Default is 1.

        Returns:
            Derivative spectrum
        """
        x, y = arraymanip.derivative(self.x, self.y, order=order)
        s    = Spectrum(x=x, y=y)

        # transfer attrs
        for attr in self.get_attrs():
            value = copy.deepcopy(self.__dict__[attr])
            s.__setattr__(attr, value)

        return s
    
    def smooth(self, n=8):
        """Returns the moving average of the data.

        Args:
            n (int, optional): number of points to average. Default is 8.

        Returns:
            spectrum of length given by (len(x)-n+1).
        """
        x = arraymanip.moving_average(self.x, n=n)
        y = arraymanip.moving_average(self.y, n=n)
        s = Spectrum(x=x, y=y)

        # special
        s._peaks = self.peaks

        # transfer attrs
        for attr in self.get_attrs():
            value = copy.deepcopy(self.__dict__[attr])
            s.__setattr__(attr, value)

        return s



        ##########################
        # reset check attributes #
        ##########################
        self._step         = None
        self._monotonicity = None

    ########################
    # calculation and info #
    ########################
    def calculate_area(self, *args, **kwargs):
        """Returns the calculated area under the curve. Wrapper for `numpy.trapz()`_.

        Usage:
            >>> s.calculate_area()  # returns the area for the whole dataset
            >>> s.calculate_area(0, 10)  # returns the area between x=0 and 10
            >>> s.calculate_area((0, 10), (90, 100))
            
        Warning:
            the last line of the `Usage section` works, but will 
            typically not return a desirable value, because the ``trapz()``
            algorithm will consider the area between 10 and 90 as a rectangle.
            The correct approach in this case would be to calculated the area
            between 0 and 10 and between 90 and 100 separately.

            >>> area = s.calculate_area(0, 10) + s.calculate_area(90, 100)


        Args:
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            number

        .. _numpy.trapz(): https://numpy.org/doc/stable/reference/generated/numpy.trapz.html
        """
        if kwargs != {} or args != ():
            s = self._extract(*args, **kwargs)
        else:
            s = self.copy()
        return np.trapz(y=s.y, x=s.x)

    def calculate_x_sum(self, *args, **kwargs):
        """Returns sum of x elements within a range.
        
        Usage:
            s.calculate_x_sum()  # returns the x sum for the whole dataset
            s.calculate_x_sum((0, 10), (90, 100))  # returns the x sum from data beteen x=0 and 10 and between x=90 and 100

        Args:
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            number
        """
        if kwargs != {} or args != ():
            s = self._extract(*args, **kwargs)
        else:
            s = self.copy()
        return sum(s.x)

    def calculate_y_sum(self, *args, **kwargs):
        """Returns sum of y elements within a range.
        
        Usage:
            s.calculate_y_sum()  # returns the y sum for the whole dataset
            s.calculate_y_sum((0, 10), (90, 100))  # returns the y sum from data beteen x=0 and 10 and between x=90 and 100

        Args:
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            number
        """
        if kwargs != {} or args != ():
            s = self._extract(*args, **kwargs)
        else:
            s = self.copy()
        return sum(s.y)

    def polyfit(self, deg, *args, **kwargs):
        """Fit data with a polynomial. Wrapper for `numpy.polyfit()`_.

        Usage:
            popt, model, R2 = polyfit(deg)

        Args:
            deg (int): degree of the fitting polynomial.
            ranges (list, ooptional): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            popt, f(x) model, R2

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        s = self._extract(*args, **kwargs)
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
    def plot(self, ax=None, offset=0, shift=0, roll=0, factor=1, calib=1, smooth=1, ranges=None, switch=False, **kwargs):
        """Plot spectrum. Wrapper for `matplotlib.pyplot.plot()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number, optional): defines a vertical offset. Default is 0.
            shift (number, optional): horizontal shift value. Default is 0.
            factor (number, optional): multiplicative factor on the y axis.
                Default is 1.
            calib (number, optional): multiplicative factor on the x axis.
                Default is 1.
            smooth (int, optional): number of points to average data (moving average).
                Default is 1.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.
            switch (bool, optional): Switch x and y axis.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.

        Returns:
            `Line2D`_

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        if ranges is not None:
            self = self._extract(ranges=ranges)
        x = self.x
        y = self.y

        if switch:
            _x = x
            x = y
            y = _x

        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figmanip.figure()

        x = (x*calib) + shift
        y = y*factor + offset

        if roll != 0:
            assert numanip.is_integer(roll), 'roll must be an integer'
            x = np.roll(x, roll)
            x, y = arraymanip.sort(x, x, y)

        # if 'label' not in kwargs and hasattr(self, 'label'):
        #     kwargs['label'] = self.label

        if smooth > 1:
            x = arraymanip.moving_average(x, int(smooth))
            y = arraymanip.moving_average(y, int(smooth))
        
        return ax.plot(x, y, **kwargs)

    def inspect(self):
        """Plots data in a inspection window.

        If the mouse pointer hovers over a datapoint, the index of this datapoint is showed.
        """
        raise NotImplementedError('This is not implemented yet.')

    def metadata(self):
        """Opens a windows with all user defined variables for inspection.
        """
        raise NotImplementedError('This is not implemented yet.')

    ###########
    # Special #
    ###########
    def find_peaks(self, prominence=None, width=4, moving_average_window=4, ranges=None):
        """Find peaks. Wrapper for `scipy.signal.find_peaks()`_.

        Sets :py:attr:`peaks` attribute.

        Args:
            prominence (number, optional): minimum prominence of peaks in percentage
                of the maximum prominence [max(y) - min(y)]. Default is 5.
            width (number, optional): minimum number of data points defining a peak.
            moving_average_window (int, optional): window size for smoothing the
                data for finding peaks. Default is 4.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            None

        .. _scipy.signal.find_peaks(): https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
        """
        if ranges is None:
            s = self
        else:
            s = self._extract(ranges)

        # check monotonicity
        if s.monotonicity is None:
            s.check_monotonicity()

        # step
        if s.step is None:
            s.check_step()

        self.peaks.find(x=s.x, y=s.y, prominence=prominence, width=width, moving_average_window=moving_average_window)

    def fit_peak(self, asymmetry=False, moving_average_window=1, method='least_squares', ranges=None):     
        """Fits one peak. Initial guess is based on the maximum y value.

        Args:
            asymmetry (bool or dict, optional): if True, fits each half of the
                with a different width. Default is False
            moving_average_window (int, optional): window size for smoothing the
                data for finding the peak. Default is 1.
            method (str, optional): Name of the fitting method to use. See methods
                available on `lmfit.minimize()`_ documentation.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            None
        """
        if ranges is None:
            x0 = self.x
            y0 = self.y
        else:
            ranges = self._validate_ranges(ranges)
            s0 = self._extract(ranges)
            x0 = s0.x
            y0 = s0.y

        # smoothing
        assert moving_average_window > 0, f'moving_average_window must be positive different than 0, not {moving_average_window}'
        if moving_average_window == 1:
            x = x0
            y = y0
        else:
            x = arraymanip.moving_average(x0, moving_average_window)
            y = arraymanip.moving_average(y0, moving_average_window)

        # guess amp and c
        amp = max(y)
        c = x[np.argmax(y)]

        # guess fwhm
        try:
            w1 = x[np.argmax(y)] - x[:np.argmax(y)][::-1][arraymanip.index(y[:np.argmax(y)][::-1], max(y)/2)]
        except ValueError:
            w1 = x[np.argmax(y):][arraymanip.index(y[np.argmax(y):], max(y)/2)] - x[np.argmax(y)]
        try:
            w2 = x[np.argmax(y):][arraymanip.index(y[np.argmax(y):], max(y)/2)] - x[np.argmax(y)]
        except ValueError:
            w2 = x[np.argmax(y)] - x[:np.argmax(y)][::-1][arraymanip.index(y[:np.argmax(y)][::-1], max(y)/2)]
        w = w1 + w2

        # x, y = derivative(moving_average(self.x, 10), moving_average(self.y, 10))
        # w = np.abs(x[np.argmin(y)] - x[np.argmax(y)])
        if w == 0:
            w = 0.1*(max(self.x)-min(self.x))

        # peaks
        self.peaks.clear()
        if asymmetry:
            self.peaks.append(amp=amp, c=c, w1=w1, w2=w2)
        else:
            self.peaks.append(amp=amp, c=c, w=w)

        # fit
        self.fit_peaks(method=method, ranges=ranges)

    def fit_peaks(self, method='least_squares', ranges=None):
        """Fit peaks. Wrapper for `lmfit.minimize()`_.

        Args:
            yerr (array, optional): data uncertainty. 
            method (str, optional): Name of the fitting method to use. See methods
                available on `lmfit.minimize()`_ documentation.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            None

        .. _lmfit.minimize(): https://lmfit.github.io/lmfit-py/fitting.html
        """
        if ranges is None:
            self.check_monotonicity()
            x = self.x
            y = self.y
        else:
            ranges = self._validate_ranges(ranges)
            s = self._extract(ranges)
            s.check_monotonicity()
            x = s.x
            y = s.y

        # check if peaks is defined
        if len(self.peaks) == 0:
            raise ValueError('No peaks to fit.\nRun Spectrum.find_peaks() or set s.peaks manually.')

        # fit
        return self.peaks.fit(x, y, method=method, update_peaks=True)


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
        step (number, read only): step size for x-coordinates.
        length (number, read only): length of x-coordinates vector.
        x (number, read only): x-coordinates.
        monotonicity (str or None, read only): monotonicity of the x-array.
        peaks (:py:attr:`Peaks`): each entry represents a parameter of a peak.

    **Computed (read-only) attributes:**
        area (list)
            area of each spectrum

        sum (:py:attr:`Spectrum`)
            sum of every spectrum

        average (:py:attr:`Spectrum`)
            average of every spectrum

        map (:py:attr:`Image`)
            2D representation of spectra

        calculated_shift (list)
            calculated x add. factor

        calculated_offset (list)
            calculated y additive factor

        calculated_roll (list)
            calculated x roll

        calculated_calib (list)
            calculated x multiplicative factor

        calculated_factor (list)
            calculated y multiplicative factor

    Write-only attributes:
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
                      'calculated_roll', 'peaks']
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
        self._data  = []
        self._folderpath = ''
        self._filepath   = ''

        # check
        self._length       = None
        self._step         = None
        self._x            = None
        self._monotonicity = None

        # calculated values
        self._calculated_calib   = None
        self._calculated_factor = None
        self._calculated_offset = None
        self._calculated_shift  = None
        self._calculated_roll   = None

        # special
        self._peaks = Peaks()

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
            self._data = [-1]*n
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

    ##############
    # attributes #
    ##############
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

    @property
    def folderpath(self):
        return self._folderpath
    @folderpath.setter
    def folderpath(self, value):
        if value is None:
            value = ''
        elif isinstance(value, str) or isinstance(value, Path):
            self._folderpath = value
        else:
            raise TypeError(r'Invalid type ' + str(type(value)) + 'for folderpath\folderpath can only be str or pathlib.Path type.')
    @folderpath.deleter
    def folderpath(self):
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
            raise TypeError(r'Invalid type ' + str(type(value)) + 'for filepath\nfilepath can only be str or pathlib.Path type.')
    @filepath.deleter
    def filepath(self):
        raise AttributeError('Cannot delete object.')  
    
    ###################################
    # computed (read-only) attributes #
    ###################################
    @property
    def area(self):
        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self[i].area
        return temp
    @area.setter
    def area(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @area.deleter
    def area(self):
        raise AttributeError('Cannot delete object.')

    @property
    def sum(self):
        return self.calculate_sum()
    @sum.setter
    def sum(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @sum.deleter
    def sum(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def average(self):
        return self.calculate_average()
    @average.setter
    def average(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @average.deleter
    def average(self):
        raise AttributeError('Attribute cannot be deleted.')

    @property
    def map(self):
        return self.calculate_map()
    @map.setter
    def map(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @map.deleter
    def map(self):
        raise AttributeError('Attribute cannot be deleted.')

    #########################
    # write-only attributes #
    #########################
    @property
    def shift(self):
        pass
    @shift.setter
    def shift(self, value):
        self.set_shift(value=value)
    @shift.deleter
    def shift(self):
        raise AttributeError('Cannot delete object.')

    @property
    def roll(self):
        pass
    @roll.setter
    def roll(self, value):
        self.set_roll(value=value)
    @roll.deleter
    def roll(self):
        raise AttributeError('Cannot delete object.')

    @property
    def calib(self):
        pass
    @calib.setter
    def calib(self, value):
        self.set_calib(value)
    @calib.deleter
    def calib(self):
        raise AttributeError('Cannot delete object.')

    @property
    def factor(self):
        pass
    @factor.setter
    def factor(self, value):
        self.set_factor(value)
    @factor.deleter
    def factor(self):
        raise AttributeError('Cannot delete object.')

    @property
    def offset(self):
        pass
    @offset.setter
    def offset(self, value):
        self.set_offset(value)
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

    def __getitem__(self, item):
        if isinstance(item, int):
            return self._data[item]
        elif isinstance(item, slice):
            ss = Spectra(self._data[item])

            # # transfer attrs
            # for attr in self._get_user_attrs():
            #     value = copy.deepcopy(self.__dict__[attr])
            #     ss.__setattr__(attr, value)

            return ss        
        else:
            raise TypeError('Index must be int or slice, not {}'.format(type(item).__name__))

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

    ###########
    # support #
    ###########
    def _default_attrs(self):
        """return list with default attrs."""
        return ['_data', '_folderpath', '_filepath', '_peaks', 
                          '_calculated_shift', '_calculated_roll',
                          '_calculated_factor', 
                          '_calculated_offset', '_calculated_calib', 
                          '_monotonicity', '_x', '_step', '_length']
    
    def _get_user_attrs(self):
        """return attrs that are user defined."""
        return [key for key in self.__dict__.keys() if key not in self._default_attrs() and key.startswith('_') == False]
    
    def _validate_ranges(self, *args, **kwargs):
        """Check if ranges is the right format.

        Possible input formats are:

        xi, xf
        (xi, xf)
        (xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3)
        ((xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3))

        Input can be passed as positional arguments or with the kwarg `ranges`.

        If xi=None, this item is replaced by the min value of the dataset.
        If xf=None, this item is replaced by the max value of the dataset.

        Returns:
            ranges in the following format:
                ((xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3))
        """
        error_message = 'Wrong input. Please check ranges input.\n' +\
                        f'args passed: {args}\n' +\
                        f'kwargs passed: {kwargs}\n\n' +\
                        'ranges should be a pair (x_init, x_final) or a list of pairs.\n' +\
                        'See examples below:\n' +\
                        '\n' +\
                        'ranges = (x_init, x_final)\n' +\
                        '\n' +\
                        'ranges = ((x_init_1, x_final_1), (x_init_2, x_final_2), ...)\n' +\
                        '\n' +\
                        'Use None to indicate the minimum or maximum x value of the data.\n'
        if kwargs != {} and args != ():
            raise AttributeError(error_message)
        if any([item not in ['ranges', ] for item in kwargs.keys()]):
            raise AttributeError(error_message)
        if 'ranges' in kwargs:
            args = kwargs['ranges']

        # get min and max range values
        if self.x is not None:
            vmin = min(self.x)
            vmax = max(self.x)
        else:
            vmin = min(min(s.x) for s in self)
            vmax = max(max(s.x) for s in self)

        ##############
        # fix format #
        ##############
        # empty input
        if args is None:
            return ((vmin, vmax),)
        elif len(args) == 0:
            return ((vmin, vmax),)
        
        # (xi, xf)
        # ((xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3))
        elif len(args) == 1:
            if isinstance(args[0][0], Iterable):
                if len(args[0][0]) == 2:
                    args = args[0]
                else:
                    raise AttributeError(error_message)

        # xi, xf
        # (xi_, xf_1), (xi_2, xf_2)
        elif len(args) == 2:
            if isinstance(args[0], Iterable):
                args = tuple(args)
            else:
                args = ((args), )

        # (xi_, xf_1), (xi_2, xf_2), (xi_3, xf_3)
        args = list(args)
        for i, pair in enumerate(args):
            pair    = list(pair)
            args[i] = pair
            if len(pair) == 2:
                if pair[0] == None: pair[0] = vmin
                if pair[1] == None: pair[1] = vmax
            else:
                raise AttributeError(error_message)        
        return args
    
    def _gather_ys(self, *args, **kwargs):
        """Return two lists, x and y's within a range.
        
        This structure speeds up some operations.
        
        Args:
            ranges (list, optional): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            x, y's
        """
        # check x
        if self.x is None:
            self.check_same_x()

        # validate ranges
        ranges = self._validate_ranges(*args, **kwargs)

        ys = np.zeros((self.length, len(self)))
        for i in range(len(self)):
            ys[:, i] = self[i].y
        if ranges is None:
            x = self[0].x
        else:
            try:
                x, ys = arraymanip.extract(self[0].x, ys, ranges=ranges)
            except RuntimeError:
                raise RuntimeError(f'It seems like all spectra has no data points within range: {ranges}.\nPlease, fix ranges so all spectra have at least one data point within range.')
        
        return x, ys
    
    #########
    # basic #
    #########
    def append(self, *args, **kwargs):
        """Append spectrum to the spectrum list.

        Usage:
            >>> ss = br.Spectra()
            >>> 
            >>> ss.append(s)
            >>> ss.append(s1, s2, s3)
            >>> ss.append([s1, s2, s3])
            >>> 
            >>> ss.append(s=s)
            >>> ss.append(s=[s1, s2, s3])

        Args:
            s (Spectrum obj or list): Spectrum object to be added or
                list of Spectrum.

        Returns:
            None

        See Also:
            :py:func:`Spectra.remove`.
        """
        ###################################
        # asserting validity of the input #
        ###################################
        error_message = 'Wrong input. Spectrum cannot be added. Please, use one ' +\
                        'of the examples below:\n' +\
                        '\n' +\
                        'ss = br.Spectra()\n' +\
                        '\n' +\
                        'ss.append(s)' +\
                        'ss.append(s1, s2, s3)' +\
                        'ss.append([s1, s2, s3])' +\
                        '\n' +\
                        'ss.append(s=s)' +\
                        'ss.append(s=[s1, s2, s3])' +\
                        '\n' +\
                        's, s1, s2, ... must be Spectrum'
        if kwargs != {} and args != ():
            raise AttributeError(error_message)
        if kwargs == {} and args == ():
            raise AttributeError(error_message)
        if len(kwargs) >= 2:
            raise AttributeError(error_message)
        
        ################
        # loading data #
        ################
        # keyword arguments
        if 's' in kwargs:
            args = kwargs['s']
        
        # positional arguments
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

        ###########################
        # reset calculated values #
        ###########################
        self._calculated_calib   = None
        self._calculated_factor = None
        self._calculated_offset = None
        self._calculated_shift  = None
        self._calculated_roll   = None

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

        ##########################
        # reset check attributes #
        ##########################
        # self._length       = None
        # self._step         = None
        # self._x            = None
        # self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        self._calculated_calib   = None
        self._calculated_factor = None
        self._calculated_offset = None
        self._calculated_shift  = None
        self._calculated_roll   = None
    
    def flip_order(self):
        """reorder spectra backwards:

        If ss.data = [s1, s2, s3], then flip_order() would make it 
        ss.data = [s3, s2, s1]

        Returns:
            None
        """
        self._data = self._data[::-1]

        ##########################
        # reset check attributes #
        ##########################
        # self._length       = None
        # self._step         = None
        # self._x            = None
        # self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        self._calculated_calib   = None
        self._calculated_factor = None
        self._calculated_offset = None
        self._calculated_shift  = None
        self._calculated_roll   = None

    def reorder_by_attr(self, attr, attrs2reorder=None, decreasing=False):
        """Reorder spectra based on a attr. attr list is also sorted.
    
        Args:
            attr (str): name of the reference attr. The attr must be a list of 
                numbers with same lenght of number of spectra.
            attrs (list of str, optional): list of other attrs that 
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

    #################
    # save and load #
    #################
    def save(self, folderpath=None, prefix='spectrum_', suffix='.dat', filenames=None, zfill=None, only_data=False, verbose=False, **kwargs):
        r"""Save spectra. Wrapper for `numpy.savetxt()`_.

        Attrs are saved as comments if only_data is False. Saving attrs to file
        is tricky because requires converting variables to string. Only attrs 
        that are of type: string, number, and list of number and strings are 
        saved somewhat correctly. Dictionaries are not saved. 

        Args:
            folderpath (string or pathlib.Path, optional): folderpath, folder handle. 
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
        
        # folderpath
        if 'filepath' in kwargs:
            raise ValueError('filepath is not a valid parameter.\nPlease, use folderpath')
        if folderpath is None:
            if folderpath == '' and self.folderpath == '':
                raise TypeError("Missing 1 required argument: 'folderpath'")
            else:
                folderpath = self.folderpath               
        folderpath = Path(folderpath)

        assert folderpath.exists(), f'dirpath does not exists.\nfolderpath: {folderpath}'
        assert folderpath.is_dir(), f'dirpath is not a directory.\nfolderpath: {folderpath}'

        # set filenames
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
                # filename = f'{prefix}' + f'{i}'.zfill(zfill) + f'{suffix}'
                filename = f'{prefix}' + f'{i_str}' + f'{suffix}'

            if verbose:  print(f'{i}/{len(self)-1}: {filename}')
            s.save(filepath=folderpath/filename, only_data=only_data, check_overwrite=False, verbose=False, **kwargs)
        
        # save folderpath
        self._folderpath = folderpath
        if verbose: print('Done!')
    
    def load(self, folderpath=None, string='*', verbose=False):
        """Load data from text files. Wrapper for `numpy.genfromtxt()`_.

        Args:
            folderpath (string or path object): folderpath. If None,
                Last used folderpath is used.
            string (str, optional): file names without this string will be ignored.
                Use '*' for matching anything. Default is '*'.
            verbose (bool, optional): turn verbose on and off. Default is `False`.

        Returns:
            None

        .. _numpy.genfromtxt(): https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html
        """
        # reset data
        self._data = []

        if folderpath is None:
            folderpath = self.folderpath
        
        # if isinstance(folderpath, str) or isinstance(folderpath, Path): 
        #     fl = filemanip.filelist(dirpath=folderpath, string=string)
        #     self.folderpath = folderpath
        # else:
        #     raise TypeError(f'folderpath is unknown format: {folderpath}')
        
        # check if folder points to a folder and exists
        folderpath = Path(folderpath)
        assert isinstance(folderpath, str) or isinstance(folderpath, Path), f'folderpath is unknown format\n{folderpath}'
        assert folderpath.exists(), f'folderpath does not exist\n{folderpath}'
        assert folderpath.is_dir(), f'folderpath must point to a file\n{folderpath}'
        
        # get filelist
        fl = filemanip.filelist(dirpath=folderpath, string=string)
        self.folderpath = folderpath
        
        # get data
        for i, filepath in enumerate(fl):
            if verbose: print(f'Loading: {folderpath}')
            if verbose: print(f'    {i+1}/{len(fl)}: {filepath.name}')
            self.append(Spectrum(filepath=filepath))

        if verbose: print('Done!')

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

    def _create_header(self, verbose=False):
        """Gather attrs to be saved to a file."""
        header = ''
        attrs = self._get_user_attrs()
        for name in attrs:
            value = self.__getattribute__(name)
            if value is None:
                header += f'{name}: None'  + '\n'
            elif isinstance(value, str):
                temp2 = str(value).replace('\n','\\n')
                header += f'{name}: \"{temp2}\"'  + '\n'
            elif isinstance(value, Iterable):
                value = list(value)
                for item in value:
                    if item is None:
                        item = 'None'
                    elif isinstance(item, str):
                        temp2 = str(item).replace('\n','\\n')
                        item = f'\"{temp2}\"'
                    elif numanip.is_number(value):
                        item = str(item)
                    else:
                        temp2 = str(item).replace('\n','\\n')
                        item = f'\"{temp2}\"'
                header += f'{name}: {value}'  + '\n'
            elif isinstance(value, dict):
                if verbose:
                    type_ = str(type(value))
                    print(r'Warning: Cannot save attr of type: ' + type_ + r'.\attr: '+ name + r'.\nTo turn off this warning, set verbose to False.')
            elif numanip.is_number(value):
                tosave = str(value)
                if tosave[-1] == '\n':
                    tosave = tosave[:-1]
                header += f'{name}: {tosave}'  + '\n'
            else:
                temp2 = str(value).replace('\n','\\n')
                header += f'{name}: \"{temp2}\"'  + '\n'
        return header[:-1]

    def save_all_single_file(self, filepath=None, only_data=False, ranges=None, check_overwrite=False, verbose=False, **kwargs):
        r"""Save all Spectra in one single file. Wrapper for `numpy.savetxt()`_.

        Args:
            filepath (string or path object, optional): filepath or file handle.
                If the filename ends in .gz, the file is automatically saved in
                compressed gzip format.
            only_data (bool, optional): If True, header and footer are ignored and
                only data is saved to the file.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                The additive factor will be calculated for spectra such the average value of
                data inside ranges is the same as for the first spectrum. Spectra
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
        # dirpath
        if filepath is None:
            if filepath == '' and self.filepath == '':
                raise TypeError("Missing 1 required argument: 'filepath'")
            else:
                filepath = self.filepath               
        filepath = Path(filepath)

        # check if filepath points to a file
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        # check x is the same
        if self.x is None:
            try:
                self.check_same_x()
            except ValueError:
                raise ValueError('Cannot save spectra in one file. x axis are different.\nMaybe try interpolating the x axis (Spectra.interp()) or use Spectra.save() to save spectra in multiple files.')
        
        # check overwrite
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if interact.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                        pass
                    else:
                        return
                else:
                    raise AttributeError('filepath not pointing to a file.')
                
        # gather ys
        x, ys = self._gather_ys(ranges=ranges)

        # kwargs
        if 'fmt' not in kwargs: # pick best format
            decimal = max([figmanip.n_decimal_places(x) for x in self.x])
            temp_decimal = max([figmanip.n_decimal_places(y) for y in arraymanip.flatten(ys)])
            if temp_decimal > decimal:
                decimal = copy.deepcopy(temp_decimal)
            kwargs['fmt'] = f'%.{decimal}f'
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'newline' not in kwargs:
            kwargs['newline'] = '\n'
        if 'comments' not in kwargs:
            kwargs['comments'] = '# '

        
        # dataset
        final = np.zeros((self.length, len(self)+1))
        final[:, 0]  = x
        final[:, 1:] = ys
        

        if only_data:
            if 'header' in kwargs:
                del kwargs['header']
            if 'footer' in kwargs:
                del kwargs['footer']
            np.savetxt(Path(filepath), final, **kwargs)
        else:
            if 'header' not in kwargs:
                kwargs['header'] = self._create_header(verbose=verbose)
            else:
                if kwargs['header'] == '':
                    kwargs['header'] = self._create_header(verbose=verbose)
                elif kwargs['header'][-1] != '\n':
                    kwargs['header'] += '\n'

        np.savetxt(Path(filepath), final, **kwargs)

    def load_from_single_file(self, filepath=None, only_data=False, verbose=False, **kwargs):
        """load multiple spectra from file.

        Args:
            filepath (str or pathlib.Path, optional): filepath. If None, the 
                last used filepath is used.
            onl
        if usecols
        """
        if 'delimiter' not in kwargs:
            kwargs['delimiter'] = ', '
        if 'comments' not in kwargs:
            kwargs['comments'] = '#'
        if 'usecols' not in kwargs:
            kwargs['usecols'] = None

        # check if filepath points to a file
        filepath = Path(filepath)
        assert filepath.is_file(), f'filepath must point to a file\n{filepath}'
        assert filepath.exists(), f'filepath does not exist\n{filepath}'

        # read data
        data = np.genfromtxt(Path(filepath), **kwargs)
        x    = data[:, 0]
        ys   = data[:, 1:]

        # artificial usecols
        if kwargs['usecols'] is None:
            kwargs['usecols'] = [i for i in range(data.shape[1])]

        # reset data
        self._data     = []

        # allocate data
        for i in range(ys.shape[1]):
            s = Spectrum(x=x, y=ys[:, i])
            self.append(s)

        # read header
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

    #################
    # check methods #
    #################
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
        # check spectra exists
        if len(self) == 0:
            raise ValueError('no spectra found')
        
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

    def check_step(self, max_error=None):
        """Check step between data points in the x-coordinates.

            If data has a well defined step size, it sets 
            :py:attr:`Spectra.step` = step.
            Otherwise, it raises an error.

            Args:
                max_error (number, optional): percentage value (in terms of the
                    average x step) of the maximum allowed error. If None, the 
                    Default value from settings will be used.

            Two checks are performed:

            1. Checks if the step between two data points is the same
            through out the x vector for each spectrum (vector uniformity).
            See Spectrum.check_step().

            (max(steps) - min(steps))/np.mean(steps) * 100 < max_error

            2. Checks if the step is the same between different spectra.

            (step[i]-step[i+1])/avg(step) * 100 < max_error

            Returns:
                None

            Raises:
                ValueError: If condition 1, or 2 are not satisfied.

            See Also:
                :py:func:`Spectra.check_length`, :py:func:`Spectra.check_same_x`.
        """
        # check spectra exists
        if len(self) == 0:
            raise ValueError('no spectra found')
        
        # max error
        if max_error is None:
            max_error = settings.MAX_ERROR_STEP_X

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

    def check_same_x(self, max_error=None):
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
        # check spectra exists
        if len(self) == 0:
            raise ValueError('no spectra found')
        
        # if only one spectra exists, then x is immediately defined
        if len(self) == 1:
            self._x = self[0].x
            self._length = len(self.x)
            return
        
        # max error
        if max_error is None:
            max_error = settings.MAX_ERROR_STEP_X

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

    ##################
    # base modifiers #
    ##################
    def set_shift(self, value=None):
        """Shift data recursively.

        if value is none, calculated values will be used.

        Args:
            value (number or list, optional): value will be added to x-coordinates.
                If None, it will look for calculated values from Spectra.calculate_shift().

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_shift()`
        """
        # if value is None, check calculated_shift
        if value is None:
            if self.calculated_shift is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_shifts().')
            value = self.calculated_shift
        else:
            # check if value is a number
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        # if all shifts are zero, does nothing
        if  all(x==0 for x in value):
            return

        # set values
        for i in range(len(self)):
            self.data[i].set_shift(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # self._length       = None
        # self._step         = None
        self._x            = None
        # self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        self._calculated_calib   = None
        # self._calculated_factor = None
        # self._calculated_offset = None
        self._calculated_shift  = None
        self._calculated_roll   = None

    def set_roll(self, value=None):
        """Roll y data recursively.

        if value is none, calculated values will be used.

        Args:
            value (number or list, optional): value to roll x-coordinates.
                If None, it will look for calculated values from Spectra.calculate_roll().

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_shift()`
        """
        # if value is None, check calculated_shift
        if value is None:
            if self.calculated_roll is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_shift().')
            value = self.calculated_roll
        else:
            # check if value is a number
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        # if all rolls are zero, does nothing
        if  all(x==0 for x in value):
            return

        # set values
        for i in range(len(self)):
            self.data[i].set_roll(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # self._length       = None
        # self._step         = None
        # self._x            = None
        # self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        self._calculated_calib   = None
        # self._calculated_factor = None
        # self._calculated_offset = None
        self._calculated_shift  = None
        self._calculated_roll   = None

    def set_calib(self, value):
        """Apply multiplicative x factor recursively.

        Args:
            value (number, list, function): calibration value. 
                If number, x-coordinates will be multiplied by this value.
                If list, polynomial coeff. (highest power first) is expected and
                a function new_x = f(x) will be set using np.polyval().
                If function, a function of type new_x = f(x) is expected.

        Returns:
            None
        """
        # #######
        # # run #
        # #######
        # if isinstance(value, Iterable):
        #     f = lambda x: np.polyval(value, x)
        #     self._x = np.array([f(x) for x in self.x])
        # elif isinstance(value, function):
        #     self._x = np.array([f(x) for x in self.x])
        # else:
        #     self._x = self.x*value


        # # check if value is a number
        # if isinstance(value, Iterable) == False:
        #     value = [value]*len(self)

        # # value must be the right length
        # assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(values)}\nnumber of spectra: {len(self)}'

        # set values ===========================================================
        for i in range(len(self)):
            self[i].set_calib(value=value)

        ##########################
        # reset check attributes #
        ##########################
        # self._length       = None
        self._step         = None
        self._x            = None
        self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        self._calculated_calib   = None
        # self._calculated_factor = None
        # self._calculated_offset = None
        self._calculated_shift  = None
        self._calculated_roll   = None

    def set_factor(self, value=None):
        """Apply multiplicative y factor recursively.

        if value is none, calculated values will be used.

        Args:
            value (number or list, optional): value will be multiplied to y-coordinates.
                If None, it will look for calculated values from Spectra.calculate_factor().

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_factor()`
        """
        if value is None:
            if self.calculated_factor is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_factor().')
            value = self.calculated_factor
        else:
            # check if value is a number
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        # set values
        for i in range(len(self)):
            self[i].set_factor(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # self._length       = None
        # self._step         = None
        # self._x            = None
        # self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        # self._calculated_calib   = None
        self._calculated_factor = None
        self._calculated_offset = None
        # self._calculated_shift  = None
        # self._calculated_roll   = None

    def set_offset(self, value=None):
        """Apply additive y factor recursively.

        if value is none, calculated values will be used.

        Args:
            value (number or list, optional): value will be added to y-coordinates.
                If None, it will look for calculated values from Spectra.calculate_offset().

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_offset()`, :py:func:`Spectra.floor()`
        """
        if value is None:
            if self.calculated_offset is None:
                raise ValueError('values not defined. Please, pass values or use Spectra.calculate_offset().')
            value = self.calculated_offset
        else:
            # check if value is a number
            if isinstance(value, Iterable) == False:
                value = [value]*len(self)

        # value must be the right length
        assert len(value) == len(self), f'value must have the same number of items as the number of spectra.\nnumber of values: {len(value)}\nnumber of spectra: {len(self)}'

        # set values
        for i in range(len(self)):
            self[i].set_offset(value=value[i])

        ##########################
        # reset check attributes #
        ##########################
        # self._length       = None
        # self._step         = None
        # self._x            = None
        # self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        # self._calculated_calib   = None
        self._calculated_factor = None
        self._calculated_offset = None
        # self._calculated_shift  = None
        # self._calculated_roll   = None

    #############
    # modifiers #
    #############
    def align(self, mode='cc', **kwargs):
        """Uses calculate_shift and set_shift.

        Args:
            same as calculate_shift() or calculate_roll()

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_shift`, :py:func:`Spectra.calculate_roll`
        """
        if mode in ['cross-correlation', 'cc']:
            self.calculate_roll(mode=mode, **kwargs)
            self.set_roll()
        else:
            self.calculate_shift(mode=mode, **kwargs)
            self.set_shift()

    def normalize(self, mode='max', **kwargs):
        """Uses Spectra.calculate_factor() and Spectra.set_factor() to normalize spectra.

        Args:
            same as calculate_factor()

        Returns:
            None

        See Also:
            :py:func:`Spectra.calculate_factor`
        """
        self.calculate_factor(mode=mode, **kwargs)
        self.set_factor()

    def floor(self, *args, **kwargs):
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
            s.floor(*args, **kwargs)

    def flip(self):
        """Flip x axis.

        Returns:
            None
        """
        for s in self:
            s.flip()
        
        ##########################
        # reset check attributes #
        ##########################
        # self._length       = None
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
            None
        """
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
                    return

        for s in self:
            s.interp(x=x, start=start, stop=stop, num=num, step=step)
        
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
        # self._calculated_calib   = None
        # self._calculated_factor = None
        # self._calculated_offset = None
        # self._calculated_shift  = None
        # self._calculated_roll   = None

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
            for s in self.data:
                s.crop(start, stop)
        elif start is not None or stop is not None:
            for s in self.data:
                s.crop(start, stop)

        ##########################
        # reset check attributes #
        ##########################
        self._length       = None
        # self._step         = None
        self._x            = None
        # self._monotonicity = None

        ###########################
        # reset calculated values #
        ###########################
        # self._calculated_calib   = None
        # self._calculated_factor = None
        # self._calculated_offset = None
        # self._calculated_shift  = None
        # self._calculated_roll   = None

    def switch(self):
        """Switch x and y axis"""
        for s in self:
            s.switch()
        
        ##########################
        # reset check attributes #
        ##########################
        # self._length       = None
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

    ##############
    # extractors #
    ##############
    def _extract(self, *args, **kwargs):
        """Same as copy(), but attributes are not copied to the new object."""
        
        # check if extract is really necessary
        if kwargs == {} and args == ():
            return copy.deepcopy(self)
        else:
            ranges = self._validate_ranges(*args, **kwargs)
        if (ranges[0][0] <= min([min(s.x) for s in self]) and ranges[0][1] >= max([max(s.x) for s in self]) ) and len(ranges) == 1:
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
        ranges = self._validate_ranges(*args, **kwargs)
        ss = Spectra(n=len(self))
        for i, s in enumerate(self.data):
            try:
                ss[i] = s.copy(ranges=ranges)
            except RuntimeError:
                raise RuntimeError(f'It seems like spectrum number {i} has no data points within range: {ranges}.\nPlease, fix ranges (or delete spectrum) so all spectra have at least one data point within range.')
        return ss

    def copy(self, *args, **kwargs):
        """Return a copy of the object with data contained in a range.

        Usage:
            >>> ss2.copy(ss1)     # ss2 is now a copy of ss1
            >>> ss2 = ss1.copy()  # ss2 is now a copy of ss1
            >>>
            >>> # ss2 will have only data between x=0 and 10 and between x=90 and 100
            >>> ss2 = ss1.copy((0, 10), (90, 100))  

        Args:
            ss (Spectra, optional): Spectra to be copied. See usage.
            ranges (list, optional): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            :py:attr:`Spectra`
        """
        #############################
        # check if input is Spectra #
        #############################
        ss = None
        if 'ss' in kwargs:
            ss = kwargs['ss']
        elif len(args) == 1:
            if isinstance(args[0], Spectra):
                ss = args[0]
        
        ######################################
        # if Spectra is passed, copy Spectra #
        ######################################
        if ss is not None:
            if isinstance(ss, Spectra):
                self._data         = ss.data
                self._folderpath   = ss.folderpath
                self._filepath     = ss.filepath
                self._length       = ss.length
                self._step         = ss.step
                self._x            = ss.x
                self._monotonicity = ss.monotonicity

                # special
                self._peaks = ss.peaks

                # user defined attrs
                for attr in self._get_user_attrs():
                    self.__delattr__(attr)
                for attr in ss._get_user_attrs():
                    self.__setattr__(attr, ss.__dict__[attr])
            else:
                raise TypeError('Only type br.Spectra can be copied to type br.Spectra')
            return
        #######################################
        # Otherwise, a type ranges is assumed #
        #######################################
        else:
            ss = self._extract(*args, **kwargs)

            # transfer attrs
            for attr in self._get_user_attrs():
                value = copy.deepcopy(self.__dict__[attr])
                ss.__setattr__(attr, value)

        return ss
    
    def concatenate(self):
        """Return spectrum of concatenate spectra.
        
        Note:
            attrs are copied to the returned spectrum, but attrs from each spectrum is lost.

        Returns:
            Spectrum
        """
        x = np.concatenate([s.x for s in self.data])
        y = np.concatenate([s.y for s in self.data])
        final = Spectrum(x=x, y=y)
        
        # transfer attrs
        default_attrs = final._default_attrs()
        for attr in self._get_user_attrs():
            if attr not in default_attrs:
                value = copy.deepcopy(self.__dict__[attr])
                final.__setattr__(attr, value)

        return final

    def calculate_sum(self, *args, **kwargs):
        """Returns Spectrum object with the sum of all spectra.

        All spectra must have the same x-coordinates. This is verified.

        attrs are copied to the final spectrum, but attrs from each spectrum is lost.

        Args:
            ranges (list, optional): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            :py:class:`Spectrum` object.
        """
        # gather ys
        x, ys = self._gather_ys(*args, **kwargs)

        # calculate sum
        y = np.zeros(len(x))
        for i in range(len(self)):
            y += ys[:, i]
        final = Spectrum(x=x, y=y)
        
        # transfer attrs
        default_attrs = final._default_attrs()
        for attr in self._get_user_attrs():
            if attr not in default_attrs:
                value = copy.deepcopy(self.__dict__[attr])
                final.__setattr__(attr, value)

        return final

    def calculate_average(self, *args, **kwargs):
        """Returns Spectrum object with the average of all spectra.

        All spectra must have the same x-coordinates. This is verified.

        attrs are copied to the final spectrum, but attrs from each spectrum is lost.

        Args:
            ranges (list, optional): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            :py:class:`Spectrum` object.
        """
        # gather ys
        x, ys = self._gather_ys(*args, **kwargs)

        # calculate sum
        y = np.zeros(len(x))
        for i in range(len(self)):
            y += ys[:, i]
        y = y/len(self)

        final = Spectrum(x=x, y=y)
        
        # transfer attrs
        default_attrs = final._default_attrs()
        for attr in self._get_user_attrs():
            if attr not in default_attrs:
                value = copy.deepcopy(self.__dict__[attr])
                final.__setattr__(attr, value)

        return final

    def calculate_map(self, axis=0, *args, **kwargs):
        """Return image representation of spectra.

        All spectra must have the same x-coordinates. This is verified.

        attrs are copied to the final image, but attrs from each spectrum is lost.

        Args:
            axis (int, optional): Image axis along which spectra will be laid out.
                If axis = 0, spectra will be placed horizontally (each spectrum 
                will be a "row of pixels"). If axis = 1, spectra will be placed
                vertically (each spectrum will be a "column"). Default is 0.
            ranges (list, optional): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            :py:class:`Image`.
        """
        # check axis
        if axis != 0 and axis != 1:
            raise ValueError('axis must be 0 or 1')

        # gather ys
        y, ys = self._gather_ys(*args, **kwargs)
        x = None

        if axis == 1:
            ys = ys.transpose()
            x = y
            y = None

        final = Image(data=ys)
        final.x_centers = x
        final.y_centers = y
        
        # transfer attrs
        default_attrs = final._default_attrs()
        for attr in self._get_user_attrs():
            if attr not in default_attrs:
                value = copy.deepcopy(self.__dict__[attr])
                final.__setattr__(attr, value)

        return final
    
    ########################
    # calculation and info #
    ########################
    def calculate_shift(self, mode='cross-correlation', **kwargs):
        """Calculate x-coord. shifts so dataset is aligned.

        Result is a list of shift values that is save in the attr:

            >>> ss.calculated_shift

        Args:
            mode (string, optional): method used for calculating shifts.
                The current options are: 

                1) 'cross-correlation' or 'cc'

                2) 'max'

                3) 'peaks' or 'peak'

        Modes may have additional parameters:

        `cross-correlation`
            ref_spectrum (int or str, optional): index of the spectrum used
                to cross-correlate with all other spectra, i.e., all spectra 
                will be aligned to ref_spectrum. If ref_spectrum = 'seq' or 
                'sequential', cross-correlation is performed between subsequent
                spectra. Default is 0.
            ranges (list, optional): a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                Use None to indicate the minimum or maximum x value of the data.

        `peaks`
            ref_spectrum (int, optional)
                index of the spectrum to which all
                other spectra will be aligned to. Default is 0.
            ref_peak (int, optional)
                peak used to calculate shifts. Default is 0.
            ref_value (number, optional)
                If not None, the center of ref_peak 
                for all spectra is set to ref_value. This overwrites ref_spectrum.
                Default is None.  

        `max`:
            ref_spectrum (int, optional)
                index of the spectrum to which all
                other spectra will be aligned to. Default is 0.
            ref_value (number, optional)
                If not None, the max y-coord. 
                for all spectra is set to ref_value. This overwrites ref_spectrum.
                Default is None.  
            ranges (list, optional)
                a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                Use None to indicate the minimum or maximum x value of the data.
                  
            
        Returns:
            None

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
        if mode in ['cross-correlation', 'cc']:
            old_calculated_roll = copy.copy(self.calculated_roll)
            self.calculate_roll(mode='cc', **kwargs)
            values = self.calculated_roll*self.step
            self._calculated_roll = old_calculated_roll

        #######
        # max #
        #######
        elif mode == 'max':
            # args
            if 'ranges' not in kwargs:
                ss = self
            elif kwargs['ranges'] is None:
                ss = self
            else:
                ss = self._extract(ranges=kwargs['ranges'])
            if 'ref_value' not in kwargs:
                ref_value = None
            else:
                ref_value = kwargs['ref_value']
            if 'ref_spectrum' not in kwargs:
                ref_spectrum = 0
            else:
                ref_spectrum = kwargs['ref_spectrum']
            
            # calculate
            if ref_value is None:
                j_ref = np.argmax(ss[ref_spectrum].y)
                ref_value = ss[ref_spectrum].x[j_ref]
            for i in range(len(ss)):
                j = np.argmax(ss[i].y)
                values[i] = -(ss[i].x[j] - ref_value)
        
        #########
        # peaks #
        #########
        elif mode in ['peaks', 'peak']:
            # check if peaks are defined
            assert len(self.peaks) > 0, 'Spectra does not have defined peaks.\nMaybe use ss.find_peaks() or ss.copy_peaks_from_spectra().'

            # args
            if 'ref_value' not in kwargs:
                ref_value = None
            else:
                ref_value = kwargs['ref_value']
            if 'ref_spectrum' not in kwargs:
                ref_spectrum = 0
            else:
                ref_spectrum = kwargs['ref_spectrum']
            if 'ref_peak' not in kwargs:
                ref_peak = 0
            else:
                ref_peak = kwargs['ref_peak']

            # calculate peak
            if ref_value is None:
                p = self.peaks._get_peaks_by_index(i1=ref_peak, i2=ref_spectrum)
                ref_value = p['c'].value
            values = self.peaks.get_values('c', ref_peak)
            values = -np.array(values)+ref_value  
        else:
            raise ValueError('mode not valid.\nValid modes: cross-correlation, max, peaks.')

        # save calculated values
        self._calculated_shift = values

    def calculate_roll(self, mode='cross-correlation', **kwargs):
        """Calculate x-coord. roll so dataset is aligned.

        Result is a list of int values that is save in the attr:

            >>> ss.calculated_roll

        Args:
            mode (string, optional): method used for calculating rolls.
                The current options are: 

                1) 'cross-correlation' or 'cc'

                2) 'max'

                3) 'peaks' or 'peak'

        Some modes may have additional parameters:

        `cross-correlation`
            ref_spectrum (int or str, optional)
                index of the spectrum used
                to cross-correlate with all other spectra, i.e., all spectra 
                will be aligned to ref_spectrum. If ref_spectrum = 'seq' or 
                'sequential', cross-correlation is performed between subsequent
                spectra. Default is 0.
            ranges (list, optional)
                a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                Use None to indicate the minimum or maximum x value of the data.

        `peaks`
            ref_spectrum (int, optional)
                index of the spectrum to which all
                other spectra will be aligned to. Default is 0.
            ref_peak (int, optional)
                peak used to calculate shifts. Default is 0.
            ref_value (number, optional)
                If not None, the center of ref_peak 
                for all spectra is set to ref_value. This overwrites ref_spectrum.
                Default is None.  

        `max`
            ref_spectrum (int, optional)
                index of the spectrum to which all
                other spectra will be aligned to. Default is 0.
            ref_value (number, optional)
                If not None, the max y-coord. 
                for all spectra is set to ref_value. This overwrites ref_spectrum.
                Default is None.  
            ranges (list, optional)
                a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                Use None to indicate the minimum or maximum x value of the data. 
            
        Returns:
            None

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
        if mode in ['cross-correlation', 'cc']:
            # args
            if 'ranges' not in kwargs:
                ss = self
            elif kwargs['ranges'] is None:
                ss = self
            else:
                ss = self._extract(ranges=kwargs['ranges'])
            if 'ref_spectrum' not in kwargs:
                ref_spectrum = 0

            if isinstance(ref_spectrum, str):
                if ref_spectrum in ['sequential', 'seq']:
                    for i in range(1, len(ss)):
                        cross_correlation = np.correlate(ss[i-1].y, ss[i].y, mode='full')
                        values[i] = np.argmax(cross_correlation) - (len(ss[i-1].y) - 1) + values[i-1]
                    ref_value = values[0]
                    values -= ref_value
                else:
                    raise ValueError('ref_spectrum must be `sequential` or a int')
            else:
                for i in range(len(self)):
                    cross_correlation = np.correlate(ss[ref_spectrum].y, ss[i].y, mode='full')
                    values[i] = np.argmax(cross_correlation)
                ref_value = values[ref_spectrum]
                values -= ref_value

        #######
        # max #
        #######
        elif mode == 'max':
            old_calculated_shift = copy.copy(self.calculated_shift)
            self.calculate_shift(mode='max', **kwargs)
            values = int(round(self.calculated_shift/self.step))
            self._calculated_shift = old_calculated_shift

        #########
        # peaks #
        #########
        elif mode in ['peaks', 'peak']:
            old_calculated_shift = copy.copy(self.calculated_shift)
            self.calculate_shift(mode='max', **kwargs)
            values = int(round(self.calculated_shift/self.step))
            self._calculated_shift = old_calculated_shift
        else:
            raise ValueError('mode not valid.\nValid modes: cross-correlation, max, peaks.')

        # save calculated values
        self._calculated_roll = values

    def calculate_factor(self, mode='fitted peaks', **kwargs):
        """Calculate mult. factor for spectra to be same height.

        Result is a list of shift values that is save in the attr:

            >>> ss.calculated_factor

        Args:
            mode (string, optional): method used to calculate the multiplicative factor.
                The current options are: 

                1) 'max'

                2) 'delta'

                3) 'area'

                4) 'peaks'

                5) 'peaks area'

        Some modes may have additional parameters:

        `max`, `delta`, `area`
            ref_spectrum (int, optional)
                index of the spectrum to be taken as 
                reference. Default is 0.
            ref_value (number, optional)
                the max, delta, or area 
                for all spectra is set to ref_value. ref_value overwrites 
                ref_spectrum. Default is None.
            ranges (list, optional)
                a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                Use None to indicate the minimum or maximum x value of the data. 

        `peaks` or `peaks area`
            ref_spectrum (int, optional)
                index of the spectrum to be used as
                reference. Default is 0.
            ref_peak (int, optional)
                peak used as reference. Default is 0.
            ref_value (number, optional)
                If not None, the amplitude (area)
                of ref_peak 
                for all spectra is set to ref_value. This overwrites ref_spectrum.
                Default is None.  

        Returns:
            None

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
            # args
            if 'ranges' not in kwargs:
                ss = self
            elif kwargs['ranges'] is None:
                ss = self
            else:
                ss = self._extract(ranges=kwargs['ranges'])
            if 'ref_value' not in kwargs:
                ref_value = None
            else:
                ref_value = kwargs['ref_value']
            if 'ref_spectrum' not in kwargs:
                ref_spectrum = 0
            else:
                ref_spectrum = kwargs['ref_spectrum']
            # calculation
            if ref_value is None:
                ref_value = max(ss.data[ref_spectrum].y)
            for i in range(len(ss)):
                values[i] = ref_value/max(ss.data[i].y)

        #########
        # delta #
        #########
        elif mode == 'delta':
            # args
            if 'ranges' not in kwargs:
                ss = self
            elif kwargs['ranges'] is None:
                ss = self
            else:
                ss = self._extract(ranges=kwargs['ranges'])
            if 'ref_value' not in kwargs:
                ref_value = None
            else:
                ref_value = kwargs['ref_value']
            if 'ref_spectrum' not in kwargs:
                ref_spectrum = 0
            else:
                ref_spectrum = kwargs['ref_spectrum']
            # calculation
            if ref_value is None:
                ref_value = max(ss.data[ref_spectrum].y) - min(ss.data[ref_spectrum].y)
            for i in range(len(ss)):
                values[i] = ref_value/(max(ss.data[i].y) - min(ss.data[i].y))
        
        ########
        # area #
        ########
        elif mode == 'area':
            # args
            if 'ranges' not in kwargs:
                ss = self
            elif kwargs['ranges'] is None:
                ss = self
            else:
                ss = self._extract(ranges=kwargs['ranges'])
            if 'ref_value' not in kwargs:
                ref_value = None
            else:
                ref_value = kwargs['ref_value']
            if 'ref_spectrum' not in kwargs:
                ref_spectrum = 0
            else:
                ref_spectrum = kwargs['ref_spectrum']
            # calculation
            if ref_value is None:
                ref_value = ss.data[ref_spectrum].area
            for i in range(len(ss)):
                values[i] = ref_value/ss.data[i].area

        #########                
        # peaks #
        #########                
        elif mode == 'peaks' or mode == 'peak':
            # check if peaks are defined
            assert len(self.peaks) > 0, 'Spectra does not have defined peaks.\nMaybe use ss.find_peaks() or ss.copy_peaks_from_spectra().'

             # args
            if 'ref_value' not in kwargs:
                ref_value = None
            else:
                ref_value = kwargs['ref_value']
            if 'ref_spectrum' not in kwargs:
                ref_spectrum = 0
            else:
                ref_spectrum = kwargs['ref_spectrum']
            if 'ref_peak' not in kwargs:
                ref_peak = 0
            else:
                ref_peak = kwargs['ref_peak']

            # calculate peak
            if ref_value is None:
                p = self.peaks._get_params_with_index(i1=ref_peak, i2=ref_spectrum)
                ref_value = p['amp'].value
            values = self.peaks.get_values('amp', ref_peak)
            values = ref_value/np.array(values)  
        
        ##############
        # peaks area #
        ##############
        elif mode == 'peaks area' or mode == 'peak area':
            raise NotImplementedError('sorry not implemented yet')
        else:
            raise ValueError('mode not valid.\nValid modes: max, delta, area, peak, peak area.')

        # save calculated values
        self._calculated_factor = values

    def calculate_offset(self, mode='average', **kwargs):
        """Calculate offsets to align data background.

        Result is a list of int values that is save in the attr:

            >>> ss.calculated_offset

        Args:
            mode (string, optional): method used to calculate the multiplicative factor.
                The current options are: 

                    1) 'average'

                    2) 'peaks'
        
        Some modes may have additional parameters:

        `average`
            ref_spectrum (int, optional)
                index of the spectrum to be taken as 
                reference. Default is 0.
            ref_value (number, optional)
                the average of data within ranges
                for all spectra is set to ref_value. ref_value overwrites 
                ref_spectrum. Default is None.  
            ranges (list, optional)
                a pair of x-coordinate values or a list of
                pairs. Each pair represents the start and stop of a data range.
                Use None to indicate the minimum or maximum x value of the data.

        Returns:
            None

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
            # args
            if 'ranges' not in kwargs:
                ss = self
            elif kwargs['ranges'] is None:
                ss = self
            else:
                ss = self._extract(ranges=kwargs['ranges'])
            if 'ref_value' not in kwargs:
                ref_value = None
            else:
                ref_value = kwargs['ref_value']
            if 'ref_spectrum' not in kwargs:
                ref_spectrum = 0
            else:
                ref_spectrum = kwargs['ref_spectrum']
            # calculation
            if ref_value is None:
                ref_value = np.mean(ss[ref_spectrum].y)
            for i in range(len(ss)):
                values[i] = ref_value-np.mean(ss[i].y)
        
        elif mode in ['peaks', 'peak']:
            raise NotImplementedError('sorry not implemented yet')
        else:
            raise ValueError('mode not valid.\nValid modes: average, peaks.')

        # final 
        self._calculated_offset = values

    def calculate_calib(self, values, mode='cross-correlation', deg=1, **kwargs):
        """Calculate calibration factor via :py:func:`Spectra.calculate_shift()`.

        The calibration factor is the shift value (x-coord) as a function of the
            values array.

        Result is a Spectrum (x=values, y=shifts) that is save in the attr:

            >>> ss.calculated_calib

        calculated_calib is a Spectrum type and has many attrs:

            ss.calculated_calib.x
                values
            ss.calculated_calib.y
                calculated shifts
            ss.calculated_calib.popt
                polynomial coeff. (highest power first.) that fit the calibration curve.
            ss.calculated_calib.model
                function f(x) of the calibration curve.
            ss.calculated_calib.R2
                R2 factor of the fitting.

        Args:
            values (list): value list. Values associated with each spectrum.
            mode (string, optional): method used to calculate the shifts. See 
                :py:func:`Spectra.calculate_shift()` for all available methods.
            deg (int, ooptional): degree of the fitting polynomial. Default is 1.
            **kwargs: kwargs are passed to ``Spectra.calculate_shift()``.

        Returns:
            None
        """
        # check number of values matches the number of spectra
        if len(self) != len(values):
            raise ValueError(f'number of values ({len(values)}) do not match the number of spectra ({len(self)})')

        # CALCULATION
        self.calculate_shift(mode=mode, **kwargs)

        # if mode == 'peaks' or mode == 'peak':
        #     centers = [c+self.calculated_shift[0] for c in self.calculated_shift]
        #     print(centers)
        # centers = -(self.calculated_shift + self.calculated_shift.ref_value)

        # save calculated values ===============================================
        self._calculated_calib = Spectrum(x=-self._calculated_shift, y=values)
        popt, model, R2 = self.calculated_calib.polyfit(deg=deg)
        self._calculated_calib.popt  = popt
        self._calculated_calib.model = model
        self._calculated_calib.R2    = R2

    def calculate_y_sum(self, *args, **kwargs):
        """Returns a list of the sum of y elements within a range for each spectra.
        
        Usage:
            >>> ss.calculate_y_sum()
            >>>
            >>> # returns the y sum from data beteen x=0 and 10 and between x=90 and 100
            >>> ss.calculate_y_sum((0, 10), (90, 100))  

        Args:
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            list
        """
        final = []
        for s in self:
            final.append(s.calculate_y_sum(*args, **kwargs))

        return final

    def polyfit(self, deg, *args, **kwargs):
        """Fit data recursively with a polynomial. Wrapper for `numpy.polyfit()`_.

        Args:
            deg (int): degree of the fitting polynomial.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.
                
        Returns:
            list with polynomial coefficients, highest power first.
            list with Model function f(x).
            list with R2.

        .. _numpy.polyfit(): https://numpy.org/doc/stable/reference/generated/numpy.polyfit.html
        """
        popt  = [0]*len(self)
        model = [0]*len(self)
        R2    = [0]*len(self)
        for i in range(len(self)):
            popt[i], model[i], R2[i] = self[i].polyfit(deg=deg, *args, **kwargs)

        return popt, model, R2

    ##########################        
    # plot and visualization #
    ##########################  
    def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, smooth=1, ranges=None, switch=False, **kwargs):
        """Plot spectra. Wrapper for `matplotlib.pyplot.plot()`_.

        Args:
            ax (matplotlib.axes, optional): axes for plotting on.
            offset (number or list, optional): defines a vertical offset. Default is 0.
            shift (number or list, optional): horizontal shift value. Default is 0.
            factor (number or list, optional): multiplicative factor on the y axis.
                Default is 1.
            calib (number or list, optional): multiplicative factor on the x axis.
                Default is 1.
            smooth (int, optional): number of points to average data (moving average).
                Default is 1.
            vertical_increment or vi (number): vertical increment between curves
                in terms of percentage of the maximum delta y.
            horizontal_increment or hi (number): horizontal increment between curves
                in terms of percentage of the x max range.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.


        Returns:
            `Line2D`_ list

        .. _matplotlib.pyplot.plot(): https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.pyplot.plot.html
        .. _Line2D: https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D
        """
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figmanip.figure()


        # percentage wise vertical increment
        if 'vi' in kwargs and 'vertical_increment' in kwargs:
            raise SyntaxError('keyword argument repeated: vertical increment and vi')
        elif 'vi' in kwargs or 'vertical_increment' in kwargs:
            if 'vi' in kwargs:
                vi = kwargs['vi']
                del kwargs['vi']
            if 'vertical_increment' in kwargs:
                vi = kwargs['vertical_increment']
                del kwargs['vertical_increment']
            temp = [0]*len(self)
            for i in range(len(self)):
                temp[i] = -(max(self.data[i].y) - min(self.data[i].y))
            vi = max(temp)*factor*vi/100
        else:
            vi = 0

        # percentage wise horizontal increment
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

        # smooth
        if isinstance(smooth, Iterable):
            if len(smooth) != len(self):
                raise ValueError(f'smooth must be a number of a list with length compatible with the number of spectra.\nnumber of smooth: {len(factor)}\nnumber of spectra: {len(self)}')
        else:
            smooth = [smooth]*len(self)

        # plot
        # temp = [0]*len(self)
        # for i in range(len(self)-1, -1, -1):
        #     temp[i] = self.data[i].plot(ax=ax, offset=offset[i], shift=shift[i], factor=factor[i], calib=calib[i], **kwargs)

        temp = [0]*len(self)
        for i in range(len(self)):
            temp[i] = self.data[i].plot(ax=ax, offset=offset[i], shift=shift[i], factor=factor[i], calib=calib[i], smooth=smooth[i], switch=switch, ranges=ranges, **kwargs)

        return temp#, offset, shift

    # experimental
    def sequential_plot(self, xlim=None, ylim=None, keep=None, update_function=None, **kwargs):
        """[EXPERIMENTAL] plot where you can use up and down keys to flip through spectra.

        Warning: 
            Only one plot can be open at a time. Also, some default 
            matplotlib keybidings are changed (should not be a problem)

        Args:
            xlim (tuple, optional): min and max ``x`` value. Default is None (full
                data range).
            ylim (tuple, optional): min and max ``y`` value. Default is None.
            keep (int, optional): index of the spectrum to keep on screen.
            update_function (function, optional): function that is called when
                the left or right arrows are pressed. update_function must be a 
                function of type:

                >>> def update_function(ss, __i):
                >>>     plt.title(__i)
                >>>     ss[__i].plot(color='black', marker='o')
            
                where ``__i`` is updated in every iteraction.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
        
        Returns:
            None
        """
        # vars
        self.__i = 0
        self.__xlim = xlim
        self.__ylim = ylim
        self.__ax   = None

        # change keybindings
        try:
            matplotlib.rcParams['keymap.back'].remove('left')
            matplotlib.rcParams['keymap.forward'].remove('right')
        except ValueError:
            pass
        
        # lims
        if self.__xlim is None:
            self.__xlim = [min(self[0].x), max(self[0].x)]
            for s in self:
                if min(s.x) < self.__xlim[0]:
                    self.__xlim[0] = min(s.x)
                if max(s.x) > self.__xlim[1]:
                    self.__xlim[1] = max(s.x)
        if self.__ylim is None:
            self.__ylim = [min(self[0].y), max(self[0].y)]
            for s in self:
                if min(s.y) < self.__ylim[0]:
                    self.__ylim[0] = min(s.y)
                if max(s.y) > self.__ylim[1]:
                    self.__ylim[1] = max(s.y)

        # keep
        if keep is not None:
            if isinstance(keep, str):
                assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
            else:
                assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

        # kwargs
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        # core update function
        if update_function is None:
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                plt.title(self.__i)
                if keep is not None:
                    if keep == 'next':
                        if self.__i+1 < len(self):
                            ss[self.__i+1].plot(color='red', alpha=0.5)
                    elif keep == 'previous':
                        if self.__i-1 >= 0:
                            ss[self.__i-1].plot(color='red', alpha=0.5)
                    else:
                        ss[keep].plot(color='red', alpha=0.5)
                ss[self.__i].plot(**kwargs)
                
                if self.__xlim is not None:
                    plt.xlim(self.__xlim)
                if self.__ylim is not None:
                    plt.ylim(self.__ylim)
        else:
            # add counter and xlim/ylim to update function
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                update_function(ss, __i=self.__i)
                
                if self.__xlim is not None:
                    plt.xlim(self.__xlim)
                if self.__ylim is not None:
                    plt.ylim(self.__ylim)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'right' or event.key == 'up':
                self.__i = self.__i + 1

                # self.__ax.cla()
                # for line in self.__ax.lines:
                #     line.remove()
                self.__ax.lines.clear()
                _update(ss)
                plt.draw()

            elif event.key == 'left' or event.key == 'down':
                self.__i = self.__i - 1

                # self.__ax.cla()
                # for line in self.__ax.lines:
                #     line.remove()
                self.__ax.lines.clear()
                _update(ss)
                plt.draw()

        # mouse events
        def mouse(event):
            """Mouse can be used with a keyboard key"""
            # print('mouse')
            # print(event.key)
            # print(event.button)
            pass

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
        fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))

        self.__ax.callbacks.connect('xlim_changed', on_xlims_change)
        self.__ax.callbacks.connect('ylim_changed', on_ylims_change)
        return

    # experimental
    def sequential_plot2(self, xlim=None, ylim=None, keep=None, update_function=None, **kwargs):
        """[EXPERIMENTAL] plot where you can use arrows to flip through spectra.

        Warning:
            Old version, zoom is reset every time. Only one plot can be open at a time.
            Also, some default matplotlib keybidings are changed (should not be a problem)

        Args:
            xlim (tuple, optional): min and max ``x`` value. Default is None (full
                data range).
            ylim (tuple, optional): min and max ``y`` value. Default is None.
            keep (int, optional): index of the spectrum to keep on screen.
            update_function (function, optional): function that is called when
                the left or right arrows are pressed. update_function must be a 
                function of type:

                >>> def update_function(ss, __i):
                >>>     plt.title(__i)
                >>>     ss[__i].plot(color='black', marker='o')
            
                where ``__i`` is updated in every iteraction.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
        
        Returns:
            None
        """
        self.__i = 0

        # change keybindings
        try:
            matplotlib.rcParams['keymap.back'].remove('left')
            matplotlib.rcParams['keymap.forward'].remove('right')
        except ValueError:
            pass

        # keep
        if keep is not None:
            if isinstance(keep, str):
                assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
            else:
                assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

        # kwargs
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        # core update function
        if update_function is None:
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                plt.title(self.__i)
                if keep is not None:
                    if keep == 'next':
                        if self.__i+1 < len(self):
                            ss[self.__i+1].plot(color='red', alpha=0.5)
                    elif keep == 'previous':
                        if self.__i-1 >= 0:
                            ss[self.__i-1].plot(color='red', alpha=0.5)
                    else:
                        ss[keep].plot(color='red', alpha=0.5)
                ss[self.__i].plot(**kwargs)
                
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)
        else:
            # add counter and xlim/ylim to update function
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                update_function(ss, __i=self.__i)
                
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'right' or event.key == 'up':
                self.__i = self.__i + 1

                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'left' or event.key == 'down':
                self.__i = self.__i - 1

                plt.cla()
                _update(ss)
                plt.draw()

        # mouse events
        def mouse(event):
            """Mouse can be used with a keyboard key"""
            # print('mouse')
            # print(event.key)
            # print(event.button)
            pass

        # plotting
        fig = figmanip.figure()

        # register callbacks
        fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
        fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))

        _update(self)
        return

    # experimental
    def shift_plot(self, xlim=None, ylim=None, step=None, vlines=None, keep=None, update_function=None, **kwargs):
        """[EXPERIMENTAL] flip through spectra (up/down keys), shift spectrum (a/d keys)

        Warning: 
            Only one plot can be open at a time. Also, some default 
            matplotlib keybidings are changed (should not be a problem)

        Args:
            xlim (tuple, optional): min and max ``x`` value. Default is None (full
                data range).
            ylim (tuple, optional): min and max ``y`` value. Default is None.
            step (number, optional): step. If None, 1% of the min datapoint 
                separation is used.
            vlines (list or number, optional): vertical dashed lines for 
                reference, default is None.
            keep (int, optional): index of the spectrum to keep on screen.
            update_function (function, optional): function that is called when
                the left or right arrows are pressed. update_function must be a 
                function of type:

                >>> def update_function(ss, __i):
                >>>     plt.title(__i)
                >>>     ss[__i].plot(shift=ss._calculated_shift[__i])
            
                where ``__i`` is updated in every iteraction.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
        
        Returns:
            None
        """
        self.__i = 0
        self.__xlim = xlim
        self.__ylim = ylim
        self.__ax   = None
        self._calculated_shift = np.array([0.0]*len(self))

        # set step
        if step is None: 
            if self.step is None:
                try:
                    self.check_step()
                    self.__step = self.step*0.5
                except:
                    self.__step = min([np.mean(np.diff(s.x)) for s in self])*0.5
            else:
                self.__step = self.step*0.5
        else:
            self.__step = step

        # change keybindings
        try:
            matplotlib.rcParams['keymap.back'].remove('left')
            matplotlib.rcParams['keymap.forward'].remove('right')
        except ValueError:
            pass

        # keep
        if keep is not None:
            if isinstance(keep, str):
                assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
            else:
                assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

        # kwargs
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        # core update function
        if update_function is None:
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                plt.title(f'{self.__i}: {self._calculated_shift[self.__i]}')
                if keep is not None:
                    if keep == 'next':
                        if self.__i+1 < len(self):
                            ss[self.__i+1].plot(shift=self._calculated_shift[self.__i+1], color='red', alpha=0.5)
                    elif keep == 'previous':
                        if self.__i-1 >= 0:
                            ss[self.__i-1].plot(shift=self._calculated_shift[self.__i-1], color='red', alpha=0.5)
                    else:
                        ss[keep].plot(shift=self._calculated_shift[keep], color='red', alpha=0.5)
                ss[self.__i].plot(shift=self._calculated_shift[self.__i], **kwargs)
                
                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)

                if self.__xlim is not None:
                    plt.xlim(self.__xlim)
                if self.__ylim is not None:
                    plt.ylim(self.__ylim)
        else:
            # add counter and xlim/ylim to update function
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                update_function(ss, __i=self.__i)

                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)
                
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'up':
                self.__i = self.__i + 1

                self.__ax.callbacks.disconnect(2)
                self.__ax.callbacks.disconnect(3)
                self.__newx = False
                self.__newy = False

                plt.cla()
                _update(ss)
                plt.draw()

                self.__ax.callbacks.connect('xlim_changed', on_xlims_change)
                self.__ax.callbacks.connect('ylim_changed', on_ylims_change)
            elif event.key == 'down':
                self.__i = self.__i - 1
                
                self.__ax.callbacks.disconnect(2)
                self.__ax.callbacks.disconnect(3)
                self.__newx = False
                self.__newy = False

                plt.cla()
                _update(ss)
                plt.draw()

                self.__ax.callbacks.connect('xlim_changed', on_xlims_change)
                self.__ax.callbacks.connect('ylim_changed', on_ylims_change)
            elif event.key == 'right':
                self._calculated_shift[self.__i] += self.__step
                
                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'left':
                self._calculated_shift[self.__i] -= self.__step
                
                plt.cla()
                _update(ss)
                plt.draw()

        # mouse events
        def mouse(event):
            """Mouse can be used with a keyboard key"""
            # print('mouse')
            # print(event.key)
            # print(event.button)
            pass

        # axis zoom changes
        def on_xlims_change(event_ax):
            if self.__newx:
                # print("updated xlims: ", event_ax.get_xlim())
                self.__xlim = event_ax.get_xlim()
            self.__newx = True

        def on_ylims_change(event_ax):
            if self.__newy:
                # print("updated ylims: ", event_ax.get_ylim())
                self.__ylim = event_ax.get_ylim()
            self.__newy = True

        # plotting
        fig, self.__ax = plt.subplots(1, 1)

        # register callbacks
        fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
        fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))
        
        _update(self)
        return
    
    # experimental
    def shift_plot2(self, xlim=None, ylim=None, step=None, vlines=None, keep=None, update_function=None, **kwargs):
        """[EXPERIMENTAL] flip through spectra (up/down keys), shift spectrum (left/right keys)

        Warning:
            Old version, zoom is reset every time. Only one plot can be open at a time.
            Also, some default matplotlib keybidings are changed (should not be a problem)

        Args:
            xlim (tuple, optional): min and max ``x`` value. Default is None (full
                data range).
            ylim (tuple, optional): min and max ``y`` value. Default is None.
            step (number, optional): step. If None, 1% of the min datapoint 
                separation is used.
            vlines (list or number, optional): vertical dashed lines for 
                reference, default is None.
            keep (int, optional): index of the spectrum to keep on screen.
            update_function (function, optional): function that is called when
                the left or right arrows are pressed. update_function must be a 
                function of type:

                >>> def update_function(ss, __i):
                >>>     plt.title(__i)
                >>>     ss[__i].plot(color='black', marker='o')
            
                where ``__i`` is updated in every iteraction.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
        
        Returns:
            None
        """
        # vars
        self.__i = 0
        self._calculated_shift = np.array([0.0]*len(self))

        # set step
        if step is None: 
            if self.step is None:
                try:
                    self.check_step()
                    self.__step = self.step*0.5
                except:
                    self.__step = min([np.mean(np.diff(s.x)) for s in self])*0.5
            else:
                self.__step = self.step*0.5
        else:
            self.__step = step

        # change keybindings
        try:
            matplotlib.rcParams['keymap.back'].remove('left')
            matplotlib.rcParams['keymap.forward'].remove('right')
        except ValueError:
            pass

        # keep
        if keep is not None:
            if isinstance(keep, str):
                assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
            else:
                assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

        # kwargs
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        # core update function
        if update_function is None:
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                plt.title(f'{self.__i}: {self._calculated_shift[self.__i]}')
                if keep is not None:
                    if keep == 'next':
                        if self.__i+1 < len(self):
                            ss[self.__i+1].plot(shift=self._calculated_shift[self.__i+1], color='red', alpha=0.5)
                    elif keep == 'previous':
                        if self.__i-1 >= 0:
                            ss[self.__i-1].plot(shift=self._calculated_shift[self.__i-1], color='red', alpha=0.5)
                    else:
                        ss[keep].plot(shift=self._calculated_shift[keep], color='red', alpha=0.5)
                ss[self.__i].plot(shift=self._calculated_shift[self.__i], **kwargs)
                
                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)

                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)
        else:
            # add counter and xlim/ylim to update function
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                update_function(ss, __i=self.__i)

                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)
                
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'up':
                self.__i = self.__i + 1

                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'down':
                self.__i = self.__i - 1
                
                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'right':
                self._calculated_shift[self.__i] += self.__step
                
                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'left':
                self._calculated_shift[self.__i] -= self.__step
                
                plt.cla()
                _update(ss)
                plt.draw()

        # mouse events
        def mouse(event):
            """Mouse can be used with a keyboard key"""
            # print('mouse')
            # print(event.key)
            # print(event.button)
            pass

        # plotting
        fig = figmanip.figure()
        fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
        fig.canvas.mpl_connect('button_press_event', lambda event: mouse(event))
        _update(self)
        return
    
    # experimental
    def roll_plot(self, xlim=None, ylim=None, vlines=None, keep=None, update_function=None, **kwargs):
        """[EXPERIMENTAL] flip through spectra (up/down keys), roll spectrum (a/d keys).

        Warning: 
            Only one plot can be open at a time. Also, some default 
            matplotlib keybidings are changed (should not be a problem)

        Args:
            xlim (tuple, optional): min and max ``x`` value. Default is None (full
                data range).
            ylim (tuple, optional): min and max ``y`` value. Default is None.
            vlines (list or number, optional): vertical dashed lines for 
                reference, default is None.
            keep (int or str, optional): index of the spectrum to keep on screen. 
                Use 'previous' or 'next' to show previous or next spectrum.
            update_function (function, optional): function that is called when
                the left or right arrows are pressed. update_function must be a 
                function of type:

                >>> def update_function(ss, __i):
                >>>     plt.title(__i)
                >>>     ss[__i].plot(color='black', marker='o')
            
                where ``__i`` is updated in every iteraction. If update_function
                is given, `keep` is disregarded.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
        
        Returns:
            None
        """
        # check step
        if self.check_step is None:
            try:
                self.check_step()
            except ValueError:
                raise ValueError('cannot apply roll because step is not the same. Use shift_plot().')

        # vars
        self.__i = 0
        self.__xlim = xlim
        self.__ylim = ylim
        self.__ax   = None
        self._calculated_roll = np.array([0]*len(self))

        # change keybindings
        try:
            matplotlib.rcParams['keymap.back'].remove('left')
            matplotlib.rcParams['keymap.forward'].remove('right')
        except ValueError:
            pass

        # lims
        if self.__xlim is None:
            self.__xlim = [min(self[0].x), max(self[0].x)]
            for s in self:
                if min(s.x) < self.__xlim[0]:
                    self.__xlim[0] = min(s.x)
                if max(s.x) > self.__xlim[1]:
                    self.__xlim[1] = max(s.x)
        if self.__ylim is None:
            self.__ylim = [min(self[0].y), max(self[0].y)]
            for s in self:
                if min(s.y) < self.__ylim[0]:
                    self.__ylim[0] = min(s.y)
                if max(s.y) > self.__ylim[1]:
                    self.__ylim[1] = max(s.y)

        # check if all data has well defined step
        for i, s in enumerate(self):
            try:
                s.check_step()
            except ValueError:
                raise ValueError(f'Spectrum number {i} has non-uniform x-coord.\nRoll can only be applied with uniform x-coord.\nMaybe used ss.interp() to fix that.')

        # keep
        if keep is not None:
            if isinstance(keep, str):
                assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
            else:
                assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

        # kwargs
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        # core update function
        if update_function is None:
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                plt.title(f'{self.__i}: {self._calculated_roll[self.__i]}')
                if keep is not None:
                    if keep == 'next':
                        if self.__i+1 < len(self):
                            ss[self.__i+1].plot(shift=self._calculated_roll[self.__i+1]*ss[self.__i+1].step, color='red', alpha=0.5)
                    elif keep == 'previous':
                        if self.__i-1 >= 0:
                            ss[self.__i-1].plot(shift=self._calculated_roll[self.__i-1]*ss[self.__i-1].step, color='red', alpha=0.5)
                    else:
                        ss[keep].plot(shift=self._calculated_roll[keep]*ss[keep].step, color='red', alpha=0.5)
                ss[self.__i].plot(shift=self._calculated_roll[self.__i]*ss[self.__i].step, **kwargs)
                
                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)

                if self.__xlim is not None:
                    plt.xlim(self.__xlim)
                if self.__ylim is not None:
                    plt.ylim(self.__ylim)
        else:
            # add counter and xlim/ylim to update function
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                update_function(ss, __i=self.__i)

                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)
                
                if self.__xlim is not None:
                    plt.xlim(self.__xlim)
                if self.__ylim is not None:
                    plt.ylim(self.__ylim)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'up':
                self.__i = self.__i + 1

                # self.__ax.cla()
                self.__ax.lines.clear()
                _update(ss)
                plt.draw()
            elif event.key == 'down':
                self.__i = self.__i - 1

                # self.__ax.cla()
                self.__ax.lines.clear()
                _update(ss)
                plt.draw()
            elif event.key == 'right':
                self._calculated_roll[self.__i] += 1

                # self.__ax.cla()
                self.__ax.lines.clear()
                _update(ss)
                plt.draw()
            elif event.key == 'left':
                self._calculated_roll[self.__i] -= 1

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

    # experimental
    def roll_plot2(self, xlim=None, ylim=None, vlines=None, keep=None, update_function=None, **kwargs):
        """[EXPERIMENTAL] flip through spectra (up/down keys), roll spectrum (left/right keys)

        Warning:
            Old version, zoom is reset every time. Only one plot can be open at a time.
            Also, some default matplotlib keybidings are changed (should not be a problem)

        Args:
            xlim (tuple, optional): min and max ``x`` value. Default is None (full
                data range).
            ylim (tuple, optional): min and max ``y`` value. Default is None.
            vlines (list or number, optional): vertical dashed lines for 
                reference, default is None.
            keep (int, optional): index of the spectrum to keep on screen.
            update_function (function, optional): function that is called when
                the left or right arrows are pressed. update_function must be a 
                function of type:

                >>> def update_function(ss, __i):
                >>>     plt.title(__i)
                >>>     ss[__i].plot(color='black', marker='o')
            
                where ``__i`` is updated in every iteraction.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
        
        Returns:
            None
        """
        # check step
        if self.check_step is None:
            try:
                self.check_step()
            except ValueError:
                raise ValueError('cannot apply roll because step is not the same. Use shift_plot().')

        # vars
        self.__i = 0
        self._calculated_roll = np.array([0]*len(self))

        # change keybindings
        try:
            matplotlib.rcParams['keymap.back'].remove('left')
            matplotlib.rcParams['keymap.forward'].remove('right')
        except ValueError:
            pass

        # check if all data has well defined step
        for i, s in enumerate(self):
            try:
                s.check_step()
            except ValueError:
                raise ValueError(f'Spectrum number {i} has non-uniform x-coord.\nRoll can only be applied with uniform x-coord.\nMaybe used ss.interp() to fix that.')

        # keep
        if keep is not None:
            if isinstance(keep, str):
                assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
            else:
                assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

        # kwargs
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        # core update function
        if update_function is None:
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                plt.title(f'{self.__i}: {self._calculated_roll[self.__i]}')
                if keep is not None:
                    if keep == 'next':
                        if self.__i+1 < len(self):
                            ss[self.__i+1].plot(shift=self._calculated_roll[self.__i+1]*ss[self.__i+1].step, color='red', alpha=0.5)
                    elif keep == 'previous':
                        if self.__i-1 >= 0:
                            ss[self.__i-1].plot(shift=self._calculated_roll[self.__i-1]*ss[self.__i-1].step, color='red', alpha=0.5)
                    else:
                        ss[keep].plot(shift=self._calculated_roll[keep]*ss[keep].step, color='red', alpha=0.5)
                ss[self.__i].plot(shift=self._calculated_roll[self.__i]*ss[self.__i].step, **kwargs)
                
                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)

                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)
        else:
            # add counter and xlim/ylim to update function
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                update_function(ss, __i=self.__i)

                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)
                
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'up':
                self.__i = self.__i + 1

                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'down':
                self.__i = self.__i - 1
                
                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'right':
                self._calculated_roll[self.__i] += 1
                
                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'left':
                self._calculated_roll[self.__i] -= 1
                
                plt.cla()
                _update(ss)
                plt.draw()


        # plotting
        fig = figmanip.figure()
        fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))

        _update(self)
        return

    # experimental
    def roll_plot3(self, xlim=None, ylim=None, vlines=None, keep=None, update_function=None, **kwargs):
        """[experimental] flip through spectra (up/down keys), roll spectrum (left/right keys).

        Warning:
            [Experimental] Use mouse to click and drag spectra to the left or right.
            Only one plot can be open at a time. Also, some default matplotlib 
            keybidings are changed (should not be a problem)

        Args:
            xlim (tuple, optional): min and max ``x`` value. Default is None (full
                data range).
            ylim (tuple, optional): min and max ``y`` value. Default is None.
            vlines (list or number, optional): vertical dashed lines for 
                reference, default is None.
            keep (int, optional): index of the spectrum to keep on screen.
            update_function (function, optional): function that is called when
                the left or right arrows are pressed. update_function must be a 
                function of type:

                >>> def update_function(ss, __i):
                >>>     plt.title(__i)
                >>>     ss[__i].plot(color='black', marker='o')
            
                where ``__i`` is updated in every iteraction.
            **kwargs: kwargs are passed to ``plt.plot()`` that plots the data.
        
        Returns:
            None
        """
        # check step
        if self.step is None:
            try:
                self.check_step()
            except ValueError:
                raise ValueError('cannot apply roll because step is not the same. Use shift_plot().')

        # raise NotImplementedError('This is not implemented yet.')
        self.__i = 0
        self._calculated_roll = np.array([0]*len(self))
        self._mouse_press   = 0

        # change keybindings
        try:
            matplotlib.rcParams['keymap.back'].remove('left')
            matplotlib.rcParams['keymap.forward'].remove('right')
        except ValueError:
            pass

        # check if all data has well defined step
        for i, s in enumerate(self):
            try:
                s.check_step()
            except ValueError:
                raise ValueError(f'Spectrum number {i} has non-uniform x-coord.\nRoll can only be applied with uniform x-coord.\nMaybe used ss.interp() to fix that.')
        
        # keep
        if keep is not None:
            if isinstance(keep, str):
                assert keep == 'previous' or keep == 'next', 'keep must be a valid spectrum index, or "previous/next"'
            else:
                assert abs(keep) < len(self), 'keep must be a valid spectrum index, or "previous/next"'

        # kwargs
        if 'color' not in kwargs:
            kwargs['color'] = 'black'
        if 'marker' not in kwargs:
            kwargs['marker'] = 'o'

        # core update function
        if update_function is None:
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                plt.title(f'{self.__i}: {self._calculated_roll[self.__i]}')
                if keep is not None:
                    if keep == 'next':
                        if self.__i+1 < len(self):
                            ss[self.__i+1].plot(shift=self._calculated_roll[self.__i+1]*ss[self.__i+1].step, color='red', alpha=0.5)
                    elif keep == 'previous':
                        if self.__i-1 >= 0:
                            ss[self.__i-1].plot(shift=self._calculated_roll[self.__i-1]*ss[self.__i-1].step, color='red', alpha=0.5)
                    else:
                        ss[keep].plot(shift=self._calculated_roll[keep]*ss[keep].step, color='red', alpha=0.5)
                ss[self.__i].plot(shift=self._calculated_roll[self.__i]*ss[self.__i].step, **kwargs)
                
                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)

                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)
        else:
            # add counter and xlim/ylim to update function
            def _update(ss):
                if self.__i >= len(self):
                    self.__i = len(self) - 1
                elif self.__i < 0:
                    self.__i = 0

                update_function(ss, __i=self.__i)

                if vlines is not None:
                    figmanip.vlines(vlines, color='red', ls='--', lw=1)
                
                if xlim is not None:
                    plt.xlim(xlim)
                if ylim is not None:
                    plt.ylim(ylim)

        # keyboard events
        def keyboard(event, ss):
            # print(event.key)
            # print('keyboard')
            # print(event.key)
            if event.key == 'up':
                self.__i = self.__i + 1

                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'down':
                self.__i = self.__i - 1
                
                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'right':
                self._calculated_roll[self.__i] += 1
                
                plt.cla()
                _update(ss)
                plt.draw()
            elif event.key == 'left':
                self._calculated_roll[self.__i] -= 1
                
                plt.cla()
                _update(ss)
                plt.draw()

        # mouse events
        def mouse_press(event):
            """Mouse can be used with a keyboard key"""
            if event.button is not None:
                if event.button == 1:
                    self._mouse_press = event.xdata

        def mouse_release(event, ss):
            """Mouse can be used with a keyboard key"""
            if event.button is not None:
                if event.button == 1:
                    self._calculated_roll[self.__i] += int(round((event.xdata - self._mouse_press)/self.step))   
                    plt.cla()
                    _update(ss)
                    plt.draw()

        # plotting
        fig = figmanip.figure()
        fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ss=self))
        fig.canvas.mpl_connect('button_press_event', lambda event: mouse_press(event))
        fig.canvas.mpl_connect('button_release_event', lambda event: mouse_release(event, ss=self))

        _update(self)
        return
    
    ###########
    # Special #
    ###########
    def copy_peaks_from_spectra(self):
        """Copy peaks from each spectrum.

        Returns:
            None
        """
        self.peaks._copy_from_spectra(self)

    def copy_peaks_to_spectra(self):
        """Copy peaks to each spectrum

        Returns:
            None
        """
        self.peaks._copy_to_spectra(self)

    def find_peaks(self, prominence=None, width=4, moving_average_window=8, ranges=None):
        """Find peaks recursively. Wrapper for `scipy.signal.find_peaks()`_.

        Args:
            prominence (number, optional): minimum prominence of peaks in percentage
                of the maximum prominence [max(y) - min(y)]. Default is 5.
            width (number, optional): minimum number of data points defining a peak.
            moving_average_window (int, optional): window size for smoothing the
                data for finding peaks. Default is 4.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            None

        .. _scipy.signal.find_peaks(): https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.find_peaks.html
        """
        for s in self:
            s.find_peaks(prominence=prominence, width=width, moving_average_window=moving_average_window, ranges=ranges)
        self.copy_peaks_from_spectra()

    def fit_peak(self, asymmetry=False, moving_average_window=1, method='least_squares', ranges=None, verbose=False):
        """Fit one peak for all spectra recursively. Wrapper for `lmfit.minimize()`_.

        Args:
            asymmetry (bool or dict, optional): if True, fits each half of the
                with a different width. Default is False
            moving_average_window (int, optional): window size for smoothing the
                data for finding the peak. Default is 1.
            method (str, optional): Name of the fitting method to use. See methods
                available on `lmfit.minimize()`_ documentation.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            None

        .. _lmfit.minimize(): https://lmfit.github.io/lmfit-py/fitting.html
        """
        for i, s in enumerate(self):
            if verbose:
                print(f'{i+1}/{len(self)}')
            s.fit_peak(asymmetry=asymmetry, moving_average_window=moving_average_window, method=method, ranges=ranges)
        self.copy_peaks_from_spectra()

        # THIS IS FOR FITTING ONE PEAK SIMULTANOUSLY
        # self.peaks.clear()

        # if ranges is None:
        #     temp = self
        # else:
        #     temp = self._extract(ranges=ranges)
        # xs = [s.x for s in temp]
        # ys = [s.y for s in temp]

        # for i2 in range(len(self)):
        #     # guess amp and c
        #     amp = max(ys[i2])
        #     c   = xs[i2][np.argmax(ys[i2])]

        #     # guess fwhm
        #     try:
        #         w1 = xs[i2][np.argmax(ys[i2])] - xs[i2][:np.argmax(ys[i2])][::-1][arraymanip.index(ys[i2][:np.argmax(ys[i2])][::-1], max(ys[i2])/2)]
        #     except ValueError:
        #         w1 = xs[i2][np.argmax(ys[i2]):][arraymanip.index(ys[i2][np.argmax(ys[i2]):], max(ys[i2])/2)] - xs[i2][np.argmax(ys[i2])]
        #     try:
        #         w2 = xs[i2][np.argmax(ys[i2]):][arraymanip.index(ys[i2][np.argmax(ys[i2]):], max(ys[i2])/2)] - xs[i2][np.argmax(ys[i2])]
        #     except ValueError:
        #         w2 = xs[i2][np.argmax(ys[i2])] - xs[i2][:np.argmax(ys[i2])][::-1][arraymanip.index(ys[i2][:np.argmax(ys[i2])][::-1], max(ys[i2])/2)]
        #     w = w1 + w2
        #     if w <= 0:
        #         w = 0.1*(max(xs[i2])-min(xs[i2]))
        #     if w <= 0:
        #         w = 1

        #     # peaks
        #     if asymmetry:
        #         self.peaks.append(i2=i2, amp=amp, c=c, w1=w1, w2=w2)
        #     else:
        #         self.peaks.append(i2=i2, amp=amp, c=c, w=w)

        # self.peaks.fit(xs=xs, ys=ys, method=method)

    def fit_peaks(self, method='least_squares', ranges=None):
        """Fit peaks for all spectra SIMULTANEOUSLY. Wrapper for `lmfit.minimize()`_.

        Args:
            method (str, optional): Name of the fitting method to use. See methods
                available on `lmfit.minimize()`_ documentation.
            ranges (list): a pair of values or a list of pairs. Each pair represents
                the start and stop of a data range from x. Use None to indicate
                the minimum or maximum x value of the data.

        Returns:
            None

        .. _lmfit.minimize(): https://lmfit.github.io/lmfit-py/fitting.html
        """
        if ranges is None:
            temp = self
        else:
            temp = self._extract(ranges=ranges)
        xs = [s.x for s in temp]
        ys = [s.y for s in temp]

        self.peaks.fit(xs=xs, ys=ys, method=method)


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
    _read_only = ['calculated_roll', 'calculated_shift']
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
        self._filepath  = ''

        self._calculated_roll  = None
        self._calculated_shift = None

        #######################
        # non-user attributes #
        #######################
        # vmin, vmax aren't user attributes. 
        # they are the default minimum and maximum value for plotting
        # vmin is set the the max of the intensity distribution (max of histogram)
        # vmax is set when intensity distribution (histogram) drops below 0.01% of the max.
        self._vmin = None
        self._vmax = None

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
                self.data = args[0]
                return

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

        ##############
        # vmin, vmax #
        ##############
        hist = self.calculate_histogram()
        x = hist.x[np.argmax(hist.y)+1:]
        y = hist.y[np.argmax(hist.y)+1:]
        filtered = np.array([[i, j] for i, j in zip(x, y) if j > abs(max(y)*0.0001)])  # clean zeros
        try:
            self._vmin = hist.x[np.argmax(hist.y)]
            self._vmax = filtered[-1, 0]
        except IndexError:  # in case the max of y is too high
            self._vmin = min([min(x) for x in self.data])
            self._vmax = max([max(x) for x in self.data])
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
        return (self.data.shape[0], self.data.shape[1])
    @shape.setter
    def shape(self, value):
        raise AttributeError('Attribute is "read only". Cannot set attribute.')
    @shape.deleter
    def shape(self):
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
        # for attr in self._get_user_attrs():
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
        # for attr in self._get_user_attrs():
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
        # for attr in self._get_user_attrs():
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
        # for attr in self._get_user_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     final.__setattr__(attr, value)

        return final
    
    def __truediv__(self, object):
        return self.__div__(object)

    ###########
    # support #
    ###########
    def _default_attrs(self):
        """return list with default attrs."""
        return ['_data', '_shape', '_calculated_shift', '_calculated_roll',
                           '_x_centers', '_y_centers']
    
    def _get_user_attrs(self):
        """return attrs that are user defined."""
        return [key for key in self.__dict__.keys() if key not in self._default_attrs() and key.startswith('_') == False]

    #################
    # save and load #
    #################
    def _create_header(self, verbose=False):
        """Gather attrs to be saved to a file."""
        header = ''
        attrs = self._get_user_attrs()
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
            try: 
                filepath = Path(self.filepath)
            except AttributeError:
                raise TypeError("Missing 1 required argument: 'filepath'")
            
        # check if filepath points to a file
        assert filepath.parent.exists(), f'filepath folder does not exists.\nfolderpath: {filepath.parent}'
        if filepath.exists():
            assert filepath.is_file(), 'filepath must point to a file'

        # check overwrite
        if check_overwrite:
            if filepath.exists() == True:
                if filepath.is_file() == True:
                    if interact.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
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
    def transpose(self):
        """Transpose image

        Returns:
            None
        """
        self._data = self.data.transpose()

        # fixing attrs
        temp = copy.deepcopy(self.c_centers)
        self._x_centers = self.y_centers
        self._y_centers = temp        
        self._calculated_roll = None

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
        
        ##############
        # vmin, vmax #
        ##############
        hist = self.calculate_histogram()
        x = hist.x[np.argmax(hist.y)+1:]
        y = hist.y[np.argmax(hist.y)+1:]
        filtered = np.array([[i, j] for i, j in zip(x, y) if j > abs(max(y)*0.0001)])  # clean zeros
        try:
            self._vmin = hist.x[np.argmax(hist.y)]
            self._vmax = filtered[-1, 0]
        except IndexError:  # in case the max of y is too high
            self._vmin = min([min(x) for x in self.data])
            self._vmax = max([max(x) for x in self.data])

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

    def crop(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Crop Image.

        Warning:
            In this version, crop overwrites the data in this object. For croping
            to a new object use im.copy(). This might change in the future.

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
        self._data      = self.data[int(y_start):int(y_stop), int(x_start):int(x_stop)]
        self._x_centers = self.x_centers[x_start:x_stop]
        self._y_centers = self.y_centers[y_start:y_stop]

        # save to a new object (old version)
        # im = Image(data=self.data[int(y_start):int(y_stop), int(x_start):int(x_stop)])
        # im.x_centers = self.x_centers[x_start:x_stop]
        # im.y_centers = self.y_centers[y_start:y_stop]

        # # transfer attrs
        # for attr in self._get_user_attrs():
        #     value = copy.deepcopy(self.__dict__[attr])
        #     im.__setattr__(attr, value)

        # return im
        return
    
    ##############
    # extractors #
    ##############
    def copy(self, *args, **kwargs):
        """Return a copy of the object.

        Usage:
            >>> im2.copy(im1)     # im2 is now a copy of im1
            >>> im2 = im1.copy()  # im2 is now a copy of im1
            >>>
            >>> # im3 will be a croped image of im1
            >>> im3 = im1.copy(x_start, x_stop, y_start, y_stop) 

        Args:
            im (Image, optional): Image to be is copied. See usage.
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
                        'im2.copy(im1)     # im2 is now a copy of im1\n'  +\
                        'im2 = im1.copy()  # im2 is now a copy of im1\n'  +\
                        'im3 = im1.copy(x_start, x_stop, y_start, y_stop)\n'
        #############################
        # check if input is Spectra #
        #############################
        im      = None
        x_start = False
        if 'im' in kwargs:
            im = kwargs['im']
        elif len(args) == 1:
            if isinstance(args[0], Spectra):
                im = args[0]
        ###################################
        # check if crop ranges are passed #
        ###################################
        elif 'x_start' in kwargs or 'x_stop' in kwargs or 'y_start' in kwargs or 'y_stop' in kwargs:
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
        elif len(args) != 4 and len(args) != 0 and len(args) != 1:
            raise ValueError(error_message)
        elif len(kwargs) != 0:
            for x in list(kwargs.keys()):
                if x not in ['im', 'x_start', 'x_stop', 'y_start', 'y_stop']:
                    raise ValueError(f'{x} is not a recognized input for copy function\n'+error_message)

        ##################################
        # if Image is passed, copy Image #
        ##################################
        if im is not None:
            if isinstance(im, Image):
                self._data      = im.data
                self._x_centers = im.x_centers
                self._y_centers = im.y_centers
                self._filepath  = im.filepath
                self._calculated_roll = im.calculated_roll

                # user defined attrs
                for attr in self._get_user_attrs():
                    self.__delattr__(attr)
                for attr in im._get_user_attrs():
                    self.__setattr__(attr, im.__dict__[attr])
            else:
                raise TypeError('Only type br.Spectra can be copied to type br.Spectra')
            return
        ##################
        # identical copy #
        ##################
        elif x_start == False:
            im = Image(data=self.data)
            im._x_centers = self.x_centers
            im._y_centers = self.y_centers
            # im._filepath  = self.filepath

            # transfer attrs
            for attr in self._get_user_attrs():
                value = copy.deepcopy(self.__dict__[attr])
                im.__setattr__(attr, value)

            return im
        #############################
        # if crop ranges are passed #
        #############################
        else:
            im = self.copy()
            im.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
            # im._filepath = self.filepath

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

        # transfer attrs
        for attr in self._get_user_attrs():
            value = copy.deepcopy(self.__dict__[attr])
            reduced.__setattr__(attr, value)

        return reduced
    
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
        final = Spectrum(x=x, y=hist)

        # transfer attrs
        default_attrs = final._default_attrs()
        for attr in self._get_user_attrs():
            if attr not in default_attrs:
                value = copy.deepcopy(self.__dict__[attr])
                final.__setattr__(attr, value)

        return final

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
            final = Spectrum(x=self.x_centers, y=np.sum(self._data, axis=0))
        elif axis == 1:
            final = Spectrum(x=self.y_centers, y=np.sum(self._data, axis=1))
        
        # transfer attrs
        default_attrs = final._default_attrs()
        for attr in self._get_user_attrs():
            if attr not in default_attrs:
                value = copy.deepcopy(self.__dict__[attr])
                final.__setattr__(attr, value)

        return final

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
            None
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
        
        # figure
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figmanip.figure()

        # kwargs
        if 'cmap' not in kwargs:
            kwargs['cmap'] = 'jet'
        if 'vmin' not in kwargs:
            kwargs['vmin'] = self._vmin
        if 'vmax' not in kwargs:
            kwargs['vmax'] = self._vmax

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
        
        # figure
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figmanip.figure()
        
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
        if 'vmin' not in kwargs:
            kwargs['vmin'] = self._vmin
        if 'vmax' not in kwargs:
            kwargs['vmax'] = self._vmax

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
        if 'vmin' not in kwargs:
            kwargs['vmin'] = self._vmin
        if 'vmax' not in kwargs:
            kwargs['vmax'] = self._vmax

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

    _read_only = ['calculated_shift']
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

    def __getitem__(self, item):
        if isinstance(item, int):
            return self.x[item], self.y[item]
        elif isinstance(item, slice):
            x = self.x[item]
            y = self.y[item]

            pe = PhotonEvents(x=x, y=y)

            # transfer attrs
            for attr in self._get_user_attrs():
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
            
    ###########
    # support #
    ###########    
    def _default_attrs(self):
        """return list with default attrs."""
        return ['_x', '_y', '_filepath', '_calculated_shift']
    
    def _get_user_attrs(self):
        """return attrs that are user defined."""
        return [key for key in self.__dict__.keys() if key not in self._default_attrs() and key.startswith('_') == False]
    
    #################
    # save and load #
    #################
    def _create_header(self, verbose=False):
        """Gather attrs to be saved to a file."""
        header = ''
        attrs = self._get_user_attrs()
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
                    if interact.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
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
                f(y_centers) for axis=1. Default is None.
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
    def transpose(self):
        """Switch x and y positions.

        Returns:
            None
        """
        x = copy.deepcopy(self.x)
        self._x = copy.deepcopy(self.y)
        self._y = x

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

    def crop(self, x_start=None, x_stop=None, y_start=None, y_stop=None):
        """Crop photon events out.

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
        self.x = list(temp[:, 0])
        self.y = list(temp[:, 1])

        return 

    def clip(self, mask):
        """Clip photon events.

        Usage:
            >>> pe.clip((0, 12, 3, 12))
            >>> pe.clip([(0, 12, 3, 12), (1, 4, 15, 18)])

        Args:
            mask (list): list with ractangles coordinates (x_start, x_stop, y_start, y_stop).

        Returns:
            None
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
        self._x = x
        self._y = y
        return

    ##############
    # extractors #
    ##############
    def _extract(self, *args, **kwargs):
        """Same as copy(), but attributes are not copied to the new object."""
        
        # check if extract is really necessary
        if kwargs == {} and args == ():
            return copy.deepcopy(self)
        else:
            ranges = self._validate_ranges(*args, **kwargs)
        if (ranges[0][0] <= min([min(s.x) for s in self]) and ranges[0][1] >= max([max(s.x) for s in self]) ) and len(ranges) == 1:
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
        ranges = self._validate_ranges(*args, **kwargs)
        ss = Spectra(n=len(self))
        for i, s in enumerate(self.data):
            try:
                ss[i] = s.copy(ranges=ranges)
            except RuntimeError:
                raise RuntimeError(f'It seems like spectrum number {i} has no data points within range: {ranges}.\nPlease, fix ranges (or delete spectrum) so all spectra have at least one data point within range.')
        return ss
    
    def copy(self, *args, **kwargs):
        """Return a copy of the object.

        Usage:
            >>> pe2.copy(pe1)     # pe2 is now a copy of ss1
            >>> pe2 = pe1.copy()  # pe2 is now a copy of ss1
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
                        'pe2.copy(pe1)     # pe2 is now a copy of pe1\n'  +\
                        'pe2 = pe1.copy()  # pe2 is now a copy of pe1\n'  +\
                        'pe3 = pe1.copy(x_start, x_stop, y_start, y_stop)\n' +\
                        'pe3 = pe1.copy([(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])'
        #############################
        # check if input is Spectra #
        #############################
        pe   = None
        mask = None
        if 'pe' in kwargs:
            pe = kwargs['pe']
        elif len(args) == 1:
            if isinstance(args[0], PhotonEvents):
                pe = args[0]
            elif isinstance(args[0], Iterable):
                mask = args[0]
        ############################
        # check if mask are passed #
        ############################
        elif 'mask' in kwargs:
            mask = kwargs['mask']
        elif len(args) == 4:
            mask = [[args[0], args[1], args[2], args[3]], ]
        elif len(args) != 4 and len(args) != 0 and len(args) != 1:
            raise ValueError(error_message)
        elif len(kwargs) != 0:
            for x in list(kwargs.keys()):
                if x not in ['im', 'mask']:
                    raise ValueError(f'{x} is not a recognized input for copy function\n'+error_message)

        ##################################
        # if Image is passed, copy Image #
        ##################################
        if pe is not None:
            if isinstance(pe, PhotonEvents):
                self._x        = pe.x
                self._y        = pe.y
                self._xlim     = pe.xlim
                self._ylim     = pe.ylim
                self._filepath = pe.filepath

                # user defined attrs
                for attr in self._get_user_attrs():
                    self.__delattr__(attr)
                for attr in pe._get_user_attrs():
                    self.__setattr__(attr, pe.__dict__[attr])
            else:
                raise TypeError('Only type br.PhotonEvents can be copied to type br.PhotonEvents')
            return
        elif mask is not None:
            pe = self.copy()
            pe.clip(mask=mask)
            return pe
        ####################################
        # Otherwise, return a copy of self #
        ####################################
        else:
            pe = PhotonEvents(x=self.x, y=self.y)

            pe._xlim     = self.xlim
            pe._ylim     = self.ylim
            pe._filepath = self.filepath
            
            # transfer attrs
            for attr in self._get_user_attrs():
                value = copy.deepcopy(self.__dict__[attr])
                pe.__setattr__(attr, value)

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
        final = Image(temp.transpose())
        final._x_centers = arraymanip.moving_average(_x_edges, n=2)
        final._y_centers = arraymanip.moving_average(_y_edges, n=2)

        # transfer attrs
        default_attrs = final._default_attrs()
        for attr in self._get_user_attrs():
            if attr not in default_attrs:
                value = copy.deepcopy(self.__dict__[attr])
                final.__setattr__(attr, value)

        return final

    ########################
    # calculation and info #
    ########################
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
        final = im.calculate_spectrum(axis=axis)

        # transfer attrs
        default_attrs = final._default_attrs()
        for attr in self._get_user_attrs():
            if attr not in default_attrs:
                value = copy.deepcopy(self.__dict__[attr])
                final.__setattr__(attr, value)

        return final

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
        # figure
        if ax is None:
            ax = plt
            if settings.FIGURE_FORCE_NEW_WINDOW:
                figmanip.figure()

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
