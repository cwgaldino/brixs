from collections.abc import Iterable
from pathlib import Path

class _settings():

    def __init__(self):
        self._METADATA = {'rixs': 0, 'xas': 0, 'linescans': 0}
        self._FOLDERPATH = ''
        self._PREFIX = 'i21-'
        self._CALIB = 1
        self._SLOPE  = None
        self._MAX_IMAGES_TO_LOAD_AT_ONCE = 200

    @property
    def METADATA(self):
        return self._METADATA
    @METADATA.setter
    def METADATA(self, value):
        _text = 'Cannot modify object\n' +\
                'Please, use `del i21.settings.METADATA[<type>][<name>]` to delete metadata\n' +\
                'Please, use `i21.settings.METADATA[<type>][<name>] = [<category>, <address>]` to delete metadata\n' +\
                'Available types are: `rixs`, `xas`, `linescan`'
        raise AttributeError(_text)
    @METADATA.deleter
    def METADATA(self):
        raise AttributeError('Cannot delete object.')

    @property
    def FOLDERPATH(self):
        return self._FOLDERPATH
    @FOLDERPATH.setter
    def FOLDERPATH(self, value):
        value = Path(value)
        assert value.exists(), 'Cannot find folderpath'
        assert value.is_dir(), 'folderpath does not point to a folder'
        self._FOLDERPATH = value
    @FOLDERPATH.deleter
    def FOLDERPATH(self):
        raise AttributeError('Cannot delete object.')

    @property
    def PREFIX(self):
        return self._PREFIX
    @PREFIX.setter
    def PREFIX(self, value):
        assert isinstance(value, Iterable), 'PREFIX must be a string'
        self._PREFIX = value
    @PREFIX.deleter
    def PREFIX(self):
        raise AttributeError('Cannot delete object.')

    @property
    def CALIB(self):
        return self._CALIB
    @CALIB.setter
    def CALIB(self, value):
        _text = 'Invalid calibration value\n' +\
                'CALIB must be a number between -100 and 100'
        assert value<100 and value>-100, _text
        self._CALIB = value
    @CALIB.deleter
    def CALIB(self):
        raise AttributeError('Cannot delete object.')

    @property
    def SLOPE(self):
        return self._SLOPE
    @SLOPE.setter
    def SLOPE(self, value):
        self._SLOPE = value
    @SLOPE.deleter
    def SLOPE(self):
        raise AttributeError('Cannot delete object.')

    @property
    def MAX_IMAGES_TO_LOAD_AT_ONCE(self):
        return self._MAX_IMAGES_TO_LOAD_AT_ONCE
    @MAX_IMAGES_TO_LOAD_AT_ONCE.setter
    def MAX_IMAGES_TO_LOAD_AT_ONCE(self, value):
        self._MAX_IMAGES_TO_LOAD_AT_ONCE = value
    @MAX_IMAGES_TO_LOAD_AT_ONCE.deleter
    def MAX_IMAGES_TO_LOAD_AT_ONCE(self):
        raise AttributeError('Cannot delete object.')
    
    def __str__(self):
        _text  = ''
        _text += f"METADATA:                   Use i21.settings.METADATA[type] with type='rixs', 'xas', 'linescans' to see all defined metadata\n" +\
                 f'FOLDERPATH:                 {self.FOLDERPATH}\n' +\
                 f'PREFIX:                     {self.PREFIX}\n' +\
                 f'CALIB:                      {self.CALIB}\n' +\
                 f'SLOPE:                      {self.SLOPE}\n' +\
                 f'MAX_IMAGES_TO_LOAD_AT_ONCE: {self.MAX_IMAGES_TO_LOAD_AT_ONCE}\n'
        return _text

    def __repr__(self):
        return self.__str__()

    def help(self):
        print(self._help)
    
    def pretty_print(self):
        print(self.__str__())
