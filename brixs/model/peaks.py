#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Peak objects"""

# standard libraries
import copy
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# specific libraries
from collections.abc import Iterable, MutableMapping
from scipy.signal import find_peaks

# BRIXS
import brixs as br

lmfitok = False
try:
    import lmfit
    lmfitok = True
except:
    pass    

# Peaks ========================================================================
if lmfitok:
    class Peaks(lmfit.Parameters):
        """Creates a ``Peaks`` object.

        Peak attrs will have two number tags, i1 and i2.

        i1 = peak number
        i2 = spectrum number

        e.g., amp_1_4 = amplitude of peak 1 in the spectrum 4.
        i2 is absent if peaks are pointing to only one spectrum.

        How to create a peaks object:
            p = br.Peaks()

            p = br.Peaks(filepath=<filepath>)
            p = br.Peaks(<filepath>)
        
        The filepath must point to a file with a JSON string.

        Attributes:
            additional (str or function): additional function f(x) to be included.
            filepath (str or pathlib.Path): filepath associated with data.

        Computed (read-only) attributes:
            spectrum (br.Spectrum): spectrum
            spectra (br.Spectra): spectra

        Write-only attributes:
            None
        """
        def __init__(self, *args, **kwargs):
            """Initialize the object instance.

            Args:
                filepath (str or pathlib.Path, optional): filepath to read.
            
            Raises:
                AttributeError: if kwargs and args cannot be read.
            """
            #######################
            # Initializing object #
            #######################
            super().__init__()
            
            ###########################
            # Initializing attributes #
            ###########################
            self._additional = ''
            self._filepath   = ''

            ###################################
            # asserting validity of the input #
            ###################################
            error_message = 'Wrong input. Peaks object cannot be created. Please, use one ' +\
                            'of the examples below to create a Peaks object:\n' +\
                            '\n' +\
                            'p = br.Peaks()\n' +\
                            '\n' +\
                            'p = br.Peaks(filepath=<filepath>)\n' +\
                            'p = br.Peaks(<filepath>)\n' +\
                            '\n' +\
                            'filepath must be a string or pathlib.Path object'
            if kwargs != {} and args != ():
                raise AttributeError(error_message)
            if any([item not in ['filepath'] for item in kwargs.keys()]):
                raise AttributeError(error_message)
            if len(args) > 1 or len(kwargs) > 1:
                raise AttributeError(error_message)

            ################
            # loading data #
            ################
            # keyword arguments
            if 'filepath' in kwargs:
                self.load(kwargs['filepath'])
                return

            # positional arguments
            if len(args) == 1:
                if isinstance(args[0], str) or isinstance(args[0], Path):
                    self.load(args[0])
                    return

        ##############
        # attributes #
        ##############
        @property
        def additional(self):
            return self._additional
        @additional.setter
        def additional(self, value):
            if value is None:
                value = ''
            elif isinstance(value, str) or isinstance(value, function):
                self._additional = value
            else:
                raise TypeError(r'Invalid type ' + str(type(value)) + 'for additional\nAdditional can only be str or function type.')  
        @additional.deleter
        def additional(self):
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

        #################
        # magic methods #
        #################
        def __setitem__(self, name, value):
            super().__setitem__(name, value)

        def __getitem__(self, name):

            # has_i2 (so one doesn't have to call it every time)
            has_i2 = self._has_i2()
            
            #######
            # int #
            #######
            if isinstance(name, int):
                if has_i2:
                    return self._get_peaks_by_index(i1=None, i2=name)
                else:
                    return self._get_peaks_by_index(i1=name, i2=None) 
            ##############
            # additional #
            ##############
            elif name == 'A':
                if has_i2:
                    raise ValueError('Peaks have multiple spectrum indexes i2.\nIndex i2 must be a integer\n"A" is not valid')
                return self._get_peaks_by_index(i1=name, i2=None) 
            ##########
            # string #
            ##########
            elif isinstance(name, str):
                # if name is written correctly, return it
                if name in self:
                    return super().__getitem__(name)

                # if `Peaks` have i2, the user has to give the full attr path with i1 and i2
                if has_i2 == True:
                    assert '_' in name, 'i1 and i2 must be defined, which peak do you want?'
                    assert len(name.split('_')) == 3, 'i1 and i2 must be defined, which peak do you want?'

                    i1 = int(name.split('_')[1])
                    i2 = int(name.split('_')[2])

                    if name.startswith('area') and self._use_area(i1=i1, i2=i2) == False:
                        # return self.get_area(i1=i1, i2=i2)
                        final        = lmfit.Parameter(name=name + ' (Pseudo)', value=self.get_area(i1=i1, i2=i2), vary=False)
                        final.stderr = self.get_area_error(i1=i1, i2=i2)
                        return final
                    elif name.startswith('amp') and self._use_area(i1=i1, i2=i2) == True:
                        # return self.get_amp(i1=i1, i2=i2)
                        final        = lmfit.Parameter(name=name + ' (Pseudo)', value=self.get_amp(i1=i1, i2=i2), vary=False)
                        final.stderr = self.get_amp_error(i1=i1, i2=i2)
                        return final
                    elif name.startswith('w') and self.is_asymmetric(i1=i1, i2=i2) == True:
                        # return self.get_w(i1=i1, i2=i2)
                        final        = lmfit.Parameter(name=name + ' (Pseudo)', value=self.get_w(i1=i1, i2=i2), vary=False)
                        final.stderr = self.get_w_error(i1=i1, i2=i2)
                        return final
                    # else:
                    #     return super().__getitem__(name)
                # if `Peaks` does not have i2 and only one peak, then the user just needs the attr name
                else:
                    if self.number_of_peaks() == 1:
                        if name in ['amp', 'area', 'c', 'w', 'w1', 'w2', 'm', 'm1', 'm2']:
                            i1   = self._get_indexes_i1()[0]
                            name = name + f'_{i1}'
                    assert '_' in name, 'i1 must be defined, which peak do you want?'     
                    i1 = int(name.split('_')[1])
                    if name in self:
                        return super().__getitem__(name)
                    elif name.startswith('area'):
                        # return self.get_area(i1=i1)
                        final        = lmfit.Parameter(name=name + ' (Pseudo)', value=self.get_area(i1=i1), vary=False)
                        final.stderr = self.get_area_error(i1=i1, i2=None)
                        return final
                    elif name.startswith('amp'):
                        # return self.get_amp(i1=i1)
                        final        = lmfit.Parameter(name=name + ' (Pseudo)', value=self.get_amp(i1=i1), vary=False)
                        final.stderr = self.get_amp_error(i1=i1, i2=None)
                        return final
                    elif name.startswith('w'):
                        # return self.get_w(i1=i1)
                        final        = lmfit.Parameter(name=name + ' (Pseudo)', value=self.get_w(i1=i1), vary=False)
                        final.stderr = self.get_w_error(i1=i1, i2=None)
                        return final
                    # else:
                    #     raise ValueError(f'Cannot find parameter: {name}')  
                
                # if the function reaches this point, it mean one could not find name
                # return super().__getitem__(name)
                raise ValueError(f'`{name}` is not a recognizable attr')

        def __delitem__(self, name):
            """TODO I want to avoid deleting single parameters from a peak.
            """
            super().__delitem__(name)
        
        ###########
        # support #
        ###########
        def is_asymmetric(self, i1=None, i2=None):
            """Return true if peak has w1 and w2 objects.
            
            Args:
                i1 (int): primary index i1. If i1='A', raises an error.
                i2 (int, optional): secondary (spectrum) index i2. Only used if 
                    peaks has secondary indexes. Default is None.
            
            Returns:
                bool
            """
            #######################
            # check primary index #
            #######################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            if i1 is None:
                    if self.number_of_peaks(i2=i2) == 1:
                        i1 = self._get_indexes_i1()[0]
                    else:
                        raise ValueError('Peaks have multiple peaks and i1 must be defined')
            assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'


            #################
            # is asymmetric #
            #################
            if self._has_i2():         
                if f'w1_{i1}_{i2}' in self and f'w2_{i1}_{i2}' in self:
                    return True
                elif f'w_{i1}_{i2}' in self:
                    return False
            else:
                if f'w1_{i1}' in self and f'w2_{i1}' in self:
                    return True
                elif f'w_{i1}' in self:
                    return False
        
        ##################
        # support hidden #
        ##################
        def _get_user_attrs(self):
            """return attrs that are user defined."""
            default_attrs =  ['_asteval', '_additional']
            return [key for key in self.__dict__.keys() if key not in default_attrs and key.startswith('_') == False]
        
        def _has_i2(self):
            """return True if parameters have secondary (spectrum) index i2.

            Raises:
                ValueError: if no parameters are found, or parameters do not have 
                    indexes, or have too many indexes.
            
            Returns:
                bool
            """
            if len(self) == 0:
                return None
            else:
                test = max([len(name.split('_')) for name in list(self.keys())])
                if test == 3:
                    return True
                elif test == 2:
                    return False
                else:
                    raise ValueError('some parameters do not have or have too many indexes')
        
        def _has_additional(self):
            """returns True if additional function is defined.

            Returns:
                bool
            """
            return self.additional != ''

        def _get_indexes_i1(self, i2=None):
            """return list of all i1's.

            Includes additional parameters.
            
            Args:
                i2 (int, optional): spectrum index. Only used if peaks has double 
                    indexes.
            
            Returns:
                list
            """
            #########################
            # check secondary index #
            #########################
            error_message1 = 'Peaks() has parameters from multiple spectra.' +\
                            'Therefore, you must indicate the spectrum index i2.'
            error_message2 = 'Peaks() does not have parameters from multiple spectra.' +\
                            f'Therefore, i2 must be set to None and not {i2}'
            has_i2 = self._has_i2()
            if has_i2 is None:
                return []
            elif has_i2:
                if i2 is None:
                    raise ValueError(error_message1)
                elif br.is_integer(i2):
                    if i2 not in self._get_indexes_i2():
                        raise ValueError(f'i2={i2} cannot be found\nvalid i2: {self._get_indexes_i2()}')
                else:
                    raise TypeError('i2 must be an integer')
            else:
                if i2 is not None:
                    raise ValueError(error_message2)
            
            # get indexes
            final = []
            for name in self:
                split = name.split('_')
                _i1 = int(split[1])
                if has_i2:
                    _i2 = int(split[2])
                    if _i2 == i2:
                        if _i1 not in final:
                            final.append(_i1) 
                elif _i1 not in final:
                    final.append(_i1) 
            return final

        def _get_indexes_i2(self):
            """return list of all secondary (spectrum) indexes i2.
            
            Returns:
                list
            """
            #########################
            # check secondary index #
            #########################
            error_message = 'Peaks() only have parameters from a single spectrum.' +\
                            'Therefore, there is no need of secondary indexes i2.'
            has_i2 = self._has_i2()
            if has_i2 is None:
                return []
            if has_i2 == False:
                raise ValueError(error_message)
            
            # get indexes
            final = []
            for name in self:
                _i2 = int(name.split('_')[2])
                if _i2 not in final:
                    final.append(_i2) 
            return final

        def _get_names_by_index(self, i1=None, i2=None):
            """return list of names with index i1 and i2.

            Args:
                i1 (int): primary index i1. Use i1='A' for additional parameters.
                i2 (int, optional): secondary (spectrum) index i2. Only used if 
                    peaks has secondary indexes. Default is None.
            
            Returns:
                list
            """           
            #########################
            # check secondary index #
            #########################
            error_message = 'Peaks() does not have parameters from multiple spectra.' +\
                            f'Therefore, i2 must be set to None and not {i2}'
            assert i2 == None or br.numanip.is_integer(i2), 'i2 must be an integer or None'
            has_i2 = self._has_i2()
            if has_i2 is None:
                return []
            elif has_i2 == True and br.numanip.is_integer(i2):
                if i2 not in self._get_indexes_i2():
                    raise ValueError(f'i2={i2} cannot be found\nvalid i2: {self._get_indexes_i2()}')
            elif has_i2 == False and i2 is not None:
                    raise ValueError(error_message)

            #######################
            # check primary index #
            #######################
            # if has_i2 == True and i1 is not None:
            if i1 is not None:
                if i1 < 0:
                    i1 = np.sort(self._get_indexes_i1(i2=i2))[i1]
                assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'

            #############
            # get names #
            #############
            final = []
            if has_i2 == True and i2 is not None:
                if i1 == 'A':
                    for name in self:
                        split = name.split('_')
                        if split[1] == 'A' and int(split[2]) == i2:
                            final.append(name)
                elif i1 == None:
                    for name in self:
                        split = name.split('_')
                        if int(split[2]) == i2:
                            final.append(name)
                else:
                    for name in self:
                        split = name.split('_')
                        if int(split[1]) == i1 and int(split[2]) == i2:
                            final.append(name)
            elif has_i2 == True and i2 is None:
                if i1 == 'A':
                    for name in self:
                        split = name.split('_')
                        if split[1] == 'A':
                            final.append(name)
                elif i1 == None:
                    return list(self.keys())
                else:
                    for name in self:
                        split = name.split('_')
                        if int(split[1]) == i1:
                            final.append(name)
            elif has_i2 == False:
                if i1 == 'A':
                    for name in self:
                        if name.split('_')[1] == 'A':
                            final.append(name)
                elif i1 is None:
                    for name in self:
                        final.append(name)
                else:
                    for name in self:
                        if int(name.split('_')[1]) == i1:
                            final.append(name)
            return final
        
        def _get_peaks_by_index(self, i1=None, i2=None):
            """return Peaks object with only selected peak.

            The returned peaks object is "connect" to the parent Peaks object, i.e.
                they point to the same object. Modifications done in the returned
                peaks are reflected in the parent peaks. (Developers note, not sure
                this will work well long term. This might change eventually).

            Args:
                i1 (int): primary index i1. Use i1='A' for additional parameters.
                i2 (int, optional): secondary (spectrum) index i2. Only used if 
                    peaks has secondary indexes. Default is None.
            
            Returns:
                Peaks
            """ 
            names = self._get_names_by_index(i1=i1, i2=i2)

            if names == []:
                return Peaks()
                # raise ValueError('Empty peaks object')
            
            peaks = Peaks()
            for name in names:
                peaks[name] = self[name]
            
            try:
                peaks.check_i2()
            except ValueError:
                peaks.remove_i2()

            return peaks

        def _use_area(self, i1=None, i2=None):
            """returns True if area is found as a peak parameter.

            Args:
                i1 (int): primary index i1. If i1='A', raises an error.
                i2 (int, optional): secondary (spectrum) index i2. Only used if 
                    peaks has secondary indexes. Default is None.
            
            Returns:
                bool
            """
            #######################
            # check primary index #
            #######################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            if i1 is None:
                    if self.number_of_peaks(i2=i2) == 1:
                        i1 = self._get_indexes_i1()[0]
                    else:
                        raise ValueError('Peaks have multiple peaks and i1 must be defined')
            assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'

            ############
            # use area #
            ############
            if self._has_i2():
                if f'area_{i1}_{i2}' in self:
                    return True
                elif f'amp_{i1}_{i2}' in self:
                    return False
            else:    
                if f'area_{i1}' in self:
                    return True
                elif f'amp_{i1}' in self:
                    return False

        def _get_peaks_iterable(self, i2=None):
            """Returns a list where each element is a Peak object with one peak.

            Developers note: For now, we are returning a list, but we might return
                a dict in the future.

            Args:
                i2 (int, optional): secondary (spectrum) index i2. Only used if 
                    peaks has secondary indexes. Default is None.
            
            Returns:
                list
            """         
            final = [self._get_peaks_by_index(i1=i1, i2=i2) for i1 in self._get_indexes_i1(i2=i2)]
            for peak in final:
                # remove i2
                try:
                    peak.check_i2()
                except ValueError:
                    peak.remove_i2()

                # user defined attrs
                for attr in peak._get_user_attrs():
                    peak.__setattr__(attr, self.__dict__[attr])
                
                # additional attr
                if 'A' in peak._get_indexes_i1():
                    peak._additional = self._additional

            return final
                
        def _get_spectrum_iterable(self, i1=None):
            """Returns a list where each element a Peak object of each spectrum.

            Developers note: For now, we are returning a list, but we might return
                a dict in the future.
            
            If some spectrum with i2 does not have i1 defined, the respective list
            element comes out empty.

            Args:
                i1 (int): primary index i1. Use i1='A' for additional parameters.      

            Returns:
                list
            """
            #########################
            # check secondary index #
            #########################
            error_message = 'Peaks() does not have parameters from multiple spectra.'
            has_i2 = self._has_i2()
            if has_i2 is None:
                return []
            elif has_i2 == False:
                raise ValueError(error_message)
                    
            # final
            final = [self._get_peaks_by_index(i1=i1, i2=i2) for i2 in self._get_indexes_i2()]
            for peaks in final:
                # remove i2
                try:
                    peaks.check_i2()
                except ValueError:
                    peaks.remove_i2()

                # user defined attrs
                for attr in peaks._get_user_attrs():
                    peaks.__setattr__(attr, self.__dict__[attr])
                # additional attr
                if 'A' in peaks._get_indexes_i1():
                    peaks._additional = self._additional

            return final

        def _find_suitable_x(self):
            """returns a array with x values for plotting.

            Returns:
                array
            """
            if self._has_i2():
                peaks_for_each_spectra = self._get_spectrum_iterable()

                vmin  = []
                vmax  = []
                step  = []
                for peaks in peaks_for_each_spectra:
                    temp = peaks._find_suitable_x()
                    vmin.append(min(temp))
                    vmax.append(max(temp))
                    step.append(abs(temp[1]-temp[0]))
                return np.arange(min(vmin), max(vmax), min(step))

            else:
                peaks = self._get_peaks_iterable()

                vmin  = []
                vmax  = []
                step  = []
                for peak in peaks:

                    # get width
                    if peak.is_asymmetric():
                        w = peak[f'w1'].value + peak[f'w2'].value
                        m = max([peak[f'm1'].value, peak[f'm2'].value])
                    else:
                        w = peak[f'w'].value
                        m = peak[f'm'].value

                    if w == 0:
                        w = peak[f'c'].value*0.1
                    if w == 0:
                        w = 1

                    # vmin and vmax
                    _vmin = peak[f'c'].value - w * (m*4 + 4)
                    _vmax = peak[f'c'].value + w * (m*4 + 4)

                    if _vmin == 0 and _vmax == 0:
                        _vmin = -10
                        _vmax = 10

                    # return np.arange(vmin, vmax, w/20)
                    vmin.append(_vmin)
                    vmax.append(_vmax)
                    step.append(w/20)
                return np.arange(min(vmin), max(vmax), min(step))


            if i1 == 'A':
                raise ValueError('i1 = "A" is not valid')
            if i1 is None and i2 is None:  # has_i2 is assumed False. Gets suitable x for full spectra
                vmin  = []
                vmax  = []
                step  = []
                for peak in self._get_peaks_iterable(i2=i2):
                    temp = peak._find_suitable_x()
                    vmin.append(min(temp))
                    vmax.append(max(temp))
                    step.append(abs(temp[1]-temp[0]))
                return np.arange(min(vmin), max(vmax), min(step))
            elif i1 is None:               # has_i2 is assumed True. Gets suitable x for full spectra
                vmin  = []
                vmax  = []
                step  = []
                for peak in self._get_peaks_iterable(i2=i2):
                    temp = peak._find_suitable_x()
                    vmin.append(min(temp))
                    vmax.append(max(temp))
                    step.append(abs(temp[1]-temp[0]))
                return np.arange(min(vmin), max(vmax), min(step))
            elif i2 is None:  # has_i2 is assumed False. Gets suitable x for peak i1
                peak = self._get_peaks_by_index(i1=i1, i2=i2)

                # get width
                if peak.is_asymmetric(i):
                    w = peak[f'w1_{i1}'].value + peak[f'w2_{i1}'].value
                    m = max([peak[f'm1_{i1}'].value, peak[f'm2_{i1}'].value])
                else:
                    w = peak[f'w_{i1}'].value
                    m = peak[f'm_{i1}'].value

                if w == 0:
                    w = peak[f'c_{i1}'].value*0.1
                if w == 0:
                    w = 1

                # vmin and vmax
                vmin = peak[f'c_{i1}'].value - w * (m*4 + 4)
                vmax = peak[f'c_{i1}'].value + w * (m*4 + 4)

                if vmin == 0 and vmax == 0:
                    vmin = -10
                    vmax = 10

                return np.arange(vmin, vmax, w/20)
            else:             # has_i2 is assumed True.  Gets suitable x for peak i1 and i2
                peak = self._get_peaks_by_index(i1=i1, i2=i2)

                # get width
                if peak.is_asymmetric(i):
                    w = peak[f'w1_{i1}_{i2}'].value + peak[f'w2_{i1}_{i2}'].value
                    m = max([peak[f'm1_{i1}_{i2}'].value, peak[f'm2_{i1}_{i2}'].value])
                else:
                    w = peak[f'w_{i1}_{i2}'].value
                    m = peak[f'm_{i1}_{i2}'].value

                if w == 0:
                    w = peak[f'c_{i1}_{i2}'].value*0.1
                if w == 0:
                    w = 1

                # vmin and vmax
                vmin = peak[f'c_{i1}_{i2}'].value - w * (m*4 + 4)
                vmax = peak[f'c_{i1}_{i2}'].value + w * (m*4 + 4)

                if vmin == 0 and vmax == 0:
                    vmin = -10
                    vmax = 10

                return np.arange(vmin, vmax, w/20)

        #########
        # basic #
        #########
        def append(self, amp=None, c=None, w=None, w1=None, w2=None, m=0, m1=0, m2=0, area=None, i1=None, i2=None):
            """
            """     
            #########################
            # check secondary index #
            #########################
            error_message = 'Peaks() has parameters from multiple spectra.' +\
                            'Therefore, you must indicate the spectrum index i2.' +\
                            'You can also use i2="all" to add this parameter to ' +\
                            'every spectra'
            has_i2 = self._has_i2()
            if has_i2 is None and i2 == 'all':
                raise ValueError('i2 cannot be "all", because Peaks does not have secondary (spectrum) i2 indexes')
            elif has_i2 == True:
                if i2 is None:
                    raise ValueError(error_message)
                assert i2 == 'all' or br.numanip.is_integer(i2), 'i2 must be an integer or "all"'
                # elif br.numanip.is_integer(i2):
                #     if i2 not in self._get_indexes_i2():
                #         raise ValueError(f'i2={i2} cannot be found\nvalid i2: {self._get_indexes_i2()}')
                # elif i2 == 'all':
                #     pass
                # else:
                #     raise TypeError('i2 must be an integer or "all"')
            
            #######################
            # check primary index #
            #######################
            # returns a _i1
            if i1 is None:
                if has_i2 == None and br.numanip.is_integer(i2):
                    _i1 = 0                
                elif has_i2 == True and i2 == 'all':
                    _i1 = -1
                    for _i2 in self._get_indexes_i2():
                        _indexes = self._get_indexes_i1(i2=_i2)
                        if _indexes != []:
                            temp = max(_indexes)
                            if _i1 < temp:
                                _i1 = temp
                    _i1 += 1
                    # for _i2 in self._get_indexes_i2():
                    #     self.append(amp=amp, c=c, w=w, w1=w1, w2=w2, m=m, m1=m1, m2=m2, area=area, i1=None, i2=_i2)
                    # return
                elif has_i2 == True and br.numanip.is_integer(i2):
                    if i2 in self._get_indexes_i2():
                        _indexes = self._get_indexes_i1(i2=i2)
                        if _indexes == []:
                            _i1 = 0
                        else:
                            _i1 = max(_indexes) + 1
                    else:
                        _i1 = 0
                else:
                    _indexes = self._get_indexes_i1()
                    if _indexes == []:
                        _i1 = 0
                    else:
                        _i1 = max(_indexes) + 1
            else:
                assert br.numanip.is_integer(i1), 'i must be an integer'
                if has_i2 == None and br.numanip.is_integer(i2):
                    _i1 = i1
                elif has_i2 and i2 == 'all':
                    for _i2 in self._get_indexes_i2():
                        if i1 in self._get_indexes_i1(i2=_i2):
                            raise ValueError(f'i1={i1} already exists in spectrum index i2={_i2}\nPlease, select a valid index i1')
                    _i1 = i1
                elif has_i2 and br.numanip.is_integer(i2):
                    if i1 in self._get_indexes_i1(i2=i2):
                        raise ValueError(f'i1={i1} already exists in spectrum index i2={i2}\nPlease, select a valid index i1')
                    _i1 = i1
                else:
                    if i1 in self._get_indexes_i1():
                        raise ValueError(f'i1={i1} already exists\nPlease, select a valid index i1')
                    _i1 = i1

            ###################################
            # asserting validity of the input #
            ###################################
            assert area is not None or amp is not None, 'amp or area not defined'
            assert c    is not None, 'c not defined'
            assert w    is not None or (w1 is not None and w2 is not None), 'Either w or (w1 and w2) must be defined'

            assert area is None or br.numanip.is_number(area), 'area must be a number'
            assert amp  is None or br.numanip.is_number(amp), 'amp must be a number'
            assert c    is None or br.numanip.is_number(c), 'c must be a number'
            assert w    is None or br.numanip.is_number(w), 'w must be a number'
            assert w1   is None or br.numanip.is_number(w1), 'w1 must be a number'
            assert w2   is None or br.numanip.is_number(w2), 'w2 must be a number'
            assert m    is None or br.numanip.is_number(m), 'm must be a number'
            assert m1   is None or br.numanip.is_number(m1), 'm1 must be a number'
            assert m2   is None or br.numanip.is_number(m2), 'm2 must be a number'

            if br.numanip.is_number(w):
                assert w >= 0, 'w must be positive'
            if br.numanip.is_number(w1):
                assert w1 >= 0, 'w1 must be positive'
            if br.numanip.is_number(w2):
                assert w2 >= 0, 'w2 must be positive'

            if br.numanip.is_number(m):
                assert m <= 1 and m >= 0, 'm must be positive less than 1'
            if br.numanip.is_number(m1):
                assert m1 <= 1 and m1 >= 0, 'm1 must be positive less than 1'
            if br.numanip.is_number(m2):
                assert m2 <= 1 and m2 >= 0, 'm2 must be positive less than 1'

            ############
            # add peak #
            ############
            if has_i2 == True and i2 == 'all':
                for _i2 in self._get_indexes_i2():
                    self.append(i1=_i1, i2=_i2, amp=amp, c=c, w=w, w1=w1, w2=w2, m=m, m1=m1, m2=m2, area=area)
                return
            elif has_i2 == True or (has_i2 == None and i2 is not None):
                if area is None:
                    super().add(f'amp_{_i1}_{i2}',  value=amp,  vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)
                else:
                    super().add(f'area_{_i1}_{i2}', value=area, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)
                
                super().add(f'c_{_i1}_{i2}',        value=c,    vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None)

                if w1 is not None and w2 is not None:
                    super().add(f'w1_{_i1}_{i2}', value=w1, vary=True,     min=0,       max=np.inf, expr=None, brute_step=None)
                    super().add(f'w2_{_i1}_{i2}', value=w2, vary=True,     min=0,       max=np.inf, expr=None, brute_step=None)
                    super().add(f'm1_{_i1}_{i2}', value=m1, vary=False,    min=0,       max=1,      expr=None, brute_step=None)
                    super().add(f'm2_{_i1}_{i2}', value=m1, vary=False,    min=0,       max=1,      expr=None, brute_step=None)
                else:
                    super().add(f'w_{_i1}_{i2}',  value=w,  vary=True,     min=0,       max=np.inf, expr=None, brute_step=None)
                    super().add(f'm_{_i1}_{i2}',  value=m,  vary=False,    min=0,       max=1,      expr=None, brute_step=None)
                return
            elif has_i2 == False or has_i2 == None:
                if area is None:
                    super().add(f'amp_{_i1}',  value=amp,  vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
                else:
                    super().add(f'area_{_i1}', value=area, vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
                super().add(f'c_{_i1}',        value=c,    vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
                if w1 is not None and w2 is not None:
                    super().add(f'w1_{_i1}',   value=w1,   vary=True,  min=0,       max=np.inf, expr=None, brute_step=None)
                    super().add(f'w2_{_i1}',   value=w2,   vary=True,  min=0,       max=np.inf, expr=None, brute_step=None)
                    super().add(f'm1_{_i1}',   value=m1,   vary=False, min=0,       max=1,      expr=None, brute_step=None)
                    super().add(f'm2_{_i1}',   value=m2,   vary=False, min=0,       max=1,      expr=None, brute_step=None)
                else:
                    super().add(f'w_{_i1}',    value=w,    vary=True,  min=0,       max=np.inf, expr=None, brute_step=None)
                    super().add(f'm_{_i1}',    value=m,    vary=False, min=0,       max=1,      expr=None, brute_step=None)
                return

        def add(self, name, value=None, vary=True, min=-np.inf, max=np.inf, expr=None, brute_step=None, i2=None):
            """
            """
            ###################################
            # asserting validity of the input #
            ###################################
            assert '_' not in name, 'name cannot have "_" (underscore)'
            assert 'amp'  == name,  'name cannot be "amp"\namp is a reserved word'
            assert 'area' == name,  'name cannot be "area"\area is a reserved word'
            assert 'c'    == name,  'name cannot be "c"\nc is a reserved word'
            assert 'w'    == name,  'name cannot be "w"\nw is a reserved word'
            assert 'w1'   == name,  'name cannot be "w1"\nw1 is a reserved word'
            assert 'w2'   == name,  'name cannot be "w2"\nw2 is a reserved word'
            assert 'm'    == name,  'name cannot be "m"\nm is a reserved word'
            assert 'm1'   == name,  'name cannot be "m1"\nm1 is a reserved word'
            assert 'm2'   == name,  'name cannot be "m2"\nm2 is a reserved word'

            #########################
            # check secondary index #
            #########################
            error_message = 'Peaks() has parameters from multiple spectra.' +\
                            'Therefore, you must indicate the spectrum index i2.' +\
                            'You can also use i2="all" to add this parameter to ' +\
                            'every spectra'
            has_i2 = self._has_i2()
            if has_i2:
                if i2 is None:
                    raise ValueError(error_message)
                elif br.numanip.is_integer(i2):
                    if i2 not in self._get_indexes_i2():
                        raise ValueError(f'i2={i2} cannot be found\nvalid i2: {self._get_indexes_i2()}')
                elif i2 == 'all':
                    pass
                else:
                    raise TypeError('i2 must be an integer or "all"')
            
            if i2 == 'all':
                for _i2 in self._get_indexes_i2():
                    self.add(name, value=value, vary=vary, min=min, max=max, expr=expr, brute_step=brute_step, i2=_i2)
            else:
                name += '_A_' + str(i2)
                super().add(name, value=value, vary=vary, min=min, max=max, expr=expr, brute_step=brute_step)

        def remove(self, i1, i2=None):
            """Remove peak.

            Args:
                i1 (int): peaks to be removed. If i1='A', removes additional parameters.
                i2 (int, optional): secondary (spectrum) index i2. Only used if 
                    peaks has secondary indexes. Default is None.

            Returns:
                None
            """
            #########################
            # check secondary index #
            #########################
            error_message1 = 'Peaks() has parameters from multiple spectra.' +\
                            'Therefore, you must indicate the spectrum index i2.'
            error_message2 = 'Peaks() does not have parameters from multiple spectra.' +\
                            f'Therefore, i2 must be set to None and not {i2}'
            has_i2 = self._has_i2()
            if has_i2 is None:
                raise ValueError('no parameters to remove')
            elif has_i2:
                if i2 is None:
                    raise ValueError(error_message1)
                elif br.numanip.is_integer(i2):
                    if i2 not in self._get_indexes_i2():
                        raise ValueError(f'i2={i2} cannot be found\nvalid i2: {self._get_indexes_i2()}')
                else:
                    raise TypeError('i2 must be an integer')
            else:
                if i2 is not None:
                    raise ValueError(error_message2)
                
            #######################
            # check primary index #
            #######################
            assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'

            #####################
            # remove parameters #
            #####################
            for name in self._get_names_by_index(i1=i1, i2=i2):
                del self[name]

        def clear(self):
            """Erase all parameters."""
            super().clear()

        #################
        # save and load #
        #################
        def _create_header(self, verbose=False):
            """Gather attrs to be saved to a file."""
            header = ''
            attrs = self._get_user_attrs()
            attrs += ['_additional']
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
                elif br.numanip.is_number(self.__dict__[name]):
                    tosave = str(self.__dict__[name])
                    if tosave[-1] == '\n':
                        tosave = tosave[:-1]
                    header += f'{name}: {tosave}'  + '\n'
                else:
                    temp2 = str(self.__dict__[name]).replace('\n','\\n')
                    header += f'{name}: \"{temp2}\"'  + '\n'
            return header[:-1]
        
        def save(self, filepath=None, check_overwrite=False):
            r"""Save peak to a text file. Wrapper for `json.dump()`_.

            Args:
                filepath (string or path object, optional): filepath or file handle.
                check_overwrite (bool, optional): if True, it will check if file exists
                    and ask if user want to overwrite file.

            Returns:
                None

            .. _json.dumps(): https://docs.python.org/3/library/json.html#json.dumps
            """
            # filepath
            if filepath is None:
                if filepath == '':
                    raise TypeError("Missing 1 required argument: 'filepath'")
                else:
                    filepath = self.filepath               
            filepath = Path(filepath)

            # check overwrite
            if check_overwrite:
                if filepath.exists() == True:
                    if filepath.is_file() == True:
                        if br.interact.query('File already exists!! Do you wish to overwrite it?', 'yes') == True:
                            pass
                        else:
                            return
                    else:
                        raise AttributeError('filepath not pointing to a file.')

            # to add header with user defined parameters we have to use self.dumps
            # to create a string, then add the header TODO

            # save
            self.dump(filepath.open('w'))

        def load(self, filepath):
            """Load peak from a text file. Wrapper for `json.load()`_.

            Args:
                filepath (string or path object, optional): filepath or file handle.
                    If the filename extension is .gz or .bz2, the file is first decompressed.

            Returns:
                None

            .. _json.load(): https://docs.python.org/3/library/json.html#json.load
            """
            print('gg')
            # filepath
            if filepath is None:
                if filepath == '':
                    raise TypeError("Missing 1 required argument: 'filepath'")
                else:
                    filepath = self.filepath               
            filepath = Path(filepath)

            # clear
            self.clear()

            # read header
            # TODO

            # load
            super().load(filepath.open('r'))

        #################
        # check methods #
        #################
        def check_feasibility(self):
            """Check if start values are within boundaries."""
            for p in self:
                value = self[p].value
                vmax  = self[p].max
                vmin  = self[p].min

                if value > vmax or value < vmin:
                    raise ValueError(f'Start value for parameter {p} is out of bounds\nbound max = {vmax}\nstart value = {value}\nbound min = {vmin}')

        def check_i2(self):
            """raises error if i2 is defined, but is unique (only one spectrum i2)."""
            error_message = 'Peaks() have secondary (spectrum) index i2 defined.' +\
                            f'However, i2 seems to be unique.' +\
                            'Therefore, i2 is not necessary. Please, use peaks.remove_i2()'
            
            if self._has_i2():
                if len(self._get_indexes_i2()) == 1:
                    raise ValueError(error_message)
            return

        def remove_i2(self):
            """removes secondary (spectrum) index i2 from parameters."""
            if self._has_i2():
                if len(self._get_indexes_i2()) == 1:
                    to_remove = []
                    for name in self._get_names_by_index():
                        to_remove.append(name)
                        super().add('_'.join(name.split('_')[:-1]), value=self[name].value, vary=self[name].vary, min=self[name].min, max=self[name].max, expr=self[name].expr, brute_step=self[name].brute_step)
                    for name in to_remove:
                        del self[name]
                else:
                    raise ValueError('Peaks have multiple i2. Can only remove i2 if only one i2 exists')
            return
        
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
            # apply
            for name in self:
                if name.split('_')[0] == 'c':
                    if isinstance(value, Iterable):
                        f = lambda x: np.polyval(value, x)
                        self[name].value = np.float(f(self[name].value))
                    elif callable(value):
                        self[name].value = np.float(f(self[name].value))
                    else:
                        if value == 0:  # calib cannot be zero
                            raise ValueError('cannot set calib = 0.0')
                        elif value == 1:   # if calib is 1, do nothing
                            return
                        self[name].value = self[name].value*value
                elif name.split('_')[0] == 'w':
                    if isinstance(value, Iterable):
                        f = lambda x: np.polyval(value, x)
                        self[name].value = abs(f(self[name].value))
                    elif callable(value):
                        self[name].value = abs(f(self[name].value))
                    else:
                        if value == 0:  # calib cannot be zero
                            raise ValueError('cannot set calib = 0.0')
                        elif value == 1:   # if calib is 1, do nothing
                            return
                        self[name].value = abs(self[name].value*value)
            
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
            
            for name in self:
                if name.split('_')[0] == 'c':
                    self[name].value += value

        def set_offset(self, value):
            """Set offset value (in fact, it does nothing)

            Applying an offset, should not change the amp or area of peaks.

            Args:
                value (value): offset value.

            Returns:
                None
            """
            pass

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

            for name in self:
                if name.split('_')[0] == 'amp' or name.split('_') == 'area':
                    self[name].value  = self[name].value * value
    
        #############
        # modifiers #
        #############
        def reorder_by_attr(self, attr):
            """Change peak indexes based on the position of the peak or other attr
            
            args:
                attr (str): peak attr ('c', 'amp', 'w', ...)

            Returns:
                None        
            """
            raise NotImplementedError('sorry, not implemented yet.')

        ##############
        # extractors #
        ##############
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

        ##############
        # get values #
        ##############
        def get_area(self, i1=None, i2=None):
            """Returns the area of a peak

            If the amplitude was used (amplitude is defined), the area is calculated.

            i1 (int): primary index i1. If i1='A', raises an error.
            i2 (int, optional): secondary (spectrum) index i2. Only used if 
                peaks has secondary indexes. Default is None.

            Returns:
                number
            """
            #######################
            # check primary index #
            #######################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            if i1 is None:
                    if self.number_of_peaks(i2=i2) == 1:
                        i1 = self._get_indexes_i1()[0]
                    else:
                        raise ValueError('Peaks have multiple peaks and i1 must be defined')
            assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'

            ##################
            # check use area #
            ##################
            if self._use_area(i1=i1, i2=i2):
                if self._has_i2():
                    return self[f'area_{i1}_{i2}'].value
                else:
                    return self[f'area_{i1}'].value
                
            #######################
            # general definitions #
            #######################
            g = 2*np.sqrt(np.log(2))/np.sqrt(np.pi)

            ######################
            # collect parameters #
            ######################
            if self.is_asymmetric(i1=i1, i2=i2):
                if self._has_i2():
                    amp = self[f'amp_{i1}_{i2}'].value
                    w1 = self[f'w1_{i1}_{i2}'].value
                    w2 = self[f'w2_{i1}_{i2}'].value
                    m1 = self[f'm1_{i1}_{i2}'].value
                    m2 = self[f'm2_{i1}_{i2}'].value
                else:
                    amp = self[f'amp_{i1}'].value
                    w1 = self[f'w1_{i1}'].value
                    w2 = self[f'w2_{i1}'].value
                    m1 = self[f'm1_{i1}'].value
                    m2 = self[f'm2_{i1}'].value
            else:
                if self._has_i2():
                    amp = self[f'amp_{i1}_{i2}'].value
                    w   = self[f'w_{i1}_{i2}'].value
                    m   = self[f'm_{i1}_{i2}'].value
                else:
                    amp = self[f'amp_{i1}'].value
                    w   = self[f'w_{i1}'].value
                    m   = self[f'm_{i1}'].value

            ###############
            # Calculation #
            ###############
            if self.is_asymmetric(i1=i1, i2=i2):
                d1 = m1/np.pi + (1-m1)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                d2 = m2/np.pi + (1-m2)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return amp * (w1*d1**-1 + w2*d2**-1)/2
            else:
                d = m/np.pi + (1-m)*2*np.sqrt(np.log(2))/np.sqrt(np.pi)
                return amp * w * d**-1

        def get_area_error(self, i1=None, i2=None):
            """Returns one standard deviation error of the fitted area of a peak

            If the amplitude was used for fitting, the error is propagated to the area

            i1 (int): primary index i1. If i1='A', raises an error.
            i2 (int, optional): secondary (spectrum) index i2. Only used if 
                peaks has secondary indexes. Default is None.

            Returns:
                number
            """
            #######################
            # check primary index #
            #######################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            if i1 is None:
                    if self.number_of_peaks(i2=i2) == 1:
                        i1 = self._get_indexes_i1()[0]
                    else:
                        raise ValueError('Peaks have multiple peaks and i1 must be defined')
            assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'

            ##################
            # check use area #
            ##################
            if self._use_area(i1=i1, i2=i2):
                if self._has_i2():
                    return self[f'area_{i1}_{i2}'].stderr
                else:
                    return self[f'area_{i1}'].stderr
                
            #######################
            # general definitions #
            #######################
            g = 2*np.sqrt(np.log(2))/np.sqrt(np.pi)

            ######################
            # collect parameters #
            ######################
            if self.is_asymmetric(i1=i1, i2=i2):
                if self._has_i2():
                    amp = self[f'amp_{i1}_{i2}'].value
                    w1  = self[f'w1_{i1}_{i2}'].value
                    w2  = self[f'w2_{i1}_{i2}'].value
                    m1  = self[f'm1_{i1}_{i2}'].value
                    m2  = self[f'm2_{i1}_{i2}'].value
                    
                    D_amp = self[f'amp_{i1}_{i2}'].stderr
                    D_w1  = self[f'w1_{i1}_{i2}'].stderr
                    D_w2  = self[f'w2_{i1}_{i2}'].stderr
                    D_m1  = self[f'm1_{i1}_{i2}'].stderr
                    D_m2  = self[f'm2_{i1}_{i2}'].stderr
                else:
                    amp = self[f'amp_{i1}'].value
                    w1  = self[f'w1_{i1}'].value
                    w2  = self[f'w2_{i1}'].value
                    m1  = self[f'm1_{i1}'].value
                    m2  = self[f'm2_{i1}'].value
                    
                    D_amp = self[f'amp_{i1}'].stderr
                    D_w1  = self[f'w1_{i1}'].stderr
                    D_w2  = self[f'w2_{i1}'].stderr
                    D_m1  = self[f'm1_{i1}'].stderr
                    D_m2  = self[f'm2_{i1}'].stderr
            else:
                if self._has_i2():
                    amp = self[f'amp_{i1}_{i2}'].value
                    w   = self[f'w_{i1}_{i2}'].value
                    m   = self[f'm_{i1}_{i2}'].value
                    
                    D_amp = self[f'amp_{i1}_{i2}'].stderr
                    D_w   = self[f'w_{i1}_{i2}'].stderr
                    D_m   = self[f'm_{i1}_{i2}'].stderr
                else:
                    amp = self[f'amp_{i1}'].value
                    w   = self[f'w_{i1}'].value
                    m   = self[f'm_{i1}'].value
                    
                    D_amp = self[f'amp_{i1}'].stderr
                    D_w   = self[f'w_{i1}'].stderr
                    D_m   = self[f'm_{i1}'].stderr

            ###############
            # calculation #
            ###############
            if self.is_asymmetric(i1=i1, i2=i2):
                d1 = m1/np.pi + (1-m1)*g
                d2 = m2/np.pi + (1-m2)*g
                j  = w1/d1 + w2/d2
                return np.sqrt( ((j/2)**2 * D_amp**2) + 
                                ((amp/2/d1)**2 * D_w1**2) + 
                                ((amp/2/d2)**2 * D_w2**2) + 
                                ((((amp*w1/2)/d1**2)*(1/np.pi - g))**2 * D_m1**2) +
                                ((((amp*w1/2)/d2**2)*(1/np.pi - g))**2 * D_m2**2)   )
            else:
                d = m/np.pi + (1-m)*g
                return np.sqrt( ((w/d)**2 * D_amp**2) + 
                                ((amp/d)**2 * D_w**2) + 
                                ((((amp*w)/d**2)*(1/np.pi - g))**2 * D_m**2)   )

        def get_amp(self, i1=None, i2=None):
            """Returns the amplitude of a peak

            If the area was used (area is defined), the amplitude is calculated.

            i1 (int): primary index i1. If i1='A', raises an error.
            i2 (int, optional): secondary (spectrum) index i2. Only used if 
                peaks has secondary indexes. Default is None.

            Returns:
                number
            """
            #######################
            # check primary index #
            #######################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            if i1 is None:
                if self.number_of_peaks(i2=i2) == 1:
                    i1 = self._get_indexes_i1()[0]
                else:
                    raise ValueError('Peaks have multiple peaks and i1 must be defined')
            assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'

            ##################
            # check use area #
            ##################
            if self._use_area(i1=i1, i2=i2) == False:
                if self._has_i2():
                    return self[f'amp_{i1}_{i2}'].value
                else:
                    return self[f'amp_{i1}'].value
                
            #######################
            # general definitions #
            #######################
            g = 2*np.sqrt(np.log(2))/np.sqrt(np.pi)

            ######################
            # collect parameters #
            ######################
            if self.is_asymmetric(i1=i1, i2=i2):
                if self._has_i2():
                    area = self[f'area_{i1}_{i2}'].value
                    w1 = self[f'w1_{i1}_{i2}'].value
                    w2 = self[f'w2_{i1}_{i2}'].value
                    m1 = self[f'm1_{i1}_{i2}'].value
                    m2 = self[f'm2_{i1}_{i2}'].value
                else:
                    area = self[f'area_{i1}'].value
                    w1 = self[f'w1_{i1}'].value
                    w2 = self[f'w2_{i1}'].value
                    m1 = self[f'm1_{i1}'].value
                    m2 = self[f'm2_{i1}'].value
            else:
                if self._has_i2():
                    area = self[f'area_{i1}_{i2}'].value
                    w = self[f'w_{i1}_{i2}'].value
                    m = self[f'm_{i1}_{i2}'].value
                else:
                    area = self[f'area_{i1}'].value
                    w = self[f'w_{i1}'].value
                    m = self[f'm_{i1}'].value

            ###############
            # calculation #
            ###############
            if self.is_asymmetric(i1=i1, i2=i2):
                d1 = m1/np.pi + (1-m1)*g
                d2 = m2/np.pi + (1-m2)*g
                return area * (d1/w1/2 + d2/w2/2)
            else:
                d = m/np.pi + (1-m)*g
                return area*d/w
        
        def get_amp_error(self, i1=None, i2=None):
            """Returns one standard deviation error of the fitted amplitude of a peak

            If the area was used for fitting, the error is propagated to the amplitude.

            i1 (int): primary index i1. If i1='A', raises an error.
            i2 (int, optional): secondary (spectrum) index i2. Only used if 
                peaks has secondary indexes. Default is None.

            Returns:
                number
            """
            #######################
            # check primary index #
            #######################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            if i1 is None:
                if self.number_of_peaks(i2=i2) == 1:
                    i1 = self._get_indexes_i1()[0]
                else:
                    raise ValueError('Peaks have multiple peaks and i1 must be defined')
            assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'

            ##################
            # check use area #
            ##################
            if self._use_area(i1=i1, i2=i2) == False:
                if self._has_i2():
                    return self[f'amp_{i1}_{i2}'].stderr
                else:
                    return self[f'amp_{i1}'].stderr
                
            #######################
            # general definitions #
            #######################
            g = 2*np.sqrt(np.log(2))/np.sqrt(np.pi)

            ######################
            # collect parameters #
            ######################
            if self.is_asymmetric(i1=i1, i2=i2):
                if self._has_i2():
                    area = self[f'area_{i1}_{i2}'].value
                    w1   = self[f'w1_{i1}_{i2}'].value
                    w2   = self[f'w2_{i1}_{i2}'].value
                    m1   = self[f'm1_{i1}_{i2}'].value
                    m2   = self[f'm2_{i1}_{i2}'].value

                    D_area = self[f'area_{i1}_{i2}'].stderr
                    D_w1   = self[f'w1_{i1}_{i2}'].stderr
                    D_w2   = self[f'w2_{i1}_{i2}'].stderr
                    D_m1   = self[f'm1_{i1}_{i2}'].stderr
                    D_m2   = self[f'm2_{i1}_{i2}'].stderr
                else:
                    area = self[f'area_{i1}'].value
                    w1 = self[f'w1_{i1}'].value
                    w2 = self[f'w2_{i1}'].value
                    m1 = self[f'm1_{i1}'].value
                    m2 = self[f'm2_{i1}'].value

                    D_area = self[f'area_{i1}'].stderr
                    D_w1   = self[f'w1_{i1}'].stderr
                    D_w2   = self[f'w2_{i1}'].stderr
                    D_m1   = self[f'm1_{i1}'].stderr
                    D_m2   = self[f'm2_{i1}'].stderr
            else:
                if self._has_i2():
                    area = self[f'area_{i1}_{i2}'].value
                    w    = self[f'w_{i1}_{i2}'].value
                    m    = self[f'm_{i1}_{i2}'].value

                    D_area = self[f'area_{i1}_{i2}'].stderr
                    D_w    = self[f'w_{i1}_{i2}'].stderr
                    D_m    = self[f'm_{i1}_{i2}'].stderr
                else:
                    area = self[f'area_{i1}'].value
                    w = self[f'w_{i1}'].value
                    m = self[f'm_{i1}'].value

                    D_area = self[f'area_{i1}'].stderr
                    D_w    = self[f'w_{i1}'].stderr
                    D_m    = self[f'm_{i1}'].stderr
            
            ###############
            # calculation #
            ###############
            if self.is_asymmetric(i1=i1, i2=i2):
                d1 = m1/np.pi + (1-m1)*g
                d2 = m2/np.pi + (1-m2)*g
                j  = d1/w1 + d2/w2

                return np.sqrt( (j*D_area/2)**2 +
                                (area*d1/2/w1**2*D_w1)**2 +
                                (area*d2/2/w2**2*D_w2)**2 +
                                (area/2/w1*(1/np.pi-g)*D_m1)**2 +
                                (area/2/w2*(1/np.pi-g)*D_m2)**2)              
            else:
                d = m/np.pi + (1-m)*g

                return np.sqrt( (j*D_area/2)**2 +
                                (area*d/2/w**2*D_w)**2 +
                                (area/2/w*(1/np.pi-g)*D_m)**2) 

        def get_w(self, i1=None, i2=None):
            """Returns the FWHM of a peak

            i1 (int): primary index i1. If i1='A', raises an error.
            i2 (int, optional): secondary (spectrum) index i2. Only used if 
                peaks has secondary indexes. Default is None.

            Returns:
                number
            """
            #######################
            # check primary index #
            #######################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            if i1 is None:
                if self.number_of_peaks(i2=i2) == 1:
                    i1 = self._get_indexes_i1()[0]
                else:
                    raise ValueError('Peaks have multiple peaks and i1 must be defined')
            assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'

            ######################
            # collect parameters #
            ######################
            if self.is_asymmetric(i1=i1, i2=i2):
                if self._has_i2():
                    w1 = self[f'w1_{i1}_{i2}'].value
                    w2 = self[f'w2_{i1}_{i2}'].value
                else:
                    w1 = self[f'w1_{i1}'].value
                    w2 = self[f'w2_{i1}'].value
            else:
                if self._has_i2():
                    w = self[f'w_{i1}_{i2}'].value
                else:
                    w = self[f'w_{i1}'].value

            ###############
            # calculation #
            ###############
            if self.is_asymmetric(i1=i1, i2=i2):
                return w1 + w2
            else:
                return w
            
        def get_w_error(self, i1=None, i2=None):
            """Returns one standard deviation error of the fitted FWHM of a peak

            i1 (int): primary index i1. If i1='A', raises an error.
            i2 (int, optional): secondary (spectrum) index i2. Only used if 
                peaks has secondary indexes. Default is None.

            Returns:
                number
            """   
            #######################
            # check primary index #
            #######################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            if i1 is None:
                    if self.number_of_peaks(i2=i2) == 1:
                        i1 = self._get_indexes_i1()[0]
                    else:
                        raise ValueError('Peaks have multiple peaks and i1 must be defined')
            assert i1 in self._get_indexes_i1(i2=i2), f'index {i1} does not exist\nPeak indexes available: {self._get_indexes_i1(i2=i2)}'

            ######################
            # collect parameters #
            ######################
            if self.is_asymmetric(i1=i1, i2=i2):
                if self._has_i2():
                    D_w1 = self[f'w1_{i1}_{i2}'].stderr
                    D_w2 = self[f'w2_{i1}_{i2}'].stderr
                else:
                    D_w1 = self[f'w1_{i1}'].stderr
                    D_w2 = self[f'w2_{i1}'].stderr
            else:
                if self._has_i2():
                    D_w = self[f'w_{i1}_{i2}'].stderr
                else:
                    D_w = self[f'w_{i1}'].stderr

            ###############
            # calculation #
            ###############
            if self.is_asymmetric(i1=i1, i2=i2):
                return np.sqrt(D_w1**2 + D_w2**2)
            else:
                return D_w
            
        def get_values_from_peak_for_each_i2(self, attr, i1=None):
            """Returns a list with values of a peak parameter for each i2.

            Example:
                `Peaks.get_values_from_peak_for_each_i2(attr='amp', i1=2)` or 
                `Peaks.get_values_from_peak_for_each_i2(attr='amp_2')`
                will return the amplitude of peak number 2 for each spectrum in a list.

            Args:
                attr (str): attribute to get (e.g. `amp`, `amp_2`, `c`, `w`, `w1`, ...).
                i1   (int): Selects the peak number. Only used if the peak number is
                    not defined in `attr`. Default is None.
            
            Returns:
                list of values
            """
            #######################################
            # check primary and secondary indexes #
            #######################################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            assert self._has_i2(), 'i2 must be defined\ne.g, peaks for multiple spectra must be defined'
            self.check_i2()

            if i1 is None:
                assert '_' in attr, 'i1 must be defined either through attr or directly by i1'
                i1   = int(attr.split('_')[1])
                attr = attr.split('_')[0]

            ############################################
            # assert at least on i2 has the defined i1 #
            ############################################
            assert any([i1 in self._get_indexes_i1(i2=i2) for i2 in self._get_indexes_i2()]), f'no spectra has peak i1 index = {i1}'

            ##################
            # collect values #
            ##################
            final = []
            for i2 in self._get_indexes_i2():
                name = attr + f"_{i1}" + f"_{i2}"

                if i1 in self._get_indexes_i1(i2=i2):
                    final.append(self[name].value)
                else:
                    final.append(None)
                # if name in self:
                #     final.append(self[name].value)
                # elif name.startswith('area'):
                #     final.append(self.get_area(i1=i1, i2=i2))
                # elif name.startswith('amp'):
                #     final.append(self.get_amp(i1=i1, i2=i2))
                # elif name.startswith('w'):
                #     final.append(self.get_w(i1=i1, i2=i2))
                # else:
                #     raise ValueError(f'Cannot find parameter: {name}')
            return final

        def get_errors_from_peak_for_each_i2(self, attr, i1=None):
            """Returns a list with one standard deviation of a fitted peak parameter for each i2.

            Example:
                `Peaks.get_errors_from_peak_for_each_i2(attr='amp', i1=2)` or 
                `Peaks.get_errors_from_peak_for_each_i2(attr='amp_2')`
                will return the fitted error of peak number 2 for each spectrum in a list.

            Args:
                attr (str): attribute to get (e.g. `amp`, `amp_2`, `c`, `w`, `w1`, ...).
                i1   (int): Selects the peak number. Only used if the peak number is
                    not defined in `attr`. Default is None.
            
            Returns:
                list of values
            """
            #######################################
            # check primary and secondary indexes #
            #######################################
            assert i1 != 'A', f'i1 must be a peak index (integer), not "A"'
            assert self._has_i2(), 'i2 must be defined\ne.g, peaks for multiple spectra must be defined'
            self.check_i2()

            if i1 is None:
                assert '_' in attr, 'i1 must be defined either through attr or directly by i1'
                i1   = int(attr.split('_')[1])
                attr = attr.split('_')[0]

            ############################################
            # assert at least on i2 has the defined i1 #
            ############################################
            assert any([i1 in self._get_indexes_i1(i2=i2) for i2 in self._get_indexes_i2()]), f'no spectra has peak i1 index = {i1}'

            ##################
            # collect values #
            ##################
            final = []
            for i2 in self._get_indexes_i2():
                name = attr + f"_{i1}" + f"_{i2}"

                if name in self:
                    final.append(self[name].stderr)
                else:
                    final.append(None)
                # if name in self:
                #     final.append(self[name].value)
                # elif name.startswith('area'):
                #     final.append(self.get_area(i1=i1, i2=i2))
                # elif name.startswith('amp'):
                #     final.append(self.get_amp(i1=i1, i2=i2))
                # elif name.startswith('w'):
                #     final.append(self.get_w(i1=i1, i2=i2))
                # else:
                #     raise ValueError(f'Cannot find parameter: {name}')
            return final


        
        ########################
        # calculation and info #
        ########################
        def number_of_peaks(self, i2=None):
            """returns the number of peaks

            Does not include additional parameters.
            
            Args:
                i2 (int, optional): spectrum index. Only used if peaks has double 
                    indexes.

            Returns:
                int
            """
            l = self._get_indexes_i1(i2=i2)
            if 'A' in l:
                return len(l)-1
            else:
                return len(l)
        
        def number_of_i2(self):
            """returns the number of secondary (spectrum) indexes i2.

            Returns:
                int
            """
            return len(self._get_indexes_i2())
        
        def number_of_additional(self, i2=None):
            """returns the number of additional parameters
            
            Args:
                i2 (int, optional): spectrum index. Only used if peaks has double 
                    indexes.

            Returns:
                int
            """
            return len(self._get_names_by_index(i1='A', i2=i2))

        def iterable(self):
            """Returns an iterable.

            If secondary index i2 is present, iterate over spectra, 
                otherwise, iterate over peaks.

            Returns:
                Iterable
            """
            if self._has_i2():
                for peaks in self._get_spectrum_iterable():
                    yield peaks
            else:
                for peaks in self._get_peaks_iterable():
                    yield peaks

        def calculate_spectrum(self, x=None):
            """Return peaks curve.

            Args:
                x (list, optional): x values to which the curve will be calculated.
                    If None, a suitable x, with at least 20 points within the peak,
                    will be constructed.

            Returns:
                :py:class:`Spectrum`.
            """
            ##########################
            # check multiple spectra #
            ##########################
            error_message = 'Peaks() has parameters from multiple spectra.' +\
                            'Therefore, spectrum cannot be calculated.'
            if self._has_i2():
                raise ValueError(error_message)
            
            ##################
            # check no peaks #
            ##################
            assert self.number_of_peaks() > 0, 'Cannot create spectrum. There are no peaks defined'
            
            if x is None:
                x = self._find_suitable_x()
            f = self.model()

            s = br.Spectrum(x=x, y=f(x))

            return s

        def calculate_spectra(self, x=None):
            """Return a Spectra object.
            
            Args:
                x (list, optional): x values to which the curve will be calculated.
                    If None, a suitable x will be constructed.

            Returns:
                :py:class:`Spectra`.
            """
            if x is None:
                x = self._find_suitable_x()

            ss = br.Spectra()
            if self._has_i2():
                for peaks in self._get_spectrum_iterable():
                    ss.append(peaks.calculate_spectrum(x=x))
            else:
                for peaks in self._get_peaks_iterable():
                    ss.append(peaks.calculate_spectrum(x=x))
            return ss

        def _model_str(self, i1=None, i2=None):
            """Returns string for building model.

            Args:
                i1 (int): primary index i1. Use i1='A' for additional parameters.
                i2 (int, optional): secondary (spectrum) index i2. Only used if 
                    peaks has secondary indexes. Default is None.

            Returns:
                function f(x) as string
            """
            if i2 is None:
                if i1 is None:
                    final = ''
                    for i1 in self._get_indexes_i1(i2=i2):
                        final += self._model_str(i1=i1, i2=i2) + ' + '
                    
                    # add additional if string
                    if self._has_additional():
                        if isinstance(self.additional, str):
                            final += self.additional  + ' + '
                    return final[:-3]

                if self.is_asymmetric(i1=i1):
                    if self._use_area(i1=i1):
                        final = f"np.heaviside(c_{i1}-x, 0)*br.voigt_area_fwhm(x, area_{i1}, c_{i1}, w1_{i1},  m1_{i1}) + np.heaviside(x-c_{i1}, 0)*br.voigt_area_fwhm(x, area_{i1}, c_{i1},  w2_{i1},  m2_{i1}) + br.dirac_delta(x, amp_{i1}, c_{i1})" 
                    else:
                        final = f"np.heaviside(c_{i1}-x, 0)*br.voigt_fwhm(x, amp_{i1}, c_{i1}, w1_{i1},  m1_{i1}) + np.heaviside(x-c_{i1}, 0)*br.voigt_fwhm(x, amp_{i1}, c_{i1},  w2_{i1},  m2_{i1}) + br.dirac_delta(x, amp_{i1}, c_{i1})" 
                else:
                    if self._use_area(i1=i1):
                        final = f"br.voigt_area_fwhm(x, area_{i1}, c_{i1}, w_{i1}, m_{i1})"
                    else:
                        final = f"br.voigt_fwhm(x, amp_{i1}, c_{i1}, w_{i1}, m_{i1})"
                
                # add additional if string
                if self._has_additional():
                    if isinstance(self.additional, str):
                        final += ' + ' + self.additional
                return final
            else:
                if i1 is None:
                    final = ''
                    for i1 in self._get_indexes_i1(i2=i2):
                        final += self._model_str(i1=i1, i2=i2) + ' + '
                    # add additional if string
                    if self._has_additional():
                        if isinstance(self.additional, str):
                            final += self.additional  + ' + '
                    return final[:-3]

                if self.is_asymmetric(i1=i1, i2=i2):
                    if self._use_area(i1=i1, i2=i2):
                        final = f"np.heaviside(c_{i1}_{i2}-x, 0)*br.voigt_area_fwhm(x, area_{i1}_{i2}, c_{i1}_{i2}, w1_{i1}_{i2},  m1_{i1}_{i2}) + np.heaviside(x-c_{i1}_{i2}, 0)*br.voigt_area_fwhm(x, area_{i1}_{i2}, c_{i1}_{i2},  w2_{i1}_{i2},  m2_{i1}_{i2}) + dirac_delta(x, amp_{i1}_{i2}, c_{i1}_{i2})" 
                    else:
                        final = f"np.heaviside(c_{i1}_{i2}-x, 0)*br.voigt_fwhm(x, amp_{i1}_{i2}, c_{i1}_{i2}, w1_{i1}_{i2},  m1_{i1}_{i2}) + np.heaviside(x-c_{i1}_{i2}, 0)*br.voigt_fwhm(x, amp_{i1}_{i2}, c_{i1}_{i2},  w2_{i1}_{i2},  m2_{i1}_{i2}) + dirac_delta(x, amp_{i1}_{i2}, c_{i1}_{i2})" 
                else:
                    if self._use_area(i1=i1, i2=i2):
                        final = f"br.voigt_area_fwhm(x, area_{i1}_{i2}, c_{i1}_{i2}, w_{i1}_{i2}, m_{i1}_{i2})"
                    else:
                        final = f"br.voigt_fwhm(x, amp_{i1}_{i2}, c_{i1}_{i2}, w_{i1}_{i2}, m_{i1}_{i2})"
                
                # add additional if string
                if self._has_additional():
                    if isinstance(self.additional, str):
                        final += ' + ' + self.additional
                return final
        
        def model(self, i1=None, i2=None):
            """Returns a function f(x) for the peak.
            
            Args:
                i1 (int): primary index i1. Use i1='A' for additional parameters.
                i2 (int, optional): secondary (spectrum) index i2. Only used if 
                    peaks has secondary indexes. Default is None.

            Returns:
                function f(x)
            """
            if i2 is None:
                var_str = ''
                if i1 is None:
                    for name in self:
                        var_str += f'{name} = self["{name}"].value' + ', '
                else:
                    for name in self._get_names_by_index(i1=i1):
                        var_str += f'{name} = self["{name}"].value' + ', '

                model_str = f'lambda x, {var_str[:-2]}: {self._model_str(i1=i1)}'
                return eval(model_str)
            else:
                var_str = ''
                for name in self._get_names_by_index(i1=i1, i2=i2):
                    var_str += f'{name} = self["{name}"].value' + ', '

                model_str = f'lambda x, {var_str[:-2]}: {self._model_str(i1=i1, i2=i2)}'
                return eval(model_str)

        def _generate_residual_function(self):
            """Returns residual function for fitting."""
            if self._has_i2():
                def residual(params, xs, ys):
                    residual = []

                    # make residual per data set
                    for i2, x in enumerate(xs):
                        residual.append(ys[i2] - params.model(i2=i2)(x))

                    residual = np.concatenate(residual).ravel()
                    return residual

                return residual
            else:
                def residual(params, x, y):
                    # for name in params:
                    #     eval(f'{name} = self[{name}]')
                    model = params.model()
                    # pvals = params.valuesdict()
                    # model = eval(model_str)
                    return model(x)-y
                return residual
        
        def fit(self, x, y, method='least_squares', update_peaks=True):
            """fit spectrum with peaks model

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
            self.check_feasibility()
            residual = self._generate_residual_function()
            if self._has_i2():
                out = lmfit.minimize(residual, method=method, params=self, kws={'xs':x, 'ys':y})        
            else:
                out = lmfit.minimize(residual, method=method, params=self, args=(x, y))

            # return peaks2
            if update_peaks:
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
                    amp = d['prominences'][i]+max([y2[d['right_bases'][i]], y2[d['left_bases'][i]]])
                    c = x2[peaks[i]]
                    w = abs(d['widths'][i]*np.mean(np.diff(x)))

                    self.append(amp=amp, c=c, w=w)
            except IndexError:
                pass

        ##########################        
        # plot and visualization #
        ##########################  
        def plot(self, ax=None, offset=0, shift=0, factor=1, calib=1, **kwargs):
            """Place a marker at the maximum of every peak position. Wrapper for `matplotlib.pyplot.errorbar()`_.

            Does not work if peaks has multiple secondary (spectrum) indexes i2.

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

                c    = []
                amp  = []
                xerr = []
                for i1 in self._get_indexes_i1():
                    c    += [self[f'c_{i1}'].value]
                    amp  += [self.get_amp(i1=i1)]
                    xerr += [self.get_w(i1=i1)/2]

                if 'lw' not in kwargs and 'linewidth' not in kwargs:
                    kwargs['lw'] = 0
                if 'elinewidth' not in kwargs :
                    kwargs['elinewidth'] = 2
                if 'marker' not in kwargs :
                    kwargs['marker'] = 'o'
                if 'markersize' not in kwargs and 'ms' not in kwargs:
                    kwargs['markersize'] = 5

                c = (np.array(c)*calib)+shift
                amp = np.array(amp)*factor+offset
                xerr = np.array(xerr)*calib

                return ax.errorbar(c, amp, xerr=xerr, **kwargs)

        ###########
        # special #
        ###########
        def _copy_from_spectra(self, ss):
            """copy peaks from each spectra

            Args:
                ss (Spectra): br.Spectra object
            
            Returns:
                None
            """
            ss.peaks.clear()
            for i2, s in enumerate(ss):
                temp = s.peaks.copy()
                for name in temp:
                    temp[name].set(expr='')
                    self[name + f'_{i2}'] = temp[name]

        def _copy_to_spectra(self, ss):
            """copy peaks to the corresponding spectrum

            Args:
                ss (Spectra): br.Spectra object
            
            Returns:
                None
            """
            temp = self.copy()
            for par in temp:
                temp[par].set(expr='')

            for i2, s in enumerate(ss):
                s.peaks.clear()
                names = temp._get_names_by_index(i2=i2)
                for name in names:
                    _name = '_'.join(name.split('_')[:2])
                    s.peaks[_name] = temp[name]


        # is this useful? Probably obsolete
        def _replace_index(self, i1, i1_new, i2=None):
            """Change the index of a peak."""
            raise NotImplementedError('sorry, not implemented yet.')

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
            raise NotImplementedError('sorry, not implemented yet.')

            indexes = self.indexes()
            final = len(indexes)
            
            # replace indexes
            for i_new, index in enumerate(np.sort(indexes)):
                if i_new != index:
                    self.replace_index(index, i_new)

            # remove
            for index in self.indexes():
                if index >= final:
                    self.remove(index)

        def _model_str_old(self, i1=None):
            """Returns string for building peak function.

            Returns:
                function f(x) as string
            """
            if i1 is None:
                final = ''
                for i in self._get_indexes_i1():
                    final += self._model_str(i1) + ' + '
                
                return final[:-3]

            if self.is_asymmetric(i1):
                if self.use_area(i1):
                    return f"np.heaviside(c_{i1}-x, 0)*br.voigt_area_fwhm(x, area_{i1}, c_{i1}, w1_{i},  m1_{i1}) + np.heaviside(x-c_{i1}, 0)*br.voigt_area_fwhm(x, area_{i1}, c_{i1},  w2_{i1},  m2_{i1}) + dirac_delta(x, amp_{i1}, c_{i1})" 
                else:
                    return f"np.heaviside(c_{i1}-x, 0)*br.voigt_fwhm(x, amp_{i1}, c_{i1}, w1_{i1},  m1_{i1}) + np.heaviside(x-c_{i1}, 0)*br.voigt_fwhm(x, amp_{i1}, c_{i1},  w2_{i1},  m2_{i1}) + dirac_delta(x, amp_{i1}, c_{i1})" 
            else:
                if self.is_asymmetric(i):
                    return f"br.voigt_area_fwhm(x, area_{i1}, c_{i1}, w_{i1}, m_{i1})"
                else:
                    return f"br.voigt_fwhm(x, amp_{i1}, c_{i1}, w_{i1}, m_{i1})"

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
else:
    class Peaks():
        def __init__(self):
            self._step = 0
        
        def __len__(self):
            return 0

        ########################
        # do not raise a error #
        ########################
        def set_factor(self, *args, **kwargs):
            pass
            
        def set_shift(self, *args, **kwargs):
            pass

        def set_offset(self, *args, **kwargs):
            pass

        def set_calib(self, *args, **kwargs):
            pass

        def set_roll(self, *args, **kwargs):
            pass

        def clear(self, *args, **kwargs):
            pass

        ###############
        # raise error #
        ###############
        def find(self, *args, **kwargs):
            raise ModuleNotFoundError('lmfit cannot be found')
        
        def append(self, *args, **kwargs):
            raise ModuleNotFoundError('lmfit cannot be found')

        def fit(self, *args, **kwargs):
            raise ModuleNotFoundError('lmfit cannot be found')

        def _get_peaks_by_index(self, *args, **kwargs):
            raise ModuleNotFoundError('lmfit cannot be found')
        
        def _get_params_with_index(self, *args, **kwargs):
            raise ModuleNotFoundError('lmfit cannot be found')

        def _copy_from_spectra(self, *args, **kwargs):
            raise ModuleNotFoundError('lmfit cannot be found')

        def _copy_to_spectra(self, *args, **kwargs):
            raise ModuleNotFoundError('lmfit cannot be found')

        def get_values(self, *args, **kwargs):
            raise ModuleNotFoundError('lmfit cannot be found')



