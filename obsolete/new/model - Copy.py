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
import lmfit
from scipy.signal import find_peaks

# BRIXS
import brixs as br


# %% Peaks ========================================================================
import lmfit
import numpy as np
import matplotlib.pyplot as plt
import brixs as br
plt.ion()

def generate_residual_function():
    """Returns residual function for fitting."""

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

    return residual

class Model(lmfit.Parameters):

    def __init__(self, parent):
        super().__init__()
        self.parent = parent

        ##############
        # components #
        ##############
        self.peaks           = Peaks(self)
        # self.asympeaks       = AsymmetricPeaks(self)
        # self.areapeaks       = AreaPeaks(self)
        # self.asym_area_peaks = AsymmetricAreaPeaks(self)

class AsymmetricPeaks(lmfit.Parameters):
    def __init__(self, parent):
        super().__init__()
        self.identifiers = ['amp1', 'c1', 'w1', 'w2', 'm1', 'm2']  # make it imutable
        self.parent      = parent

    #########
    # basic #
    #########
    def add(self, amp=None, c=None, w1=None, w2=None, m1=0, m2=0, i1=None, i2='all'):  
        self.parent.add(f'amp1_{i1}_{i2}', value=amp, vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
        self.parent.add(f'c1_{i1}_{i2}',   value=c,   vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
        self.parent.add(f'w1_{i1}_{i2}',   value=w1,   vary=True,  min=0,       max=np.inf, expr=None, brute_step=None)
        self.parent.add(f'w2_{i1}_{i2}',   value=w2,   vary=True,  min=0,       max=np.inf, expr=None, brute_step=None)
        self.parent.add(f'm1_{i1}_{i2}',   value=m1,   vary=False, min=0,       max=1,      expr=None, brute_step=None)
        self.parent.add(f'm2_{i1}_{i2}',   value=m2,   vary=False, min=0,       max=1,      expr=None, brute_step=None)
        return

    def _get_all_names(self):
        return [name for name in self.parent if name.split('_')[0] in self.identifiers]

    def _get_indexes_i2(self):
        """return list of all secondary (spectrum number) indexes i2
        
        Returns:
            list
        """
        return np.unique([name.split('_')[2] for name in self._get_all_names()])

class Peaks(lmfit.Parameters):
    def __init__(self, parent):
        super().__init__()

        ##############
        # EDIT THESE #
        ##############
        self.identifiers = ['amp', 'c', 'w', 'm']  # make it immutable
        def function(i1, i2):
            return f"br.voigt_fwhm(x, amp_{i1}_{i2}, c_{i1}_{i2}, w_{i1}_{i2}, m_{i1}_{i2})"
        
        self.base_modifiers = {'shift':  ['c', ],
                               'offset': [],
                               'calib':  ['c', 'w'],
                               'factor': ['amp', ]}   

        #########
        # final #
        #########
        self.parent      = parent
        self.function    = function


    ####################
    # must be modified #
    ####################
    def add(self, amp=None, c=None, w=None, m=0, i1=None, i2='all'):  
        """ok"""
        self.parent.add(f'amp_{i1}_{i2}', value=amp, vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
        self.parent.add(f'c_{i1}_{i2}',   value=c,   vary=True,  min=-np.inf, max=np.inf, expr=None, brute_step=None)
        self.parent.add(f'w_{i1}_{i2}',   value=w,   vary=True,  min=0,       max=np.inf, expr=None, brute_step=None)
        self.parent.add(f'm_{i1}_{i2}',   value=m,   vary=False, min=0,       max=1,      expr=None, brute_step=None)
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
        for name in self._get_all_names():
            if name.split('_')[0] in self.base_modifiers['calib']:
                if isinstance(value, Iterable):
                    f = lambda x: np.polyval(value, x)
                    new_value = np.float(f(self[name].value))
                elif callable(value):
                    new_value = np.float(f(self[name].value))
                else:
                    new_value = self[name].value*value
                
                self[name].value = new_value

            # EXTRA
            if name.split('_')[0] == 'w':
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
        for name in self._get_all_names():
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
        for name in self._get_all_names():
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
        for name in self._get_all_names():
            if name.split('_')[0] in self.base_modifiers['factor']:
                self[name].value += value

    #########
    # basic #
    #########
    def _get_all_names(self):
        """ok"""
        return [name for name in self.parent if name.split('_')[0] in self.identifiers]

    def _get_names_by_index(self, i1='all', i2='all'):
        """return list of names with index i1 and i2.

        Args:
            i1 (int, optional): primary index i1
            i2 (int, optional): secondary (spectrum) index i2
        
        Returns:
            list
        """           

        # transform in list
        if i2 is None:
            i2 = self._get_indexes_i2()
            if i1 is None:
                i1 = [self._get_indexes_i1(i2=_i2) for _i2 in i2]
            else:
                # certify that all i2 have i1
                i1_exists_for_all_i2 = {f'i2={_i2}': i1 in self._get_indexes_i1(i2=_i2) for _i2 in i2}
                for _i2 in i1_exists_for_all_i2:
                    if i1_exists_for_all_i2[_i2] == False: 
                        raise ValueError(f'i1 does not exists for all i2\n{i1_exists_for_all_i2}')
                # if passes the test, then:
                i1 = [[i1, ] for _ in i2]
        else:
            if i1 is None:
                i1 = [self._get_indexes_i1(i2=i2), ]
            else:
                i1 = [[i1], ]
            i2 = [i2, ]
      
        # gather attrs
        final = []
        for i, _i2 in enumerate(i2):
            for _i1 in i1[i]:
                for name in self:
                    split = name.split('_')
                    raw  = split[0]
                    if raw != 'A':
                        __i1 = int(split[1])
                        __i2 = int(split[2])
                        if __i1 == _i1 and __i2 == _i2:
                            final.append(name)
        return final
    
    def _get_indexes_i2(self):
        """return list of all secondary (spectrum number) indexes i2
        
        Returns:
            list
        """
        return list(np.unique([int(name.split('_')[2]) for name in self._get_all_names()]))

    def _get_indexes(self):
        names = self._get_all_names()
        _i2 = np.unique([int(name.split('_')[2]) for name in names])
        return {i2:list(np.unique([int(name.split('_')[1]) for name in names if int(name.split('_')[2]) == i2])) for i2 in _i2}
        
    #########
    # model #
    #########
    def check_feasibility(self):
        """Check if start values are within boundaries."""
        final = []
        for name in self._get_all_names():
            value = self.parent[name].value
            vmax  = self.parent[name].max
            vmin  = self.parent[name].min

            if value > vmax or value < vmin:
                final.append(f'parameter `{name}`\n    bound max = {vmax}\n    start value = {value}\n    bound min = {vmin}\n')
        if final != []:
            raise ValueError(f'Feasibility error. Start values outside boundaries:\n\n' + '\n'.join(final))
    
    def get_model_str(self, i1='all', i2='all'):
        """Returns string for building model.

        Args:
            i1 (int): primary index i1
            i2 (int, optional): secondary (spectrum) index i2

        Returns:
            function f(x) as string
        """
        # check i1 and i2
        # TODO
        
        
        # Initialalization
        indexes  = self._get_indexes()
        final    = ''
        variables = ''

        # model
        if i1 == 'all' and i2 == 'all':
            for _i2 in indexes:
                for _i1 in indexes[_i2]:
                    _variables, _final  = self.get_model_str(i1=_i1, i2=_i2)
                    final     +=  _final     + ' + '
                    variables +=  _variables + ', '
        elif i1 != 'all' and i2 == 'all':
            for _i2 in indexes:
                _variables, _final  = self.get_model_str(i1=i1, i2=_i2)
                final     +=  _final     + ' + '
                variables +=  _variables + ', '
        elif i1 == 'all' and i2 != 'all':
            for _i1 in indexes[_i2]:
                _variables, _final  = self.get_model_str(i1=_i1, i2=i2)
                final     +=  _final     + ' + '
                variables +=  _variables + ', '
        else:
            # final = f"br.voigt_fwhm(x, amp_{i1}_{i2}, c_{i1}_{i2}, w_{i1}_{i2}, m_{i1}_{i2}) + "
            final = self.function(i1=i1, i2=i2) + ' + '
            for name in self.identifiers:
                variables += f'{name}_{i1}_{i2} = self.parent["{name}_{i1}_{i2}"].value, '
        return variables[:-2], final[:-3]

    def get_model(self, i1='all', i2='all'):
        """Returns a function f(x) for the peak.
        
        Args:
            i1 (int): primary index i1
            i2 (int, optional): secondary (spectrum) index i2

        Returns:
            function f(x)
        """
        variables, _model = self.get_model_str(i1=i1, i2=i2)
        model = f'lambda x, {variables}: {_model}'
        return eval(model)





    
    def _get_indexes_i1(self, i2=None):
        """return list of all primary (peak number) indexes i1's
        
        Args:
            i2 (int, optional): spectrum index
        
        Returns:
            list
        """
        return np.unique([int(name.split('_')[2]) for name in self._get_all_names()])

        # get indexes
        final = []
        for name in self:
            if 'A' not in name:
                _i1   = int(name.split('_')[1])
                if _i1 not in final:
                    final.append(_i1) 
        return final

# %%
    #################
    # magic methods #
    #################
    def __setitem__(self, name, value):
        super().__setitem__(name, value)

    def __getitem__(self, name):
        
        #######
        # int #
        #######
        if isinstance(name, int):
            if self._unique_i2(): 
                i1 = name      
                i2 = self._get_indexes_i2()[0]             
            else:
                i2 = name      
                i1 = None
            return self._get_peaks_by_index(i1=i1, i2=i2) 
        ##############
        # additional #
        ##############
        elif name.startswith('A'):
            raise NotImplementedError('')
            # if self._unique_i2(): 

            #     raise ValueError('Index i2 must be a integer\n"A" is not valid')
            # return self._get_peaks_by_index(i1=name, i2=None) 
        ##########
        # string #
        ##########
        elif isinstance(name, str):
            # (1) `Peaks` have only one peak, and user just needs the attr name without i1 and i2
            # (2) `Peaks` have multiple peaks for one i2. The user must indicate i1 only
            # (3) `Peaks` have multiple peaks for multiple spectra. i1 and i2 must be indicated
            
            # if name is written exactly as in the dict, return it
            if name in self:
                return super().__getitem__(name)

            # if name cannot be find directly, try and see if i1 and i2 is defined in name
            i1 = None
            i2 = None
            if len(name.split('_')) == 2:    # i1 is defined
                i1 = name.split('_')[1]
            elif len(name.split('_')) == 3:  # i1 and i2 are defined
                i1 = name.split('_')[1]
                i2 = name.split('_')[2]
            elif len(name.split('_')) > 3:   # error
                raise ValueError(f'`{name}` is not a recognizable attr\nindexes are separated by `_`, like attr_i1_i2\n`{name}` appears to have more than two underscores')

            # if i1 and i2 are found in the name, assert that they are a number
            if i1 is not None:
                if br.numanip.is_integer(i1): i1 = int(i1)
                else: raise ValueError(f'`{name}` is not a recognizable attr')
            if i2 is not None:
                if br.numanip.is_integer(i2): i2 = int(i2)
                else: raise ValueError(f'`{name}` is not a recognizable attr')

            # if i2 is defined, check if it exists in `Peaks
            if i2 is not None:
                assert i2 in self._get_indexes_i2(), f'i2={i2} cannot be found\navailable i2={self._get_indexes_i2()}'
            
            # if i2 is not defined, check if i2 is unique and assign
            if i2 is None and self._unique_i2():
                i2 = self._get_indexes_i2()[0]
            else:
                raise ValueError(f'`{name}` is not a recognizable attr\n`{name}` does not specify i2')

            # if i1 is defined, check if it exists
            if i1 is not None:
                assert i1 in self._get_indexes_i1(i2=i2), f'i1={i1} cannot be found for i2={i2}\navailable i1={self._get_indexes_i1(i2=i2)}'

            # if i1 is not defined, check if i1 is unique and assign
            if i1 is None and self._unique_i1(i2=i2):
                i1 = self._get_indexes_i1(i2=i2)[0]
            else:
                raise ValueError(f'`{name}` is not a recognizable attr\n`{name}` does not specify i1')

            # at this point, i1 and i2 have been defined
            # rewrite name
            raw = name.split('_')[0]
            if i2 is not None:
                name = raw + '_' + str(i1) + '_' + str(i2)
            else:
                name = raw + '_' + str(i1)

            # try getting the attr again
            if name in self:
                return super().__getitem__(name)

            # if this doesn't work, maybe the user is trying to access a pseudo attr
            if raw == 'area':
                final        = lmfit.Parameter(name=name + ' (Pseudo)', value=self.get_area(i1=i1, i2=i2), vary=False)
                final.stderr = self.get_area_error(i1=i1, i2=i2)
                return final
            elif raw == 'amp':
                final        = lmfit.Parameter(name=name + ' (Pseudo)', value=self.get_amp(i1=i1, i2=i2), vary=False)
                final.stderr = self.get_amp_error(i1=i1, i2=i2)
                return final
            elif raw == 'w':
                final        = lmfit.Parameter(name=name + ' (Pseudo)', value=self.get_w(i1=i1, i2=i2), vary=False)
                final.stderr = self.get_w_error(i1=i1, i2=i2)
                return final 
            
            # if the function reaches this point, it means it could not find name
            raise ValueError(f'`{name}` is not a recognizable attr')

    def __delitem__(self, name):
        """TODO I want to avoid deleting single parameters from a peak.
        """
        super().__delitem__(name)
    
    ###########
    # support #
    ###########

    #########################
    # additional parameters #
    #########################

    ##################
    # support hidden #
    ##################
    def _get_user_attrs(self):
        """return attrs that are user defined."""
        default_attrs =  ['_asteval', '_additional']
        return [key for key in self.__dict__.keys() if key not in default_attrs and key.startswith('_') == False]
    
    def _unique_i2(self):
        """returns True if i2 is unique"""
        if len(self._get_indexes_i2()) == 1: return True
        return False

    def _unique_i1(self, i2=None):
        """returns True if i1 is unique for a certain i2
        
        it raises an error if i2=None and i2 is not unique"""
        ###################
        # check unique i2 #
        ###################
        if i2 is None:
            if self._unique_i2():
                i2 = self._get_indexes_i2()[0]
            else:
                raise ValueError(f'i2 must be defined\navailable i2={self._get_indexes_i2()}')
        else:
            assert i2 in self._get_indexes_i2(), f'i2={i2} cannot be found\navailable i2={self._get_indexes_i2()}'

        if len(self._get_indexes_i1(i2=i2)) == 1: return True
        return False
    
    def _check_i1_i2(self, i1, i2, i1_req=False, i2_req=False):
        """Check if i1 and i2 exist and are valid

        Does not consider additional parameters A

        check if they are numbers: no
        check if they are integer: no
        check if they exist: yes
        if none, assign to unique value if possible: yes
        
        Args:
            i1, i2 (number or None): indexes
            i1_req, i2_req (bool, optional): if True, raise error if i1, i2 aren't
                unique and not defined.

        Returns:
            i1, i2
        """
        #########################
        # check secondary index #
        #########################
        if i2 is None:
            # unique i2?
            if len(self._get_indexes_i2()) == 1: unique_i2 = True
            unique_i2 = False

            # check i2
            if self._unique_i2():
                i2 = self._get_indexes_i2()[0]
            else:
                if i2_req:
                    raise ValueError(f'i2 must be defined\navailable i2={self._get_indexes_i2()}')
        else:
            if i2 < 0:
                i2 = np.sort(self._get_indexes_i2())[i2]
            assert i2 in self._get_indexes_i2(), f'i2={i2} cannot be found\navailable i2={self._get_indexes_i2()}'

        #######################
        # check primary index #
        #######################
        # if i2 is not unique and not defined, i1 cannot be checked
        if i2 is not None:

            if i1 is None:
                # unique i1?
                if len(self._get_indexes_i1(i2=i2)) == 1: unique_i1 = True
                unique_i1 = False

                # check i1
                if unique_i1:
                    i1 = self._get_indexes_i1(i1=i2)[0]
                else:
                    if i1_req:
                        raise ValueError(f'i1 must be defined for i2={i2}\navailable i1 (for i2={i2})={self._get_indexes_i1(i2=i2)}')
            else:
                if i1 < 0:
                    i1 = np.sort(self._get_indexes_i1(i2=i2))[i1]
                assert i1 in self._get_indexes_i1(i2=i2), f'i1={i1} cannot be found\navailable i1 (for i2={i2})={self._get_indexes_i1(i2=i2)}'

        return i1, i2

    
    def _get_peaks_by_index(self, i1=None, i2=None):
        """return Peaks object with only selected peak.

        Does not return additional 'A' parameters

        The returned peaks object is "connect" to the parent Peaks object, i.e.
            they point to the same object. Modifications done in the returned
            peaks are reflected in the parent peaks. (Developers note, not sure
            this will work well long term. This might change eventually).

        Args:
            i1 (int): primary index i1. Use i1='A' for additional parameters.
            i2 (int, optional): secondary (spectrum) index i2. Only us
            
        Raises:
            ValueError if i2 is None, and i1 is not None, then i1 is required to exists
            for all secondary indexes i2. Otherwise it raises an error.
        
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

        return peaks




    
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
        
    
    


