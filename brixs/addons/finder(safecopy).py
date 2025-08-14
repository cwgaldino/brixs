#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Use finder to avoid running functions multiple types with same input parameters

Let's say you have a function processing_function(a, b, c) that returns a 
br.Spectrum type based on parameters `a`, `b`, and `c`. 

You can use the finder functionality to avoid running processing_function() 
multiple times with the same set of `a`, `b`, and `c` parameters.

How it works: The first time you run processing_function(a, b, c), 
it will process/calculate the desired spectrum based on `a`, `b`, and `c` 
parameters. It will then save the spectrum somewhere. If you try to run
processing_function(a, b, c) with the same parameters set, it will load and
return the already calculated spectrum. 







USAGE:

There are two ways to user the finder functionality



1) via decorators (recommended)

**1a) function that saves one Spectrum**

# define path where files will be saved

>>> br.finder.folderpath = '<some-path>'

# verbose for finder can be set (default is True)

>>> br.finder.verbose = True

# apply the decorator to you function

>>> @br.finder.track
>>> def processing_function(a, b, c):
>>>    s = <does something with a, b and c and returns s>
>>>    return s


**1b) function that saves one multiple Spectrum**

# define path where files will be saved

>>> br.finder.folderpath = '<some-path>'

# verbose for finder can be set (default is True)

>>> br.finder.verbose = True

# apply the decorator to you function

>>> @br.finder.track2
>>> def processing_function(a, b, c):
>>>    s1 = <does something with a, b and c and returns s1>
>>>    s2 = <does something with a, b and c and returns s2>
>>>    s3 = <does something with a, b and c and returns s3>
>>>    return s1, s2, s3


**1c) function that saves one Spectra**

# define path where files will be saved

>>> br.finder.folderpath = '<some-path>'

# verbose for finder can be set (default is True)

>>> br.finder.verbose = True

# apply the decorator to you function

>>> @br.finder.trackss
>>> def processing_function(a, b, c):
>>>    ss = <does something with a, b and c and returns ss>
>>>    return ss


    
2) manually setting up the finder function, where files will be save in `folderpath`

# verbose for finder can be set (default is True)
br.finder.verbose = True

# your processing function must be defined like this, where parameters is a dictionary
# with the calculation parameters

>>> def processing_function(parameters, folderpath):
>>> 
>>>     # try and find if spectrum has already been calculated
>>>     s = br.finder.search(parameters, folderpath=br.finder.folderpath)
>>>     if s is not None:
>>>         return s
>>>     
>>>     #### process something
>>>     <include-code-here>
>>> 
>>>     # save spectra so it is not needed to run it again
>>>     save(s=s, parameters=parameters, folderpath=folderpath)
>>> 
>>>     return s

Developers note: in the future, maybe we can make a better way to get the 
last file number in function _save()
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np

# %% ------------------------------- brixs -------------------------------- %% #
import brixs as br
# %%

# %% ============================== finder ================================ %% #
############################
# run finder via decorator #
############################
folderpath = ''
verbose    = True

import inspect
def track(func):
    """track one spectrum (should work for spectra too, if saved in one file)"""
    def inner(*args, **kwargs):
        
        ######################################################
        # run function directly if folderpath is not defined #
        ######################################################
        if folderpath == '':
            if verbose:
                print('cannot check if data was already processed, because finder folder is not defined')
            return func(*args, **kwargs)
        
        ##########################
        # get function arguments #
        ##########################
        attr_names, default_values = _get_function_args_and_default_values(func)
        kwargs = _args2kwargs(attr_names, default_values, args, kwargs)

        ########################################################
        # try and find if spectrum has already been calculated #
        ########################################################
        s = search(parameters=kwargs, folderpath=folderpath)
        if s is not None:
            return s

        ###############
        # calculation #
        ###############
        s = func(**kwargs)

        ######################
        # save args as attrs #
        ######################
        # for key in kwargs:
        #     s.__setattr__(key, kwargs[key])
        
        ####################################################
        # save spectra so it is not needed to run it again #
        ####################################################
        save(s=s, parameters=kwargs, folderpath=folderpath)

        # returning the value to the original frame
        return s
    return inner

def track2(func):
    """track 2 returns a list of spectra saved individually. function must return list"""
    def inner(*args, **kwargs):
        
        ######################################################
        # run function directly if folderpath is not defined #
        ######################################################
        if folderpath == '':
            if verbose:
                print('cannot check if data was already processed, because finder folderpath is not defined')
            return func(*args, **kwargs)
        
        ##########################
        # get function arguments #
        ##########################
        attr_names, default_values = _get_function_args_and_default_values(func)
        kwargs = _args2kwargs(attr_names, default_values, args, kwargs)

        ########################################################
        # try and find if spectrum has already been calculated #
        ########################################################
        if 'i1212___' not in kwargs:
            kwargs['i1212___'] = 0
        else:
            raise ValueError('Cannot run finder because one of the parameters is a reserved word "i1212___"')
        s = search(parameters=kwargs, folderpath=folderpath)
        if s is not None:
            total = s.i1212___total
            del s.i1212___total

            mult    = [None]*total
            mult[0] = s
            for i in range(1, total):
                kwargs['i1212___'] = i
                mult[i] = search(parameters=kwargs, folderpath=folderpath)
                del mult[i].i1212___total
            # return br.Spectra(data=mult)
            return mult
        _ = kwargs.pop('i1212___')
        
        ###############
        # calculation #
        ###############
        mult = func(**kwargs)
        
        ####################################################
        # save spectra so it is not needed to run it again #
        ####################################################
        for i, s in enumerate(mult):
            kwargs['i1212___'] = i
            s.i1212___total    = len(mult)
            save(s=s, parameters=kwargs, folderpath=folderpath)
            del s.i1212___total

        # returning the value to the original frame
        return mult
    return inner

def trackss(func):
    """function must return ss. ss is saved in one file
    """
    def inner(*args, **kwargs):
        
        ######################################################
        # run function directly if folderpath is not defined #
        ######################################################
        if folderpath == '':
            if verbose:
                print('cannot check if data was already processed, because finder folder is not defined')
            return func(*args, **kwargs)
        
        ##########################
        # get function arguments #
        ##########################
        attr_names, default_values = _get_function_args_and_default_values(func)
        kwargs = _args2kwargs(attr_names, default_values, args, kwargs)

        ########################################################
        # try and find if spectrum has already been calculated #
        ########################################################
        ss = searchss(parameters=kwargs, folderpath=folderpath)
        if ss is not None:
            return ss

        ###############
        # calculation #
        ###############
        ss = func(**kwargs)

        ######################
        # save args as attrs #
        ######################
        # for key in kwargs:
        #     s.__setattr__(key, kwargs[key])
        
        ####################################################
        # save spectra so it is not needed to run it again #
        ####################################################
        save(s=ss, parameters=kwargs, folderpath=folderpath)

        # returning the value to the original frame
        return ss
    return inner

def args2attrs(func):
    def inner(*args, **kwargs):
               
        ##########################
        # get function arguments #
        ##########################
        attr_names, default_values = _get_function_args_and_default_values(func)
        kwargs = _args2kwargs(attr_names, default_values, args, kwargs)

        ###############
        # calculation #
        ###############
        s = func(**kwargs)
        
        ######################
        # save args as attrs #
        ######################
        for key in kwargs:
            s.__setattr__(key, kwargs[key])

        # returning the value to the original frame
        return s
    return inner

######################################
# functions for manual finder set up #
######################################
def reset(folderpath_=None):
    """delete all spectrum and restart finder file"""
    if folderpath_ is None:
        if folderpath != '':
            folderpath_ = folderpath
        else:
            raise ValueError('invalid folderpath\nfolderpath argument missing')
    br.rmdir(folderpath_, not_found_ok=True)
    br.mkdir(folderpath_)
    f = open(folderpath_/'finder.txt', 'w')
    f.close()

def search(parameters, folderpath):
    """Check if data has already been calculated/processed with similar parameters

    Args:
        string (str): string to search inside filepath. String must be 
            representative of all parameters used to process the spectrum
        folderpath (string or Path): folderpath to save spectra. This folder path
            must contain a file named 'finder.txt'. If not, one will be created.


    Raises:
        keyerror: if finder_values have keys that do not match with finder_tags.

    Returns 
        spectrum if data is found, or None.   
    """
    # get names in alphabetical order
    names = np.sort(list(parameters.keys()))

    # search string
    search_string = ''
    for name in names:
        search_string += name + str(parameters[name]) + '_'

    search_result = _search(folderpath=folderpath, string=search_string)
    if isinstance(search_result, Path):
        if verbose:
            print(f'Loading data already processed: {search_result.name}')
        return br.Spectrum(filepath=search_result)
    else:
        return None

def searchss(parameters, folderpath):
    """Check if data has already been calculated/processed with similar parameters

    Args:
        string (str): string to search inside filepath. String must be 
            representative of all parameters used to process the spectrum
        folderpath (string or Path): folderpath to save spectra. This folder path
            must contain a file named 'finder.txt'. If not, one will be created.

    Raises:
        keyerror: if finder_values have keys that do not match with finder_tags.

    Returns 
        spectra if data is found, or None.   
    """
    # get names in alphabetical order
    names = np.sort(list(parameters.keys()))

    # search string
    search_string = ''
    for name in names:
        search_string += name + str(parameters[name]) + '_'

    search_result = _search(folderpath=folderpath, string=search_string)
    if isinstance(search_result, Path):
        if verbose:
            print(f'Loading data already processed: {search_result.name}')
        return br.Spectra(filepath=search_result)
    else:
        return None
    
def save(parameters, s, folderpath):
    """saves processed/calculated spectrum so one does not have to process it again

    Args:
        folderpath (string or Path): folderpath to save spectra.
        string (string): string to save in the finder file
        s (br.spectrum): spectrum to be saved
        folderpath (string or Path): folderpath to save spectra. This folder path
            must contain a file named 'finder.txt'. If not, one will be created.

    Returns:
        None
    """
    # get names in alphabetical order
    names = np.sort(list(parameters.keys()))

    # search string
    search_string = ''
    for name in names:
        search_string += name + str(parameters[name]) + '_'

    # save spectrum and put entry in the finder file
    _save(s=s, string=search_string, folderpath=folderpath)

##################
# Core functions #
##################
def _search(folderpath, string):
    """Returns spectrum filepath or filename.

    Args:
        string (str): string to search inside filepath. String must be 
            representative of all parameters used to process the spectrum
        folderpath (string or Path): folderpath to save spectra. This folder path
            must contain a file named 'finder.txt'. If not, one will be created.

    Returns:
        spectrum filepath, or None
    """
    # check folderpath
    folderpath = Path(folderpath)
    assert folderpath.exists(), f'folderpath does not exist\n{folderpath}'
    
    # check finder file
    if (folderpath/'finder.txt').exists() == False:
        f = open(folderpath/'finder.txt', 'w')
        f.close()
        return None

    # read finder file
    f = open(folderpath/'finder.txt', 'r')
    lines = f.read().splitlines()

    # search for string
    for i, line in enumerate(lines):
        if i%2 != 0:
            pass
        elif line == string:
            f.close()
            return Path(lines[i+1])
    return None#int(len(lines))

def _save(s, string, folderpath):
    """core function that saves processed/calculated spectrum
    
    Args:
        s (br.spectrum): spectrum to be saved
        string (string): string to save in the finder file
        folderpath (string or Path): folderpath to save spectra. This folder path
            must contain a file named 'finder.txt'. If not, one will be created.

    Returns:
        None
    """
    # check folderpath
    folderpath = Path(folderpath)
    assert folderpath.exists(), f'folderpath does not exist\n{folderpath}'
    
    # check finder file
    if (folderpath/'finder.txt').exists() == False:
        f = open(folderpath/'finder.txt', 'w')
        f.close()

    # get number of files to create filename
    # print(br.filelist(folderpath, string='finderfile*.dat'))
    # print('aaaaa')

    fl = [int(filename.name.split('.')[0].split('_')[1]) for filename in br.filelist(folderpath, string='finderfile*.dat')]
    if len(fl) == 0:
        fl = (-1, )
    # print(fl)
    filename = 'finderfile_' + str(max(fl) + 1) + '.dat'
    filepath2save = folderpath/filename
    # print(filepath2save)

    # save spectrum
    if isinstance(s, br.Spectrum):
        s.save(filepath2save)
    elif isinstance(s, br.Spectra):
        s.save_all_single_file(filepath2save)

    # save string to finder file
    f = open(folderpath/'finder.txt', 'a')
    f.write(string + '\n' + str(filepath2save) + '\n')
    f.close() 

#####################
# Support functions #
#####################
def _get_function_args_and_default_values(func):
    """returns function attrs and default values"""
    arguments      = inspect.signature(func).parameters
    attr_names     = list(arguments.keys())
    default_values = [_.default for _ in list(arguments.values())]
    return attr_names, default_values

def _args2kwargs(attr_names, default_values, args=[], kwargs={}):
    """convert (*args, **kwargs) input type to **kwargs and includes default values.

    Args:
        attr_name (list): list with attr names
        default_values (list): list with default values for the attrs
        args, kwargs (list and dict, optional): *args and **kwargs from function call.

    Returns:
        dict {attr_name: value}
    """
    # update values passed via args
    for i, value in enumerate(args):
        default_values[i] = value

    # define parameters
    _kwargs = {name:default_values[i] for i, name in enumerate(attr_names)}

    # update parameters via kwargs
    for key in kwargs:
        _kwargs[key] = kwargs[key]

    return _kwargs


