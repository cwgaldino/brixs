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


#########
# NOTES #
#########
As for right now, finder does not work with br.Image() type, but implementation 
should be straight-forward.

Spectra is saved in multiple files (one for each spectrum). Therefore, metadata
from the spectra object itself is not saved (only metadata for each spectrum).



#########
# Usage #
#########

There are two ways to user the finder functionality:


1) via decorator

>>> # finder settings
>>> br.finder.folderpath = <folder-to-save-temporary-files>
>>> br.finder.verbose    = if True, a text is printed every time an object is loaded via finder 
>>>
>>> @br.finder.track
>>> def processing_function(a, b, c):
>>>    s = <does something with a, b and c and returns s>
>>>    return s

Using finder via decorator will work in most simple cases. For more complex implementations
one needs to use method 2. See below.
    
2) manually setting up the finder function

>>> # finder settings
>>> br.finder.folderpath = <folder-to-save-temporary-files>
>>> br.finder.verbose    = if True, a text is printed every time an object is loaded via finder 
>>>
>>> def processing_function(a, b, c):
>>>
>>>     br.finder.kwargs = vars()
>>>     s = br.finder.search()
>>>     if s != False:
>>>         return s
>>>
>>>    s = <does something with a, b and c and returns s>
>>>
>>>    br.finder.save(s)
>>>    
>>>    return s

Method 2 has the advantage that one can have different folderpaths for different
functions by explicitly passing the folderpath. Also, the arguments can be 
modified manually:

>>> def processing_function(a, b, c):
>>>
>>>     s = br.finder.search(kwargs=<kwargs>, folderpath=<folderpath>)
>>>     if s != False:
>>>         return s
>>>
>>>    s = <does something with a, b and c and returns s>
>>>
>>>    br.finder.save(s, folderpath=<folderpath>)
>>>    
>>>    return s


############
# Warnings #
############

1. Watch out for functions that use variables defined outside of the scope of 
the function as they need to be included manually

>>> c = 100
>>>
>>> def processing_function(a, b):
>>>
>>>     # ---------> manually include c in kwargs!
>>>     br.finder.kwargs = vars()
>>>     br.finder.kwargs.update({'c': c})
>>>     s = br.finder.search()
>>>     if s != False:
>>>         return s
>>>
>>>    s = <uses a, b and c, but c which is defined outside of the function>
>>>
>>>    br.finder.save(s)
>>>    
>>>    return s


2. watch out for function that change dicts or lists inside the function

>>> def processing_function(a, b):
>>>
>>>     br.finder.kwargs = vars()
>>>     s = br.finder.search()
>>>     if s != False:
>>>         return s
>>>
>>>    s = <do something using a, b>
>>>
>>>    # changes b
>>>    b[2] = 300  # b is now [1, 2, 300]
>>>
>>>    br.finder.save(s)
>>>    
>>>    return s
>>> 
>>> # setting b
>>> b = [1, 2, 3]
>>>
>>> # if you run the function the first time, s will be calculated
>>> s = processing_function(a=0, b=b)  
>>>
>>> # if you run the function a second time, 
>>> s = processing_function(a=0, b=b)
>>>
>>> s will also be calculated because b changed and is now different. 
>>> This creates the impression that the input parameters are the same, when in 
>>> fact they are not, therefore, the function runs again.


3. watch out for functions that calls finder multiple times

>>> def process1(a, b):
>>>  
>>>     br.finder.kwargs = vars()
>>>     s = br.finder.search()
>>>     if s != False:
>>>         return s
>>>
>>>     s = <uses a, b>
>>>    
>>>     br.finder.save(s)
>>>    
>>>     return s
>>>
>>> def process2(a, b, c, d):
>>>
>>>     br.finder.kwargs = vars()
>>>     s = br.finder.search()
>>>     if s != False:
>>>         return s
>>>
>>>     s = process1(a, b)
>>>     s = <do something else with s, c, and d>
>>> 
>>>     br.finder.save(s)
>>>     
>>>     return s

because when you run process2(), the function br.finder.search() will set up the 
variable br.finder._search_string. However, when you run process1(), _search_string 
will be modified by the br.finder.search() inside process1(). This way, the finder
will work well for process1(), but at the end of process1() it will save the file
and set br.finder._search_string back to None. There will be an error when you 
try to run br.finder.save() in process2() because the search_string is not existing
anymore. To fix this, you have to manually save the search string obtained in 
process2() and again manually feed it to br.finder.save() in process2(). See below

>>> def process2(a, b, c, d):
>>>
>>>     br.finder.kwargs = vars()
>>>     s = br.finder.search()
>>>     _search_string = br.finder._search_string
>>>     if s != False:
>>>         return s
>>>
>>>     s = process1(a, b)
>>>     s = <do something else with s, c, and d>
>>>
>>>     br.finder._search_string = _search_string
>>>     br.finder.save(s)
>>>    
>>>     return s

###################
# Developers note #
###################
1. maybe we can make a better way to get the last file number in function _save().
I am thinking about using a variable as a counter. Finder would have the check 
the number of files only once when finder is restarted or when folderpath changes.

"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np
import inspect
# %%

# %% ------------------------------- brixs -------------------------------- %% #
import brixs as br
# %%

# --------------------------------- Settings ------------------------------ %% #
folderpath = ''
verbose    = True
kwargs     = ''
_search_string = None

# %% ============================= Support ================================ %% #
def reset(folderpath=None):
    """delete all spectrum/spectra and restart finder file"""
    ##############
    # folderpath #
    ##############
    if folderpath is None:
        folderpath = br.finder.folderpath
    folderpath = Path(folderpath)
    assert folderpath.exists(), f'folderpath does not exist\n{folderpath}'

    ################
    # remove files #
    ################
    br.rmdir(folderpath, not_found_ok=True)
    br.mkdir(folderpath)

    ######################
    # create finder file #
    ######################
    f = open(folderpath/'finder.txt', 'w')
    f.close()

    ##################
    # create counter #
    ##################
    f = open(folderpath/'_counter.txt', 'w')
    f.write('-1')
    f.close()

    return
# %%

# %% ============================== finder ================================ %% #
def _search(folderpath=None):
    """try to find _search_string inside finder.txt file

    Args:
        folderpath (string or Path, optional): folderpath to save spectra. 
            If None, filepath will be taken from `br.finder.folderpath`. 
            This folderpath must contain a file named 'finder.txt'. 
            If not, one will be created.

    Returns:
        object filepath, or None
    """   
    ##############
    # folderpath #
    ##############
    if folderpath is None:
        folderpath = br.finder.folderpath
    folderpath = Path(folderpath)
    assert folderpath.exists(), f'folderpath does not exist\n{folderpath}'

    #####################
    # check finder file #
    #####################
    if (folderpath/'finder.txt').exists() == False:
        f = open(folderpath/'finder.txt', 'w')
        f.close()
        return None

    ####################
    # read finder file #
    ####################
    f = open(folderpath/'finder.txt', 'r')
    lines = f.read().splitlines()

    #####################
    # search for string #
    #####################
    for i, line in enumerate(lines):
        if i%2 != 0:
            pass
        elif line == br.finder._search_string:
            f.close()
            br.finder._search_string = None
            return Path(lines[i+1])
    return None

def search(kwargs=None, folderpath=None):
    """Check if data has already been calculated/processed with similar parameters

    Args:
        kwargs (dict, optional): arguments (parameters) used. If None, it will 
            try to get kwargs from `br.finder.kwargs`.
        folderpath (string or Path, optional): folderpath to save spectra. 
            If None, filepath will be taken from `br.finder.folderpath`. 
            This folderpath must contain a file named 'finder.txt'. 
            If not, one will be created.

    Returns 
        False, or spectrum/spectra if data is found
    """
    ##############
    # folderpath #
    ##############
    if folderpath is None:
        folderpath = br.finder.folderpath
    folderpath = Path(folderpath)
    assert folderpath.exists(), f'folderpath does not exist\n{folderpath}'

    ##########
    # kwargs #
    ##########
    if kwargs is None:
        kwargs = br.finder.kwargs
        if kwargs == '':
            raise ValueError('kwargs must be defined, e.g., `br.finder.kwargs = vars()`')
    assert isinstance(kwargs, dict), f'kwargs must be a dict, not type {type(kwargs)}'

    ##################################
    # get vars in alphabetical order #
    ##################################
    names = np.sort(list(kwargs.keys()))

    ########################
    # create search string #
    ########################
    search_string = ''
    for name in names:
        search_string += name + str(kwargs[name]).replace('\n', '') + '_'
                         
    ######################
    # save search_string #
    ######################
    br.finder._search_string = search_string

    #################
    # search string #
    #################
    search_result = _search(folderpath=folderpath)

    ###########################
    # load obj or return None #
    ###########################
    if isinstance(search_result, Path):
        if verbose:
            print(f'Loading data already processed: {search_result.name}')

        if 'Spectrum' in search_result.name:
            return br.Spectrum().load(filepath=folderpath/search_result)
        elif 'PhotonEvents' in search_result.name:
            return br.PhotonEvents().load(filepath=folderpath/search_result)
        elif 'Spectra' in search_result.name:
            split = search_result.name.split('_')
            start = int(split[2])
            stop  = int(split[3])
            ss = br.Spectra()
            for i in range(start, stop+1):
                filename = f'finderfile_Spectra_{i}_{stop}'
                ss.append(br.Spectrum().load(filepath=folderpath/filename))
            
            if hasattr(ss[0], 'spectra_attrs_123_finder_copy'):
                for _attr in ss[0].spectra_attrs_123_finder_copy:
                    ss.__setattr__(_attr, ss[0].spectra_attrs_123_finder_copy[_attr])
            del ss[0].spectra_attrs_123_finder_copy

            return ss
        else:
            pass
    else:
        return False

def save(obj, folderpath=None):
    """saves processed/calculated spectrum so one does not have to process it again

    Args:
        obj (Spectrum, Spectra, or PhotonEvents): object to be saved. Does not 
            work with Image() yet.
        folderpath (string or Path, optional): folderpath to save spectra. 
            If None, filepath will be taken from `br.finder.folderpath`. 
            This folderpath must contain a file named 'finder.txt'. 
            If not, one will be created.
        
    Returns:
        None
    """
    ##############
    # folderpath #
    ##############
    if folderpath is None:
        folderpath = br.finder.folderpath
    folderpath = Path(folderpath)
    assert folderpath.exists(), f'folderpath does not exist\n{folderpath}'

    #####################
    # check finder file #
    #####################
    if (folderpath/'finder.txt').exists() == False:
        f = open(folderpath/'finder.txt', 'w')
        f.close()

    #################
    # check counter #
    #################
    if (folderpath/'_counter.txt').exists() == False:
        f = open(folderpath/'_counter.txt', 'w')
        f.write('-1')
        f.close()

    ###################
    # get next number #
    ###################
    f = open(folderpath/'_counter.txt', 'r')
    next_file_number = f.read()
    f.close() 
    next_file_number = int(next_file_number) + 1

    ###############
    # save object #
    ###############
    if isinstance(obj, br.Spectrum):
        filename = f'finderfile_Spectrum_{next_file_number}'
        obj.save(folderpath/filename)
    if isinstance(obj, br.PhotonEvents):
        filename = f'finderfile_PhotonEvents_{next_file_number}'
        obj.save(folderpath/filename)
    elif isinstance(obj, br.Spectra):
        filename = f'finderfile_Spectra_{next_file_number}_{next_file_number+len(obj)-1}'
        if len(obj) > 0:
            if hasattr(obj[0], 'spectra_attrs_123_finder_copy') == False:
                _attrsdict = {}
                for _attr in obj.get_attrs():
                    _attrsdict[_attr] = obj.__getattribute__(_attr)
                obj[0].spectra_attrs_123_finder_copy = _attrsdict
            obj[0].save(filepath=folderpath/filename)
            del obj[0].spectra_attrs_123_finder_copy

        start = next_file_number
        for i, _s in enumerate(obj[1:]):
            next_file_number += 1
            _filename = f'finderfile_Spectra_{next_file_number}_{start+len(obj)-1}'
            _s.save(filepath=folderpath/_filename)


    ###########################################
    # save string and filepath to finder file #
    ###########################################
    f = open(folderpath/'finder.txt', 'a')
    f.write(br.finder._search_string + '\n' + str(filename) + '\n')
    f.close() 

    ################
    # tick counter #
    ################
    f = open(folderpath/'_counter.txt', 'w')
    f.write(str(next_file_number))
    f.close() 

    #######################
    # reset search string #
    #######################
    br.finder._search_string = None

    return
# %%

# %% ============================ decorator =============================== %% #
def track(func):
    """track spectrum or spectra"""
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
        arguments      = inspect.signature(func).parameters
        attr_names     = list(arguments.keys())
        default_values = [_.default for _ in list(arguments.values())]

        # update values passed via args
        for i, value in enumerate(args):
            default_values[i] = value

        # define parameters
        _kwargs = {name:default_values[i] for i, name in enumerate(attr_names)}

        # update parameters via kwargs
        for key in kwargs:
            _kwargs[key] = kwargs[key]

        # save kwargs
        br.finder.kwargs = _kwargs

        ########################################################
        # try and find if spectrum has already been calculated #
        ########################################################
        s = search()
        if s != False:
            return s

        ###############
        # calculation #
        ###############
        obj = func(**_kwargs)
        
        ####################################################
        # save spectra so it is not needed to run it again #
        ####################################################
        save(obj)

        #######################
        # returning the value #
        #######################
        return obj
    return inner
# %%
