#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core functions for I21 beamline

See the advanced.py file for higher level functions

Data is saved in hdf and nxs files, but we only need to read nxs files because 
they contain the data as well the metadata
"""

# %% ========================== Standard Imports ========================= %% #
from pathlib import Path
import numpy as np
import datetime

# %% =============================== brixs =============================== %% #
import brixs as br

# %% ========================== Special Imports ========================== %% #
try:
    import h5py
except:
    pass
# %%

# %% ============================= support =============================== %% #
def _str2datetime(string):
    """convert I21 date string pattern to date --> '2022-07-20T21:08:36.921'

    Args:
        string (str): string with I21 date string

    Returns:
        datetime.datetime object
    """
    #########
    # split #
    #########
    date, time = string.split('T')
    
    ########
    # date #
    ########
    year, month,  day = (int(_) for _ in date.split('-'))

    ########
    # time #
    ########
    hour, minute, seconds = time.split(':')

    hour    = int(hour)
    minute  = int(minute)
    seconds = int(round(float(seconds)))

    if seconds >= 60:
        minute  = minute + 1
        seconds = 0
    if minute >= 60:
        hour   = hour + 1
        minute = 0

    ############
    # datetime #
    ############
    return datetime.datetime(year=year, month=month, day=day, hour=hour, minute=minute, second=seconds)

def scanlist(folderpath):
    """Returns list with all available scan numbers inside folderpath
    
    Args:
        folderpath (str, Path): filepath where .nxs files can be found.

    Returns:
        list
    """
    # return [int(_.name.split('.')[0].split('-')[1]) for _ in br.parsed_filelist(folderpath, string='*.nxs', ref=1)]
    return br.parsed_filelist(folderpath, string='*.nxs', ref=1, return_type='dict')
# %%

# %% ========================== read functions =========================== %% #
def _read(filepath, verbose=False):
    """Read rixs scan file from I21 beamline at Diamond light source.
    
    Each scan should be composed by one or multiple images.

    Usage:
        im, ims = read(<filepath>)

    Args:
        filepath (str or path): path to .nxs file
        verbose (bool, optional): if True, print error message when metadata 
            cannot be retrieved. Default is False.

    Returns
        im, ims
        summed up final br.Image object and list with individual images 
    """
    #################
    # metadata list #
    #################
    # attrs: same for all images
    # attrs2: different for each image
    attrs = {'string':   {'script_name':      'entry/current_script_name',
                          'command':          'entry/diamond_scan/scan_command',
                          'facility_user_id': 'entry/user01/facility_user_id',
                          'username':         'entry/user01/name',
                          'beamline':         'entry/instrument/beamline',
                          'pol':              'entry/instrument/id/polarisation'},
             'number':   {'sample_x':         'entry/instrument/manipulator/x',
                          'sample_y':         'entry/instrument/manipulator/y',
                          'sample_z':         'entry/instrument/manipulator/z'},
             'round':    {'E':                'entry/sample/beam/incident_energy',
                          'exit_slit':        'entry/instrument/s5/v1_gap',
                          'T':                'entry/instrument/lakeshore336/sample',                                                      
                          'T_setpoint':       'entry/instrument/lakeshore336/demand',                                                    
                          'tth':              'entry/instrument/spectrometer/armtth',
                          'chi':              'entry/instrument/manipulator/chi',                                                                       
                          'phi':              'entry/instrument/manipulator/phi',                                                                         
                          'th':               'entry/instrument/manipulator/th',                             },
             'datetime': {'start_time':       'entry/start_time',
                          'end_time':         'entry/end_time'}, 
             'bool':     {'finished':         'entry/diamond_scan/scan_finished'}}
    attrs2 = {'round':   {'exposure_time':    'entry/instrument/andor/count_time'},
              'bool':    {'checkbeam' :       'entry/andor/checkbeam'}}

    ##########################
    # open file and get data #
    ##########################
    ind = br.Dummy()
    with h5py.File(Path(filepath), 'r') as f:

        ########
        # data #
        ########
        for data in f['entry']['andor']['data']:
            ind.append(br.Image(data=data))
        # number_of_images = f['entry']['scan_shape'][()]

        ###################
        # common metadata #
        ###################
        for _type in attrs:
            for name in attrs[_type]:
                address = attrs[_type][name]
                try: 
                    if _type == 'string':   attrs[_type][name] = f[address][()].decode("utf-8")
                    if _type == 'number':   attrs[_type][name] = f[address][()]
                    if _type == 'round':    attrs[_type][name] = round(f[address][()], 2)
                    if _type == 'datetime': attrs[_type][name] = _str2datetime(f[address][()].decode("utf-8"))
                    if _type == 'bool':     attrs[_type][name] = f[address][()][0] == 1
                except Exception as e:
                    if verbose: print(e)
                    attrs[_type][name] = None

        #####################
        # specific metadata #
        #####################
        for _type in attrs2:
            for name in attrs2[_type]:
                address = attrs2[_type][name]
                try: 
                    if _type == 'round':    attrs2[_type][name] = [round(_, 2) for _ in f[address][()]]
                    if _type == 'bool':     attrs2[_type][name] = [_ == 1 for _ in f[address][()]]
                except Exception as e:
                    if verbose: print(e)
                    attrs2[_type][name] = None

    ###############
    # final image #
    ###############
    # sum images
    data = np.zeros(ind[0].shape)
    for _im in ind:
        data += _im.data
    # create object
    im = br.Image(data)

    ##############
    # save attrs #
    ##############
    # print(attrs)
    for _type in attrs:
        for name in attrs[_type]:
            im.__setattr__(name, attrs[_type][name])
            ind.__setattr__(name, attrs[_type][name])
            for _im in ind:
                _im.__setattr__(name, attrs[_type][name])
    for _type in attrs2:
        for name in attrs2[_type]:
            im.__setattr__(name, attrs2[_type][name])
            ind.__setattr__(name, attrs2[_type][name])
            for i, _im in enumerate(ind):
                if attrs2[_type][name] is None:
                    _im.__setattr__(name, None)
                else:
                    _im.__setattr__(name, attrs2[_type][name][i])

    im.type = 'RIXS'
    return im, ind

def _read_xas(filepath, verbose=False):
    """Read xas scan file from I21 beamline at Diamond light source.
    
    Usage:
        diff1, TEY, TFY, I0 = read_xas(<filepath>)

    Args:
        filepath (str or path): path to .nxs file
        verbose (bool, optional): if True, print error message when metadata 
            cannot be retrieved. Default is False.

    Returns
        diff1, TEY, TFY, I0
    """
    #################
    # metadata list #
    #################
    attrs = {'string':   {#'script_name':      'entry/current_script_name',
                          'command':          'entry/scan_command',
                          'facility_user_id': 'entry/user01/facility_user_id',
                          'username':         'entry/user01/name',
                          'beamline':         'entry/instrument/beamline',
                          'pol':              'entry/instrument/id/polarisation'},
             'number':   {'sample_x':         'entry/instrument/manipulator/x',
                          'sample_y':         'entry/instrument/manipulator/y',
                          'sample_z':         'entry/instrument/manipulator/z'},
             'round':    {#'E':                'entry/sample/beam/incident_energy',
                          'exit_slit':        'entry/instrument/s5/v1_gap',
                          'T':                'entry/instrument/lakeshore336/sample',                                                      
                          'T_setpoint':       'entry/instrument/lakeshore336/demand',                                                    
                          'tth':              'entry/instrument/spectrometer/armtth',
                          'chi':              'entry/instrument/manipulator/chi',                                                                       
                          'phi':              'entry/instrument/manipulator/phi',                                                                         
                          'th':               'entry/instrument/manipulator/th',                             },
             'datetime': {'start_time':       'entry/start_time',
                          'end_time':         'entry/end_time'}, 
             'bool':     {'finished':         'entry/diamond_scan/scan_finished'}}
    attrs2 = {'split':   {'exposure_time':    'entry/scan_command'}}

    ##########################
    # open file and get data #
    ##########################
    ind    = []
    with h5py.File(Path(filepath), 'r') as f:

        ########
        # data #
        ########
        try: 
            # name = list(f['entry/diff1_c'].keys())
            # name.remove('diff1_c')
            # name = name[0]
            diff1 = br.Spectrum(x=f[f'entry/diff1_c/energy'][()], y=f['entry/diff1_c/data'][()])
            diff1.ok = True
        except Exception as e:
            if verbose: print(e)
            diff1 = br.Spectrum()
            diff1.ok = False
        diff1.mode = 'diff1'
        try: 
            # name = list(f['entry/fy2_c'].keys())
            # name.remove('fy2_i')
            # name = name[0]
            TFY   = br.Spectrum(x=f['entry/fy2_c/energy'][()], y=f['entry/fy2_c/data'][()])
            TFY.ok = True
        except Exception as e:
            if verbose: print(e)
            TFY = br.Spectrum()
            TFY.ok = False
        TFY.mode = 'TFY'
        try:            
            # name = list(f['entry/draincurrent_c'].keys())
            # name.remove('draincurrent_c')
            # name = name[0] 
            TEY   = br.Spectrum(x=f['entry/draincurrent_c/energy'][()], y=f['entry/draincurrent_c/data'][()])
            TEY.ok = True
        except Exception as e:
            if verbose: print(e)
            TEY = br.Spectrum()
            TEY.ok = False
        TEY.mode = 'TEY'
        try: 
            I0 = br.Spectrum(x=f['entry/m4c1_c/energy'][()], y=f['entry/m4c1_c/data'][()])
            I0.ok = True
        except Exception as e:
            if verbose: print(e)
            I0 = br.Spectrum()
            I0.ok = False
        I0.mode = 'I0'

        ###################
        # common metadata #
        ###################
        for _type in attrs:
            for name in attrs[_type]:
                address = attrs[_type][name]
                try: 
                    if _type == 'string':   attrs[_type][name] = f[address][()].decode("utf-8")
                    if _type == 'number':   attrs[_type][name] = f[address][()]
                    if _type == 'round':    attrs[_type][name] = round(f[address][()], 2)
                    if _type == 'datetime': attrs[_type][name] = _str2datetime(f[address][()].decode("utf-8"))
                    if _type == 'bool':     attrs[_type][name] = f[address][()][0] == 1
                except Exception as e:
                    if verbose: print(e)
                    attrs[_type][name] = None

        #####################
        # specific metadata #
        #####################
        for _type in attrs2:
            for name in attrs2[_type]:
                address = attrs2[_type][name]
                try: 
                    if _type == 'split':  
                        temp = f[address][()].decode("utf-8").split(' ')
                        time = []  
                        for name2 in ('diff1_c', 'fy2_c', 'draincurrent_c', 'm4c1_c'):
                            try:
                                i = temp.index(name2)
                                time.append(temp[i+1])
                            except:
                                time.append('err')
                        attrs2[_type][name] =  str(time).replace("'", '').replace(',', ';')
                    if _type == 'bool':     attrs2[_type][name] = [_ == 1 for _ in f[address][()]]
                except Exception as e:
                    if verbose: print(e)
                    attrs2[_type][name] = None

    ###############
    # final image #
    ###############
    # # sum images
    # data = np.zeros(ind[0].shape)
    # for _im in ind:
    #     data += _im.data
    # # create object
    # im = br.Image(data)

    ##############
    # save attrs #
    ##############
    # print(attrs)
    for _type in attrs:
        for name in attrs[_type]:
            for s in (diff1, TEY, TFY, I0):
                s.__setattr__(name, attrs[_type][name])
    for s in (diff1, TEY, TFY, I0):
        s.__setattr__('exposure_time', attrs2['split']['exposure_time'])
        s.type = 'XAS'
    # for _type in attrs2:
    #     for name in attrs2[_type]:
    #         for s in (diff1, TEY, TFY, I0):
    #             if attrs2[_type][name] is None:
    #                 s.__setattr__(name, None)
    #             else:
    #                 s.__setattr__(name, attrs2[_type][name][i])
    #         # s.__setattr__(name, attrs2[_type][name])


    return diff1, TEY, TFY, I0

def _read_line(filepath, verbose=False):
    """Read xas scan file from I21 beamline at Diamond light source.
    
    Usage:
        im, ims = read(<filepath>)

    Args:
        filepath (str or path): path to .nxs file
        verbose (bool, optional): if True, print error message when metadata 
            cannot be retrieved. Default is False.

    Returns
        im, ims
        summed up final br.Image object and list with individual images 
    """
    #################
    # metadata list #
    #################
    attrs = {'string':   {#'script_name':      'entry/current_script_name',
                          'command':          'entry/scan_command',
                          'facility_user_id': 'entry/user01/facility_user_id',
                          'username':         'entry/user01/name',
                          'beamline':         'entry/instrument/beamline',
                          'pol':              'entry/instrument/id/polarisation'},
             'number':   {'sample_x':         'entry/instrument/manipulator/x',
                          'sample_y':         'entry/instrument/manipulator/y',
                          'sample_z':         'entry/instrument/manipulator/z'},
             'round':    {'E':                'entry/sample/beam/incident_energy',
                          'exit_slit':        'entry/instrument/s5/v1_gap',
                          'T':                'entry/instrument/lakeshore336/sample',                                                      
                          'T_setpoint':       'entry/instrument/lakeshore336/demand',                                                    
                          'tth':              'entry/instrument/spectrometer/armtth',
                          'chi':              'entry/instrument/manipulator/chi',                                                                       
                          'phi':              'entry/instrument/manipulator/phi',                                                                         
                          'th':               'entry/instrument/manipulator/th',                             },
             'datetime': {'start_time':       'entry/start_time',
                          'end_time':         'entry/end_time'}, 
             'bool':     {'finished':         'entry/diamond_scan/scan_finished'}}
    attrs2 = {'split':   {'exposure_time':    'entry/scan_command'}}

    ##########################
    # open file and get data #
    ##########################
    ind    = []
    with h5py.File(Path(filepath), 'r') as f:

        ########
        # data #
        ########
        try: 
            name = list(f['entry/diff1_i'].keys())
            name.remove('diff1_i')
            name = name[0]
            diff1 = br.Spectrum(x=f[f'entry/diff1_i/{name}'][()], y=f['entry/diff1_i/diff1_i'][()])
            diff1.ok = True
        except Exception as e:
            if verbose: print(e)
            diff1 = br.Spectrum()
            diff1.ok = False
        diff1.mode = 'diff1'
        try: 
            name = list(f['entry/fy2_i'].keys())
            name.remove('fy2_i')
            name = name[0]
            TFY   = br.Spectrum(x=f[f'entry/fy2_i/{name}'][()], y=f['entry/fy2_i/fy2_i'][()])
            TFY.ok = True
        except Exception as e:
            if verbose: print(e)
            TFY = br.Spectrum()
            TFY.ok = False
        TFY.mode = 'TFY'
        try:           
            name = list(f['entry/draincurrent_i'].keys())
            name.remove('draincurrent_i')
            name = name[0]  
            TEY   = br.Spectrum(x=f[f'entry/draincurrent_i/{name}'][()], y=f['entry/draincurrent_i/draincurrent_i'][()])
            TEY.ok = True
        except Exception as e:
            if verbose: print(e)
            TEY = br.Spectrum()
            TEY.ok = False
        TEY.mode = 'TEY'
        try: 
            name = list(f['entry/m4c1_c'].keys())
            name.remove('m4c1_c')
            name = name[0]  
            I0 = br.Spectrum(x=f[f'entry/m4c1_c/{name}'][()], y=f['entry/m4c1_c/m4c1_c'][()])
            I0.ok = True
        except Exception as e:
            if verbose: print(e)
            I0 = br.Spectrum()
            I0.ok = False
        I0.mode = 'I0'

        ###################
        # common metadata #
        ###################
        for _type in attrs:
            for name in attrs[_type]:
                address = attrs[_type][name]
                try: 
                    if _type == 'string':   attrs[_type][name] = f[address][()].decode("utf-8")
                    if _type == 'number':   attrs[_type][name] = f[address][()]
                    if _type == 'round':    attrs[_type][name] = round(f[address][()], 2)
                    if _type == 'datetime': attrs[_type][name] = _str2datetime(f[address][()].decode("utf-8"))
                    if _type == 'bool':     attrs[_type][name] = f[address][()][0] == 1
                except Exception as e:
                    if verbose: print(e)
                    attrs[_type][name] = None

        #####################
        # specific metadata #
        #####################
        for _type in attrs2:
            for name in attrs2[_type]:
                address = attrs2[_type][name]
                try: 
                    if _type == 'split':  
                        temp = f[address][()].decode("utf-8").split(' ')
                        time = []  
                        for name2 in ('diff1_c', 'fy2_c', 'draincurrent_c', 'm4c1_c'):
                            try:
                                i = temp.index(name2)
                                time.append(temp[i+1])
                            except:
                                time.append('err')
                        attrs2[_type][name] =  str(time).replace("'", '').replace(',', ';')
                    if _type == 'bool':     attrs2[_type][name] = [_ == 1 for _ in f[address][()]]
                except Exception as e:
                    if verbose: print(e)
                    attrs2[_type][name] = None

    ###############
    # final image #
    ###############
    # # sum images
    # data = np.zeros(ind[0].shape)
    # for _im in ind:
    #     data += _im.data
    # # create object
    # im = br.Image(data)

    ##############
    # save attrs #
    ##############
    # print(attrs)
    for _type in attrs:
        for name in attrs[_type]:
            for s in (diff1, TEY, TFY, I0):
                s.__setattr__(name, attrs[_type][name])
    for s in (diff1, TEY):
        s.__setattr__('exposure_time', attrs2['split']['exposure_time'])
    # for _type in attrs2:
    #     for name in attrs2[_type]:
    #         for s in (diff1, TEY, TFY, I0):
    #             if attrs2[_type][name] is None:
    #                 s.__setattr__(name, None)
    #             else:
    #                 s.__setattr__(name, attrs2[_type][name][i])
    #         # s.__setattr__(name, attrs2[_type][name])

    return diff1, TEY, TFY, I0
# %%