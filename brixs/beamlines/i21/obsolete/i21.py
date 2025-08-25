#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""base functions for I21 beamline

Data is saved in hdf and nxs files, but we only need to read nxs files because 
they contain the data as well the metadata



"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np
import datetime

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br

# %% ------------------------- Special Imports ---------------------------- %% #
try:
    import h5py
except:
    pass

# %% ============================= support ================================ %% #
def _str2datetime(string):
    """convert I21 date string pattern to date --> '2022-07-20T21:08:36.921'

    Args:
        string (str): string with I21 date string

    Return
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

# temp = r'C:\Users\galdin_c\Documents\current\2024_10_01-CuSbO6\apollo\2022_07_I21_MM30806-1\raw\mm30806-1'
# temp = r'C:\Users\galdin_c\Documents\current\2025_01_29_DLS_Zhijia\fixtures'
# diff1, TEY, TFY, I0 = i21.read_xas(filepath = Path(temp)/'i21-243868.nxs')

# try:
#     I0 = i21.read(filepath = Path(temp)/'i21-243890.nxs')
def scanlist(folderpath):
    # return [int(_.name.split('.')[0].split('-')[1]) for _ in br.parsed_filelist(folderpath, string='*.nxs', ref=1)]
    return br.parsed_filelist(folderpath, string='*.nxs', ref=1, return_type='dict')

# %% ================================ read ================================ %% #
def read(filepath, verbose=False):
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

def read_xas(filepath, verbose=False):
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

def read_line(filepath, verbose=False):
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

# %% ============================== obsolete =============================== %% #
def _h52dict(f):
    """Converts h5 to dict"""
    data = {}
    if type(f) == h5py._hl.files.File:
        data = _h52dict(f['entry'])
    elif type(f) == h5py._hl.group.Group:
        for key in f.keys():
            data[key] = _h52dict(f[key])
    elif type(f) == h5py._hl.dataset.Dataset:
        if f.shape == ():
            return None
        else:
            return f[:]
    else:
        print(f'Cannot recognize type {type(f)}. Skipping {f.name}')
    return data

def _read(filepath):
    """Read files from I21 beamline at Diamond light source.

    Args:
        filepath (str or Path object): filepath.

    Returns:
        :py:class:`brixs.Image`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    filepath = Path(filepath)

    f = h5py.File(Path(filepath), 'r')
    f = _h52dict(f)

    data = np.zeros(f['andor']['data'].shape[1:])
    for temp in f['andor']['data']:
        data += temp
    _ = f['andor'].pop('data')
    _ = f['instrument']['andor'].pop('data')

    im = br.Image(data=data)
    im.nd = f

    return im

def _read1(filepath):
    """Read file from I21 beamline at Diamond light source and return a list.

    Args:
        filepath (str or Path object): filepath.

    Returns:
        list with :py:class:`brixs.Image`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    filepath = Path(filepath)

    f = h5py.File(Path(filepath), 'r')
    f = _h52dict(f)

    data = []
    for temp in f['andor']['data']:
        temp = np.zeros(f['andor']['data'].shape[1:])
        data.append(br.Image(data=temp))

    _ = f['andor'].pop('data')
    _ = f['instrument']['andor'].pop('data')
    for im in data:
        im.nd = f

    return data

# %% HDF5 attributes list for RIXS scans ---------------------------------------
"""
# script
f['entry']['current_script_name']  # COPYING
# scan
f['entry']['scan_command'][()]  # COPYING
f['entry']['scan_fields'][()]   # not needed
# user data
f['entry']['user01']['facility_user_id'][()]  # COPYING
f['entry']['user01']['name'][()]              # COPYING
# date and time
f['entry']['start_time'][()]  # COPYING
f['entry']['end_time'][()]    # COPYING
f['entry']['duration'][()]    # not needed
# detector
f['entry']['andor'].keys()
f['entry']['andor']['checkbeam'][()]  # COPYING
f['entry']['andor']['data'][()]       # DATA
f['entry']['andor']['ds'][()]         # not needed
# scan2
f['entry']['diamond_scan'].keys() 
f['entry']['diamond_scan']['scan_finished'][()]            # COPYING
f['entry']['diamond_scan']['scan_command'][()]             # repeated
f['entry']['diamond_scan']['point_start_times'][()]        # not needed
f['entry']['diamond_scan']['scan_dead_time'][()]           # not needed
f['entry']['diamond_scan']['scan_dead_time_percent'][()]   # not needed
f['entry']['diamond_scan']['scan_estimated_duration'][()]  # not needed
f['entry']['diamond_scan']['scan_fields'][()]              # repeated
f['entry']['diamond_scan']['scan_rank'][()]                # not needed
f['entry']['diamond_scan']['scan_shape'][()]               # repeated
f['entry']['diamond_scan']['start_time'][()]               # repeated
# whatever 1
f['entry']['m4c1'].keys()  # ['checkbeam', 'ds', 'm4c1']
f['entry']['m4c1']['checkbeam'][()]  # repeated
f['entry']['m4c1']['ds'][()]         # repeated
f['entry']['m4c1']['m4c1'][()]       # not needed???
# whatever 2
f['entry']['m4c1_gain'].keys()  # ['checkbeam', 'ds', 'gain']
f['entry']['m4c1_gain']['checkbeam'][()]  # repeated
f['entry']['m4c1_gain']['ds'][()]         # repeated
f['entry']['m4c1_gain']['gain'][()]       # not needed???
# whatever 3
f['entry']['sample'].keys()  # ['beam']
f['entry']['sample']['beam'].keys()
f['entry']['sample']['beam']['incident_energy'][()]           # COPYING
f['entry']['sample']['beam']['incident_polarization'][()]     # not needed
f['entry']['sample']['beam']['beamExtentScannableName'][()]              # not needed
f['entry']['sample']['beam']['distance'][()]                             # not needed
f['entry']['sample']['beam']['fluxScannableName'][()]                    # not needed
f['entry']['sample']['beam']['incidentBeamDivergenceScannableName'][()]  # not needed
# andor
andor = f['entry']['instrument']['andor']
andor['count_time'][()]    # COPYING    
andor['data'][()]          # repeated
andor['description'][()]   # not needed
andor['local_name'][()]    # not needed
andor['manufacturer'][()]  # not needed
andor['model'][()]         # not needed
# andor settings
andor['andor_settings'].keys()
andor['andor_settings']['accumulation_period'][()]  # not needed
andor['andor_settings']['adc_speed'][()]            # not needed
andor['andor_settings']['binning_x'][()]            # not needed
andor['andor_settings']['binning_y'][()]            # not needed
andor['andor_settings']['cooler_control'][()]       # not needed
andor['andor_settings']['cooler_status'][()]        # not needed
andor['andor_settings']['cooler_temperature'][()]   # not needed
andor['andor_settings']['em_ccd_gain'][()]          # not needed
andor['andor_settings']['image_size_x'][()]         # not needed
andor['andor_settings']['image_size_y'][()]         # not needed
andor['andor_settings']['preamp_gain'][()]          # not needed
andor['andor_settings']['region_size_x'][()]        # not needed
andor['andor_settings']['region_size_y'][()]        # not needed
andor['andor_settings']['region_start_x'][()]       # not needed
andor['andor_settings']['region_start_y'][()]       # not needed
andor['andor_settings']['sensor_size_x'][()]        # not needed
andor['andor_settings']['sensor_size_y'][()]        # not needed
andor['andor_settings']['shutter_ext_TTL'][()]      # not needed
andor['andor_settings']['shutter_mode'][()]         # not needed
andor['andor_settings']['temperature_actual'][()]   # not needed
andor['andor_settings']['vertical_shift_amplitude'][()]  # not needed
andor['andor_settings']['vertical_shift_speed'][()]      # not needed
# poliriser
f['entry']['instrument']['polariser'].keys()
f['entry']['instrument']['polariser']['gamma'][()]  # not needed
f['entry']['instrument']['polariser']['stick'][()]  # not needed
# beamline
f['entry']['instrument']['beamline'][()]  # COPYING
f['entry']['instrument']['name'][()]  # repeated
# ring
f['entry']['instrument']['source'].keys()
f['entry']['instrument']['source']['name'][()]
f['entry']['instrument']['source']['current'][()]
f['entry']['instrument']['source']['energy'][()]
f['entry']['instrument']['source']['probe'][()]
f['entry']['instrument']['source']['type'][()]
# check beam
f['entry']['instrument']['checkbeam'].keys()
f['entry']['instrument']['checkbeam']['checkrc_beamok'][()]          # not needed
f['entry']['instrument']['checkbeam']['checktopup_time_beamok'][()]  # not needed
f['entry']['instrument']['checkbeam']['name'][()]                    # not needed
# whatever
f['entry']['instrument']['ds'].keys()
f['entry']['instrument']['ds']['name'][()]
f['entry']['instrument']['ds']['soft_limit_max'][()]
f['entry']['instrument']['ds']['soft_limit_min'][()]
f['entry']['instrument']['ds']['value'][()]
# fast_shutter
f['entry']['instrument']['fast_shutter'].keys()
f['entry']['instrument']['fast_shutter']['x'][()]
# undulator
f['entry']['instrument']['id'].keys()
f['entry']['instrument']['id']['bottomInner'][()]
f['entry']['instrument']['id']['bottomOuter'][()]
f['entry']['instrument']['id']['enabled'][()]
f['entry']['instrument']['id']['gap'][()]
f['entry']['instrument']['id']['harmonic'][()]
f['entry']['instrument']['id']['mode'][()]
f['entry']['instrument']['id']['polarisation'][()]  # COPYING
f['entry']['instrument']['id']['rowPhase'][()]
f['entry']['instrument']['id']['taper'][()]
f['entry']['instrument']['id']['topInner'][()]
f['entry']['instrument']['id']['topOuter'][()]
f['entry']['instrument']['id']['type'][()]
# mirrors motors
f['entry']['instrument']['m1'].keys()  # 'feedback', 'fine_pitch', 'fpsetpoint', 'height', 'pitch', 'roll', 'x', 'yaw'
f['entry']['instrument']['m2'].keys()  # 'feedback', 'fine_pitch', 'fpsetpoint', 'height', 'pitch', 'roll', 'x', 'yaw'
f['entry']['instrument']['m4'].keys()  # 'femto1', 'femto2', 'longy', 'rx', 'ry', 'rz', 'x', 'y', 'z'
f['entry']['instrument']['m4c1'].keys()  # 'count_time', 'coupling', 'gain', 'local_name', 'm4c1', 'mode'
f['entry']['instrument']['m5'].keys()  # 'hqrx', 'hqry', 'hqrz', 'hqx', 'hqy', 'hqz', 'longy', 'lqrx', 'lqry', 'lqrz', 'lqx', 'lqy', 'lqz', 'tth'
f['entry']['instrument']['pgm'].keys()  # 'b2_shadow', 'cff', 'energy', 'grating_pitch', 'grating_select', 'mirror_pitch', 'mirror_select'
f['entry']['instrument']['sgm'].keys()  # 'grating_select', 'h', 'pitch', 'r1', 'roll', 'wedge_nearside', 'wedge_offside', 'x'
# slit motors (I think s5 is the exit slit)
f['entry']['instrument']['s1'].keys()  # 'x_gap', 'x_pos', 'y_gap', 'y_pos'
f['entry']['instrument']['s2'].keys()  # 'x_gap', 'x_pos', 'y_gap', 'y_pos'
f['entry']['instrument']['s3'].keys()  # 'x_gap', 'x_pos', 'y_gap', 'y_pos'
f['entry']['instrument']['s4'].keys()  # 'lower', 'nearside', 'offside', 'upper', 'x_gap', 'x_pos', 'y_gap', 'y_pos'
f['entry']['instrument']['s5'].keys()  # 'h_gap', 'hdso', 'sut', 'v1_gap', 'v2_gap', 'vdso1', 'vdso2'
f['entry']['instrument']['s6'].keys()  # 'x_gap', 'x_pos', 'y_gap', 'y_pos'
# temperature
f['entry']['instrument']['lakeshore336'].keys()
f['entry']['instrument']['lakeshore336']['cryostat'][()]
f['entry']['instrument']['lakeshore336']['demand'][()]        # COPYING
f['entry']['instrument']['lakeshore336']['heater'][()]
f['entry']['instrument']['lakeshore336']['heater_range'][()]
f['entry']['instrument']['lakeshore336']['sample'][()]        # COPYING
f['entry']['instrument']['lakeshore336']['shield'][()]
# arm
f['entry']['instrument']['spectrometer'].keys()
f['entry']['instrument']['spectrometer']['armtth'][()]     # COPYING
f['entry']['instrument']['spectrometer']['specgamma'][()]
f['entry']['instrument']['spectrometer']['spech'][()]
f['entry']['instrument']['spectrometer']['specl'][()]
# manipulator
f['entry']['instrument']['manipulator'].keys() 
f['entry']['instrument']['manipulator']['chi'][()]  # COPYING
f['entry']['instrument']['manipulator']['phi'][()]  # COPYING
f['entry']['instrument']['manipulator']['th'][()]   # COPYING
f['entry']['instrument']['manipulator']['x'][()]    # COPYING
f['entry']['instrument']['manipulator']['y'][()]    # COPYING
f['entry']['instrument']['manipulator']['z'][()]    # COPYING
f['entry']['instrument']['manipulator']['difftth'][()]
f['entry']['instrument']['manipulator']['draincurrent'][()]
f['entry']['instrument']['manipulator']['sapara'][()]
f['entry']['instrument']['manipulator']['saperp'][()]
# momentum transfer
f['entry']['Q'].keys()
f['entry']['Q']['H'][()]
"""

# %% HDF5 attributes list for XAS scans ----------------------------------------
"""
f = h5py.File(TOP/'i21-341699.nxs', 'r')
f['entry'].keys()  # 'diamond_scan', 'diff1_c', 'draincurrent_c', 'duration', 'end_time', 
# 'experiment_identifier', 'fy2_c', 'instrument', 'm4c1_c', 'program_name', 'sample', 
# 'scan_command', 'scan_fields', 'scan_shape', 'start_time', 'user01', 'user_input'

# f['entry/diff1_c'].keys()   # 'data', 'energy'
# f['entry/draincurrent_c'].keys()  # 'data', 'energy'
# f['entry/fy2_c'].keys()    # 'data', 'energy'
# f['entry/m4c1_c'].keys()   # 'data', 'energy'

# f['entry/duration'][()] 
# f['entry/end_time'][()] 
# f['entry/start_time'][()]

f['entry/diamond_scan/duration'][()]

f['entry/diamond_scan'].keys()  # 'duration', 'end_time', 'keys', 'point_end_times', 'point_start_times', 'scan_command', 'scan_dead_time', 'scan_dead_time_percent', 'scan_estimated_duration', 'scan_fields', 'scan_finished', 'scan_rank', 'scan_shape', 'start_time'
f['entry/instrument'].keys()    # 'beamline', 'diff1_c', 'draincurrent_c', 'energy', 'fast_shutter', 'fy2_c', 'id', 'lakeshore336', 'm1', 'm2', 'm4', 'm4c1_c', 'm5', 'manipulator', 'name', 'pgm', 'polariser', 's1', 's2', 's3', 's4', 's5', 's6', 'sgm', 'source', 'spectrometer'
f['entry/sample/beam'].keys()  # 'distance', 'extent', 'flux', 'incident_beam_divergence', 'incident_energy', 'incident_polarization'

# f['entry/instrument/diff1_c/type'][()]#.keys()   

# f['entry/experiment_identifier'][()] 
# f['entry/scan_command'][()]
# f['entry/scan_fields'][()]
# f['entry/scan_shape'][()]
# f['entry/program_name'][()]  # b'GDA 9.32.0'

# f['entry/user01'].keys()  # 'facility_user_id', 'name'
f['entry/user_input/command'][()]  # b'cvscan energy 700 730 0.05 diff1_c 0.1 fy2_c 0.1 draincurrent_c 0.1 m4c1_c 0.1 checkbeam_cv'
"""