#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from VERITAS beamline of MAX-IV.

######################
# Mesh scan profiles #
######################
         ┌──── b       /\        a
         |            /  \      /
    ┌────┘           /    \    /
    |               /      \  /
────┘              /        \/

#########
# Usage #
#########

>>> import brixs as br
>>> import brixs.beamlines.VERITAS as VERITAS
>>> 
>>> # filepaths
>>> RIXS = <filepath-rixs-hdf5>
>>> XAS  = <filepath-xas-hdf5>
>>> 
>>> # scanlist
>>> VERITAS.scanlist(RIXS)
>>> VERITAS.scanlist(XAS)
>>> 
>>> # read rixs
>>> pe = VERITAS.read(RIXS, 200)
>>> 
>>> # read xas
>>> ss = VERITAS.read(XAS, 142)
>>> TEY, TFY, EXF, RMU = VERITAS.read(XAS, 164)
>>> 
>>> # metadata tree
>>> tree = VERITAS.tree(RIXS, 200)
>>> print(tree)
>>> br.save_text(tree, 'tree.txt', encoding='utf-8')
>>>
>>> # get metadata value
>>> value = VERITAS.get_metadata(RIXS, 200, 'startmetadata/Sample')
>>>
>>> # Advanced: add metadata to read function
>>> # read() function extracts metadata from h5 file based on two dictionaries
>>> # VERITAS.rixs_attrs and VERITAS.xas_attrs
>>> VERITAS.rixs_attrs['raw']['<new_name>'] = <hdf5-entry-address>
>>> VERITAS.xas_attrs['raw']['<new_name>']  = <hdf5-entry-address>
>>>

Last edited: Carlos Galdino 2024-05-05
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np
import datetime
import warnings
import copy

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br

# %% ------------------------- Special Imports ---------------------------- %% #
import h5py
# %%

# %% ========================= useful functions =========================== %% #
def get_metadata(filepath, scan, entry):
    """returns metadata from hdf5 file for inspection (VERITAS)

    Args:
        filepath (str, Path): hdf5 filepath
        scan (number): scan number
        entry (str): hdf5 entry address, e.g., 'External/a_mp1_x/position'

    Returns:
        metadata for inspection
    """
    with h5py.File(Path(filepath), 'r') as f:
        prefix = ''.join([i for i in list(f.keys())[0] if not i.isdigit()])
    
        final = f[prefix + str(scan)][entry]
        try:
            final = final[()]
        except TypeError:
            raise TypeError(f"No data to return. Object is group with keys {final.keys()}")
    return final

def tree(filepath, scan):
    """Returns a text with the structure of a hdf5 file

    Args:
        filepath (str or path): file to hdf5 file
        scan (int): scan number

    Returns:
        str
    """
    with h5py.File(Path(filepath), 'r') as f:
        prefix = ''.join([i for i in list(f.keys())[0] if not i.isdigit()])
        text = _tree(f[prefix + str(scan)])
    return text

def scanlist(filepath):
    """Return list of scans available in filepath"""
    with h5py.File(Path(filepath), 'r') as f:
        prefix = ''.join([i for i in list(f.keys())[0] if not i.isdigit()])
        return [int(_) for _ in np.sort([int(scan.split(prefix)[1]) for scan in list(f.keys())])]

# %% ======================== support functions =========================== %% #
def _tree(f, pre=''):
    """Returns a text with the structure of a hdf5 file

    Usage:
        >>> import h5py
        >>> 
        >>> with h5py.File(<filepath>, 'r') as f:
        >>>     text = tree(f)
        >>>
        >>> with h5py.File(<filepath>, 'r') as f:
        >>>     text = tree(f[entry])
    
    Args:
        f (h5 object): h5 file object
        pre (str): For internal use only. This function needs to run recursively

    Returns:
        str
    """
    text = ''
    items = len(f)
    for key, f in f.items():
        items -= 1
        if items == 0:
            if type(f) == h5py._hl.group.Group:
                text += pre + '└── ' + key + '\n'
                text += _tree(f, pre + '    ')
            else:
                try:
                    text += pre + '└── ' + key + ' (%d)' % len(f) + '\n'
                except TypeError:
                    text += pre + '├── ' + key + ' (scalar)' + '\n'
        else:
            if type(f) == h5py._hl.group.Group:
                text += pre + '├── ' + key + '\n'
                text += _tree(f, pre+'│   ')
            else:
                try:
                    text += pre + '├── ' + key + ' (%d)' % len(f) + '\n'
                except TypeError:
                    text += pre + '├── ' + key + ' (scalar)' + '\n'
    return text

def _str2datetime(string):
    """convert VERITAS date/time string pattern to date --> '2022-07-20 21:08:36.921'

    Args:
        string (str): string with I21 date string

    Return
        datetime.datetime object
    """
    #########
    # split #
    #########
    if 'T' in string:
        date, time = string.split('T')
    else:
        date, time = string.split(' ')
    
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

def _sort_metadata(attrs, f, verbose=True):
    """get metadata from hdf5 file and pre-process it

    Usage: 
        >> attrs = {'raw':   {'A': 'External/b316a-o01/dia/tco-02/Temperature',
        >>                    'B': 'External/beamline_energy/position'},
        >>         'string': {'C': 'startmetadata/Sample',
        >>                    'D': 'endmetadata/Sample}}
        >> 
        >> pe = br.PhotonEvents()
        >> 
        >> with h5py.File(Path(filepath), 'r') as f:
        >>     metadata = _sort_metadata(attrs=attrs, f=f['entry1'], verbose=True)
        >>     for attr in metadata:
        >>         setattr(pe, attr, metadata[attr])

    attrs (dict): dict with dicts with attr names and hdf5 addresses. First dict
        entry will indicate the type o pre-processing (`raw`, `string`, `mean`, 
            `mean1`, `mean_round2`, `round2`, `datetime`, `bool`). To learn about 
            what they do, read the function declaration.  
    f (hdf5 object): hdf5 object, e.g., use `with h5py.File(Path(filepath), 'r') as f:`
    verbose (bool, optional): if True, if attr value cannot be extracted from 
        hdf5 file it will print a warning. Default is True.

    Returns:
        dict[attr] : value
    """
    # check unique names
    values = {}
    for _type in attrs:
        for name in attrs[_type]:
            if name in values:
                raise KeyError(f"name `{name}` is duplicated in VERITAS attrs list. Names must be unique")
            values[name] = None

    # get attr values
    for _type in attrs:
        for name in attrs[_type]:
            address = attrs[_type][name]
            try: 
                if _type == 'raw':   
                    values[name] = f[address][()]

                elif _type == 'string':   
                    values[name] = f[address][()].decode("utf-8")

                elif _type == 'mean':   
                    values[name] = np.mean(f[address][()])

                elif _type == 'mean1':   
                    values[name] = np.mean(f[address][()][:, 1])

                elif _type == 'mean1_round2':    
                    values[name] = round(f[address][()][:, 1], 2)

                elif _type == 'round2':    
                    values[name] = round(f[address][()], 2)

                elif _type == 'datetime': 
                    values[name] = _str2datetime(f[address][()].decode("utf-8"))

                elif _type == 'bool':     
                    values[name] = f[address][()][0] == 1

            except Exception as e:
                if verbose: 
                    print(f'get attr `{name}` error: {e}')
                # attrs[_type][name] = None

    return values
            
def _pol(gap, phase):
    if copy.deepcopy(phase) > 23:
        return 'LV'
    else:
        return 'LH'
# %%

# %% ========================== metadata list ============================= %% #
rixs_attrs = {'raw':    {'temperature':      'External/b316a-o01/dia/tco-02/Temperature',
                        'beamline_energy':  'External/beamline_energy/position',
                        # 'startmetadata_mono_energy': 'startmetadata/cffreal/Sample',
                        # 'endmetadata_mono_energy':   'endmetadata/cffreal/Sample',
                        # 'startmetadata_temperature': 'startmetadata/b316a-o01/dia/tco-02/Temperature',
                        # 'endmetadata_temperature':   'endmetadata/b316a-o01/dia/tco-02/Temperature'
                        },
             'string': {'startmetadata_sample': 'startmetadata/Sample',
                        'endmetadata_sample':   'endmetadata/Sample'
                        },
             'mean':   {'exposure_time':    'Instrument/DLD8080/exposure_time'
                        },
             'mean1':  {'a_mp1_x':           'External/a_mp1_x/position',
                        'a_mp1_y':           'External/a_mp1_y/position',
                        'a_mp1_z':           'External/a_mp1_z/position',
                        'a_mp1_yaw':         'External/a_mp1_yaw/position',
                        'sample_x':          'External/a_mp1_x/position',
                        'sample_y':          'External/a_mp1_y/position',
                        'sample_z':          'External/a_mp1_z/position',
                        'epu_r3_316_gap':    'External/epu_r3_316_gap/position',
                        'epu_r3_316_phase':  'External/epu_r3_316_phase/position',
                        'mono_energy_calib': 'External/mono_energy_calib/position',
                        'veritasarm_energy': 'External/veritasarm_energy/position'
                        }, 
             'mean_round2':  {'E':  'External/beamline_energy/position',
                              'T':  'External/b316a-o01/dia/tco-02/Temperature',
                              'th': 'External/a_mp1_yaw/position',
                              },                         
             'datetime':     {'start_time':       'start_time',
                              'end_time':         'end_time',
                             }}

xas_attrs = {'string': {'definition':       'definition',
                        'entry_identifier': 'entry_identifier',
                        'program_name':     'program_name',
                        'command':          'title',
                        'user':             'user/name',
                        },
             'raw':  {'a_mp1_x':       'measurement/pre_scan_snapshot/a_mp1_x',
                      'a_mp1_y':       'measurement/pre_scan_snapshot/a_mp1_y',
                      'a_mp1_z':       'measurement/pre_scan_snapshot/a_mp1_z',
                      'a_mp1_yaw':     'measurement/pre_scan_snapshot/a_mp1_yaw',
                      'sample_x':      'measurement/pre_scan_snapshot/a_mp1_x',
                      'sample_y':      'measurement/pre_scan_snapshot/a_mp1_y',
                      'sample_z':      'measurement/pre_scan_snapshot/a_mp1_z',
                      'a_slit1_hz':    'measurement/pre_scan_snapshot/a_slit1_hz',
                      'exit_slit ':    'measurement/pre_scan_snapshot/a_slit1_hz',
                      'a_m4_lateral':  'measurement/pre_scan_snapshot/a_m4_lateral',
                      'a_m4_pitch':    'measurement/pre_scan_snapshot/a_m4_pitch',
                      'a_m4_roll':     'measurement/pre_scan_snapshot/a_m4_roll',
                      'a_m4_vertical': 'measurement/pre_scan_snapshot/a_m4_vertical',
                      'a_m4_yaw':      'measurement/pre_scan_snapshot/a_m4_yaw',
                      'm1_yaw':        'measurement/pre_scan_snapshot/m1_yaw',
                      'm3_yaw':        'measurement/pre_scan_snapshot/m3_yaw',
                      }, 
             'round2':  {'th': 'measurement/pre_scan_snapshot/a_mp1_yaw',
                        },                         
             'datetime':     {'start_time':       'start_time',
                              'end_time':         'end_time',
                              'macro_start_time': 'macro_start_time',
                             }}
# %%

# %% ============================= read =================================== %% #
def read(filepath, scan, verbose=True):
    """Return data from filepath
    
    RIXS        --> returns PhotonEvent 
    XAS         --> returns Spectra (TEY, TFY, EXF, RMU) 
    sample scan --> returns Spectra (TEY, TFY, EXF, RMU) 
    mesh scan   --> returns Image (TEY, TFY, EXF, RMU) 

    No extra processing is done on RIXS spectra
    For XAS, sample scan, and mesh scan, datapoints are divided by the exposure

    Args:
        filepath (str or path): file to hdf5 file
        scan (int): scan number

    Returns:
        PhotonEvents for RIXS
        Spectra (TEY, TFY, EXF, RMU) for XAS and sample scans
        Image list (TEY, TFY, EXF, RMU) for mesh scans
    """

    with h5py.File(Path(filepath), 'r') as f:

        # find prefix
        prefix = ''.join([i for i in list(f.keys())[0] if not i.isdigit()])

        # scan list
        sl = [int(scan.split(prefix)[1]) for scan in list(f.keys())]
        assert scan in sl, f'scan {scan} not found in file {filepath}\n\nScan available: {sl}'

        ###########
        # verbose #
        ###########
        if verbose: 
            print(f'loading scan: {scan}')

        ########
        # RIXS #
        ########
        if prefix == 'acq':
            ##################
            # get scan entry #
            ##################
            h =  f[prefix+str(scan)]
    
            #############
            # read data #
            #############
            points = h['Instrument']['DLD8080']['data']['points'][:]
            x      = points[:, 0]  # detector x position of a photon hit
            y      = points[:, 1]  # detector y position of a photon hit
            t      = points[:, 2]  # time of a photon hit relative to bunch period
            t2     = points[:, 3]  # absolute time of a photon hit

            #############################
            # Create PhotonEvent object #
            #############################
            pe = br.PhotonEvents(x=x, y=y)
            try:
                pe.time_bunch    = t - min(t)  # make it start from zero
                pe.time_absolute = t2
                pe.exposure_time_calculated = max(t2)-min(t2)
            except ValueError:
                pe.time_bunch    = t
                pe.time_absolute = t2
                pe.exposure_time_calculated = 0
            
            #########
            # attrs #
            #########
            metadata = _sort_metadata(attrs=rixs_attrs, f=h, verbose=verbose)
            for attr in metadata:
                setattr(pe, attr, metadata[attr])

            # scan
            pe.scan = scan

            # undulator
            if 'epu_r3_316_phase' in pe and 'epu_r3_316_gap' in pe:
                pe.pol = _pol(pe.epu_r3_316_phase, pe.epu_r3_316_gap)


            # post processing
            # pe.cutoff       = None
            # pe.time         = None
            # pe.tmask        = None
            # pe.folded_time  = None
            # pe.folded_tmask = None
            # pe.nbunchs      = None
            # pe.period       = None
            # pe.mask         = None

            ###########
            # verbose #
            ###########
            if verbose:
                print(f'scan type: {scan} - RIXS')
            
            return pe
    
        #######
        # XAS #
        #######
        if prefix == 'entry': 
            ##################
            # get scan entry #
            ##################
            h  = f[prefix+str(scan)]
            h2 = h['measurement']
            
            ######################
            # get scanned motors #
            ######################
            motors = []
            for key in h2:
                if key.startswith('aemexp2_') == False and key not in ('Pt_No', 'dt', 'ux_ct', 'pre_scan_snapshot'):
                    motors.append(key)

            ############################
            # get x axis and scan type #
            ############################
            if len(motors) == 0:
                warnings.warn(f'cannot find motors for this scan. Using generic motor positions')
                # x axis
                x = h2['Pt_No'][()]
                # scan type
                scan_type = 'None'
            elif len(motors) == 1:
                # x axis
                x = h2[motors[0]][()]
                # scan type
                if motors[0] == 'beamline_energy':
                    scan_type = 'XAS'
                elif motors[0] == 'a_mp1_x':
                    scan_type = 'sample scan x'
                elif motors[0] == 'a_mp1_y':
                    scan_type = 'sample scan y'
                elif motors[0] == 'a_mp1_z':
                    scan_type = 'sample scan z'
                else:
                    scan_type = motors[0]
            elif len(motors) == 2:
                # x axis
                x = h2['Pt_No'][()]
                # scan type
                scan_type = f'mesh scan {motors}'
            else:
                command = h['title'][()].decode("utf-8")
                if command.startswith('ascan'):
                    motors = [command.split(' ')[1], ]

                    # x axis
                    x = h2[motors[0]][()]
                    # scan type
                    if motors[0] == 'beamline_energy':
                        scan_type = 'XAS'
                    elif motors[0] == 'a_mp1_x':
                        scan_type = 'sample scan x'
                    elif motors[0] == 'a_mp1_y':
                        scan_type = 'sample scan y'
                    elif motors[0] == 'a_mp1_z':
                        scan_type = 'sample scan z'
                    else:
                        scan_type = motors[0]
                else:
                    # x axis
                    x = h2['Pt_No'][()]
                    # scan type
                    scan_type = 'None'
                    warnings.warn(f'cannot recognize command for scan {scan}: {command}')

            #########################
            # create Spectra object #
            #########################
            RMU = br.Spectrum(x=x, y=h2['aemexp2_ch1'][:] / h2['aemexp2_timer'][:])
            TFY = br.Spectrum(x=x, y=h2['aemexp2_ch2'][:] / h2['aemexp2_timer'][:])
            EXF = br.Spectrum(x=x, y=h2['aemexp2_ch3'][:] / h2['aemexp2_timer'][:])
            TEY = br.Spectrum(x=x, y=h2['aemexp2_ch4'][:] / h2['aemexp2_timer'][:])   
            ss  = br.Spectra(TEY, TFY, EXF, RMU)

            #####################
            # reshape mesh scan #
            #####################
            if len(motors) == 2:
                # get which motor was scanned first
                # motor a is frozen while b is scanned
                command = h['title'][()].decode("utf-8")
                if command.find(motors[1]) > command.find(motors[0]):
                    motors = motors[0], motors[1]
                else:
                    motors = motors[1], motors[0]
                a, b = h2[motors[0]], h2[motors[1]]

                # find `b` motor points
                # `b` motor is frozen while `a` is moving. Therefore, the position of `b` will be like a 
                # 'stair' (see drawing at the beginning of this file), since the position won't 
                # change for multiple collected data points. We can then calculate a `threshold`
                # which above that we consider that `b` is at it's next step. With that we can 
                # infer all the `b` points used
                b1     = np.diff(b)
                bfinal = []
                if b[-1] > b[0]: # if `b` is going up
                    threshold = max(b1)/2
                    temp = []
                    for _b, _b1 in zip(b, b1):
                        temp.append(_b)
                        if _b1 > threshold:
                            bfinal.append(np.mean(temp))
                            temp = []
                    bfinal.append(np.mean(temp))
                else:  # if `b` is going down
                    threshold = min(b1)/2
                    temp = []
                    for _b, _b1 in zip(b, b1):
                        temp.append(_b)
                        if _b1 < threshold:
                            bfinal.append(np.mean(temp))
                            temp = []
                    bfinal.append(np.mean(temp))
                
                # find `a` motor points (method 2)
                try:
                    afinal = np.mean([_ if i%2 ==0 else _[::-1] for i, _ in enumerate(np.split(a, len(bfinal)))], axis=0)
                except ValueError:
                    raise ValueError(f'reshape failed for mesh {scan} ({command}) [{motors[1]}={len(bfinal)}, len({motors[0]})={len(a)}, len({motors[0]})/len({motors[1]})={len(a)/len(bfinal)}]')

                # reshape data
                ss = br.Dummy(TEY, TFY, EXF, RMU)
                for j, s in enumerate(ss):
                    ss[j] = br.Image(data=[row if i%2 == 0 else row[::-1] for i, row in enumerate(s.y.reshape(len(bfinal), len(afinal)))])
                    ss[j].x_centers = afinal
                    ss[j].y_centers = bfinal

            #############
            # raw attrs #
            #############   
            metadata = _sort_metadata(attrs=xas_attrs, f=h, verbose=verbose)
            for s in [_ for _ in ss] + [ss]:
                for attr in metadata:
                    setattr(s, attr, metadata[attr])

                # # other measurement attrs
                # if 'measurement' in h:
                #     for key in h['measurement']:
                #         if key != 'pre_scan_snapshot' and key.startswith('aemexp2_') == False:
                #             setattr(s, 'measurement_' + key, h['measurement'][key][()])

            ################
            # pretty attrs #
            ################
            for s in [_ for _ in ss] + [ss]:
                s.scan           = scan
                s.scanned_motors = motors
                s.scan_type      = scan_type

                # Not in metadata
                s.E   = None
                s.T   = None
                s.Pol = None

                # motors
                if hasattr(s, 'pre_scan_a_mp1_x'):   s.sample_x = s.pre_scan_a_mp1_x
                if hasattr(s, 'pre_scan_a_mp1_y'):   s.sample_y = s.pre_scan_a_mp1_y
                if hasattr(s, 'pre_scan_a_mp1_z'):   s.sample_z = s.pre_scan_a_mp1_z
                if hasattr(s, 'pre_scan_a_mp1_yaw'): s.th       = round(s.pre_scan_a_mp1_yaw, 2)
            
            ss.header = ('TEY', 'TFY', 'EXF', 'RMU')
            RMU.label = 'RMU'
            TFY.label = 'TFY'
            TEY.label = 'TEY'
            EXF.label = 'EXF'

            ###########
            # verbose #
            ###########
            if verbose:
                print(f'scan type: {scan} - {scan_type}')
            
            return ss



# %%

# %% script for mesh
# find `a` motor points - method 1 - independent of `b` - WORKS
# a1        = np.diff(a)
# # find length
# alength = 0
# if a1[0] > 0:
#     for _a1 in a1:
#         alength += 1
#         if _a1 < 0:
#             alength -= 1
#             break
#
# afinal    = np.zeros(alength)
# temp      = []
# going_up = True
# number_of_zig_zags = 0
# for _a, _a1 in zip(a, a1):
#     if going_up:
#         threshold = max(a1)/2
#         if _a1 > threshold:
#             temp.append(_a)
#         else:
#             temp.append(_a)
#             afinal += np.array(temp)
#             going_up = False
#             number_of_zig_zags += 1
#             temp = []
#     else:
#         threshold = min(a1)/2
#         if _a1 < threshold:
#             temp.append(_a)
#         else:
#             temp.append(_a)
#             afinal += np.array(temp[::-1])
#             going_up = True
#             number_of_zig_zags += 1
#             temp = []
# afinal = afinal/number_of_zig_zags