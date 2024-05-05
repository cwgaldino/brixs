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


RIXS metadata
{'External': <Closed HDF5 group>,
 'Instrument': <Closed HDF5 group>,

 'data': <Closed HDF5 group>,

 'endmetadata': <Closed HDF5 group>,
 'startmetadata': <Closed HDF5 group>}

 
 'start_time': <Closed HDF5 dataset>,
 'end_time': <Closed HDF5 dataset>,


Last edited: Carlos Galdino 2024-04-27
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from collections.abc import MutableSequence
from collections.abc import Iterable
from pathlib import Path
import numpy as np
import datetime
import warnings
import copy

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br
# import matplotlib.pyplot as plt
# %% ------------------------- Special Imports ---------------------------- %% #
import h5py
# %%

# %% ============================ support ================================= %% #
def tree(f, pre=''):
    """Returns a text with the structure of a hdf5 file

    Usage:
        >>> import brixs as br
        >>> import h5py
        >>> 
        >>> with h5py.File(<filepath>, 'r') as f:
        >>>     text = tree(f)
        >>>
        >>> with h5py.File(<filepath>, 'r') as f:
        >>>     text = tree(f[entry])
        >>> 
        >>> print(text)
        >>>
        >>> br.save_text(text, 'text.txt')
    
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
                text += tree(f, pre + '    ')
            else:
                try:
                    text += pre + '└── ' + key + ' (%d)' % len(f) + '\n'
                except TypeError:
                    text += pre + '├── ' + key + ' (scalar)' + '\n'
        else:
            if type(f) == h5py._hl.group.Group:
                text += pre + '├── ' + key + '\n'
                text += tree(f, pre+'│   ')
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

def scanlist(filepath):
    """Return list of scans available in filepath"""
    with h5py.File(Path(filepath), 'r') as f:
        prefix = ''.join([i for i in list(f.keys())[0] if not i.isdigit()])
        return [int(_) for _ in np.sort([int(scan.split(prefix)[1]) for scan in list(f.keys())])]

def metadata(attrs, f, verbose=True):
    for _type in attrs:
        for name in attrs[_type]:
            address = attrs[_type][name]
            try: 
                if _type == 'string':   
                    attrs[_type][name] = f[address][()].decode("utf-8")
                elif _type == 'number':   
                    attrs[_type][name] = f[address][()]
                elif _type == 'round':    
                    attrs[_type][name] = round(f[address][()], 2)
                elif _type == 'datetime': 
                    attrs[_type][name] = _str2datetime(f[address][()].decode("utf-8"))
                elif _type == 'bool':     
                    attrs[_type][name] = f[address][()][0] == 1
            except Exception as e:
                if verbose: print(e)
                attrs[_type][name] = None

class Images():
    
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

        ########
        # RIXS #
        ########
        if prefix == 'acq':

            #################
            # metadata list #
            #################
            attrs = {'string':   {'script_name':      'entry/current_script_name',
                                  'command':          'entry/diamond_scan/scan_command',
                                  'facility_user_id': 'entry/user01/facility_user_id',
                                  'username':         'entry/user01/name',
                                  'beamline':         'entry/instrument/beamline',
                                  'pol':              'entry/instrument/id/polarisation'},
                    'number':   {'sample_x':         'entry/instrument/manipulator/x',
                                'sample_y':         'entry/instrument/manipulator/y',
                                'sample_z':         'entry/instrument/manipulator/z'},
                    'mean':     {},
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
        
            #############
            # raw attrs #
            #############
            # external
            pe.a_mp1_x           = np.mean(h['External']['a_mp1_x']['position'][:][:, 1])
            pe.a_mp1_y           = np.mean(h['External']['a_mp1_y']['position'][:][:, 1])
            pe.a_mp1_z           = np.mean(h['External']['a_mp1_z']['position'][:][:, 1])
            pe.a_mp1_yaw         = np.mean(h['External']['a_mp1_yaw']['position'][:][:, 1])
            pe.temperature       = np.mean(h['External']['b316a-o01']['dia']['tco-02']['Temperature'][:][:, 1])
            pe.beamline_energy   = h['External']['beamline_energy']['position'][:][:, 1]
            pe.epu_r3_316_gap    = np.mean(h['External']['epu_r3_316_gap']['position'][:][:, 1])
            pe.epu_r3_316_phase  = np.mean(h['External']['epu_r3_316_phase']['position'][:][:, 1])
            pe.mono_energy_calib = np.mean(h['External']['mono_energy_calib']['position'][:][:, 1])
            pe.veritasarm_energy = np.mean(h['External']['veritasarm_energy']['position'][:][:, 1])
        
            # time
            pe.exposure_time    = np.mean(h['Instrument']['DLD8080']['exposure_time'][()])          
            pe.start_time       = _str2datetime(h['start_time'][()].decode("utf-8"))
            pe.end_time         = _str2datetime(h['end_time'][()].decode("utf-8"))

            # start/end metadata
            
            for name in ['startmetadata', 'endmetadata']:
                # print(name)
                h2 =  h[name]
                for key in h2:
                    # print(key)
                    new_attr_name = name + '_' + key
                    if key == 'Sample':
                        setattr(pe, new_attr_name, h2[key][()].decode("utf-8"))
                    elif key == 'mono_energy':
                        setattr(pe, new_attr_name, h2[key]['cffreal'][()])
                    elif key == 'b316a-o01':
                        new_attr_name = name + '_' + 'temperature'
                        setattr(pe, new_attr_name, h2['b316a-o01/dia/tco-02/Temperature'][()])
                    elif key == 'veritasarm_energy':
                        pass
                    else:
                        setattr(pe, new_attr_name, h2[key]['position'][()])

            ################
            # pretty attrs #
            ################
            # scan
            pe.scan     = scan

            # motors
            pe.sample_x = copy.deepcopy(pe.a_mp1_x)
            pe.sample_y = copy.deepcopy(pe.a_mp1_y)
            pe.sample_z = copy.deepcopy(pe.a_mp1_z)
            pe.th       = round(pe.a_mp1_yaw - 261, 2)

            # temperature
            pe.T        = round(pe.temperature + 273.15, 2)

            # undulator
            pe.E        = round(np.mean(pe.beamline_energy), 2)
            if copy.deepcopy(pe.epu_r3_316_phase) > 23:
                pe.pol = 'LV'
            else:
                pe.pol = 'LH'

            # time
            pe.time_bunch    = t - min(t)  # make it start from zero
            pe.time_absolute = t2
            pe.exposure_time_calc = max(t2)-min(t2)

            # post processing
            pe.cutoff       = None
            pe.time         = None
            pe.tmask        = None
            pe.folded_time  = None
            pe.folded_tmask = None
            pe.nbunchs      = None
            pe.period       = None
            pe.mask         = None

            ###########
            # verbose #
            ###########
            if verbose:
                print(f'loading scan: {scan} - RIXS')
            
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
                ss = Images(TEY, TFY, EXF, RMU)
                for j, s in enumerate(ss):
                    ss[j] = br.Image(data=[row if i%2 == 0 else row[::-1] for i, row in enumerate(s.y.reshape(len(bfinal), len(afinal)))])
                    ss[j].x_centers = afinal
                    ss[j].y_centers = bfinal

            #############
            # raw attrs #
            #############   
            for s in [_ for _ in ss] + [ss]:
                # bytes
                if 'definition' in h:       s.definition       = h['definition'][()].decode("utf-8")
                if 'entry_identifier' in h: s.entry_identifier = h['entry_identifier'][()].decode("utf-8")
                if 'program_name' in h:     s.program_name     = h['program_name'][()].decode("utf-8")
                if 'title' in h:            s.command          = h['title'][()].decode("utf-8")
                if 'user' in h: 
                    if 'name' in h['user']:   s.user             = h['user/name'][()].decode("utf-8")

                # other measurement attrs
                if 'measurement' in h:
                    for key in h['measurement']:
                        if key != 'pre_scan_snapshot' and key.startswith('aemexp2_') == False:
                            setattr(s, 'measurement_' + key, h['measurement'][key][()])

                # measurement snapshot
                h2 =  h['measurement/pre_scan_snapshot']
                for key in h2:
                    new_attr_name = 'pre_scan_' + key
                    setattr(s, new_attr_name, h2[key][()])

                # plot 2 (whatever this is)
                if 'plot_2' in h:
                    for key in h['plot_2']:
                        new_attr_name = 'plot_2_' + key
                        s.plot_2 = h['plot_2'][key][()]
                
                # plot 1 
                if 'plot_1' in h:
                    for key in h['plot_1']:
                        new_attr_name = 'plot_1_' + key
                        setattr(s, new_attr_name, h['plot_1'][key][()])

                # time
                if 'start_time' in h:       s.start_time       = _str2datetime(h['start_time'][()].decode("utf-8"))
                if 'end_time' in h:         s.end_time         = _str2datetime(h['end_time'][()].decode("utf-8"))
                if 'macro_start_time' in h: s.macro_start_time = _str2datetime(h['macro_start_time'][()].decode("utf-8"))

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
                print(f'loading scan: {scan} - {scan_type}')
            
            return ss



# %%

# script for mesh
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