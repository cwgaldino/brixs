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
import matplotlib.pyplot as plt
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
    """for Cr edge"""
    if gap > 20 and gap < 30 and phase < 2 and phase > -2:
        return 'LH'
    elif gap > 20 and gap < 30 and phase < 30 and phase > 20:
        return 'LV'
    elif gap > 12 and gap < 20 and phase > 14 and phase < 20:
        return 'CR'
    elif gap > 12 and gap < 20 and phase < -14 and phase > -20:
        return 'CN'
    else:
        return 'do not know'
        
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
                              'tth': 'endmetadata/q_angle/position'
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
                      'exit_slit':    'measurement/pre_scan_snapshot/a_slit1_hz',
                      'temperature1': 'measurement/pre_scan_snapshot/lakeshore_1',
                      'temperature2': 'measurement/pre_scan_snapshot/lakeshore_2',
                      'energy':       'measurement/pre_scan_snapshot/beamline_energy',
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
        Spectra (TEY, MCP, TFY, RMU) for XAS and sample scans
        Image list (TEY, MCP, TFY, RMU) for mesh scans
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
                pe.pol_calculated = _pol(pe.epu_r3_316_phase, pe.epu_r3_316_gap)


            # post processing
            pe.cutoff       = None
            pe.time_bunch_histogram         = None
            pe.time_bunch_histogram_nbins        = None
            pe.folded_time_bunch  = None
            pe.folded_tmask = None
            pe.tmask = None
            pe.nbunchs      = None
            pe.period       = None
            pe.mask         = None

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

            ##################################
            # get command and scanned motors #
            ##################################
            command = h['title'][()].decode("utf-8")
            if command.startswith('ascan '):
                motors = [command.split(' ')[1], ]
            elif command.startswith('mesh '):
                motors = [command.split(' ')[1], command.split(' ')[5],]
            else: 
                raise ValueError(f'can only read `ascan` and `mesh`. command not recognized: `{command}`')

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
                # command = h['title'][()].decode("utf-8")
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
                bfinal = np.linspace(float(command.split(' ')[6]), float(command.split(' ')[7]), int(command.split(' ')[8]) + 1)
                
                # find `a` motor points (method 2)
                afinal = np.linspace(float(command.split(' ')[2]), float(command.split(' ')[3]), int(command.split(' ')[4]) + 1)

                # reshape data
                ss = br.Dummy(TEY, TFY, EXF, RMU)
                for j, s in enumerate(ss):
                    # ss[j] = br.Image(data=[row if i%2 == 0 else row[::-1] for i, row in enumerate(s.y.reshape(len(afinal), len(bfinal)))])
                    # print(len(s.y))
                    # ss[j] = br.Image(data=[row if i%2 == 0 else row[::-1] for i, row in enumerate(s.y.reshape(int(command.split(' ')[4])+1, int(command.split(' ')[8])+1))])
                    ss[j] = br.Image(data=[row if i%2 == 0 else row[::] for i, row in enumerate(s.y.reshape(int(command.split(' ')[4])+1, int(command.split(' ')[8])+1))])
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



# %% support
def _clip(self, mask):
    """Clip photon events.

    Usage:
        >>> pe.clip((xmin, xmax, ymin, ymax))
        >>> pe.clip((0, 12, 3, 12))
        >>> pe.clip([(0, 12, 3, 12), (1, 4, 15, 18)])

    Args:
        mask (list): list with rectangles coordinates (x_start, x_stop, y_start, y_stop).

    Returns:
        PhotonEvents
    """
    # assert mask is the right format
    assert isinstance(mask, br.Iterable), 'mask must be iterable'
    if len(mask) == 4:
        if isinstance(mask[0], br.Iterable) == False:
            mask = [mask, ]

    ########
    # clip #
    ########
    x = []
    y = []
    time_bunch = []
    time_absolute = []
    for r in mask:
        temp = np.array([(x, y, t, t2) for x, y, t, t2 in zip(self.x, self.y, self.time_bunch, self.time_absolute) if ((x > r[0] and x < r[1]) and (y > r[2] and y < r[3]))])
        x += list(temp[:, 0])
        y += list(temp[:, 1])
        time_bunch    += list(temp[:, 2])
        time_absolute += list(temp[:, 3])


    #########
    # final #
    #########
    pe = br.PhotonEvents(x=x, y=y)

    ##################
    # transfer attrs #
    ##################
    pe.copy_lims_from(self)
    pe.copy_attrs_from(self)

    pe.time_absolute = time_absolute
    pe.time_bunch = time_bunch

    return pe
br.PhotonEvents.clip2 = _clip

def _find_empty(self, mask, axis=0):
    """Find pixel cols or rows that are outside of mask
    
    axis = 0 : looking for empty rows
    axis = 1 : looking for empty cols
    
    Returns:
        list
    """
    if axis == 1:
        x = self.x_centers
    if axis == 0:
        x = self.y_centers
    else:
        raise ValueError('axis must be 0 or 1')
        
    empty = [True]*len(x)
    for i, x in enumerate(x):
        for r in mask:
            if (x > r[2] and x < r[3]):
                empty[i] = False
                
    return empty
br.Image.find_empty = _find_empty

def _get_filtered_spectra(self, mask, axis=0):
    empty = self.find_empty(mask=mask, axis=axis)
    
    assert axis == 1 or axis == 0, 'axis must be 0 or 1'
    if axis == 0:
        ss = self.rows
        x  = self.y_centers
    else:
        ss = self.columns
        x  = self.x_centers
    
    filtered = br.Spectra()
    for i, is_empty in enumerate(empty):
        if not is_empty:
            s          = ss[i]
            s.x_center = x[i]
            filtered.append(s)
        
    filtered._x_centers = [s.x_center for s in filtered]
    return filtered
br.Image.get_filtered_spectra = _get_filtered_spectra


# %% ========================== Time support ============================== %% #
def _apply_time_cutoff(self, value=3e7):
    # cutoff time
    indexes = np.array([i for i in range(len(self.time_bunch)) if self.time_bunch[i] > value]) # remove this one last point
    final  = [True]*len(self.time_bunch)
    for i in indexes:
        final[i] = False
    _cutoff = [not elem for elem in final]
    
    # saving
    cutoff = br.PhotonEvents(x=self.x[_cutoff], y=self.y[_cutoff])
    cutoff.copy_attrs_from(self)
    cutoff.time_absolute = self.time_absolute[_cutoff]
    
    good = br.PhotonEvents(x=self.x[final], y=self.y[final])
    good.copy_attrs_from(self)
    good.time_absolute = self.time_absolute[final]
    good.time_bunch    = self.time_bunch[final]


    return good, cutoff
br.PhotonEvents.apply_time_cutoff = _apply_time_cutoff    

def _calculate_time_bunch_histogram(self, nbins=10000):
    # calculate time histogram
    hist, bin_edges = np.histogram(self.time_bunch, nbins)
    time       = br.Spectrum(x=br.moving_average(bin_edges, 2), y=hist)
    # self.time.shift = -min(self.time.x)
    
    self.time_bunch_histogram       = time
    self.time_bunch_histogram_nbins = nbins
    return
br.PhotonEvents.calculate_time_bunch_histogram = _calculate_time_bunch_histogram        

def _calculate_folded_time_bunch(self, period=1459, offset=None):
    
    assert self.time_bunch_histogram is not None, 'time is not binned. Please use calculate_time_bunch_histogram'
    
    if offset is None:
        offset = min(self.time_bunch_histogram.x)
    x      = (self.time_bunch_histogram.x-offset)%period
    
    hist, bin_edges  = np.histogram(x, bins=int(period/10), weights=self.time_bunch_histogram.y)
    folded_time_bunch = br.Spectrum(x=br.moving_average(bin_edges, 2), y=hist)
    # folded_t        = (self.t-min(self.t))%period  # min(self.t) should not be necessary because this is already set to zero during read()
    
    self.folded_time_bunch = folded_time_bunch
    self.nbunchs = int(round(max(self.time_bunch_histogram.x)/period))
    self.period  = period
    
    return
br.PhotonEvents.calculate_folded_time_bunch = _calculate_folded_time_bunch


def _calculate_tmask(self, w=300, c='gauss'):
    """calculate tmask from folded time."""
    
    assert self.folded_time_bunch is not None, 'folded time must be defined. Please use pe.calculate_folded_time_bunch()'
    
    # get folded_time peak center
    if c == 'gauss':
        fit, popt, err, f = self.folded_time_bunch.fit_peak()
        c = popt[1]
    elif c == 'max':
        c = self.folded_time_bunch.x[np.argmax(self.folded_time_bunch.y)]
    
    # tmask
    right = max(self.folded_time_bunch.x)
    left  = min(self.folded_time_bunch.x)
    
    assert c >= left and c <= right, 'c outside range'
    
    if c+w/2 > right: 
        rest = w - (right - (c-w/2))
        folded_tmask = [(left, left+rest), (c-w/2, right)]
    elif c-w/2 < left: 
        rest = w - ((c+w/2) - left)
        folded_tmask = [(left, c+w/2), (right-rest, right)]
    else:
        folded_tmask = [(c-w/2, c+w/2), ]
    tmask        = [(c1-w/2+c, c1+w/2+c) for c1 in np.arange(min(self.time_bunch_histogram.x), max(self.time_bunch_histogram.x), self.period)][:-1]
    
    self.folded_tmask = folded_tmask
    self.tmask        = tmask
    
    return
br.PhotonEvents.calculate_tmask = _calculate_tmask   

def _plot_folded_tmask(self, ax=None, **kwargs):
    
    if ax is None:
        ax = plt.gca()
    
    if 'edgecolor' not in kwargs:
        kwargs['edgecolor'] = 'red'
    if 'fill' not in kwargs:
        kwargs['fill'] = False
    if 'lw' not in kwargs and 'linewidth' not in kwargs:
        kwargs['lw'] = 3
    # if 'height' not in kwargs:
    #     kwargs['height'] = max(self.folded_time_bunch.y)
    
    # folded time
    assert self.folded_tmask is not None, 'cannot find folded tmask'
    for m in self.folded_tmask:
        rect = br.rectangle((m[0], m[1]), (0, max(self.folded_time_bunch)), ax=ax, **kwargs)
        # rect = Rectangle((m[0], 0), m[1]-m[0], **kwargs)
        # ax.add_patch(rect)  
    
    return
br.PhotonEvents.plot_folded_tmask = _plot_folded_tmask   

def _plot_tmask(self, ax=None, **kwargs):
    
    if ax is None:
        ax = plt.gca()
        
    if 'edgecolor' not in kwargs:
        kwargs['edgecolor'] = 'red'
    if 'fill' not in kwargs:
        kwargs['fill'] = False
    if 'lw' not in kwargs and 'linewidth' not in kwargs:
        kwargs['lw'] = 3
    # if 'height' not in kwargs:
    #     kwargs['height'] = max(self.time.y)/3
    
    # folded time
    assert self.tmask is not None, 'cannot find tmask'

    for m in self.tmask:
        rect = br.rectangle((m[0], m[1]), (0, max(self.time_bunch_histogram)), ax=ax, **kwargs)
        # rect = Rectangle((m[0], 0), m[1]-m[0], **kwargs)
        # ax.add_patch(rect)  
    
    return
br.PhotonEvents.plot_tmask = _plot_tmask   


def _get_count_per_window(self):
    centers = []
    counts = []
    for m in self.tmask:
        centers += [m[0] + (m[1] - m[0])/2]
        counts += [sum([1 for t in self.time_bunch if t > m[0] and t < m[1]])]            
    return br.Spectrum(centers, counts)
br.PhotonEvents.get_count_per_window = _get_count_per_window

def _apply_tmask(self):
    """tmask is applied on the regular time (not folded time)
    
    folded_time doesn't even need to be defined
    """
    
    indexes = []
    for m in self.tmask:
        indexes += [i for i, _t in enumerate(self.time_bunch) if _t>=m[0] and _t<=m[1]]
    # print(len(indexes))
    final = [False]*len(self.time_bunch)

    for i in indexes:
        final[i] = True
    opposite = [not elem for elem in final]

    good = br.PhotonEvents(x=self.x[final],    y=self.y[final])
    bad  = br.PhotonEvents(x=self.x[opposite], y=self.y[opposite])

    good.percentage = round(len(good.x)/len(self.x)*100, 2)
    bad.percentage  = round(len(bad.x)/len(self.x)*100, 2)
    
    # transfer attrs
    for attr in self.get_attrs():
        value = copy.deepcopy(self.__dict__[attr])
        good.__setattr__(attr, value)
        bad.__setattr__(attr, value)

    
    # good.folded_t = self.folded_t[final]
    # bad.folded_t  = self.folded_t[opposite]
    good.time_bunch = np.array(self.time_bunch)[final]    #- min(self.t[final])    # set min value to zero
    bad.time_bunch  = np.array(self.time_bunch)[opposite] #- min(self.t[opposite]) # set min value to zero
                    
    for pe in (good, bad):
        pe.cutoff  = None
        pe.time    = None
        pe.tmask   = None
        pe.folded_time  = None
        pe.folded_tmask = None
        pe.nbunchs      = self.nbunchs
    
    return good, bad
br.PhotonEvents.apply_tmask = _apply_tmask
# %%

# %% verify and process
def verify(filepath, scan, mask, tcutoff=3e7, tnbins=10000, period=1458, twidth=320, tcenter='max', curv='self', curv_nbins=(20, 1000), sbins=1200, calib=None, **kwargs):
    pe = read(filepath, scan)

    # pe = copy.deepcopy(pe_safe)
    if tcutoff is not None:
        pe, cutoff = pe.apply_time_cutoff(value=tcutoff)

    # mask ===========================================
    fig, axes = br.subplots(2, 4, figsize=(12, 4))
    plt.subplots_adjust(left=0.1, right=0.99)
    # fig.suptitle('Physical mask')

    # axes 0
    axes[0].set_title('1: Raw (pe1)')
    pe.plot(axes[0], color='black', **kwargs)
    if tcutoff is not None:
        cutoff.plot(axes[0], color='magenta', s=.6)

    # plot mask
    for m in mask:
        br.rectangle(m[:2], m[2:], ax=axes[0], lw=2, edgecolor='red')

    # axes 1
    axes[1].set_title('2: After mask (pe2)')
    pe2 = pe.clip2(mask)
    pe2.mask = copy.deepcopy(mask)
    pe2.plot(axes[1], color='firebrick', **kwargs)

    for ax in axes:
        ax.set_xlabel('x (mm)')
        ax.set_ylabel('y (mm)')

    # time ============================================
    # folded time
    axes[4].set_title('3a: Folded Time')
    pe2.calculate_time_bunch_histogram(nbins=tnbins)
    pe2.calculate_folded_time_bunch(period=period)
    try:
        pe2.calculate_tmask(w=twidth, c=tcenter)
    except AssertionError:
        print('error finding suitable c for time window')
        pe2.calculate_tmask(w=twidth, c='max')
    pe2.folded_time_bunch.plot(axes[4])
    pe2.plot_folded_tmask(axes[4], lw=3, edgecolor='red')
    axes[4].set_xlabel('Folded time (ps)')
    axes[4].set_ylabel('Photon count')

    # unfolded time
    for ax in (axes[5], axes[6]):
        ax.set_title('3b: Unfolded Time')
        pe2.time_bunch_histogram.plot(ax, color='black', lw=.2, marker='o', ms=1)
        pe2.plot_tmask(ax)
        ax.set_xlabel('Time (ps)')
        ax.set_ylabel('Photon count')
    br.zoom(min(pe2.time_bunch_histogram.x), period*5, axes[5])

    # photons per time window
    axes[7].set_title('3b: Photons per time window')
    s = pe2.get_count_per_window()
    s.plot(axes[7], color='black', marker='o', ms=2)
    axes[7].set_xlabel('Time (ps)')
    axes[7].set_ylabel('Photon count per window')

    # time rejected photons
    axes[2].set_title('4: Rejected photons')
    pe3, bad = pe2.apply_tmask()
    bad.plot(axes[2], color='green', **kwargs)
    axes[2].set_xlabel('x (mm)')
    axes[2].set_ylabel('y (mm)')

    # curvature
    if curv is not None and curv != False:
        axes[3].set_title('5: curvature (pe3)')
        if curv == 'self':
            im = pe3.binning(nbins=curv_nbins)
            ss = im.get_filtered_spectra(mask=pe3.mask, axis=0)
            ss.calculate_shift()
            s_ = br.Spectrum(x=ss._x_centers, y=ss.calculated_shift)
            curv, _, _ = s_.polyfit(2)
        pe4 = copy.deepcopy(pe3)
        pe4.set_shift(p=curv, axis=1)
        pe4.plot(axes[3], color='red', **kwargs)
    else:
        pe4 = None
        
    # spectrum
    if curv is not None and isinstance(curv, br.Iterable):
        print(pe.exposure_time)
        s = pe4.calculate_spectrum(axis=0, nbins=sbins)
        s.set_factor(1/sum([m[3]-m[2] for m in mask]))
        s.set_factor(1/s.exposure_time)
        s.set_factor(sbins)
        # s.t = None

        s.scan = scan
        # # intensity calibration correction ====================
        # if calib_intensity is not None:
        #     model = lambda x: np.polyval(calib_intensity, x)
        #     newy = np.zeros(len(s))
        #     for i in range(len(s)):
        #         newy[i] =  s.y[i]/model(s.x[i])
        #     s.y = newy
    else:
        s = None
        
    # calib ==============================================
    if calib is not None and s is not None:
        s.calib = calib
        s.shift = -s.E
    return pe, pe2, pe3, pe4, curv, s

@br.finder.track
def process(filepath, scan, mask, tcutoff=3e7, tnbins=10000, period=1458, twidth=320, tcenter='max', curv='self', curv_nbins=(20, 1000), sbins=1200, calib=None, **kwargs):
    pe = read(filepath, scan)

    # pe = copy.deepcopy(pe_safe)
    if tcutoff is not None:
        pe, cutoff = pe.apply_time_cutoff(value=tcutoff)

    # mask
    pe2 = pe.clip2(mask)
    pe2.mask = copy.deepcopy(mask)


    # folded time
    pe2.calculate_time_bunch_histogram(nbins=tnbins)
    pe2.calculate_folded_time_bunch(period=period)
    try:
        pe2.calculate_tmask(w=twidth, c=tcenter)
    except AssertionError:
        print('error finding suitable c for time window')
        pe2.calculate_tmask(w=twidth, c='max')

    # time rejected photons
    pe3, bad = pe2.apply_tmask()

    # curvature
    if curv is not None and curv != False:
        if curv == 'self':
            im = pe3.binning(nbins=curv_nbins)
            ss = im.get_filtered_spectra(mask=pe3.mask, axis=0)
            ss.calculate_shift()
            s_ = br.Spectrum(x=ss._x_centers, y=ss.calculated_shift)
            curv, _, _ = s_.polyfit(2)
        pe4 = copy.deepcopy(pe3)
        pe4.set_shift(p=curv, axis=1)
    else:
        pe4 = None
        
    # spectrum
    if curv is not None and isinstance(curv, br.Iterable):
        print(pe.exposure_time)
        s = pe4.calculate_spectrum(axis=0, nbins=sbins)
        s.set_factor(1/sum([m[3]-m[2] for m in mask]))
        s.set_factor(1/s.exposure_time)
        s.set_factor(sbins)
        # s.t = None

        s.scan = scan
    else:
        s = None
        
    # calib ==============================================
    if calib is not None and s is not None:
        s.calib = calib
        s.shift = -s.E
    return s
# %%


# %%
# %% alternative script for mesh
# Find b - method 1
# b1     = np.diff(b)
# bfinal = []
# if b[-1] > b[0]: # if `b` is going up
#     threshold = max(b1)/2
#     temp = []
#     for _b, _b1 in zip(b, b1):
#         temp.append(_b)
#         if _b1 > threshold:
#             bfinal.append(np.mean(temp))
#             temp = []
#     bfinal.append(np.mean(temp))
# else:  # if `b` is going down
#     threshold = min(b1)/2
#     temp = []
#     for _b, _b1 in zip(b, b1):
#         temp.append(_b)
#         if _b1 < threshold:
#             bfinal.append(np.mean(temp))
#             temp = []
#     bfinal.append(np.mean(temp))
                
# find `a` motor points (method 2)
# try:
#     afinal = np.mean([_ if i%2 ==0 else _[::-1] for i, _ in enumerate(np.split(a, len(bfinal)))], axis=0)
# except ValueError:
#     raise ValueError(f'reshape failed for mesh {scan} ({command}) [{motors[1]}={len(bfinal)}, len({motors[0]})={len(a)}, len({motors[0]})/len({motors[1]})={len(a)/len(bfinal)}]')

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