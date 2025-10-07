#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from VERITAS beamline of MAX-IV."""

# %% ------------------------- Standard Imports --------------------------- %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import datetime
import copy

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br
import brixs.addons.h5 as h5

# %% ------------------------- Special Imports ---------------------------- %% #
try:
    import h5py
except:
    pass
# %%

# %% ========================= useful functions =========================== %% #
def get_metadata(scan, filepath, entry):
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

def tree(scan, filepath):
    """Returns a text with the structure of a hdf5 file

    Args:
        filepath (str or path): file to hdf5 file
        scan (int): scan number

    Returns:
        str
    """
    with h5py.File(Path(filepath), 'r') as f:
        prefix = ''.join([i for i in list(f.keys())[0] if not i.isdigit()])
        text = h5._tree(f[prefix + str(scan)])
    return text

def scanlist(filepath):
    """Return list of scans available in filepath"""
    with h5py.File(Path(filepath), 'r') as f:
        prefix = ''.join([i for i in list(f.keys())[0] if not i.isdigit()])
        return [int(_) for _ in np.sort([int(scan.split(prefix)[1]) for scan in list(f.keys())])]
# %%

# %% =================== metadata support functions ======================= %% #
def _str2datetime(string):
    """convert VERITAS date/time string pattern to date --> '2022-07-20 21:08:36.921'

    Args:
        string (str): string with VERITAS date string

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
       
def _pol(gap, phase):
    """for Cr edge"""
    if gap > 20 and gap < 30 and phase < 2 and phase > -2:
        return 'LH'
    elif gap > 20 and gap < 30 and phase < 30 and phase > 20:
        return 'LV'
    elif gap > 12 and gap < 20 and phase > 14 and phase < 20:
        return 'C+'
    elif gap > 12 and gap < 20 and phase < -14 and phase > -20:
        return 'C-'
    else:
        return 'do not know' 
# %%

# %% ========================== rixs metadata ============================= %% #
rixs_attrs = {'ignore':{}, 'raw':{}, 'string':{}, 'round2':{},
              'mean':{}, 'min':{}, 'max':{}, 'sigma':{},
              'second_column_mean':{}, 'second_column_min':{}, 'second_column_max':{}, 'second_column_sigma':{}}

h = rixs_attrs['ignore']  # attrs that are added in the post processing -> keep them here because scanlist uses attrs_dict as reference
h['scan']                     = ''
h['pol']                      = ''
h['th']                       = ''
h['exposure_time_calculated'] = ''
h['number_of_photons']        = ''
h['exposure_time_calculated_via_photon_time'] = ''

h = rixs_attrs['raw']
h['exposure_time'] = 'Instrument/DLD8080/exposure_time'

h = rixs_attrs['string']
h['start_time']    = 'start_time'
h['end_time']      = 'end_time'

for _string in ('mean', 'min', 'max', 'sigma'):
    h = rixs_attrs[f'second_column_{_string}']
    if _string == 'mean': _string = ''
    else: _string = '_' + _string
    h['sample_x' + _string]    = 'External/a_mp1_x/position'
    h['sample_y' + _string]    = 'External/a_mp1_y/position'
    h['sample_z' + _string]    = 'External/a_mp1_z/position'
    h['gap' + _string]         = 'External/epu_r3_316_gap/position'
    h['phase' + _string]       = 'External/epu_r3_316_phase/position'
    h['mono_energy' + _string] = 'External/mono_energy_calib/position'
    h['arm_energy' + _string]  = 'External/veritasarm_energy/position'
    h['E' + _string]           = 'External/beamline_energy/position'
    h['T' + _string]           = 'External/b316a-o01/dia/tco-02/Temperature'
    h['th_veritas' + _string]  = 'External/a_mp1_yaw/position'

h = rixs_attrs['round2']
h['exit_slit'] = 'startmetadata/a_slit1_v/position'
h['tth']       = 'endmetadata/q_angle/position'                        
# %%

# %% =========================== xas metadata ============================= %% #
xas_attrs = {'ignore': {}, 'raw':{}, 'string':{}, 'round2':{}}

h = xas_attrs['ignore']
h['scan']           = ''
h['pol']            = ''
h['th']             = ''
h['scanned_motors'] = ''
h['scan_type']      = ''
h['elapsed_time']   = ''

h = xas_attrs['raw']
h['sample_x']      = 'measurement/pre_scan_snapshot/a_mp1_x'
h['sample_y']      = 'measurement/pre_scan_snapshot/a_mp1_y'
h['sample_z']      = 'measurement/pre_scan_snapshot/a_mp1_z'
h['gap']           = 'measurement/pre_scan_snapshot/EPU_R3_316_GAP'
h['phase']         = 'measurement/pre_scan_snapshot/EPU_R3_316_PHASE'
# h['a_m4_lateral']  = 'measurement/pre_scan_snapshot/a_m4_lateral'
# h['a_m4_pitch']    = 'measurement/pre_scan_snapshot/a_m4_pitch'
# h['a_m4_roll']     = 'measurement/pre_scan_snapshot/a_m4_roll'
# h['a_m4_vertical'] = 'measurement/pre_scan_snapshot/a_m4_vertical'
# h['a_m4_yaw']      = 'measurement/pre_scan_snapshot/a_m4_yaw'
# h['m1_yaw']        = 'measurement/pre_scan_snapshot/m1_yaw'
# h['m3_yaw']        = 'measurement/pre_scan_snapshot/m3_yaw'

h = xas_attrs['string']
# h['definition']       = 'definition'
h['entry_identifier'] = 'entry_identifier'
# h['program_name']     = 'program_name'
h['command']          = 'title'
# h['user']             = 'user/name'
h['start_time']       = 'start_time'
h['end_time']         = 'end_time'

h = xas_attrs['round2']
h['exit_slit']     = 'measurement/pre_scan_snapshot/a_slit1_v'
h['E']             = 'measurement/pre_scan_snapshot/beamline_energy'
h['temperature1']  = 'measurement/pre_scan_snapshot/lakeshore_1'
h['temperature2']  = 'measurement/pre_scan_snapshot/lakeshore_2'
h['th_veritas']    = 'measurement/pre_scan_snapshot/a_mp1_yaw'                        
# %%

# %% ============================= read =================================== %% #
def read(scan, filepath, verbose=False):
    """Return data from filepath
    
    RIXS        --> returns PhotonEvent 
    XAS         --> returns Spectra (TEY, MCP, TFY, RMU) 
    sample scan --> returns Spectra (TEY, MCP, TFY, RMU) 
    mesh scan   --> returns Image (TEY, MCP, TFY, RMU) 

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
            h =  f[prefix + str(scan)]

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
            pe.xlim = (0, 8200)
            pe.ylim = (0, 8200)
            try:
                pe.time_bunch    = t - min(t)  # make it start from zero
                pe.time_absolute = t2
                pe.exposure_time_calculated_via_photon_time = int(round(max(t2) - min(t2)))
            except ValueError:
                if verbose: print('cannot zero time-bunch')
                pe.time_bunch    = t
                pe.time_absolute = t2
                pe.exposure_time_calculated_via_photon_time = None
            if pe.x is None:
                pe.number_of_photons = 0
            else:
                pe.number_of_photons = len(pe)

            #########
            # attrs #
            #########
            metadata = h5.sort_metadata(f=h, attrs_dict=rixs_attrs, verbose=verbose)
            for attr in metadata:
                setattr(pe, attr, metadata[attr])

            # scan
            pe.scan = scan

            # undulator
            if 'phase' in pe.get_attrs() and 'gap' in pe.get_attrs():
                if pe.phase is not None and pe.gap is not None:
                    pe.pol = _pol(phase=pe.phase, gap=pe.gap)
                else:
                    pe.pol = None
            else:
                pe.pol = None

            # round 2
            for _attr in ('th_veritas', 'tth', 'E', 'T', 'gap', 'arm_energy', 'mono_energy', 'phase'):
                for suffix in ('', '_max', '_min', '_sigma'):
                    attr = _attr + suffix
                    if hasattr(pe, attr):
                        if pe.__getattribute__(attr) is not None:
                            pe.__setattr__(attr, round(pe.__getattribute__(attr), 2))

            # round 4
            for _attr in ('sample_x', 'sample_y', 'sample_z'):
                for suffix in ('', '_max', '_min', '_sigma'):
                    attr = _attr + suffix
                    if hasattr(pe, attr):
                        if pe.__getattribute__(attr) is not None:
                            pe.__setattr__(attr, round(pe.__getattribute__(attr), 4))

            # datetime
            for attr in ('start_time', 'end_time'):
                if hasattr(pe, attr):
                    if pe.__getattribute__(attr) is not None:
                        pe.__setattr__(attr, _str2datetime(pe.__getattribute__(attr)))

            # calculated time
            if hasattr(pe, 'start_time') and hasattr(pe, 'end_time'):
                if isinstance(pe.start_time, datetime.datetime) and isinstance(pe.end_time, datetime.datetime):
                    try:
                        pe.exposure_time_calculated = (pe.end_time - pe.start_time).total_seconds()
                    except:
                        pe.exposure_time_calculated = None
                else:
                    pe.exposure_time_calculated = None
            else:
                pe.exposure_time_calculated = None

            # miliseconds to seconds
            for attr in ('exposure_time', ):
                if hasattr(pe, attr):
                    if pe.__getattribute__(attr) is not None:
                        pe.__setattr__(attr, pe.__getattribute__(attr)/1000)

            # real th
            # if hasattr(pe, 'th_veritas'):
            #     if pe.__getattribute__('th_veritas') is not None:
            #         pe.__setattr__('th', pe.__getattribute__('th_veritas')-285+90)
            # if hasattr(pe, 'th_veritas_max'):
            #     if pe.__getattribute__('th_veritas_max') is not None:
            #         pe.__setattr__('th_max', pe.__getattribute__('th_veritas_max')-285+90)
            # if hasattr(pe, 'th_veritas_min'):
            #     if pe.__getattribute__('th_veritas_min') is not None:
            #         pe.__setattr__('th_min', pe.__getattribute__('th_veritas_min')-285+90)
            # if hasattr(pe, 'th_veritas_sigma'):
            #     if pe.__getattribute__('th_veritas_sigma') is not None:
            #         pe.__setattr__('th_sigma', pe.__getattribute__('th_veritas_sigma'))
                           

            # # post processing
            # pe.cutoff                     = None
            # pe.time_bunch_histogram       = None
            # pe.time_bunch_histogram_nbins = None
            # pe.folded_time_bunch          = None
            # pe.folded_tmask               = None
            # pe.tmask                      = None
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
            h  = f[prefix + str(scan)]
            
            ######################################################
            # get command, scanned motors, x axis, and scan type #
            ######################################################
            assert 'title' in h, 'cannot find command that launched this scan'
            assert 'measurement' in h, 'cannot find data inside file'
            command = h['title'][()].decode("utf-8")
            h2      = h['measurement']

            #####################
            # find kind of scan #
            #####################
            if command.startswith('ascan '):
                motors = [command.split(' ')[1], ]
                x      = h2[motors[0]][()]
                if motors[0] == 'beamline_energy':
                    scan_type = 'XAS'
                elif motors[0] == 'a_mp1_x':
                    scan_type = 'ascan x'
                elif motors[0] == 'a_mp1_y':
                    scan_type = 'ascan y'
                elif motors[0] == 'a_mp1_z':
                    scan_type = 'ascan z'
                else:
                    scan_type = 'ascan ' + motors[0]
            elif command.startswith('mesh '):
                motors    = [command.split(' ')[1], command.split(' ')[5],]
                x         = h2['Pt_No'][()]
                scan_type = f'mesh {motors}'
            elif command.startswith('a2scan '):
                motors    = [command.split(' ')[1], command.split(' ')[4],]
                x         = h2['Pt_No'][()]
                scan_type = f'a2scan {motors}'
            else: 
                raise ValueError(f'command not recognized: `{command}`')

            #########################
            # create Spectra object #
            #########################
            try: RMU = br.Spectrum(x=x, y=h2['aemexp2_ch1'][:] / h2['aemexp2_timer'][:])
            except: RMU = br.Spectrum()
            try: MCP = br.Spectrum(x=x, y=h2['aemexp2_ch2'][:] / h2['aemexp2_timer'][:])
            except: MCP = br.Spectrum()
            try: TFY = br.Spectrum(x=x, y=h2['aemexp2_ch3'][:] / h2['aemexp2_timer'][:])
            except: TFY = br.Spectrum()
            try: TEY = br.Spectrum(x=x, y=h2['aemexp2_ch4'][:] / h2['aemexp2_timer'][:])   
            except: TEY = br.Spectrum()
            ss  = br.Spectra([TEY, MCP, TFY, RMU])

            ###############
            # check empty #
            ###############
            if TEY.y is None:
                for s in [_ for _ in ss] + [ss]:
                    s.error = 'empty scan'
                if verbose: print(f'Warning: empty scan ({scan})')

            #####################
            # reshape mesh scan #
            #####################
            if scan_type.startswith('mesh') and TEY.y is not None:
                # get which motor was scanned first
                # motor a is frozen while b is scanned
                # command = h['title'][()].decode("utf-8")
                # if command.find(motors[1]) > command.find(motors[0]):
                #     motors = motors[0], motors[1]
                # else:
                #     motors = motors[1], motors[0]
                # a, b = h2[motors[0]], h2[motors[1]]

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
                ss = br.Dummy([TEY, MCP, TFY, RMU])
                for j, s in enumerate(ss):
                    y = copy.deepcopy(s.y)
                    if len(y) < (len(afinal) * len(bfinal)):
                        y = np.array(list(y) + [y[-1]]*((len(afinal) * len(bfinal))-len(y)))
                        ss.error = 'mesh interrupted early'
                    if command.split(' ')[-1] == 'True':
                        ss[j] = br.Image(data=[row if i%2 == 0 else row[::-1] for i, row in enumerate(y.reshape(int(command.split(' ')[8])+1, int(command.split(' ')[4])+1))])
                    else:
                        ss[j] = br.Image(data=[row if i%2 == 0 else row[::] for i, row in enumerate(y.reshape(int(command.split(' ')[8])+1, int(command.split(' ')[4])+1))])
                    ss[j].x_centers = afinal
                    ss[j].y_centers = bfinal

            #######################
            # reshape a2scan scan #
            #######################
            elif scan_type.startswith('a2scan') and TEY.y is not None:
                for s in ss:
                    s.a2scan_x      = [list(s.x), list(h2[motors[0]][:]), list(h2[motors[1]][:])]
                    s.a2scan_labels = ['Pt_No', ] + motors

            #############
            # raw attrs #
            #############   
            metadata = h5.sort_metadata(f=h, attrs_dict=xas_attrs, verbose=verbose)
            for s in [_ for _ in ss] + [ss]:
                for attr in metadata:
                    setattr(s, attr, metadata[attr])

            ################
            # pretty attrs #
            ################
            for s in [_ for _ in ss] + [ss]:
                # basic
                s.scan           = scan
                s.scanned_motors = motors
                s.scan_type      = scan_type

                # undulator
                if 'phase' in s.get_attrs() and 'gap' in s.get_attrs():
                    if s.phase is not None and s.gap is not None:
                        s.pol = _pol(phase=s.phase, gap=s.gap)
                    else:
                        s.pol = None
                else:
                    s.pol = None

                # datetime
                for attr in ('start_time', 'end_time'):
                    if hasattr(s, attr):
                        if s.__getattribute__(attr) is not None:
                            s.__setattr__(attr, _str2datetime(s.__getattribute__(attr)))

                # elapsed time
                if hasattr(s, 'start_time') and hasattr(s, 'end_time'):
                    if isinstance(s.start_time, datetime.datetime) and isinstance(s.end_time, datetime.datetime):
                        try:
                            s.elapsed_time = (s.end_time - s.start_time).total_seconds()
                        except:
                            s.elapsed_time = None
                    else:
                        s.elapsed_time = None
                else:
                    s.elapsed_time = None

                # # real th
                # if hasattr(s, 'th_veritas'):
                #     if s.__getattribute__('th_veritas') is not None:
                #         s.__setattr__('th', s.__getattribute__('th_veritas')-285+90)
     
            ss.header = ('TEY', 'MCP', 'TFY', 'RMU')
            TEY.label = 'TEY'
            MCP.label = 'MCP'
            TFY.label = 'TFY'
            RMU.label = 'RMU'

            ###########
            # verbose #
            ###########
            if verbose:
                print(f'scan type: {scan} - {scan_type}')
            
            return ss
    return
# %%
