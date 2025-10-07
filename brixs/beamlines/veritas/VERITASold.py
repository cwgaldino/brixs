#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from VERITAS beamline of MAX-IV.

>>>>> br.finder.folderpath = <temp-file>

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
>>> TEY, MCP, TFY, RMU = VERITAS.read(XAS, 164)
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
from collections.abc import Iterable
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import datetime
import warnings
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
def read(filepath, scan, verbose=True):
    """Return data from filepath
    
    RIXS        --> returns PhotonEvent 
    XAS         --> returns Spectra (TEY, MCP, TFY, RMU) 
    sample scan --> returns Spectra (TEY, MCP, TFY, RMU) 
    mesh scan   --> returns Image (TEY, MCP, TFY, RMU) 

    No extra processing is done on RIXS spectra
    For XAS, sample scan, and mesh scan, datapoints are divided by the exposure

    # in seconds
    exposure_time
    exposure_time_calculated
    exposure_time_calculated_via_photon_time


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
            assert 'measurement' in h, 'cannot find scan data inside file'
            command = h['title'][()].decode("utf-8")
            h2      = h['measurement']
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
            RMU = br.Spectrum(x=x, y=h2['aemexp2_ch1'][:] / h2['aemexp2_timer'][:])
            MCP = br.Spectrum(x=x, y=h2['aemexp2_ch2'][:] / h2['aemexp2_timer'][:])
            TFY = br.Spectrum(x=x, y=h2['aemexp2_ch3'][:] / h2['aemexp2_timer'][:])
            TEY = br.Spectrum(x=x, y=h2['aemexp2_ch4'][:] / h2['aemexp2_timer'][:])   
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

# %% ============================= RIXS =================================== %% #
def _process(filepath, scan, mask=None, tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max', curv='self', curv_nbins=(20, 1000), sbins=1200, calib=None):
    """
        s = s.set_factor(1/sum([m[3]-m[2] for m in mask]))
        s = s.set_factor(1/s.exposure_time)
        s = s.set_factor(sbins)
        s = s.set_factor(1000)

        future: put an argument like, time_start, time_finish to control if we
            will read all photons or just photos at the beginning or end of a scan.
             Should be easy if we convert from ps to s. We can do just like 
             "time cutoff".

    """
    #############
    # read file #
    #############
    pe = read(filepath, scan)

    ###############
    # time cutoff #
    ###############
    if tcutoff is not None:
        pe, cutoff = pe.apply_time_cutoff(value=tcutoff)
    else:
        cutoff = None

    ########
    # mask #
    ########
    if mask is None:
        mask = [[pe.xlim[0], pe.xlim[1], pe.ylim[0], pe.ylim[1]], ]
    else:
        mask = pe._check_mask(mask)
    pe2      = pe.clip2(mask)
    pe2.mask = copy.deepcopy(mask)

    # get separate PhotonEvents'
    pes = br.Dummy()
    for m in mask:
        pes.append(pe.clip2(m))

    ###############
    # folded time #
    ###############
    time_bunch_histogram = pe2.calculate_time_bunch_histogram(nbins=tnbins)  # returned spectrum has attr: nbins
    folded_time = time_bunch_histogram.calculate_folded_time(period=period, offset=offset)  # returned spectrum has attrs: nbunchs and period

    ##########
    # t mask #
    ##########
    try:
        folded_tmask, tmask = calculate_folded_tmask_and_tmask(folded_time, time_bunch_histogram, w=twidth, c=tcenter)
    except AssertionError:
        print('error finding suitable c for time window')
        folded_tmask, tmask = calculate_folded_tmask_and_tmask(folded_time, time_bunch_histogram, w=twidth, c='max')

    #########################
    # time rejected photons #
    #########################
    pe3, bad = pe2.apply_tmask(tmask=tmask)

    # curvature
    if curv is not None:
        if type(curv) == str:
            if curv == 'self':
                # curvature correction with broken intervals
                fit, curv = pes.curvature_correction_with_broken_intervals(curv_nbins=curv_nbins)
                
                # curvature correction for solid intervals (works fine)
                # im = pe3.binning(ncols=curv_nbins[1], nrows=curv_nbins[0])
                # limits = [[m[2], m[3]] for m in mask]
                # _s, fit, popt, R2, model = im.calculate_horizontal_shift_curvature(deg=2, mode='cc', limits=limits)
                # curv = popt

        pe4 = pe3.set_horizontal_shift_via_polyval(p=curv)
    else:
        pe4 = None
    
    # spectrum
    if curv is not None and isinstance(curv, br.Iterable):
        s = pe4.integrated_columns_vs_x_centers(ncols=sbins)
        s.scan = scan

        # normalization
        s = s.set_factor(1/sum([m[3]-m[2] for m in mask]))
        s = s.set_factor(1/s.exposure_time)
        s = s.set_factor(sbins)
        s = s.set_factor(1000)
    else:
        s = None
        
    # calib ==============================================
    if calib is not None and s is not None:
        s.calib = calib
        s.shift = -s.E

    return {'pe':pe, 'cutoff':cutoff, 'pe2':pe2, 
            'time_bunch_histogram':time_bunch_histogram, 
            'folded_time':folded_time, 
            'tmask':tmask, 
            'folded_tmask':folded_tmask, 
            'bad':bad, 
            'pe3':pe3, 
            'pe4':pe4, 
            'curv':curv, 
            'curv_fit':fit,
            's':s}

def verify(filepath, scan, mask, tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max', curv='self', curv_nbins=(20, 1000), sbins=1200, calib=None, maximize=True, left=0.04, right=0.99, top=0.97, bottom=0.05):
    """open a figure with step-by-step rixs data reduction

    Args:

    
    Returns:
        pe, pe2, pe3, pe4, curv, s
    """
    ##### process ######
    d = _process(filepath=filepath, scan=scan, mask=mask, tcutoff=tcutoff, tnbins=tnbins, period=period, offset=offset, twidth=twidth, tcenter=tcenter, curv=curv, curv_nbins=curv_nbins, sbins=sbins, calib=calib)
    pe                   = d['pe']
    cutoff               = d['cutoff']
    pe2                  = d['pe2']
    time_bunch_histogram = d['time_bunch_histogram']
    folded_time          = d['folded_time']
    tmask                = d['tmask']
    folded_tmask         = d['folded_tmask']
    bad                  = d['bad']
    pe3                  = d['pe3']
    pe4                  = d['pe4']
    curv_fit             = d['curv_fit']
    s                    = d['s']


    ###### figure ######
    fig, axes = br.subplots(2, 5, figsize=(12, 4))
    plt.subplots_adjust(left=left, right=right, top=top, bottom=bottom)
    if maximize:
        br.maximize()

    ###### axes 0 ######
    ax = axes[0]
    ax.set_title('1: Raw (pe1)')
    # plot pe
    pe.plot(ax, color='black')  
    if tcutoff is not None:
        cutoff.plot(ax, color='magenta', s=.6)
    # plot mask
    for m in pe2.mask:
        br.rectangle(m[:2], m[2:], ax=ax, lw=2, edgecolor='red')
    
    ###### axes 1 ######
    ax = axes[1]
    ax.set_title('2: After mask (pe2)')
    # plot pe2
    pe2.plot(ax, color='dodgerblue')   
    
    ###### axes 2 ######
    ax = axes[2]
    ax.set_title('4a: Time-rejected photons (bad)')
    bad.plot(ax, color='lightcoral')

    ###### axes 3 ######
    ax = axes[3]
    ax.set_title('4b: Time-accepted photons (pe3)')
    pe3.plot(ax, color='green')
    # plot fit
    if pe4 is not None:
        temp  = pe3.binning(nrows=20, ncols=1000).rows[0]
        shift = temp.x[np.argmax(temp.y)]
        curv_fit.switch_xy().flip_x().plot(ax, offset=0, shift=shift, color='red')
 

    ###### axes 4 ######
    ax = axes[4]
    if pe4 is not None:
        ax.set_title('5: curvature corrected (pe4)')
        pe4.plot(ax, color='red')

    ###### axes 5 ######
    ax = axes[5]
    ax.set_title('3a: Folded Time (folded_time)')
    folded_time.plot(ax)
    for m in folded_tmask:
        br.rectangle(m, (0, max(folded_time.y)), ax=ax, lw=2, edgecolor='red')
    
    ###### axes 6 ######
    ax = axes[6]
    ax.set_title('3b: Begining of unfolded Time')
    time_bunch_histogram.plot(ax, color='black', lw=.2, marker='o', ms=1)
    for m in tmask:
        i = time_bunch_histogram.index(m[0])
        f = time_bunch_histogram.index(m[1])
        br.rectangle(m, (0, max(time_bunch_histogram.y[i:f])), ax=ax, lw=2, edgecolor='red')
    if offset is None: offset = 0
    br.zoom(min(time_bunch_histogram.x)-min(time_bunch_histogram.x)*-0.1, period*5+offset, ax)        

    # ###### axes 7 ######
    # ax = axes[7]
    # ax.set_title('3b: End of unfolded Time (time_bunch_histogram)')
    # time_bunch_histogram.plot(ax, color='black', lw=.2, marker='o', ms=1)
    # for m in tmask:
    #     i = time_bunch_histogram.index(m[0])
    #     f = time_bunch_histogram.index(m[1])
    #     br.rectangle(m, (0, max(time_bunch_histogram.y[i:f])), ax=ax, lw=2, edgecolor='red')
    # if offset is None: offset = 0
    # br.zoom(max(time_bunch_histogram.x) - period*5+offset, max(time_bunch_histogram.x)+max(time_bunch_histogram.x)*0.1, ax)  

    ###### axes 7 ######
    ax = axes[7]
    ax.set_title('3b: Full unfolded Time (time_bunch_histogram)')
    time_bunch_histogram.plot(ax, color='black', lw=.2, marker='o', ms=1)
    for m in tmask:
        i = time_bunch_histogram.index(m[0])
        f = time_bunch_histogram.index(m[1])
        br.rectangle(m, (0, max(time_bunch_histogram.y[i:f])), ax=ax, lw=2, edgecolor='red')

    ###### axes 8 ######
    ax = axes[8]
    ax.set_title('3b: Photons per time window')
    get_count_per_time_window(pe2.time_bunch, tmask).plot(ax, color='black', marker='o', ms=2)

    ###### axes 9 ######
    ax = axes[9]
    ax.set_title('6: Final spectrum')
    if s is not None:
        s.plot(ax, color='black', marker='o', ms=2)
    
    ###### axis labels ######
    for i in (0, 1, 2, 3, 4):
        axes[i].set_xlabel('x (mm)')
        axes[i].set_ylabel('y (mm)')
    axes[5].set_xlabel('Folded time (ps)')
    axes[5].set_ylabel('Photon count')
    for i in (6, 7):
        axes[i].set_xlabel('Time (ps)')
        axes[i].set_ylabel('Photon count')
    axes[8].set_xlabel('Time (ps)')
    axes[8].set_ylabel('Photon count per window')
    if calib is None:
        axes[9].set_xlabel('x (mm)')
    else:
        axes[9].set_xlabel('Energy loss (eV)')
    axes[9].set_ylabel('Intensity (arb. units)')

    return {'pe':pe, 'cutoff':cutoff, 'pe2':pe2, 'time_bunch_histogram':time_bunch_histogram, 'folded_time':folded_time, 'tmask':tmask, 'folded_tmask':folded_tmask, 'bad':bad, 'pe3':pe3, 'pe4':pe4, 'curv':curv, 's':s}

@br.finder.track
def process(filepath, scan, mask, tcutoff=3e7, tnbins=10000, period=1458, offset=None, twidth=320, tcenter='max', curv='self', curv_nbins=(20, 1000), sbins=1200, calib=None):
    """
    """
    d = _process(filepath=filepath, scan=scan, mask=mask, tcutoff=tcutoff, tnbins=tnbins, period=period, offset=offset, twidth=twidth, tcenter=tcenter, curv=curv, curv_nbins=curv_nbins, sbins=sbins, calib=calib)
    return d['s']
# %%


# %%

# %% ========================== Time support ============================== %% #
def _apply_time_cutoff(self, value=3e7):
    """return good and cutoff photon events

    Args:
        Value (number, optional): cutoff time. Photon events with time higher then
            this value will be removed
    
    Returns:
        good, cutoff
    """
    # cutoff time
    indexes = np.array([i for i in range(len(self.time_bunch)) if self.time_bunch[i] > value]) # remove this one last point
    final  = [True]*len(self.time_bunch)
    for i in indexes:
        final[i] = False
    _cutoff = [not elem for elem in final]
    
    # saving
    cutoff = br.PhotonEvents(x=self.x[_cutoff], y=self.y[_cutoff])
    cutoff.copy_attrs_from(self)
    cutoff.xlim = self.xlim
    cutoff.ylim = self.ylim
    cutoff.time_absolute = self.time_absolute[_cutoff]
    cutoff.time_bunch    = self.time_bunch[_cutoff]
    
    good = br.PhotonEvents(x=self.x[final], y=self.y[final])
    good.copy_attrs_from(self)
    good.xlim = self.xlim
    good.ylim = self.ylim
    good.time_absolute = self.time_absolute[final]
    good.time_bunch    = self.time_bunch[final]

    return good, cutoff
br.PhotonEvents.apply_time_cutoff = _apply_time_cutoff    

def _calculate_time_bunch_histogram(self, nbins=10000):
    """returns the histogram of time_bunch

    Args:
        nbins (int, optional): number of bins 

    Return:
        Spectrum
        returned spectrum has attr: nbins
    """
    # calculate time histogram
    hist, bin_edges = np.histogram(self.time_bunch, nbins)
    time            = br.Spectrum(x=br.moving_average(bin_edges, 2), y=hist)
    
    # self.time.shift = -min(self.time.x)
    time.nbins = nbins
    # self.time_bunch_histogram       = time
    # self.time_bunch_histogram_nbins = nbins
    return time
br.PhotonEvents.calculate_time_bunch_histogram = _calculate_time_bunch_histogram        

def _calculate_folded_time(self, period=1459, offset=None):
    """return folded time

    Args:
        period (int, optional): time period to fold. After this value, time is supposed to repeat itself
        offset (number, optional): use `offset` to shift time before applying period, 
            in case time does not start from 0
    
    Return:
        Spectrum
        returned spectrum has attrs: nbunchs and period
    """    
    if offset is None:
        offset = min(self.x)
    x = (self.x - offset)%period
    
    hist, bin_edges = np.histogram(x, bins=int(period/10), weights=self.y)
    folded_time     = br.Spectrum(x=br.moving_average(bin_edges, 2), y=hist)
    # folded_t        = (self.t-min(self.t))%period  # min(self.t) should not be necessary because this is already set to zero during read()
    folded_time.nbunchs = int(round((max(self.x) - offset)/period))
    folded_time.period  = period
    return folded_time
br.Spectrum.calculate_folded_time = _calculate_folded_time

def calculate_folded_tmask_and_tmask(folded_time, time_bunch_histogram, w=300, c='max'):
    """return a suitable time mask from folded time (a peak -bunch- is expected)

    Args:
        folded_time (br.Spectrum): folded time calculated via br.Spectrum.calculate_folded_time() function. 
            Spectrum must have attr `period` used to calculate the folding
        time_bunch_histogram (br.Spectrum): time bunch histogram calculated via
            br.PhotonEvents.calculate_time_bunch_histogram() function
        w (number, optional): size of the mask
        c (number or str, optional): time center of the bunch or method for guessing the center. 
            options for guessing are 'gauss' and 'max'.
        
    Return:
        folded_tmask, tmask
    """
    # w must be less or equal the full time range
    assert w <= max(folded_time.x) - min(folded_time.x), f'w ({w}) must be equal or less than full time range ({max(self.x) - min(self.x)})'

    # get folded_time peak center
    if c == 'gauss':
        # check if fitting was imported
        if hasattr(folded_time, 'fit_peak') == False and callable(folded_time.fit_peak) == False:
            raise ValueError('cannot calculate shifts via `peaks` because fitting functions are not imported\nPlease import fitting function via `import brixs.addons.fitting`')
        fit, popt, err, f = folded_time.fit_peak()
        c = popt[1]
    elif c == 'max':
        c = folded_time.x[np.argmax(folded_time.y)]
    
    # tmask
    right = max(folded_time.x)
    left  = min(folded_time.x)
    assert c >= left and c <= right, f'c={c} outside range'
    
    # find if mask is outside
    if c+w/2 > right: 
        rest = w - (right - (c-w/2))
        folded_tmask = [(left, left+rest), (c-w/2, right)]
    elif c-w/2 < left: 
        rest = w - ((c+w/2) - left)
        folded_tmask = [(left, c+w/2), (right-rest, right)]
    else:
        folded_tmask = [(c-w/2, c+w/2), ]

    # calculate tmask
    tmask = [(c1-w/2+c, c1+w/2+c) for c1 in np.arange(min(time_bunch_histogram.x), max(time_bunch_histogram.x), folded_time.period)][:-1]

    return folded_tmask, tmask

def _apply_tmask(self, tmask):
    """Return photon events with point only within tmask

    Args:
        tmask (list): time mask with the following format: ((t1_start, t1_stop), (t2_start, t2_stop), ...)

    Note:
        tmask is applied on the regular time (not folded time)

    Return:
        good, bad
        PhotonEvent's have attr: percentage
    """
    indexes = []
    for m in tmask:
        indexes += [i for i, _t in enumerate(self.time_bunch) if _t>=m[0] and _t<=m[1]]
    final = [False]*len(self.time_bunch)

    for i in indexes:
        final[i] = True
    opposite = [not elem for elem in final]

    # final
    good = br.PhotonEvents(x=self.x[final], y=self.y[final])
    good.time_bunch    = np.array(self.time_bunch)[final]    #- min(self.t[final])    # set min value to zero
    good.time_absolute = np.array(self.time_absolute)[final]    #- min(self.t[final])    # set min value to zero
    good.copy_attrs_from(self)
    good.xlim = self.xlim
    good.ylim = self.ylim
    good.percentage = round(len(good.x)/len(self.x)*100, 2)

    bad  = br.PhotonEvents(x=self.x[opposite], y=self.y[opposite])
    bad.time_bunch  = np.array(self.time_bunch)[opposite] #- min(self.t[opposite]) # set min value to zero
    good.time_absolute   = np.array(self.time_absolute)[final]    #- min(self.t[final])    # set min value to zero
    good.copy_attrs_from(self)
    good.xlim = self.xlim
    good.ylim = self.ylim
    bad.percentage  = round(len(bad.x)/len(self.x)*100, 2)
    
    return good, bad
br.PhotonEvents.apply_tmask = _apply_tmask

def get_count_per_time_window(time_bunch, tmask):
    centers = []
    counts  = []
    for m in tmask:
        centers += [m[0] + (m[1] - m[0])/2]
        counts += [sum([1 for t in time_bunch if t > m[0] and t < m[1]])]            
    return br.Spectrum(centers, counts)
# %%

# %% ======================== spacial support ============================= %% #
def _clip(self, mask):
    """Return a masked copy of the object (VERITAS).

    Note:
        This is different from pe.clip() because this also clips the time.

    Args:
        mask (list): list with rectangular coordinates `(x_start, x_stop, y_start, y_stop)`
            or a list with multiple rectangular coordinates, i.e., `[(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])`

    Returns:
        :py:attr:`PhotonEvents`
    """
    return self.clip(mask=mask, attrs2clip=['time_bunch', 'time_absolute'])
br.PhotonEvents.clip2 = _clip

def _plot_a2scan(self):
    assert self.command.startswith('a2scan '), 'this is not a a2scan'

    # save initial x
    x0 = self.x

    # figure
    br.figure()

    # plot with x = 1 (ax1)
    self.select_x_for_a2scan(1, verbose=False)
    ax1 = self.plot().axes

    # plot with x = 2 (ax2)
    ax2 = ax1.twiny()
    self.select_x_for_a2scan(x=2, verbose=False)
    self.plot(ax=ax2, color='red')

    # linear functions for converting to motor 1 to motor 2
    # x1 = (TEY.a2scan_x[1][-1] - TEY.a2scan_x[1][0])/(TEY.a2scan_x[0][-1] - TEY.a2scan_x[0][0])*pt + TEY.a2scan_x[1][0]
    # x2 = (TEY.a2scan_x[2][-1] - TEY.a2scan_x[2][0])/(TEY.a2scan_x[0][-1] - TEY.a2scan_x[0][0])*pt + TEY.a2scan_x[2][0]
    def x1_to_pt(x):
        return (x - self.a2scan_x[1][0])*(self.a2scan_x[0][-1] - self.a2scan_x[0][0])/(self.a2scan_x[1][-1] - self.a2scan_x[1][0])
    def pt_to_x2(pt):
        return (self.a2scan_x[2][-1] - self.a2scan_x[2][0])/(self.a2scan_x[0][-1] - self.a2scan_x[0][0])*pt + self.a2scan_x[2][0]
    def x1_to_x2(x):
        return pt_to_x2(x1_to_pt(x))

    # set ax2 ticks
    ax2.set_xlim([x1_to_x2(x) for x in ax1.get_xlim()])
    ax2.spines['top'].set_color('red') 

    # labels
    ax1.set_xlabel(self.a2scan_labels[1])
    ax2.set_xlabel(self.a2scan_labels[2])

    return ax1, ax2
br.Spectrum.plot_a2scan = _plot_a2scan

def _select_x_for_a2scan(self, x, verbose=True):
    """
    x = int or motor name
    """
    if isinstance(x, str):
        assert x in self.a2scan_labels, f'x=`{x}` is not a valid motor name. Available options: {self._a2scan_x_labels}'
        x = self.a2scan_labels.index(x)
    elif isinstance(x, int):
        assert x >= 0 and x <= 2, f'x must be a int between 0 and 2 or a motor name, not `{x}`'
    else:
        raise ValueError(f'x must be a int or a motor name, not type `{type(x)}`')
    self.x = self.a2scan_x[x]
    if verbose: print(f'x selected: {self.a2scan_labels[x]}')
br.Spectrum.select_x_for_a2scan = _select_x_for_a2scan
# %%

# %% ====================== curvature correction ========================== %% #
def _curvature_correction_with_broken_intervals(self, curv_nbins):
    # binning PhotonEvents'
    ims = br.Dummy()
    for i in range(len(self)):
        ims.append(self[i].binning(ncols=curv_nbins[1], nrows=int(curv_nbins[0]/len(self))))

    # collect rows
    ss = br.Spectra()
    y_centers = []
    for i in range(len(self)):
        for y in ims[i].y_centers:
            y_centers.append(y)
        for s in ims[i].rows:
            ss.append(s)

    # calculate shift
    try:
        ss.check_same_x()
    except ValueError:
        ss = ss.interp()
    values = ss.calculate_shift(mode='cc')
    _s = br.Spectrum(y_centers, values)
    fit, popt, R2, model = _s.polyfit(deg=2)
    curv = popt

    return fit, curv
br.Dummy.curvature_correction_with_broken_intervals = _curvature_correction_with_broken_intervals
# %%

# %% ============================ obsolete ================================ %% #
def _clip(self, mask):
    """Return a masked copy of the object (VERITAS).

    Note:
        This is different from pe.clip() because this also clips the time.

    Args:
        mask (list): list with rectangular coordinates `(x_start, x_stop, y_start, y_stop)`
            or a list with multiple rectangular coordinates, i.e., `[(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])`

    Returns:
        :py:attr:`PhotonEvents`
    """
    ######################
    # assert mask format #
    ######################
    assert isinstance(mask, Iterable), 'mask must be iterable'
    if len(mask) == 4:
        if isinstance(mask[0], Iterable) == False:
            mask = [mask, ] 
    for m in mask:
        assert len(m) == 4, 'mask must have the format: [(x1_start, x1_stop, y1_start, y1_stop), (x2_start, x2_stop, y2_start, y2_stop), ...])'
        
    ########
    # clip #
    ########
    x = []
    y = []
    time_bunch = []
    time_absolute = []
    xstart = None
    xstop  = None
    ystart = None
    ystop  = None
    for r in mask:
        temp = np.array([(x, y, t, t2) for x, y, t, t2 in zip(self.x, self.y, self.time_bunch, self.time_absolute) if ((x > r[0] and x < r[1]) and (y > r[2] and y < r[3]))])
        x += list(temp[:, 0])
        y += list(temp[:, 1])
        time_bunch    += list(temp[:, 2])
        time_absolute += list(temp[:, 3])
        if xstart is None: xstart = r[0]
        if xstop is None:  xstop  = r[1]
        if ystart is None: ystart = r[2]
        if ystop is None:  ystop  = r[3]
        if r[0] < xstart: xstart = r[0]
        if r[1] > xstop:  xstop  = r[1]
        if r[2] < ystart: ystart = r[2]
        if r[3] > ystop:  ystop  = r[3]


    #########
    # final #
    #########
    pe = br.PhotonEvents(x=x, y=y)
    pe.copy_attrs_from(self)
    pe.xlim = [xstart, xstop]
    pe.ylim = [ystart, ystop]
    pe.time_absolute = time_absolute
    pe.time_bunch    = time_bunch

    return pe
# br.PhotonEvents.clip2 = _clip

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
# br.Image.find_empty = _find_empty

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
# br.Image.get_filtered_spectra = _get_filtered_spectra
# %%


# %% ============================== other ================================= %% #

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
# %%

