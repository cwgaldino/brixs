#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from IPE beamline - Sirius.

Last edited: Felipe Custódio and Carlos Galdino 2024-08-15
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from collections.abc import Iterable
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import matplotlib
import datetime
import warnings
import copy

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br
import brixs.addons.h5 as h5
br.Dummy.create_attr_from_spectra = br.Spectra.create_attr_from_spectra
br.Dummy.get_by_attr = br.Spectra.get_by_attr

# %% ------------------------- Special Imports ---------------------------- %% #
import h5py
# %%

# %% ========================= useful functions =========================== %% #
def scanlist(folderpath):
    """Return list of scans available in folderpath"""
    folderpath = Path(folderpath)
    assert folderpath.exists(), f'fpath does not exist ({folderpath})'
    return br.parsed_filelist(folderpath, string='*', ref=0, return_type='dict')
# %%

# %% =================== metadata support functions ======================= %% #
def _str2datetime(string):
    """convert IPE date/time string pattern to date
    
    Example:
        '2022/07/20 21:08:36.0'
        '2022/07/20T21:08:36.0'
        '2022-07-20 21:08:36.0'
        '2022-07-20T21:08:36.0'

    Args:
        string (str): string with IPE date string

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
    if '/' in date:
        year, month, day = (int(_) for _ in date.split('/'))
    elif '-' in date:
        year, month, day = (int(_) for _ in date.split('-'))
    else:
        raise ValueError(f'cannot split date: {date}')

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
# %%

# %% =========================== metadata ============================= %% #
_attrs = {}

# %% ========================== rixs metadata ============================= %% #
_attrs['rixs'] = {'ignore': {}, 'raw':{}}

h = _attrs['rixs']['ignore']
h['modified_date'] = ''
h['scan']          = ''
h['error']         = ''

h = _attrs['rixs']['raw']
# h['Energy']              = 'entry/instrument/NDAttributes/Energy'
# h['Energy_SP']           = 'entry/instrument/NDAttributes/Energy_SP'
# h['NDArrayEpicsTSSec']   = 'entry/instrument/NDAttributes/NDArrayEpicsTSSec'
# h['NDArrayEpicsTSnSec']  = 'entry/instrument/NDAttributes/NDArrayEpicsTSnSec'
# h['NDArrayTimeStamp']    = 'entry/instrument/NDAttributes/NDArrayTimeStamp'
# h['NDArrayUniqueId']     = 'entry/instrument/NDAttributes/NDArrayUniqueId'
# h['PGM_Cff']             = 'entry/instrument/NDAttributes/PGM_Cff'
# h['PGM_GR']              = 'entry/instrument/NDAttributes/PGM_GR'
# h['PGM_GT']              = 'entry/instrument/NDAttributes/PGM_GT'
# h['PGM_MR']              = 'entry/instrument/NDAttributes/PGM_MR'
# h['PGM_MT']              = 'entry/instrument/NDAttributes/PGM_MT'
# h['RIXSCam_ActualImage'] = 'entry/instrument/NDAttributes/RIXSCam_ActualImage'
# h['RIXSCam_NumImages']   = 'entry/instrument/NDAttributes/RIXSCam_NumImages'
# h['RIXSCam_exposure']    = 'entry/instrument/NDAttributes/RIXSCam_exposure'
# h['RIXS_Ry']             = 'entry/instrument/NDAttributes/RIXS_Ry'
# h['RIXS_X']              = 'entry/instrument/NDAttributes/RIXS_X'
# h['RIXS_Y']              = 'entry/instrument/NDAttributes/RIXS_Y'
# h['RIXS_Z']              = 'entry/instrument/NDAttributes/RIXS_Z'
# h['RIXS_GX']             = 'entry/instrument/NDAttributes/RIXS_GX'
# h['RIXS_GY']             = 'entry/instrument/NDAttributes/RIXS_GY'
# h['RIXS_GZ']             = 'entry/instrument/NDAttributes/RIXS_GZ'
# h['RIXS_GRx1']           = 'entry/instrument/NDAttributes/RIXS_GRx1'
# h['RIXS_GRz1']           = 'entry/instrument/NDAttributes/RIXS_GRz1'
# h['RIXS_DZ']             = 'entry/instrument/NDAttributes/RIXS_DZ'
# h['RIXS_DY']             = 'entry/instrument/NDAttributes/RIXS_DY'
# h['TempA']               = 'entry/instrument/NDAttributes/TempA'
# h['TempB']               = 'entry/instrument/NDAttributes/TempB'
# h['Undulator']           = 'entry/instrument/NDAttributes/Undulator'  


# %% ========================== xas metadata ============================ %% #
_attrs['xas'] = {'ignore': {}, 'raw':{}, 'string': {}, 'bool': {}}

h = _attrs['xas']['ignore']
h['modified_date']    = ''
h['motors']           = ''
h['error']            = ''

h = _attrs['xas']['string']
h['motors']           = 'entry/instrument/bluesky/metadata/motors'
h['detectors']        = 'entry/instrument/bluesky/metadata/detectors'
h['scan_type']        = 'entry/instrument/bluesky/metadata/scan_type'
h['title']            = 'entry/title'
h['entry_identifier'] = 'entry/entry_identifier'
h['start_time']       = 'entry/start_time'
h['end_time']         = 'entry/end_time'

h = _attrs['xas']['raw']
h['duration']         = 'entry/duration'
h['num_points']       = 'entry/instrument/bluesky/metadata/num_points'
h['proposal']         = 'entry/instrument/bluesky/metadata/proposal'
h['scan']             = 'entry/instrument/bluesky/metadata/scan'
h['exposure']         = 'entry/instrument/bluesky/metadata/exposure'

# %% ========================== ascan metadata ============================ %% #
_attrs['ascan'] = {'ignore': {}, 'raw':{}, 'string': {}, 'bool': {}}

h = _attrs['ascan']['ignore']
h['modified_date']    = ''
h['motors']           = ''
h['error']            = ''

h = _attrs['ascan']['string']
h['motors']           = 'entry/instrument/bluesky/metadata/motors'
h['detectors']        = 'entry/instrument/bluesky/metadata/detectors'
h['scan_type']        = 'entry/instrument/bluesky/metadata/scan_type'
h['title']            = 'entry/title'
h['entry_identifier'] = 'entry/entry_identifier'
h['start_time']       = 'entry/start_time'
h['end_time']         = 'entry/end_time'

h = _attrs['ascan']['raw']
h['duration']         = 'entry/duration'
h['num_points']       = 'entry/instrument/bluesky/metadata/num_points'
h['proposal']         = 'entry/instrument/bluesky/metadata/proposal'
h['scan']             = 'entry/instrument/bluesky/metadata/scan'
h['exposure']         = 'entry/instrument/bluesky/metadata/exposure'

# %% ========================== mesh metadata ============================ %% #
_attrs['mesh'] = {'ignore': {}, 'raw':{}, 'string': {}, 'bool': {}}

h = _attrs['mesh']['ignore']
h['modified_date']    = ''
h['motors']           = ''
h['error']            = ''
h['main_motor']       = 'entry/instrument/bluesky/metadata/main_motor'

h = _attrs['mesh']['string']
h['motors']           = 'entry/instrument/bluesky/metadata/motors'
h['detectors']        = 'entry/instrument/bluesky/metadata/detectors'
h['snaking']          = 'entry/instrument/bluesky/metadata/snaking'
h['motor_x']          = 'entry/instrument/bluesky/metadata/motor_x'
h['motor_y']          = 'entry/instrument/bluesky/metadata/motor_y'
h['scan_type']        = 'entry/instrument/bluesky/metadata/scan_type'
h['title']            = 'entry/title'
h['entry_identifier'] = 'entry/entry_identifier'
h['start_time']       = 'entry/start_time'
h['end_time']         = 'entry/end_time'


h = _attrs['mesh']['raw']
h['duration']         = 'entry/duration'
h['num_points']       = 'entry/instrument/bluesky/metadata/num_points'
h['proposal']         = 'entry/instrument/bluesky/metadata/proposal'
h['nstep_x']          = 'entry/instrument/bluesky/metadata/nstep_x'
h['nstep_y']          = 'entry/instrument/bluesky/metadata/nstep_y'
h['start_x']          = 'entry/instrument/bluesky/metadata/start_x'
h['start_y']          = 'entry/instrument/bluesky/metadata/start_y'
h['stop_x']           = 'entry/instrument/bluesky/metadata/stop_x'
h['stop_y']           = 'entry/instrument/bluesky/metadata/stop_y'
h['scan']             = 'entry/instrument/bluesky/metadata/scan'
h['exposure']         = 'entry/instrument/bluesky/metadata/exposure'


h = _attrs['mesh']['bool']
h['snake']            = 'entry/instrument/bluesky/metadata/snake'

# %% ========================= read (RIXS and XAS) ======================== %% #
def _read_rixs(filepath, curv=True, verbose=True):
    """return PhotonEvents list

    Args:
        filepath (str or path): filepath
        curv (bool, optional): if True, returns the curvature corrected photon events
        verbose (bool, optional): if True, warns when metadata cannot be read

    Returns:
        pe1, pe2 (photon events for each ccd)
    """
    ##################
    # check filepath #
    ##################
    filepath = Path(filepath)
    assert filepath.exists(), f'filepath does not exist ({filepath})'

    #############
    # open file #
    #############
    with h5py.File(Path(filepath), 'r') as f:
        #############
        # read data #
        #############
        data = (f['entry/data/data'][:])
        x = data[:, 2]
        if curv:
            y = data[:, 4]
        else:
            y = data[:, 3]
            
        #############################
        # Create PhotonEvent object #
        #############################
        pe1 = br.PhotonEvents(x=x, y=y, xlim=(18, 1650),   ylim=(0, 1608)).crop(18,   1650, None, None)
        pe2 = br.PhotonEvents(x=x, y=y, xlim=(1668, 3300), ylim=(0, 1608)).crop(1668, 3300, None, None)
        
        #########
        # attrs #
        #########
        metadata = h5.sort_metadata(f=f, attrs_dict=_attrs['rixs'], verbose=verbose)
        for attr in metadata:
            setattr(pe1, attr, metadata[attr][0])
            setattr(pe1, attr, metadata[attr][0])

        _group = 'entry/instrument/NDAttributes'
        for attr in list(f[_group].keys()):
            setattr(pe1, attr, f[f'{_group}/{attr}'][0])
            setattr(pe1, attr, f[f'{_group}/{attr}'][0])

        # date
        temp = br.get_modified_date(filepath)
        setattr(pe1, 'modified_date', temp)
        setattr(pe2, 'modified_date', temp)

        # filename
        setattr(pe1, 'filename', filepath.name)
        setattr(pe2, 'filename', filepath.name)
        
        setattr(pe1, 'scan_type', 'rixs')
        setattr(pe2, 'scan_type', 'rixs')

        # scan (this must be romoved when scan number and image number are included as metadata)
        name = str(filepath.name)
        if len(name.split('_')) == 2:
            scan         = name.split('_')[0]
            image_number = name.split('_')[1].split('.')[0]
            setattr(pe1, 'scan', scan)
            setattr(pe2, 'scan', scan)
            setattr(pe1, 'image_number', image_number)
            setattr(pe2, 'image_number', image_number)
        elif len(name.split('_')) == 3:
            scan         = name.split('_')[1]
            image_number = name.split('_')[2].split('.')[0]
            setattr(pe1, 'scan', scan)
            setattr(pe2, 'scan', scan)
            setattr(pe1, 'image_number', image_number)
            setattr(pe2, 'image_number', image_number)

        # ccd
        setattr(pe1, 'ccd', 1)
        setattr(pe2, 'ccd', 2)

    return pe1, pe2

def read(fpath, verbose=True, start=0, stop=None, skip=[], curv=True):
    """Return data from folderpath (RIXS) or filepath (XAS, ascan, mesh) for IPE beamline

    Args:
        fpath (filepath or folderpath): filepath for xas and folderpath for RIXS
        verbose (bool, optional): Verbose, default is True.
        start, stop, skip (list, optional): For RIXS only. Start and stop are 
            the indexes for the first and last image to sum (inclusive). 
            Default start is 0 and the default for stop is the None, which
            will get up to the last image available. skip should be a list with
            image number indexes to not read (skip). Default is an empty list [].
        
    For RIXS, this function returns the PhotonEvents of all images summed up (pe1
    and pe2), as well as a list with the individual PhotonEvents for each image (
    pe1's and pe2's).

    Returns:
        pe1, pe2, pe1's, pe2's for RIXS
        TEY, TFY, I0, PD for XAS
        TEY, TFY, I0, PD for ascan
        (not implemented)for as2can
        (not implemented)for as3can
        TEY, TFY, I0, PD for mesh
    """
    ####################
    # check folderpath #
    ####################
    fpath = Path(fpath)
    assert fpath.exists(), f'fpath does not exist ({fpath})'

    ########
    # RIXS #
    ########
    if fpath.is_dir():
        ##################
        # get each image #
        ##################
        # filelist = br.parsed_filelist(dirpath=fpath, string='.h5', ref=5)
        filelist = br.parsed_filelist(dirpath=fpath, string='.h5', ref=1)
        assert len(filelist) > 0, f'no h5 files found in folderpath: {fpath}'

        # set stop
        if stop is None:
            stop = len(filelist) - 1

        # assert stop is number TODO
        # assert start is number TODO
        assert stop >= 0 and start >= 0, f'start ({start}) and stop ({stop}) must be equal or higher than zero' 
        assert stop >= start, f'stop ({stop}) must be equal or bigger than start ({start})'
        assert stop <= len(filelist) -1, f'stop ({stop}) must be equal or smaller than number of images indexes ({len(filelist)-1})'

        # check skip
        for i in skip:
            assert i >= start or i <= stop, f'skip index ({i}) outside of start ({start}) stop ({stop}) image indexes'

        # filter filelist
        filelist = filelist[start:stop + 1]

        # collect data
        x1 = list()
        y1 = list()
        x2 = list()
        y2 = list()
        dummy1 = br.Dummy()
        dummy2 = br.Dummy()
        for j, filepath in enumerate(filelist):
            if j not in skip:
                try:
                    _pe1, _pe2 = _read_rixs(filepath, curv, verbose)
                    x1.extend(_pe1.x)
                    y1.extend(_pe1.y)
                    x2.extend(_pe2.x)
                    y2.extend(_pe2.y)
                    dummy1.append(_pe1)
                    dummy2.append(_pe2)
                except Exception as e:
                    print(f' === ERROR! Image {j} cannot be loaded: {e} ===\n{filepath}')
        
        pe1 = br.PhotonEvents(x=x1, y=y1, xlim=_pe1.xlim, ylim=_pe1.ylim)
        pe2 = br.PhotonEvents(x=x2, y=y2, xlim=_pe2.xlim, ylim=_pe2.ylim)

        # attrs pe1
        for attr in _pe1.get_attrs():
            temp = [getattr(pe, attr) for pe in dummy1]
            setattr(dummy1, attr, temp)
            if attr not in ('modified_date', 'ccd', 'filename'):
                try:
                    setattr(pe1, attr, np.mean(temp))
                    setattr(pe1, attr + '_min', np.min(temp))
                    setattr(pe1, attr + '_max', np.max(temp))
                    setattr(pe1, attr + '_sigma', np.std(temp))
                except:
                    setattr(pe1, attr, None)
        pe1.ccd = 1
        pe1.modified_date = _pe1.modified_date
        pe1.number_of_images = len(filelist)
        try:
            pe1.scan = int(fpath.name)
        except:
            pass

        # attrs pe2
        for attr in _pe2.get_attrs():
            temp = [getattr(pe, attr) for pe in dummy1]
            setattr(dummy2, attr, temp)
            if attr not in ('modified_date', 'ccd', 'filename'):
                try:
                    setattr(pe2, attr, np.mean(temp))
                    setattr(pe2, attr + '_min', np.min(temp))
                    setattr(pe2, attr + '_max', np.max(temp))
                    setattr(pe2, attr + '_sigma', np.std(temp))
                except:
                    setattr(pe2, attr, None)
        pe2.ccd = 2
        pe2.modified_date = _pe2.modified_date
        pe2.number_of_images = len(filelist) - len(skip)
        try:
            pe2.scan = int(fpath.name)
        except:
            pass

        return pe1, pe2, dummy1, dummy2

    #########
    # scans #
    #########
    elif fpath.suffix == '.nxs':
        with h5py.File(Path(fpath), 'r') as f:
            #################
            # get scan type #
            #################
            scan_type = f['entry/instrument/bluesky/metadata/scan_type'][()].decode("utf-8")
            plan_name = f['entry/instrument/bluesky/plan_name'][()].decode("utf-8")
            
            #########
            # attrs #
            #########
            if plan_name == 'adaptive_scan' and 'num_points' in _attrs[scan_type]['raw']:
                del _attrs[scan_type]['raw']['num_points']
            
            metadata = h5.sort_metadata(f=f, attrs_dict=_attrs[scan_type], verbose=verbose)
            try:
                metadata['modified_date'] = br.get_modified_date(fpath)
                metadata['start_time']    = _str2datetime(metadata['start_time'])
                metadata['end_time']      = _str2datetime(metadata['end_time'])
            except:
                pass
            
            ##############
            # get motors #
            ##############
            metadata['motors'] = [m.strip().lstrip('-').strip() for m in metadata['motors'].split('\n') if m.strip().startswith('-')]
            motors = metadata['motors']

            #################
            # get detectors #
            #################
            metadata['detectors'] = [d.strip().lstrip('-').strip() for d in metadata['detectors'].split('\n') if d.strip().startswith('-')]
            detectors = metadata['detectors'] 

            ###################
            # addtional attrs #
            ###################
            try:
                metadata_ = {}
                _group = 'entry/instrument/bluesky/streams/baseline'
                for _key in list(f[_group].keys()):
                    metadata[_key] = f[f'{_group}/{_key}/value'][()]
            except:
                pass
            
            #############
            # read data #
            #############
            _motors = {}
            for motor in motors:
                _motors[motor] = f[f'entry/data/{motor}'][()]
            
            _detectors = {}
            for detector in detectors:
                _detectors[detector] = f[f'entry/data/{detector}'][()]

            ss  = br.Spectra()

            #########
            # xas #
            #########
            if scan_type == 'xas':
                x = f[f'entry/data/energy'][()]
                for _key in _detectors.keys():
                    _s = br.Spectrum(x=x, y=_detectors[_key])
                    _s.detector=_key
                    ss.append(_s)

            #########
            # ascan #
            #########
            if scan_type == 'ascan':
                assert len(motors) == 1, f'number of motors ({len(motors)}) is not compatible with scan type (ascan)'
                if 'main_motor' in metadata:
                    x = _motors[metadata['main_motor']]
                else:
                    x = _motors[list(_motors.keys())[0]]
                for _key in _detectors.keys():
                    _s = br.Spectrum(x=x, y=_detectors[_key])
                    _s.detector=_key
                    ss.append(_s)

            ##########
            # a2scan #
            ##########
            if scan_type == 'a2scan':
                raise NotImplmentedError('read a2scan not implemented')

            ##########
            # a3scan #
            ##########
            elif scan_type == 'a3scan':
                raise NotImplmentedError('read a2scan not implemented')

            ########
            # mesh #
            ########
            elif scan_type == 'mesh':
                assert len(motors) == 2, f'number of motors ({len(motors)}) is not compatible with scan type (mesh)'
                ss = br.Dummy()
                # getting snaking attr
                snake = metadata['snaking'].split('\n')[1:-1]
                metadata['snake'] = [m[2:]=='true' for m in snake]
                metadata['snaking'] = [m[2:]=='true' for m in snake]

                # find `x` motor points
                afinal = np.linspace(metadata['start_x'], metadata['stop_x'], metadata['nstep_x'])

                # find `y` motor points (frozen motor --> always motor_y)
                bfinal = np.linspace(metadata['start_y'], metadata['stop_y'], metadata['nstep_y'])

                # reshape data
                for j, _key in  enumerate(_detectors.keys()):
                    # get intensities
                    y = _detectors[_key]

                    # check if scan finished
                    if len(y) < (len(afinal) * len(bfinal)):
                        # y = np.array(list(y) + [y[-1]]*((len(afinal) * len(bfinal))-len(y)))
                        y = np.array(list(y) + [None]*((len(afinal)*len(bfinal)) - len(y)))
                        metadata['error']  = 'mesh interrupted early'

                    # reashape (check snake)
                    if metadata['snake'][motors.index(metadata['motor_x'])] == True:
                        _s = br.Image(data=[row if i%2 == 0 else row[::-1] for i, row in enumerate(y.reshape(metadata['nstep_y'], metadata['nstep_x']))])
                    else:
                        _s = br.Image(data=[row if i%2 == 0 else row[::] for i, row in enumerate(y.reshape(metadata['nstep_y'], metadata['nstep_x']))])
                    _s.x_centers = afinal
                    _s.y_centers = bfinal
                    _s.detector=_key
                    ss.append(_s)

            metadata.update(metadata_)
                    
            data_group = f['entry/data']
            for s in [_ for _ in ss] + [ss]:
                s.EPOCH = data_group['EPOCH'][()]
                for motor in motors:
                    key = motor + '_user_setpoint'
                    if key in data_group:
                        s.__setattr__('SETPOINT_' + motor, data_group[key][()])
                        
                    key = motor + '_setpoint'
                    if key in data_group:
                        s.__setattr__('SETPOINT_' + motor, data_group[key][()])

                for attr in metadata:
                    setattr(s, attr, metadata[attr])
            ss.create_attr_from_spectra('detector', '_detectors')

        return ss

    ###########
    # XAS old #
    ###########
    else:
        ################
        # get metadata #
        ################
        comments = br.load_comments(filepath=fpath, comment_flag='#', stop_flag='#')
        comments = {line.split(':')[0][2:]:line.split(':')[1][1:-1] for line in comments if ':' in line}
        assert 'scan_type' in comments, 'scan type not found. File corrupted'

        #############
        # XAS old 2 #
        #############
        if comments['scan_type'] == 'xas old':
            #############
            # load file #
            #############
            try:
                data = br.load_data(filepath=fpath, force_array=True)
            except IndexError:
                raise ValueError(f'Error loading file {filepath}')
            
            #############
            # sort data #
            #############
            TEY = br.Spectrum(x=data[:, 0], y=data[:, 2])
            TFY = br.Spectrum(x=data[:, 0], y=data[:, 3])
            I0  = br.Spectrum(x=data[:, 0], y=data[:, 4])
            PD  = br.Spectrum(x=data[:, 0], y=data[:, 5])
            ss  = br.Spectra([TEY, TFY, PD, I0])
            for s in [_ for _ in ss] + [ss]:
                s.PHASE        = data[:, 1]
                s.DVF          = data[:, 6]
                s.SP_ENERGY    = data[:, 7]
                s.SP_PHASE     = data[:, 8]
                s.RING_CURRENT = data[:, 9]
                s.TIMESTAMP    = data[:, 10]

            #########
            # attrs #
            #########
            for line in br.load_comments(fpath):
                if line.startswith('#C # '):
                    if '=' in line:
                        temp = line.split('#C # ')[1].split('=')
                        name  = temp[0].strip()
                        try:
                            value = float(temp[-1].strip())
                        except ValueError:
                            value = temp[-1].strip()
                        for s in [_ for _ in ss] + [ss]:
                            s.__setattr__(name, value)
                elif line.startswith('#S '):
                    scan = line.split('#S ')[1].split()[0].strip()
                    for s in [_ for _ in ss] + [ss]:
                        s.__setattr__('scan', scan)

                    command = ''.join(line.split('#S ')[1].split()[1:]).strip()
                    for s in [_ for _ in ss] + [ss]:
                        s.__setattr__('command', command)
                elif line.startswith('#D '):
                    start_time = _str2datetime(line.split('#D ')[1].strip())
                    for s in [_ for _ in ss] + [ss]:
                        s.__setattr__('start_time', start_time)
            TEY.mode = 'TEY'
            TFY.mode = 'TFY'
            I0.mode  = 'I0'
            PD.mode  = 'PD'
            return ss
        
        ###############
        # XAS Galdino #
        ###############
        elif comments['scan_type'] == 'xas':
            #############
            # load file #
            #############
            try:
                data = br.load_data(filepath=fpath, force_array=True, delimiter=',')
            except IndexError:
                raise ValueError(f'Error loading file: {fpath}')
            
            #############
            # sort data #
            #############
            TEY = br.Spectrum(x=data[:, 0], y=data[:, 2])
            TFY = br.Spectrum(x=data[:, 0], y=data[:, 3])
            I0  = br.Spectrum(x=data[:, 0], y=data[:, 4])
            PD  = br.Spectrum(x=data[:, 0], y=data[:, 5])
            ss  = br.Spectra([TEY, TFY, PD, I0])
            for s in [_ for _ in ss] + [ss]:
                s.PHASE        = data[:, 1]
                s.DVF          = data[:, 6]
                s.SP_ENERGY    = data[:, 7]
                s.SP_PHASE     = data[:, 8]
                s.RING_CURRENT = data[:, 9]
                s.TIMESTAMP    = data[:, 10]

            ############
            # metadata #
            ############
            ss.copy_attrs_from(br.Spectrum(filepath=fpath, usecols=[0, 1]))
            for s in ss:
                s.copy_attrs_from(ss)
            TEY.mode = 'TEY'
            TFY.mode = 'TFY'
            I0.mode  = 'I0'
            PD.mode  = 'PD'
            return ss
        
        #########
        # ERROR #
        #########
        else:
            raise ValueError('not able to identify scan type. File corrupted')

    return
# %%

# %% ============================= RIXS =================================== %% #
def _process(folderpath, sbins, calib=None, norm=True, start=0, stop=None, skip=[]):
    """Returns a dict with objects from each step of the rixs data processing
    
    PhtonEvents for each Image are summed up. The summed up PhotonEvents is 
    turned into one spectrum. The spectrum for each CCD is then aligned and summed. 
    
    Args:
        folderpath(str or path): folderpath with rixs images.
        sbins (int): number of bins for converting photon events to spectrum (
            number of points in the spectrum).
        calib (number, optional): if not None, the x axis is multipled by calib.
            You can give two numbers (linear and constant terms), like 
            calib=[calib, shift].
        norm (bool, optional): if True, spectrum is divided by the exposure time,
            number of images, and number of bins (sbins).
        start, stop, skip (list, optional): For RIXS only. Start and stop are 
            the indexes for the first and last image to sum (inclusive). 
            Default start is 0 and the default for stop is the None, which
            will get up to the last image available. skip should be a list with
            image number indexes to not read (skip). Default is an empty list [].

    Returns:
        dict {'pe1':pe1, 'pe2':pe2, 'pes1':pes1, 'pes2':pes2, 'ss1':ss1, 'ss2':ss2, 's':s}
    """
    #############
    # read file #
    #############
    pe1, pe2, pes1, pes2 = read(folderpath, start=start, stop=stop, skip=skip)
    
    ################
    # ccd spectrum #
    ################
    s1  = pe1.integrated_rows_vs_y_centers(nrows=sbins)
    s2  = pe2.integrated_rows_vs_y_centers(nrows=sbins)
    ss1 = br.Spectra([s1, s2])
    
    #########
    # calib #
    #########
    if isinstance(calib, Iterable):
        ss1.calib = calib[0]
        ss1.shift = -(pe1.Energy - calib[1])
    elif calib:
        ss1.calib = calib
    
    ################
    # sum spectrum #
    ################
    ss2 = ss1.interp().align()
    s   = ss2.calculate_sum()
    
    #########
    # attrs #
    #########
    s.copy_attrs_from(pe1)
    del s.ccd
    setattr(s, 'scan_type', 'rixs')
    
    #################
    # normalization #
    #################
    if norm:
        # s = s.set_factor(1/sum([m[3]-m[2] for m in mask]))
        s = s.set_factor(1/s.RIXSCam_exposure)
        s = s.set_factor(1/s.RIXSCam_NumImages)
        s = s.set_factor(sbins)
        # s = s.set_factor(1000)
    # else:
    #     s = None
    
    return {'pe1':pe1, 'pe2':pe2, 'pes1':pes1, 'pes2':pes2, 'ss1':ss1, 'ss2':ss2, 's':s}

def verify(folderpath, sbins, calib=None, norm=True, **kwargs):
    """open a figure with step-by-step rixs data processing

    Args:
        folderpath(str or path): folderpath with rixs images.
        sbins (int): number of bins for converting photon events to spectrum (
            number of points in the spectrum).
        calib (number, optional): if not None, the x axis is multipled by calib.
            You can give two numbers (linear and constant terms), like 
            calib=[calib, shift].
        norm (bool, optional): if True, spectrum is divided by the exposure time,
            number of images, and number of bins (sbins).
        **kwargs are passed to the scatter plot that plots photon events.
        
    Note:
        Use the argument s=10, to increase the marker size of photon evets plots.
    
    Returns:
        dict {'pe1':pe1, 'pe2':pe2, 'pes1':pes1, 'pes2':pes2, 'ss1':ss1, 'ss2':ss2, 's':s}
    """
    ################
    # process data #
    ################
    temp = _process(folderpath=folderpath, sbins=sbins, calib=calib, norm=norm)
    s    = temp['s']
    pe1  = temp['pe1']
    pe2  = temp['pe2']
    pes1 = temp['pes1']
    pes2 = temp['pes2']

    #######################
    # initial definitions #
    #######################
    pes1.__i  = 0

    ######################
    # change keybindings #
    ######################
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    ###################
    # keyboard events #
    ###################
    def keyboard(event, pes1, pes2, axes):
        if event.key == 'right':
            # increase i
            pes1.__i = pes1.__i + 1
            if pes1.__i >= len(pes1):
                pes1.__i = len(pes1) - 1
    
        elif event.key == 'left':# or event.key == 'down':
            # decrease i
            pes1.__i = pes1.__i - 1
            if pes1.__i < 0:
                pes1.__i = 0
        else:
            return
            
        # clear axis
        axes[0].cla()
        axes[1].cla()
        
        # set labels
        axes[0].set_xlabel('x (pixel)')
        axes[0].set_ylabel('y (pixel)')
        axes[1].set_xlabel('counts/bin')
        
        # change title
        axes[0].set_title('Use left/right keyboard keys to flip through images: ' + str(pes1.__i) + '/' + str(len(pes1)-1), fontsize='small')

        # plot axes 0
        pes1[pes1.__i].plot(ax=axes[0], show_limits=True, s=0.5, **kwargs)
        pes2[pes1.__i].plot(ax=axes[0], show_limits=True, s=0.5, **kwargs)

        # plot axes 1
        pes1[pes1.__i].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])
        pes2[pes1.__i].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])
    
        plt.draw()

    #######################
    # initializing figure #
    #######################
    fig, axes = br.subplots(4, 2, width_ratios=[4, 1], height_ratios=[1, 1, 1, 2], wspace=0.1, hspace=0.8, figsize=(18, 26))
    axes[1].remove_yticklabels()
    axes[3].remove_yticklabels()
    
    ##############
    # share axis #
    ##############
    br.sharey([axes[0], axes[1]])
    br.sharey([axes[2], axes[3]])
    

    ##################
    # error messages #
    ##################
    if pe1.RIXSCam_NumImages != len(pes1):
        fig.suptitle(f'WARNING: # of images ({len(pes1)}) inside folder is different from # of acquired images ({int(pe1.RIXSCam_NumImages)})', color='red')

    ######################
    # set initial titles #
    ######################
    axes[0].set_title('Use left/right keyboard keys to flip through images: ' + str(0) + '/' + str(len(pes1)-1), fontsize='small')
    axes[1].set_title(f'nbins = {sbins}', fontsize='small')
    axes[2].set_title('Summed photon events for each CCD', fontsize='small')
    axes[4].set_title('Number of photons per image', fontsize='small')
    axes[6].set_title('Final spectrum', fontsize='small')
    

    ########
    # plot #
    ########
    # plot initial photon events (axes 0)
    pes1[0].plot(ax=axes[0], show_limits=True, s=0.5, **kwargs)
    pes2[0].plot(ax=axes[0], show_limits=True, s=0.5, **kwargs)
    
    # Inverter o eixo Y no axes[0]
    axes[0].invert_yaxis()
    axes[2].invert_yaxis()
    axes[4].invert_yaxis()

    # plot initial spectra (axes 1)
    pes1[0].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])
    pes2[0].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])

    # plot photon events summed (axes 2)
    pe1.plot(ax=axes[2], show_limits=True, s=0.2, **kwargs)
    pe2.plot(ax=axes[2], show_limits=True, s=0.2, **kwargs)

    # plot spectra summed (axes 3)
    pe1.integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[3])
    pe2.integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[3])

    # plot number of photons per image (axes 4)
    for pes in (pes1, pes2):
        number_of_photons_ccd = [len(_pe) for _pe in pes]
        axes[4].plot(np.arange(0, len(number_of_photons_ccd)), number_of_photons_ccd, marker='o', lw=1)

    # plot spectrum (axes 6)
    s.plot(ax=axes[6], color='black')

    ##############
    # set labels #
    ##############
    for i in (0, 2):
        axes[i].set_xlabel('x (pixel)')
        axes[i].set_ylabel('y (pixel)')
    for i in (1, 3):
        axes[i].set_xlabel('counts/bin')
    axes[4].set_xlabel('Image number')
    axes[4].set_ylabel('Number of photons')

    if calib is None:
        axes[6].set_xlabel('y (pixel)')
    else:
        axes[6].set_xlabel('Energy (eV)')
    if norm:
        axes[6].set_ylabel('Norm. intensity (arb. units)')
    else:
        axes[6].set_ylabel('Photon count per bin')

    ######################
    # register callbacks #
    ######################
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, pes1=pes1, pes2=pes2, axes=axes))
    return temp

# @br.finder.track
def process(folderpath, sbins, calib=None, norm=True, start=0, stop=None, skip=[]):
    """Returns spectrum
    
    PhotonEvents for each Image are summed up. The summed up PhotonEvents is 
    turned into one spectrum. The spectrum for each CCD is then aligned and summed. 
    
    Args:
        folderpath(str or path): folderpath with rixs images.
        sbins (int): number of bins for converting photon events to spectrum (
            number of points in the spectrum).
        calib (number, optional): if not None, the x axis is multipled by calib.
            You can give two numbers (linear and constant terms), like 
            calib=[calib, shift].
        norm (bool, optional): if True, spectrum is divided by the exposure time,
            number of images, and number of bins (sbins).
        start, stop, skip (list, optional): For RIXS only. Start and stop are 
            the indexes for the first and last image to sum (inclusive). 
            Default start is 0 and the default for stop is the None, which
            will get up to the last image available. skip should be a list with
            image number indexes to not read (skip). Default is an empty list [].
        
    Returns:
        Spectrum
    """
    # parameters
    parameters = {}
    parameters['folderpath'] = folderpath
    parameters['sbins']      = sbins
    parameters['calib']      = calib
    parameters['norm']       = norm
    parameters['start']      = start
    parameters['stop']       = stop
    parameters['skip']       = skip
    parameters['ni']         = len(br.filelist(folderpath)) # number of added images

    # # try and find if spectrum has already been calculated
    # s = br.finder.search(parameters=parameters, folderpath=br.finder.folderpath)
    # if s is not None:
    #     return s

    # PROCESS
    d = _process(folderpath, sbins=sbins, calib=calib, norm=norm, start=start, stop=stop, skip=skip)

    # # save spectra so it is not needed to run it again
    # br.finder.save(s=d['s'], parameters=parameters, folderpath=br.finder.folderpath)

    return d['s']

def sequence(folderpath, scans, sbins, calib=True, norm=True):
    """return a list with rixs spectra

    Example:

    >>> # ss1 will have 3 spectra
    >>> ss1 = sequence([100, 101, 102])
    >>>
    >>> # ss2 will have 3 spectra. The middle one will a sum of 2 spectra
    >>> ss2 = sequence([100, [101, 102], 103])

    Args:
        scans (list): list of rixs scan number. Replace a scan number for a list to 
            sum spectra inside list.
        sbins (int): number of bins for converting photon events to spectrum (number 
            of points in the spectrum).
        calib (number, optional): if not None, the x axis is multipled by calib
        norm (bool, optional): if True, spectrum is divided by the exposure time,
            number of images, and number of bins (sbins).
        
    Returns:
        Spectra
    """
    folderpath = Path(folderpath)
    ss = br.Spectra()
    for scan in scans:
        if isinstance(scan, Iterable):
            _ss = br.Spectra()
            for _scan in scan:
                _ss.append(process(folderpath/str(_scan).zfill(4), sbins=sbins, calib=calib, norm=norm))
            # fix attrs
            for attr in _ss[0].get_attrs():
                try:
                    _ss.create_attr_from_spectra(attr)        

                    if attr.endswith('_min'):
                        _ss.__setattr__(attr, min(_ss.__getattr__(attr)))
                    elif attr.endswith('_max') or attr.endswith('_sigma'):
                        _ss.__setattr__(attr, max(_ss.__getattr__(attr)))
                    else:
                        _ss.__setattr__(attr, np.mean(_ss.__getattr__(attr)))
                except:
                    pass
            ss.append(_ss.align().interp().calculate_sum())

        else:
            ss.append(process(folderpath/str(scan).zfill(4), sbins=sbins, calib=calib, norm=norm))

    #########
    # attrs #
    #########
    for attr in ss[0].get_attrs():
        try:
            ss.create_attr_from_spectra(attr)
        except:
            pass

    return ss

# %% =========================== alignment plot =============================== %% #
def alignment(folderpath, scans, atype, sbins=2000, calib=None, norm=False, start=0, stop=None, skip=[], limits=None, external_list=None, **kwargs):
    """Plot alignment analisys
    
    Use ss.create_attr_from_spectra(attr) to substitute atype fot attr
    
    Args:
        folderpath (str or path): folderpath with rixs images.
        scans (list): scans IDs. Ex [1, 534, 458].
        atype (string): select a alignment type attribute from metadata or 'external' if it 
            is external data. There are also nicknames for some  alignemnt types as 'r1', 'r2', 'z' and 'time'.
        sbins (int): number of bins for converting photon events to spectrum (
            number of points in the spectrum).
        calib (number or list, optional): if not None, the x axis is multipled by calib.
            You can give two numbers (linear and constant terms), like calib=[calib, shift].
        norm (bool, optional): if True, spectrum is divided by the exposure time,
            number of images, and number of bins (sbins).
        start, stop, skip (list, optional): For RIXS only. Start and stop are 
            the indexes for the first and last image to sum (inclusive). 
            Default start is 0 and the default for stop is the None, which
            will get up to the last image available. skip should be a list with
            image number indexes to not read (skip). Default is an empty list [].
        limits (list, optional): start and stop to set the limits of the x axis of the plot.
        external_list (list): list of values to user if the atype is \'external\'. 
        
    Returns:
        dict {
            'ss1':    {'ss1':ss_ccd1,'ss2':ss_ccd2, 'fit1':fit_ccd1, 'fit2':fit_ccd2},
            'fwhm':   {'x1':x1,'x2':x2, 'y1':y_ccd1, 'y2':y_ccd2},
            'amp':    {'x1':x1,'x2':x2, 'y1':y_ccd1, 'y2':y_ccd2},
            'center': {'x1':x1,'x2':x2, 'y1':y_ccd1, 'y2':y_ccd2},
            'offset': {'x1':x1,'x2':x2, 'y1':y_ccd1, 'y2':y_ccd2},
            'popt1':  popt_ccd1,
            'popt2':  popt_ccd2,
            'x1':     x1,
            'x2':     x2
        }
        
        x1 is the attr
        x2 exists in special casesof atype = r1 or r2
    """
    attr = [None,None]
    align_list = [list(),list()]
    if atype == 'external':
        assert external_list is not None, 'If atype=\'motor\' you must put a external_list=[pos1, pos2]'
        assert len(external_list) == len(scans), 'external_list and scans must have the same len'
    elif atype == 'r1':
        attr[0] = 'RIXS_GZ'
        attr[1] = 'RIXS_DZ'
    elif atype == 'r2':
        attr[0] = 'RIXS_DZ'
        attr[1] = 'RIXS_DY'
    elif atype == 'z':
        attr[0] = 'RIXS_MZ'
        attr[1] = attr[0]
    elif atype == 'time':
        attr[0] = 'modified_date'
        attr[1] = attr[0]
    else:
        attr[0] = atype
        attr[1] = attr[0]
    
    basepath = Path(folderpath)
    
    ss   = br.Spectra()
    popt = list()
    ss1  = br.Spectra()
    ss2  = br.Spectra()
    
    si    = br.Spectra([br.Spectrum()]*2)
    ssi   = [ss1, ss2]
    fiti  = [br.Spectra(),br.Spectra()]
    popti = [[],[]]
    
    ss1.__i = 0
    ss2.__i = 0
    
    ######################
    # change keybindings #
    ######################
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    ###################
    # keyboard events #
    ###################
    def keyboard(event, ssi, fiti, popti, ax, x, limits=None):
        for i in range(2):
            if event.key == 'right':
                # increase i
                ssi[i].__i = ssi[i].__i + 1
                if ssi[i].__i >= len(ssi[i]):
                    ssi[i].__i = len(ssi[i]) - 1
        
            elif event.key == 'left':# or event.key == 'down':
                # decrease i
                ssi[i].__i = ssi[i].__i - 1
                if ssi[i].__i < 0:
                    ssi[i].__i = 0
            else:
                return
                    
        # clear axis
        ax.cla()
        
        # set labels
        if calib is None:
            ax.set_xlabel('pixel')
        else:
            ax.set_xlabel('Energy (eV)')
        ax.set_ylabel('counts/bin')
        
        # change title
        ax.set_title('Use left/right keyboard keys to flip through images: ' + str(ssi[0].__i) + '/' + str(len(ssi[0])-1), fontsize='small')

        # plot axes 0        
        ax1.plot([], [], marker='', linestyle='', label='ccd  scan   pos', color='black')
        ssi[0][ssi[i].__i].plot(ax=ax1, marker='.', label=f'  {1}   {str(ssi[0][ssi[i].__i].scan).zfill(4)}  {x1[ssi[i].__i]}')
        ssi[1][ssi[i].__i].plot(ax=ax1, marker='.', label=f'  {2}   {str(ssi[1][ssi[i].__i].scan).zfill(4)}  {x1[ssi[i].__i]}')
        
        fiti[0][ssi[i].__i].plot(ax=ax1, linestyle='--', color='darkblue')
        fiti[1][ssi[i].__i].plot(ax=ax1, linestyle='--', color='red')
        
        c = 7.5*max([popti[0][ssi[i].__i][2], popti[1][ssi[i].__i][2]])
        if limits is None or limits == []:
            ax.set_xlim(min([popti[0][ssi[i].__i][1], popti[1][ssi[i].__i][1]])-c, max([popti[0][ssi[i].__i][1], popti[1][ssi[i].__i][1]])+c)
        else:
            ax.set_xlim(limits[0], limits[1])
        ax.legend()
        
        plt.draw()
    
    if limits is None or limits == []:
        limits=np.array([500, 1300])
        if isinstance(calib, Iterable):
            limits = (limits-calib[1])*calib[0]
        elif calib:
            limits = limits*calib  
        
    for scan in scans:
        folderpath = basepath/str(scan).zfill(4)
        assert folderpath.exists(), f"Folderpath does not exist, {folderpath}"

        tmp = _process(folderpath, sbins=sbins, calib=calib, norm=norm, start=start, stop=stop, skip=skip)
        s = tmp['s']
        
        for i in range(2):
            si[i] = tmp['ss1'][i]
            ssi[i].append(si[i])
            fit, popt, R2, model = si[i].fit_peak(fixed_m=0, asymmetry=False, limits=limits)
            fiti[i].append(fit)
            popti[i].append(popt)
        if atype != 'external':
            assert hasattr(s, attr[0]), f"\'{atype}\' is not in {s.get_attrs()}"
            assert hasattr(s, attr[1]), f"\'{atype}\' is not in {s.get_attrs()}"
            align_list[0].append(getattr(s, attr[0]))
            align_list[1].append(getattr(s, attr[1]))
            
        ss.append(s)

    if atype != 'external':
        x1 = align_list[0]
        x2 = align_list[1]
        x1_label = attr[0]
        x2_label = attr[1]

    else:
        x1 = external_list
        x2 = external_list
        x1_label = 'external'
        x2_label = 'external'
    
    # Ordenar x, y e z juntos
    _sorted = sorted(zip(x1, x2, popti[0], popti[1]), key=lambda item: item[0])

    # Separar x, y e z_sorted
    x1_sorted, x2_sorted, popt0_sorted, popt1_sorted = zip(*_sorted)

    # Separar popt_sorted em 4 colunas
    popt0_new = [list(column) for column in zip(*popt0_sorted)]
    popt1_new = [list(column) for column in zip(*popt1_sorted)]

    fig, axs = plt.subplots(3, 2, figsize=(10, 12))
    fig.suptitle(f'{atype.upper()} - Aligment', fontweight="bold")
    
    ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((3, 2), (1, 0))
    ax3 = plt.subplot2grid((3, 2), (1, 1))
    ax4 = plt.subplot2grid((3, 2), (2, 0))
    ax5 = plt.subplot2grid((3, 2), (2, 1))
    
    
    ax1.set_title('Use left/right keyboard keys to flip through images: 0/' + str(len(scans)-1), fontsize='small')
    ax1.plot([], [], marker='', linestyle='', label='ccd  scan   pos', color='black')
    ssi[0][0].plot(ax=ax1, marker='.', label=f'  {1}   {str(ssi[0][0].scan).zfill(4)}  {x1[0]}')
    ssi[1][0].plot(ax=ax1, marker='.', label=f'  {2}   {str(ssi[1][0].scan).zfill(4)}  {x1[0]}')
    fiti[0][0].plot(ax=ax1, linestyle='--', color='darkblue')
    fiti[1][0].plot(ax=ax1, linestyle='--', color='red')
    ax1.legend()
    c = 7.5*max([popti[0][0][2], popti[1][0][2]])
    if limits is None or limits == []:
        ax1.set_xlim(min([popti[0][0][1], popti[1][0][1]])-c, max([popti[0][0][1], popti[1][0][1]])+c)
    else:
        ax1.set_xlim(limits[0], limits[1])
    ax1.set_ylabel('counts/bin')
    if calib is None:
        ax1.set_xlabel('pixel')
    else:
        ax1.set_xlabel('Energy (eV)')

    ax2.set_title('fwhm', fontweight="bold")
    ax2.plot(x1_sorted, popt0_new[2], 'o-', label='ccd 1')
    ax2.plot(x1_sorted, popt1_new[2], 'o-', label='ccd 2')
    ax2.set_xlabel(x1_label)
    if attr[1]!=attr[0]:
        ax22 = ax2.twiny()  # Criando um eixo secundário
        ax22.plot(x2_sorted, popt0_new[2], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
        ax22.plot(x2_sorted, popt1_new[2], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
        ax22.set_xlabel(x2_label)
        ax22.xaxis.set_label_coords(0.5, 0.90)
    ax2.legend()
    if calib is None:
        ax2.set_ylabel('pixel')
    else:
        ax2.set_ylabel('Energy (eV)')
    
    ax3.set_title('center', fontweight="bold")
    ax3.plot(x1_sorted, popt0_new[1], 'o-', label='ccd 1')
    ax3.plot(x1_sorted, popt1_new[1], 'o-', label='ccd 2')
    ax3.set_xlabel(x1_label)
    if attr[1]!=attr[0]:
        ax32 = ax3.twiny()  # Criando um eixo secundário
        ax32.plot(x2_sorted, popt0_new[1], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
        ax32.plot(x2_sorted, popt1_new[1], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
        ax32.set_xlabel(x2_label)
        ax32.xaxis.set_label_coords(0.5, 0.90)
    if calib is None:
        ax3.set_ylabel('pixel')
    else:
        ax3.set_ylabel('Energy (eV)')
    ax3.legend()
    
    ax4.set_title('amplitude', fontweight="bold")
    ax4.plot(x1_sorted, popt0_new[0], 'o-', label='ccd 1')
    ax4.plot(x1_sorted, popt1_new[0], 'o-', label='ccd 2')
    ax4.set_xlabel(x1_label)
    if attr[1]!=attr[0]:
        ax42 = ax4.twiny()  # Criando um eixo secundário
        ax42.plot(x2_sorted, popt0_new[0], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
        ax42.plot(x2_sorted, popt1_new[0], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
        ax42.set_xlabel(x2_label)
        ax42.xaxis.set_label_coords(0.5, 0.90)
    ax4.set_ylabel('counts/bin')
    ax4.legend()
    
    ax5.set_title('offset', fontweight="bold")
    ax5.plot(x1_sorted, popt0_new[3], 'o-', label='ccd 1')
    ax5.plot(x1_sorted, popt1_new[3], 'o-', label='ccd 2')
    ax5.set_xlabel(x1_label)
    if attr[1]!=attr[0]:
        ax52 = ax5.twiny()  # Criando um eixo secundário
        ax52.plot(x2_sorted, popt0_new[3], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
        ax52.plot(x2_sorted, popt1_new[3], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
        ax52.set_xlabel(x2_label)
        ax52.xaxis.set_label_coords(0.5, 0.90)
    ax5.set_ylabel('counts/bin')
    ax5.legend()
    
    # Formatação do eixo de tempo
    if atype == 'time':
        import matplotlib.dates as mdates
        # Função para formatar o eixo com datas
        def format_date_axis(ax):
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))  # Formato de hora
            ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=5))  # Intervalos de 5 minutos
            fig.autofmt_xdate()  # Ajusta as datas

        # Plotando nos outros eixos (com data)
        for ax in [ax2, ax3, ax4, ax5]:
            format_date_axis(ax)
    
    plt.subplots_adjust(hspace=0.7, wspace=0.4)
    
    ######################
    # register callbacks #
    ######################
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ssi=ssi, fiti=fiti, popti=popti, ax=ax1, x=x1, limits=limits))
    
    output = {}
    output['ssi']    = {'ss1':ssi[0],'ss2':ssi[1], 'fit1':fiti[0], 'fit2':fiti[1]}
    output['amp']    = {'x1':x1_sorted,'x2':x2_sorted, 'y1':popt0_new[0], 'y2':popt1_new[0]}
    output['center'] = {'x1':x1_sorted,'x2':x2_sorted, 'y1':popt0_new[1], 'y2':popt1_new[1]}
    output['fwhm']   = {'x1':x1_sorted,'x2':x2_sorted, 'y1':popt0_new[2], 'y2':popt1_new[2]}
    output['oddset'] = {'x1':x1_sorted,'x2':x2_sorted, 'y1':popt0_new[3], 'y2':popt1_new[3]}
    output['popt1']  = popt0_new
    output['popt2']  = popt1_new
    output['x1']     = x1_sorted
    output['x2']     = x2_sorted
    
    return output

# %% =========================== curvature parameters =============================== %% #
def curvature(folderpath, ccd, ncols=10, nrows=1000, deg=2, ylimits=None, xlimits=None, popt=None, offset=None, figsize=(50, 10)):
    """Calculate curvature

    Args:
        folderpath
        ccd (int): ccd number (1 or 2)
        ncols, nrows (int, optional): horizontal and vertical number of bins. 
            Default is ncols=10 and nrows=1000
        deg (int, optional): polynomial degree for fiting the curvature. Default is 2
        xlimits, ylimits (tuple): start and stop values for calculating shifts
        popt (list, optional): if not None, the optimal parameters for curvature
            will not be calculated and popt will be used instead.
        offset (number, optional): vertical offset for ploting the fited curvature
            on the third panel (red curve). If None, it will be calculated.
        figsize (tuple, optional): figure size for ploting.

    Returns:
        popt for ccd1 and ccd2
    """
    #####################
    # initialize figure #
    #####################
    fig, axes = br.subplots(1, 5, sharey='row', figsize=figsize)
    fig.subplots_adjust(top=0.99, bottom=0.05, left=0.05, right=0.99)

    #############
    # read file #
    #############
    pe1, pe2, pes1, pes2 = read(fpath=folderpath, curv=False)
    if ccd == 1:
        pe = pe1
    elif ccd == 2:
        pe = pe2
    else:
        raise ValueError('wrong ccd input. valid ccd is 1 or 2')

    ########
    # crop #
    ########
    if xlimits is not None:
        pe = pe.crop(x_start=xlimits[0], x_stop=xlimits[1])
        pe.xlim = xlimits

    ###########
    # binning #
    ###########
    im = pe.binning(ncols=ncols, nrows=nrows)

    ####################
    # calculate shifts #
    ####################
    if popt is None:
        s, fit, popt, R2, model = pe.calculate_vertical_shift_curvature(ncols=ncols, nrows=nrows, deg=deg, mode='cc', ylimits=ylimits, limit_size=1000)
        if xlimits is not None:
            fit = fit.crop(xlimits[0], xlimits[1])
        
    ########
    # plot #
    ########

    # plot photon events (raw)
    #pe.plot(axes[0], color='black')

    # plot reduced image
    pos = im.plot(axes[1])

    # plot photon events (raw) with vertical bins
    pe.plot(axes[0], color='black')
    br.axvlines(ax=axes[0], x=pos.x_edges, colors='red', lw=.5)

    # plot horizontal integration of each vertical bin
    cols = im.columns.switch_xy()#.flip_y()
    cols.plot(axes[2])

    # get max and min y (makes plot nicer)
    cols2 = im.columns
    arr100, _popt, err, model = cols2[0].fit_peak()
    ymin = _popt[1] - _popt[2]*12

    arr100, _popt, err, model = cols2[-1].fit_peak()
    ymax = _popt[1] + _popt[2]*12

    # plot fitting
    offset = cols[0].y[np.argmax(cols[0].x)]
    pe.plot(axes[3], color='black')
    fit.plot(axes[3], factor=-1, offset=offset, color='red')

    # set shifts
    pe.plot(axes[4], color='black')

    # fix curvature
    pe = pe.set_vertical_shift_via_polyval(p=popt)
    pe.plot(axes[4], color='red')

    # set y lim
    for i in range(5):
        axes[i].set_ylim(ymin, ymax)
        axes[i].set_xlabel(f'x subpixel')

    axes[0].set_ylabel(f'ccd {ccd}: y subpixel')

    return popt

# %%

# %% =========================== find calib =============================== %% #
def find_calib(folderpath, scans, sbins=2000, limits=None):
    from lmfit.models import LinearModel
    var = alignment(folderpath, scans, atype='Energy', sbins=sbins, calib=None, limits=limits)['center']
    model = LinearModel()
    params = {}
    result = {}
    x = {}
    y = {}
    _calib = {}
    _pixel = {}
    
    for i in ('1','2'):
        x[i] = np.array(var[f'x{i}'])
        y[i] = np.array(var[f'y{i}'])
        params[i] = model.make_params(slope=1, intercept=0)
        result[i] = model.fit(y[i], params[i], x=x[i])

        slope = result[i].params['slope'].value
        intercept = result[i].params['intercept'].value

        _calib[i] = 1 / slope
        _pixel[i] = -intercept / slope
    
    _calib = [np.average(list(_calib.values())), np.average(list(_pixel.values()))]
    print(f'calib = {_calib}')
    return _calib


# %% EXPERIMENTAL EXPERIMENTAL EXPERIMENTAL EXPERIMENTAL EXPERIMENTAL EXPERIMENTAL 

# %% =========================== fancy plot =============================== %% #
def _plotall(filepath, scan, axes=None, set_window_size=True):
    TEY, TFY, I0, PD = read(filepath, scan) 

    if axes is None:
        fig, axes = br.subplots(2, 3)
    if set_window_size:
        br.set_window_size((514, 1091)) 
        plt.subplots_adjust(left=0.06, right=0.99, top=0.95, bottom=0.1)

    seq = [TEY, MCP, TFY, TEY/RMU, MCP/RMU, TFY/RMU]
    seq[3]. label = 'TEY/RMU'
    seq[4]. label = 'MCP/RMU'
    seq[5]. label = 'TFY/RMU'
    # seq[6]. label = 'norm TEY/RMU'
    # seq[7]. label = 'norm MCP/RMU'
    # seq[8]. label = 'norm TFY'
    for i in range(2*3):
        ax = axes[i]
        seq[i].plot(ax=ax, label=f'{seq[i].label}, #{seq[i].scan}, x={round(seq[i].sample_x, 4)}, y={round(seq[i].sample_y, 4)}')
        ax.labels_xas()
        br.leg(ax=ax, fontsize='xx-small')
    return axes

def _sequential(*args, **kwargs):

    ################
    # process data #
    ################
    temp = _process(folderpath=folderpath, sbins=sbins, calib=calib, norm=norm)
    s    = temp['s']
    pe1  = temp['pe1']
    pe2  = temp['pe2']
    pes1 = temp['pes1']
    pes2 = temp['pes2']

    #######################
    # initial definitions #
    #######################
    pes1.__i  = 0

    ######################
    # change keybindings #
    ######################
    try:
        matplotlib.rcParams['keymap.back'].remove('left')
        matplotlib.rcParams['keymap.forward'].remove('right')
    except ValueError:
        pass

    ###################
    # keyboard events #
    ###################
    def keyboard(event, pes1, pes2, axes):
        if event.key == 'right':
            # increase i
            pes1.__i = pes1.__i + 1
            if pes1.__i >= len(pes1):
                pes1.__i = len(pes1) - 1
    
        elif event.key == 'left':# or event.key == 'down':
            # decrease i
            pes1.__i = pes1.__i - 1
            if pes1.__i < 0:
                pes1.__i = 0
        else:
            return
            
        # clear axis
        axes[0].cla()
        axes[1].cla()
        
        # set labels
        axes[0].set_xlabel('x (pixel)')
        axes[0].set_ylabel('y (pixel)')
        axes[1].set_xlabel('counts/bin')
        
        # change title
        axes[0].set_title('Use left/right keyboard keys to flip through images: ' + str(pes1.__i) + '/' + str(len(pes1)-1), fontsize='small')

        # plot axes 0
        pes1[pes1.__i].plot(ax=axes[0], show_limits=True, **kwargs)
        pes2[pes1.__i].plot(ax=axes[0], show_limits=True, **kwargs)

        # plot axes 1
        pes1[pes1.__i].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])
        pes2[pes1.__i].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])
    
        plt.draw()

    #######################
    # initializing figure #
    #######################
    fig, axes = br.subplots(4, 2, width_ratios=[4, 1], height_ratios=[1, 1, 1, 2], wspace=0.1, hspace=0.8, figsize=(18, 26))
    axes[1].remove_yticklabels()
    axes[3].remove_yticklabels()
    
    ##############
    # share axis #
    ##############
    br.sharey([axes[0], axes[1]])
    br.sharey([axes[2], axes[3]])
    

    ##################
    # error messages #
    ##################
    if pe1.RIXSCam_NumImages != len(pes1):
        fig.suptitle(f'ERROR: # of images ({len(pes1)}) inside folder is different from # of acquired images ({int(pe1.RIXSCam_NumImages)})', color='red')

    ######################
    # set initial titles #
    ######################
    axes[0].set_title('Use left/right keyboard keys to flip through images: ' + str(0) + '/' + str(len(pes1)-1), fontsize='small')
    axes[1].set_title(f'nbins = {sbins}', fontsize='small')
    axes[2].set_title('Summed photon events for each CCD', fontsize='small')
    axes[4].set_title('Number of photons per image', fontsize='small')
    axes[6].set_title('Final spectrum', fontsize='small')
    

    ########
    # plot #
    ########
    # plot initial photon events (axes 0)
    pes1[0].plot(ax=axes[0], show_limits=True, **kwargs)
    pes2[0].plot(ax=axes[0], show_limits=True, **kwargs)

    # plot initial spectra (axes 1)
    pes1[0].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])
    pes2[0].integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[1])

    # plot photon events summed (axes 2)
    pe1.plot(ax=axes[2], show_limits=True, **kwargs)
    pe2.plot(ax=axes[2], show_limits=True, **kwargs)

    # plot spectra summed (axes 3)
    pe1.integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[3])
    pe2.integrated_rows_vs_y_centers(nrows=sbins).switch_xy().plot(ax=axes[3])

    # plot number of photons per image (axes 4)
    for pes in (pes1, pes2):
        number_of_photons_ccd = [len(_pe) for _pe in pes]
        axes[4].plot(np.arange(0, len(number_of_photons_ccd)), number_of_photons_ccd, marker='o', lw=1)

    # plot spectrum (axes 6)
    s.plot(ax=axes[6], color='black')

    ##############
    # set labels #
    ##############
    for i in (0, 2):
        axes[i].set_xlabel('x (pixel)')
        axes[i].set_ylabel('y (pixel)')
    for i in (1, 3):
        axes[i].set_xlabel('counts/bin')
    axes[4].set_xlabel('Image number')
    axes[4].set_ylabel('Number of photons')

    if calib is None:
        axes[6].set_xlabel('y (pixel)')
    else:
        axes[6].set_xlabel('Energy (eV)')
    if norm:
        axes[6].set_ylabel('Norm. intensity (arb. units)')
    else:
        axes[6].set_ylabel('Photon count per bin')

    ######################
    # register callbacks #
    ######################
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, pes1=pes1, pes2=pes2, axes=axes))
    return temp

# from brixs.addons.png2clipboard import png2clipboard
# def figure2clipboard():
#     plt.savefig(TMP/'temp.png', dpi=1000)
#     png2clipboard(TMP/'temp.png')
