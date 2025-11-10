#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core functions for IPE beamline - Sirius.

See the advanced.py file for higher level functions

Last edited: Carlos Galdino 2025-08-09
"""

# %% ========================== Standard Imports ========================= %% #
from pathlib import Path
import numpy as np
import datetime

# %% =============================== brixs =============================== %% #
import brixs as br
import brixs.addons.h5 as h5

# %% ========================== Special Imports ========================== %% #
try:
    import h5py
except:
    pass
# %%

# %% ============================= support =============================== %% #
def _str2datetime(string):
    """convert IPE date/time string pattern to date
    
    Example:
        '2022/07/20 21:08:36.0'
        '2022/07/20T21:08:36.0'
        '2022-07-20 21:08:36.0'
        '2022-07-20T21:08:36.0'

    Args:
        string (str): string with IPE date string

    Returns:
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

def scanlist(folderpath):
    """Return list of scans available in folderpath"""
    folderpath = Path(folderpath)
    assert folderpath.exists(), f'fpath does not exist ({folderpath})'
    return br.parsed_filelist(folderpath, string='*', ref=0, return_type='dict')
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
# %%

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
            ss._copy_attrs_from(br.Spectrum(filepath=fpath, usecols=[0, 1]))
            for s in ss:
                s._copy_attrs_from(ss)
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
