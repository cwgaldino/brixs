#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Core functions for ID32 beamline at ESRF"""

# %% ------------------------- Standard Imports --------------------------- %% #
from collections.abc import Iterable
from pathlib import Path
import numpy as np
import datetime

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br

# %% ------------------------- Special Imports ---------------------------- %% #
try:
    import h5py
    from silx.io.specfile import SpecFile
except:
    pass
# %%

# %% ========================= metadata table ============================ %% #
metadata = {'rixs': {'E': 'energy', 
                     'H': 'H', 
                     'K': 'K', 
                     'L': 'L',
                     'th': 'th', 
                     'tth': 'tth', 
                     'phi': 'phi',
                     'chi': 'chi', 
                     'slit': 'exagap',
                     'motor_x': 'xsam', 
                     'motor_y': 'ysam', 
                     'motor_z': 'zsam',
                     'sample_x': 'beamx', 
                     'sample_y': 'ysam',
                     'sample_z': 'beamz'},
            'xmcd': {'E': 'energy', 
                     'magnetic_field': 'magnet', 
                     'slit': 'exbgap',
                     'th': 'srot', 
                     'sample_y': 'sy', 
                     'sample_z': 'sz'}}   

additional = {'rixs': {'TemperatureA': 'lak332A_rixs',
                       'TemperatureB': 'lak332B_rixs'},
              'xmcd': {}
              }
# %%

# %% ============================= support =============================== %% #
def _str2datetime(string):
    """convert ID32 date string pattern to date --> '2025-04-10T18:39:13.331583+02:00'

    Args:
        string (str): string with ID32 date string

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
    year, month, day = (int(_) for _ in date.split('-'))

    ########
    # time #
    ########
    hour, minute, seconds = time.split('+')[0].split(':')

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

def _polarization(value):
    """returns polarization string based on undulator gap"""
    pol = value

    if isinstance(value, Iterable):
        value = np.mean(value)
        
    if value < -5:
        pol = 'CM'
    elif -5 <= value < 5:
        pol = 'LH'
    elif 5 <= value < 30:
        pol = 'CP'
    elif value >= 31:
        pol = 'LV'  
    
    return pol

def list_samples(TOP):
    """returns a list with available sample names
    
    Args:
        TOP (str or path): top/main experiment folderpath where folder such as
            RAW_DATA and PROCESSED_DATA can be found.
    
    Returns:
        list
    """
    TOP = Path(TOP)
    assert TOP.exists(), f'Top directory does not exist TOP="{TOP}"'
    assert 'RAW_DATA' in [_.name for _ in br.filelist(TOP)], f'folder RAW_DATA cannot be found in TOP="{TOP}". TOP must be the top/main folderpath of the experiment'
    return [_.name for _ in br.filelist(TOP/'RAW_DATA') if _.is_dir()]

def list_datasets(TOP, sample):
    """returns a list with available datasets for a sample
    
    Args:
        TOP (str or path): top/main experiment folderpath where folder such as
            RAW_DATA and PROCESSED_DATA can be found.
        sample (str): sample name.
    
    Returns:
        list
    """
    TOP = Path(TOP)
    available_samples = list_samples(TOP)
    assert sample in available_samples, f'sample "{sample}" cannot be found. Available options are: {available_samples}'
    return [_.name for _ in br.filelist(TOP/'RAW_DATA'/sample) if _.is_dir()]

def list_scans(TOP=None, sample=None, dataset=None, filepath=None):
    """"Returns list of available scans for a given sample

    Usage:
        import brixs.beamlines.id32 as id32
        scanlist = id32.list_scans(TOP=<TOP>, sample=<sample>)
        scanlist = id32.list_scans(<TOP>, <sample>)
        scanlist = id32.list_scans(TOP=<TOP>, sample=<sample>, dataset=<dataset>)
        scanlist = id32.list_scans(<TOP>, <sample>, <dataset>)
        scanlist = id32.list_scans(filepath=<filepath_to_h5>)
    
    Args:
        TOP (str or path): top/main experiment folderpath where folder such as
            RAW_DATA and PROCESSED_DATA can be found.
        sample (str): sample name.
        dataset (str or None, optional): dataset name. If None, a dictionary
            will be returned with a list of scans for each dataset.
        filepath (str or path, optional): this overwrites all other arguments.
            Directly reads an .h5 file. This function expects the .h5 file from
            a dataset, e.g., TOP/RAW_DATA/<sample>/<dataset>/dataset.h5
    
    Returns:
        list or dict (if dataset=None)
    """   
    # read filepath directly
    if filepath is not None:
        with h5py.File(Path(filepath)) as sf:
            out = br.remove_duplicates([int(_.split('.')[0]) for _ in sf.keys()])
        out = np.array(br.sort(out, out))
        return out

    # search the appropriate file inside the TOP directory
    assert any([TOP is None, sample is None])==False or any([filepath is None, ])==False, f'Wrong input. Arguments are `TOP`, `sample`, `dataset`, `filepath`. At least `TOP` and `sample` must be given or a `filepath`. Please see help(id32.list_scans)'
    TOP = Path(TOP)
    available_datasets = list_datasets(TOP=TOP, sample=sample)

    if dataset is None:
        out = {}
        for dataset in available_datasets:
            out[dataset] = list_scans(TOP=TOP, sample=sample, dataset=dataset)
    else:
        assert dataset in available_datasets, f'dataset "{dataset}" cannot be found. Available options are: {available_datasets}'

        # filepath
        filepath = [_ for _ in br.filelist(TOP/'RAW_DATA'/sample/dataset, string='.h5') if _.is_file()][0]

        # read filepath
        with h5py.File(str(filepath)) as sf:
            out = br.remove_duplicates([int(_.split('.')[0]) for _ in sf.keys()])
        out = np.array(br.sort(out, out))
    return out

def list_available_detectors(TOP=None, sample=None, dataset=None, scan=None, filepath=None):
    """return a list of available detectors for a scan (hdf5 address: <scan>/measurement)

    Usage:
        import brixs.beamlines.id32 as id32
        detectors = id32.list_available_detectors(TOP=<TOP>, sample=<sample>, dataset=<dataset>, scan=<scan>)
        detectors = id32.list_available_detectors(<TOP>, <sample>, <dataset>, <scan>)
        detectors = id32.list_available_detectors(filepath=<filepath_to_h5>, scan=<scan>)

    Args:
        TOP (str or path): top/main experiment folderpath where folder such as
            RAW_DATA and PROCESSED_DATA can be found.
        sample (str): sample name.
        dataset (str or None, optional): dataset name. If None, a dictionary
            will be returned with a list of scans for each dataset.
        scan (int): scan number.
        filepath (str or path, optional): this overwrites all other arguments.
            Directly reads an .h5 file. This function expects the .h5 file from
            a dataset, e.g., TOP/RAW_DATA/<sample>/<dataset>/dataset.h5
    
    Returns:
        list
    """
    assert scan is not None, '`scan` was not be given'
    assert any([TOP is None, sample is None, dataset is None, scan is None])==False or any([filepath is None, scan is None])==False, f'Wrong input. Arguments are `TOP`, `sample`, `dataset`, `filepath`. Please see help(id32.list_available_detectors)'

    # get available scans
    available_scans = list_scans(TOP=TOP, sample=sample, dataset=dataset, filepath=filepath)
    assert scan in available_scans, f'scan {scan} cannot be found. Available scans for {sample}/{dataset}: {available_scans}'
    
    # get filepath
    if filepath is None:
        filepath = [_ for _ in br.filelist(TOP/'RAW_DATA'/sample/dataset, string='.h5') if _.is_file()][0]
    
    # open file
    with h5py.File(str(filepath)) as sf:
        # get initial metadata
        initial    = sf[str(scan+0.1)]

        detectors = list(initial['measurement'].keys())
    return detectors

def list_available_metadata(TOP=None, sample=None, dataset=None, scan=None, filepath=None):
    """return a list of available metadata for a scan (hdf5 address: <scan>/instrument/positioners)

    Usage:
        import brixs.beamlines.id32 as id32
        detectors = id32.list_available_metadata(TOP=<TOP>, sample=<sample>, dataset=<dataset>, scan=<scan>)
        detectors = id32.list_available_metadata(<TOP>, <sample>, <dataset>, <scan>)
        detectors = id32.list_available_metadata(filepath=<filepath_to_h5>, scan=<scan>)

    Args:
        TOP (str or path): top/main experiment folderpath where folder such as
            RAW_DATA and PROCESSED_DATA can be found.
        sample (str): sample name.
        dataset (str or None, optional): dataset name. If None, a dictionary
            will be returned with a list of scans for each dataset.
        scan (int): scan number.
        filepath (str or path, optional): this overwrites all other arguments.
            Directly reads an .h5 file. This function expects the .h5 file from
            a dataset, e.g., TOP/RAW_DATA/<sample>/<dataset>/dataset.h5
    
    
    Returns:
        list
    """
    assert scan is not None, '`scan` was not be given'
    assert any([TOP is None, sample is None, dataset is None, scan is None])==False or any([filepath is None, scan is None])==False, f'Wrong input. Arguments are `TOP`, `sample`, `dataset`, `filepath`. Please see help(id32.list_available_detectors)'

    # get available scans
    available_scans = list_scans(TOP=TOP, sample=sample, dataset=dataset, filepath=filepath)
    assert scan in available_scans, f'scan {scan} cannot be found. Available scans for {sample}/{dataset}: {available_scans}'
    
    # get filepath
    if filepath is None:
        filepath = [_ for _ in br.filelist(TOP/'RAW_DATA'/sample/dataset, string='.h5') if _.is_file()][0]
    
    # open file
    with h5py.File(str(filepath)) as sf:
        initial = sf[str(scan+0.1)]
        out = list(initial['instrument']['positioners'].keys())
    return out
# %%

# %% =============================== read ================================== %% #
def read(TOP, sample, dataset, scan, branch=None, detectors_rixs_branch=['dbig_n', 'sam_n', 'mir_rixs'], detectors_xmcd_branch=['ifluo_n', 'it_n', 'i0_n_xmcd'], processed_rixs=True, subscan='auto'):
    """return data from ID32 beamline

    Usage:
        import brixs.beamlines.id32 as id32
        ss = id32.read(TOP, sample, dataset, scan)

    Args:
        TOP (str or path): top/main experiment folderpath where folder such as
            RAW_DATA and PROCESSED_DATA can be found.
        sample (str): sample name.
        dataset (str or None, optional): dataset name. If None, a dictionary
            will be returned with a list of scans for each dataset.
        scan (int): scan number.
        branch (str or None, optional): brach where the data was collected. If
            None, it will be automatically detected. Automatic detection should
            work fine. Default is None. 
        detectors_rixs_branch, detectors_xmcd_branch (list, optional): list of
            detector data to be pulled out from the file. The output will follow
            the same order as the order of the detectors in the list. Detector
            data is pulled out from hdf5 address <scan>/measurement.
        processed_rixs (bool, optional): only for rixs scans. If True, spectrum
            will be pulled out from the processed data folder. If False, this
            function returns the detector image pulled out from the RAW folder.
            Default is True.
        subscan (int or str, optional): if 'auto', it will select the appropriate
            subscan for processed RIXS scan (subscan = 1 or 2). Subscan='auto'
            corrects the bug where first scan in the processed data is saved in 
            the subscan 2. For ascan, trigscan and amesh, subscan is always set 
            to 1. For trigscan, additional metadata is also loaded from subscan
            2. Default is 'auto'.

    Return:
        Returns a list of br.Spectrum objects (br.Spectra) if scan is `ascan` 
        or `xas` (trigscan). Each spectrum carry data from one defined detector
        following the order given in the arguments detectors_rixs_branch and 
        detectors_xmcd_branch.

        Returns a list of br.Image objects if scan is `amesh`.Each image carry 
        data from one defined detector following the order given in the 
        arguments detectors_rixs_branch and detectors_xmcd_branch.
        
        Returns a br.Spectrum object if scan is `rixs` (loopscan) and 
        processed_rixs=True.
        
        Returns a list of br.Image objects if scan is `rixs` (loopscan) and 
        processed_rixs=False. Returns all images from that scan
    """
    # assert branch
    assert branch in ['rixs', 'xmcd', None], f'branch `{branch}` is invalid. Valid options are "rixs", "xmcd", or None. If None, it will be automatically detected'

    # get available scans
    available_scans = list_scans(TOP, sample, dataset)
    assert scan in available_scans, f'scan {scan} cannot be found. Available scans for {sample}/{dataset}: {available_scans}'
    
    # get filepath
    filepath = [_ for _ in br.filelist(TOP/'RAW_DATA'/sample/dataset, string='.h5') if _.is_file()][0]
    
    # open file
    with h5py.File(str(filepath)) as sf:
        # get initial metadata
        initial    = sf[str(scan+0.1)]
        command    = initial['title'][()].decode()
        end_reason = initial['end_reason'][()].decode()
        end_time   = _str2datetime(initial['end_time'][()].decode())
        start_time = _str2datetime(initial['start_time'][()].decode())

        pol = _polarization(initial['instrument']['positioners']['hu70ap'][()])

        # unnecessary metadata
        # print(initial['writer']['date'][()])     # b'2025-04-10T18:45:34.245749+02:00'
        # print(initial['writer']['status'][()])   # b'SUCCEEDED'
        # print(sf[str(scan+0.1)]['sample'].name)  # /16.1/sample
        # print(sf[str(scan+0.1)]['plotselect']['dbig_rixs'][()])     # [-123.53896312 -125.53896312 -118.53896312 -123.53896312 -123.53896312]
        # print(sf[str(scan+0.1)]['plotselect']['elapsed_time'][()])  # [0.         2.15263224 4.32394409 6.49652314 8.67149973]

        # get branch
        if branch is None:
            branch = 'rixs' if 'dbig_n' in initial['measurement'] else 'xmcd'

        # metadata
        if branch == 'rixs':
            detectors = detectors_rixs_branch
        else:
            detectors = detectors_xmcd_branch

        # get data
        if command.startswith('loopscan'):
            if processed_rixs:
                # get rixs file from processed folder
                filepath2 = Path(min([str(_) for _ in br.filelist(TOP/'PROCESSED_DATA'/sample/dataset, string='.spec') if _.is_file()], key=len))
                sf2 = SpecFile(str(filepath2))
                if subscan == 'auto':
                    if str(scan + 0.2) in sf2:
                        initial2 = sf2[str(scan + 0.2)]
                    else:
                        initial2 = sf2[str(scan + 0.1)]
                else:
                    assert br.is_integer(subscan) and subscan > 0, f'subscan must be a integer number higher than zero, not subscan={subscan}'
                    initial2 = sf2[str(scan) + '.' + str(subscan)]

                # spectrum
                s = br.Spectrum(x=initial2.data_column_by_name("Energy (auto)"), 
                                y=initial2.data_column_by_name("SPC"))
                
                # other relevant metadata
                s.pixel   = initial2.data_column_by_name("Pixel")
                s.photons = br.Spectrum(x=initial2.data_column_by_name("Energy (auto)"), 
                                        y=initial2.data_column_by_name("Photons"))
                s.acquisition_time = initial2.data_column_by_name("Acquisition time")[0]
                s.header = [line.split('#C')[-1].strip() for line in initial2.scan_header if line.startswith('#C')]
                s.acquisition_time_per_image = float([line for line in s.header if line.startswith('Acq time')][0].split()[-2])

                # metadata
                s.scan    = scan
                s.sample  = sample
                s.dataset = dataset
                s.branch  = branch

                s.command    = command
                s.start_time = start_time
                s.end_time   = end_time
                s.end_reason = end_reason

                s.pol = pol

                s.count_time = initial['scan_parameters']['count_time'][()]

                s.metadata = {}
                for key in metadata[branch]:
                    try:    
                        s.metadata[key] = initial['instrument']['positioners'][metadata[branch][key]][()]
                    except KeyError:
                        pass
                
                s.additional = {}
                for key in additional[branch]:
                    try:    
                        s.additional[key] = initial['measurement'][additional[branch][key]][()]
                    except KeyError:
                        pass
                    
                return s
            else:
                if branch == 'rixs': # ref indicates how to find the image number via the filename of the image
                    ref = 1
                else: 
                    ref = 2
                out = br.Dummy()
                filepaths = br.parsed_filelist(TOP/'RAW_DATA'/sample/dataset/('scan'+str(scan).zfill(4)), string='.h5', ref=ref, return_type='list')
                for filepath in filepaths:
                    with h5py.File(filepath) as sf3:
                        _im = br.Image(data=sf3[list(sf3.keys())[0]]['measurement/data'][0])
                        if branch == 'rixs':
                            _im.aquisition_time = float(sf3[list(sf3.keys())[0]]['ESRF-ID32/andor1/acquisition/exposure_time'][()])
                        else:
                            _im.aquisition_time = float(sf3[list(sf3.keys())[0]]['ESRF-ID32/dhyana95v2/acquisition/exposure_time'][()])
                        out.append(_im)
        elif command.startswith('amesh'):
            motor1 = command.split(' ')[1]
            motor2 = command.split(' ')[5]
            n1     = int(command.split(' ')[4])
            n2     = int(command.split(' ')[8])
            
            centers = initial['measurement'][motor2][::n1+1]
            
            out = br.Dummy()
            for detector in detectors:
                ss = br.Spectra()
                x = initial['measurement'][motor1][0:n1+1]
                for i in range(n2+1):
                    # x = initial['measurement'][motor1][(n1+1)*i:(n1+1)*(i+1)]
                    y = initial['measurement'][detector][(n1+1)*i:(n1+1)*(i+1)]
                    _s = br.Spectrum(x=x, y=y)
                    ss.append(_s)
                _im = ss.stack_spectra_as_columns()
                _im.x_centers = centers
                _im.detector = detector
                _im.moving_motor = [motor1, motor2]
                out.append(_im)
        elif command.startswith('ascan'):
            motor = command.split(' ')[1]
            x = initial['measurement'][motor][()]

            # create spectra
            out = br.Spectra()
            for detector in detectors:
                _s = br.Spectrum(x=x, y=initial['measurement'][detector][()])
                _s.detector = detector
                _s.moving_motor = motor
                out.append(_s)
        elif command.startswith('trigscan'):
            motor = 'energy_enc'
            x = initial['measurement'][motor][()]

            # create spectra
            out = br.Spectra()
            for detector in detectors:
                _s = br.Spectrum(x=x, y=initial['measurement'][detector][()])
                _s.detector = detector
                _s.moving_motor = motor
                out.append(_s)
        else:
            raise ValueError(f'scan command not valid `{command}`. This function can only read ascan, trigscan, amesh, and loopscan (rixs)')

        # metadata
        _metadata = {}
        for key in metadata[branch]:
            try:    
                _metadata[key] = initial['instrument']['positioners'][metadata[branch][key]][()]
            except KeyError:
                pass

        _additional = {}
        for key in additional[branch]:
            try:    
                _additional[key] = initial['measurement'][additional[branch][key]][()]
            except KeyError:
                if command.startswith('trigscan'):
                    try:
                        _additional[key] = sf[str(scan+0.2)]['measurement'][additional[branch][key]][()]
                    except KeyError:
                        pass

        for _s in out:
            _s.scan    = scan
            _s.sample  = sample
            _s.dataset = dataset
            _s.branch  = branch

            _s.command    = command
            _s.start_time = start_time
            _s.end_time   = end_time
            _s.end_reason = end_reason
            # _s.count_time = count_time

            _s.pol = pol

            _s.metadata = _metadata
            _s.additional = _additional

        return out
# %%

# %% ==================== read directly from filepath ==================== %% #
def readxas(filepath, scan, branch=None, detectors_rixs_branch=['dbig_n', 'sam_n', 'mir_rixs'], detectors_xmcd_branch=['ifluo_n', 'it_n', 'i0_n_xmcd']):
    """return data from ID32 beamline

    Usage:
        import brixs.beamlines.id32 as id32
        ss = id32.read(filepath, scan)

    Args:
        filepath (str or path): filepath is expected to point to
            the .h5 file from a dataset, 
            e.g., TOP/RAW_DATA/<sample>/<dataset>/dataset.h5
        scan (int): scan number.
        branch (str or None, optional): brach where the data was collected. If
            None, it will be automatically detected. Automatic detection should
            work fine. Default is None. 
        detectors_rixs_branch, detectors_xmcd_branch (list, optional): list of
            detector data to be pulled out from the file. The output will follow
            the same order as the order of the detectors in the list. Detector
            data is pulled out from hdf5 address <scan>/measurement.

    Return:
        Returns a list of br.Spectrum objects (br.Spectra) if scan is `ascan` 
        or `xas` (trigscan). Each spectrum carry data from one defined detector
        following the order given in the arguments detectors_rixs_branch and 
        detectors_xmcd_branch.

        Returns a list of br.Image objects if scan is `amesh`.Each iamge carry 
        data from one defined detector following the order given in the 
        arguments detectors_rixs_branch and detectors_xmcd_branch.
    """
    # assert branch
    assert branch in ['rixs', 'xmcd', None], f'branch `{branch}` is invalid. Valid options are "rixs", "xmcd", or None. If None, it will be automatically detected'

    # get available scans
    available_scans = list_scans(filepath=filepath)
    assert scan in available_scans, f'scan {scan} cannot be found. Available scans: {available_scans}'
    
    # get filepath
    filepath = Path(filepath)
    
    # open file
    with h5py.File(str(filepath)) as sf:
        # get initial metadata
        initial    = sf[str(scan+0.1)]
        command    = initial['title'][()].decode()
        end_reason = initial['end_reason'][()].decode()
        end_time   = _str2datetime(initial['end_time'][()].decode())
        start_time = _str2datetime(initial['start_time'][()].decode())

        pol = _polarization(initial['instrument']['positioners']['hu70ap'][()])

        # get branch
        if branch is None:
            branch = 'rixs' if 'dbig_n' in initial['measurement'] else 'xmcd'

        # metadata
        if branch == 'rixs':
            detectors = detectors_rixs_branch
        else:
            detectors = detectors_xmcd_branch

        # get data
        if command.startswith('loopscan'):
            raise ValueError(f'scan {scan} seems to be a rixs scan. This function is for reading xas, ascan, and mesh')
        elif command.startswith('amesh'):
            motor1 = command.split(' ')[1]
            motor2 = command.split(' ')[5]
            n1     = int(command.split(' ')[4])
            n2     = int(command.split(' ')[8])
            
            centers = initial['measurement'][motor2][::n1+1]
            
            out = br.Dummy()
            for detector in detectors:
                ss = br.Spectra()
                x = initial['measurement'][motor1][0:n1+1]
                for i in range(n2+1):
                    # x = initial['measurement'][motor1][(n1+1)*i:(n1+1)*(i+1)]
                    y = initial['measurement'][detector][(n1+1)*i:(n1+1)*(i+1)]
                    _s = br.Spectrum(x=x, y=y)
                    ss.append(_s)
                _im = ss.stack_spectra_as_columns()
                _im.x_centers = centers
                _im.detector = detector
                _im.moving_motor = [motor1, motor2]
                out.append(_im)
        elif command.startswith('ascan'):
            motor = command.split(' ')[1]
            x = initial['measurement'][motor][()]

            # create spectra
            out = br.Spectra()
            for detector in detectors:
                _s = br.Spectrum(x=x, y=initial['measurement'][detector][()])
                _s.detector = detector
                _s.moving_motor = motor
                out.append(_s)
        elif command.startswith('trigscan'):
            motor = 'energy_enc'
            x = initial['measurement'][motor][()]

            # create spectra
            out = br.Spectra()
            for detector in detectors:
                _s = br.Spectrum(x=x, y=initial['measurement'][detector][()])
                _s.detector = detector
                _s.moving_motor = motor
                out.append(_s)
        else:
            raise ValueError(f'scan command not valid `{command}`. This function can only read ascan, trigscan, amesh, and loopscan (rixs)')

        # metadata
        _metadata = {}
        for key in metadata[branch]:
            try:    
                _metadata[key] = initial['instrument']['positioners'][metadata[branch][key]][()]
            except KeyError:
                pass

        _additional = {}
        for key in additional[branch]:
            try:    
                _additional[key] = initial['measurement'][additional[branch][key]][()]
            except KeyError:
                if command.startswith('trigscan'):
                    try:
                        _additional[key] = sf[str(scan+0.2)]['measurement'][additional[branch][key]][()]
                    except KeyError:
                        pass

        for _s in out:
            _s.scan    = scan
            _s.branch  = branch

            _s.command    = command
            _s.start_time = start_time
            _s.end_time   = end_time
            _s.end_reason = end_reason

            _s.pol = pol

            _s.metadata = _metadata

        return out

def readrixs(filepath, scan):
    """return processed spectrum from ID32 beamline

    Usage:
        import brixs.beamlines.id32 as id32
        ss = id32.readrixs(filepath, scan)

    Warning:
        Metadata is not loaded like in id32.read() and id32.readxas()

    Args:
        filepath (str or path): filepath is expected to point to
            the .spec file from a dataset, 
            e.g., TOP/PROCESSED_DATA/<sample>/<dataset>/Online_analysis_..._andor1.SPEC
        scan (int): scan number.

    Return:
        br.Spectrum object
    """
    # get filepath
    filepath2 = Path(filepath)
    
    # open file
    sf2 = SpecFile(str(filepath2))
    initial2 = sf2[str(scan+0.1)]

    # spectrum
    s = br.Spectrum(x=initial2.data_column_by_name("Energy (auto)"), 
                    y=initial2.data_column_by_name("SPC"))
    
    # other relevant metadata
    s.pixel   = initial2.data_column_by_name("Pixel")
    s.photons = br.Spectrum(x=initial2.data_column_by_name("Energy (auto)"), 
                            y=initial2.data_column_by_name("Photons"))
    s.acquisition_time = initial2.data_column_by_name("Acquisition time")[0]
    s.header = [line.split('#C')[-1].strip() for line in initial2.scan_header if line.startswith('#C')]
    s.acquisition_time_per_image = float([line for line in s.header if line.startswith('Acq time')][0].split()[-2])

    # metadata
    s.scan = scan

    return s        

def readraw(filepath):
    """return raw image from ID32 beamline

    Usage:
        import brixs.beamlines.id32 as id32
        ss = id32.readraw(filepath, scan)
    
    Warning:
        Metadata is not loaded like in id32.read() and id32.readxas()

    Args:
        filepath (str or path): filepath is expected to point to
            the .h5 file from the detector, e.g., andor1_0000.h5

    Return:
        br.Image object
    """
    filepath = Path(filepath)
    with h5py.File(filepath) as sf3:
        im = br.Image(data=sf3[list(sf3.keys())[0]]['measurement/data'][0])
    return im
# %%

