#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""The brixs.beamlines.id32 module offers a python solution for data processing
and data analysis for data collected at the ID32 beamline of ESRF.

Besides brixs requiriments, this module also requires h5py and silx.


Onsite, the folderpath for the data is 

/data/visitor/<proposalnumber>/id32/<fulldate>/

for example

/data/visitor/hc6098/id32/20250408/

The folder structure is:

GALLERY
jupyter
NOBACKUP
NOTES
PROCESSED_DATA
RAW_DATA
SCRIPTS

Data is stored in RAW_DATA and PROCESSED_DATA. Inside these two folders,
each dataset will get a folder, i.e., these two folders have the same folder
structure.

In the experiment, 




"""

# %%
cd C:\Users\galdin_c\Documents\work\CrSBr\ares\2025_09_10_ESRF_EDA
ipython
# %%

# %% ========================== standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np


# %% =========================== brixs imports =========================== %% #
import brixs as br
import brixs.beamlines.id32 as id32



# %% ============================== settings ============================= %% #
# matplotlib
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()

# brixs
br.settings.FIGURE_POSITION = (283, 567)
# br.get_window_position()

# %% ============================= folderpaths ============================ %% #
TOP = Path(r'C:\Users\galdin_c\Documents\work\CrSBr\ares\2025_09_10_ESRF_EDA')

OUT    = TOP/'out'   # output
TMP    = TOP/'tmp'   # temporary files (can be deleted)
STORE  = TOP/'store' # storage (not output file, cannot be deleted)


# RAW  = Path(r'/data/visitor/hc6098/id32/20250408')
TOP = Path(r'D:\esrf_CrTe2_data')
# XAS  = TOP/'RAW_DATA'
# RIXS = TOP/'PROCESSED_DATA'
# SAMPLE_ALIGN = RIXS/'align/align_0001/Online_analysis_align_0001_dhyana95v2.spec'
# XAS1 = XAS/'CrTe2/CrTe2_0001/CrTe2_0001.h5'
# %%

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
                     'sample_z': 'sz'}
            }   



# %%
import h5py
import datetime
from silx.io.specfile import SpecFile
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

def list_scans(TOP, sample, dataset=None):
    """"Returns list of available scans for a given sample
    
    Args:
        TOP (str or path): top/main experiment folderpath where folder such as
            RAW_DATA and PROCESSED_DATA can be found.
        sample (str): sample name.
        dataset (str or None, optional): dataset name. If None, a dictionary
            will be returned with a list of scans for each dataset.
    
    Returns:
        list or dict (if dataset=None)
    """   
    TOP = Path(TOP)
    available_datasets = list_datasets(TOP=TOP, sample=sample)

    if dataset is None:
        out = {}
        for dataset in available_datasets:
            out[dataset] = list_scans(TOP=TOP, sample=sample, dataset=dataset)
    else:
        assert dataset in available_datasets, f'dataset "{dataset}" cannot be found. Available options are: {available_datasets}'

        # filepath
        specfile = [_ for _ in br.filelist(TOP/'RAW_DATA'/sample/dataset, string='.h5') if _.is_file()][0]

        # read specfile
        with h5py.File(str(specfile)) as sf:
            out = br.remove_duplicates([int(_.split('.')[0]) for _ in sf.keys()])
        out = np.array(br.sort(out, out))
    return out

def list_available_detectors(TOP, sample, dataset, scan):
    """return a list of available detectors for a scan

    Args:
        TOP (str or path): top/main experiment folderpath where folder such as
            RAW_DATA and PROCESSED_DATA can be found.
        sample (str): sample name.
        dataset (str or None, optional): dataset name. If None, a dictionary
            will be returned with a list of scans for each dataset.
        scan (int): scan number.

    
    Returns:
        list
    """

    initial['measurement'] 
    return

def read(TOP, sample, dataset, scan, branch=None, detectors_rixs_branch=['dbig_n', 'sam_n', 'mir_rixs'], detectors_xmcd_branch=['ifluo_n', 'it_n', 'i0_n_xmcd'], processed_rixs=True):
    """
    """
    # assert branch
    assert branch in ['rixs', 'xmcd', None], f'branch `{branch}` is invalid. Valid options are "rixs", "xmcd", or None. If None, it will be automatically detected'

    # get available scans
    available_scans = list_scans(TOP, sample, dataset)
    assert scan in available_scans, f'scan {scan} cannot be found. Available scans for {sample}/{dataset}: {available_scans}'
    
    # get specfile
    specfile = [_ for _ in br.filelist(TOP/'RAW_DATA'/sample/dataset, string='.h5') if _.is_file()][0]
    
    # open file
    with h5py.File(str(specfile)) as sf:
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

        # RIXS
        if command.startswith('loopscan'):
            # count_time = initial['scan_parameters']['count_time'][()]

            if processed_rixs:
                # get rixs file from processed folder
                specfile2 = Path(min([str(_) for _ in br.filelist(TOP/'PROCESSED_DATA'/sample/dataset, string='.spec') if _.is_file()], key=len))
                sf2 = SpecFile(str(specfile2))
                initial2 = sf2[str(scan+0.1)]

                # spectrum
                s = br.Spectrum(x=initial2.data_column_by_name("Energy (auto)"), 
                                y=initial2.data_column_by_name("SPC"))
                
                # other relevant metadata
                s.pixel   = initial2.data_column_by_name("Pixel")
                s.photons = br.Spectrum(x=initial2.data_column_by_name("Energy (auto)"), 
                                        y=initial2.data_column_by_name("Photons"))
                s.acquisition_time = initial2.data_column_by_name("Acquisition time")[0]
                header = [line.split('#C')[-1].strip() for line in sf2.scan_header(scan+0.1) if line.startswith('#C')]
                s.acquisition_time_per_image = float([line for line in header if line.startswith('Acq time')][0].split()[-2])

                return s
            else:
                ims = br.Dummy()
                filepaths = br.parsed_filelist(TOP/'RAW_DATA'/sample/dataset/('scan'+str(scan).zfill(4)), string='.h5', ref=1, return_type='list')
                for filepath in filepaths:
                    with h5py.File(filepath) as sf3:
                        _im = br.Image(data=sf3[list(sf3.keys())[0]]['measurement/data'][0])
                        _im.aquisition_time = float(sf3[list(sf3.keys())[0]]['ESRF-ID32/andor1/acquisition/exposure_time'][()])
                        ims.append(_im)
                return ims
        elif command.startswith('amesh'):
            motor1 = command.split(' ')[1]
            motor2 = command.split(' ')[5]
            n1     = int(command.split(' ')[4])
            n2     = int(command.split(' ')[8])
            
            centers = initial['measurement'][motor2][::n1+1]
            
            ims = br.Dummy()
            for detector in detectors:
                ss = br.Spectra()
                for i in range(n2+1):
                    x = initial['measurement'][motor1][(n1+1)*i:(n1+1)*(i+1)]
                    y = initial['measurement'][detector][(n1+1)*i:(n1+1)*(i+1)]
                    _s = br.Spectrum(x=x, y=y)
                    ss.append(_s)
                _im = ss.stack_spectra_as_columns()
                _im.x_centers = centers
                _im.detector = detector
                _im.moving_motor = [motor1, motor2]
                ims.append(_im)
            return ims
        elif command.startswith('ascan'):
            motor = command.split(' ')[1]
            x = initial['measurement'][motor][()]
        elif command.startswith('trigscan'):
            motor = 'mtraj'
            x = initial['measurement'][motor][()]
        else:
            raise ValueError(f'scan command not valid `{command}`. This function can only read ascan, trigscan, amesh, and loopscan (rixs)')

        # create spectra
        ss = br.Spectra()
        for detector in detectors:
            _s = br.Spectrum(x=x, y=initial['measurement'][detector][()])
            _s.detector = detector
            _s.moving_motor = motor
            ss.append(_s)

        # metadata
        for _s in ss:
            _s.scan    = scan
            _s.sample  = sample
            _s.dataset = dataset
            _s.branch  = branch

            _s.command    = command
            _s.start_time = start_time
            _s.end_time   = end_time
            _s.end_reason = end_reason
            _s.count_time = count_time

            _s.pol = pol

            _s.metadata = {}
            for key in metadata:
                try:
                    _s.metadata[key] = initial['instrument']['positioners'][metadata[key]][()]
                except KeyError:
                    pass

        return ss
# %%


# %% xas
TOP = Path(r'C:\Users\galdin_c\Documents\work\CrSBr\ares\2025_09_10_ESRF_EDA\id32')

sample  = 'sample1_rixs'
dataset = 'sample1_rixs_0001'
scan    = 18
ss = read(TOP, sample, dataset, scan)

# this sample/dataset was measured at the rixs branch
# therefore, we know the name of the detectors
# by default, read() will data from the following detectors
# detectors_rixs_branch=['dbig_n', 'sam_n', 'mir_rixs']
# 


# %%








# %%
filelist = br.filelist(TOP/'PROCESSED_DATA'/sample/dataset/'workflows')
br.rename_files(filelist, '{}_{}_{}_{}.json', 'sample1_{1}_{2}_{3}.json', ask=True)
y
TOP/'PROCESSED_DATA'/sample/dataset/'workflows'
C:\Users\galdin_c\Documents\work\CrSBr\ares\2025_09_10_ESRF_EDA\id32\PROCESSED_DATA\sample1_rixs\sample1_rixs_0001\workflows


# %% rixs raw image
sample  = 'CrSBr_rixs'
dataset = 'CrSBr_rixs_0001'
scans = [19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30]

ss = br.Spectra()
for scan in scans:
    _s = read(TOP, sample, dataset, scan)
    ss.append(_s)

br.figure()
_ = ss.interp().align().plot()
ss.interp().align().calculate_average().plot(color='black')

scan = 19
ims = read(TOP, sample, dataset, scan, processed_rixs=False)

br.figure()
ims[0].plot()

ss = br.Spectra()
for scan in scans:
    _s = read(TOP, sample, dataset, scan, processed_rixs=False)
    ss.append(_s)

# %%


list_samples(TOP)
sample = 'alignrixs'
list_datasets(TOP, sample)
dataset = 'alignrixs_0002'
list_scans(TOP, sample, dataset)
scan = 3
branch = None

list_samples(TOP)
sample = 'CrTe2'
list_datasets(TOP, sample)
dataset = 'CrTe2_0002'
list_scans(TOP, sample, dataset)
scan = 16
branch = None
ims = read(TOP, sample, dataset, scan)

br.figure()
ims[0].rotate('counterclockwise').plot()
# %%

# %%




# %%
# print(sf[str(scan+0.1)]['measurement'].keys())
# print(sf[str(scan+0.1)]['instrument'].keys())
polv = 


# %% =============================== XAS ================================== %% #
def _readxas(initial):

    
            try:
                if branch == 'XMCD':
                    _s.label = f'#{scan} {_s.type}, pol={_s.pol}, M={round(_s.magnetic_field, 1)} T, th={round(_s.th, 2)}, y={round(_s.sample_y, 4)}, z={round(_s.sample_z, 4)}, {start_time}'
                elif branch == 'RIXS':
                    _s.label = f'#{scan}, pol={_s.pol}, th={round(_s.th, 2)}, x={round(_s.sample_x, 4)}, y={round(_s.sample_y, 4)}, z={round(_s.sample_z, 4)}, {start_time}'
            except: pass

    # get branch
    with h5py.File(filepath, 'r') as f:

        # branch detectors
        if branch == 'XMCD':
            metadata = 


        # metadata
        for _s in (TEY, TFY, I0):
            for attr in metadata:
                try:
                    _s.__setattr__(attr, f[index]['instrument']['positioners'][metadata[attr]][()])
                except KeyError:
                    if verbose: print(f'Cannot load metadata: {attr}')


# 