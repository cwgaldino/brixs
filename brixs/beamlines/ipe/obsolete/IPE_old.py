<<<<<<< HEAD
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from IPE beamline - Sirius.

Last edited: Felipe Custódio 08-2023
"""

# %% ------------------------- Standard Imports --------------------------- %% #
import numpy as np
from collections.abc import Iterable
import glob
import matplotlib as mpl
import matplotlib.pyplot as plt

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br

# %% ------------------------- Special Imports ---------------------------- %% #
try:
    import h5py
except:
    pass


# %% --------------------- matplotlib Configurations ------------------------ %% #
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['lines.linestyle'] = '-'
mpl.rcParams['lines.markersize'] = 4
mpl.rcParams['lines.marker'] = ''
mpl.rcParams['axes.titlesize'] = 20
mpl.rcParams['axes.labelsize'] = 16
mpl.rcParams['axes.grid'] = True
mpl.rcParams['xtick.labelsize'] = 14
mpl.rcParams['ytick.labelsize'] = 14
mpl.rcParams['grid.color'] = 'b0b0b0'
mpl.rcParams['grid.linestyle'] = '-'
mpl.rcParams['grid.linewidth'] = 0.5
mpl.rcParams['grid.alpha'] = 0.8
mpl.rcParams['legend.fontsize'] = 14

# %% IPE beamline - SIRIUS - Brazil ============================================
def zero(self, peak = 1000, ranges = 10):
    s = self.copy()
    s.crop(peak-ranges/2, peak+ranges/2)
    s.fit_peak(asymmetry=False, moving_average_window=1, ranges=[(peak-ranges/3, peak+ranges/3)])
    self.set_shift(value=-s.peaks[0]['c'].value)
    return s.peaks[0]['c'].value

def calc_spectrum(self:br.PhotonEvents, n:float):
    return self.calculate_spectrum(nbins = int((max(self.y)-min(self.y))*n))

def set_xlabel(self:br.Spectrum):
    if self.calib and self.calib != 1:
        self.xlabel = "Energy [eV]"
    else:
        self.xlabel = "Pixel"   

br.PhotonEvents.calc_spectrum = calc_spectrum
br.Spectrum.zero = zero
br.Spectrum.set_xlabel = set_xlabel
br.Spectra.set_xlabel = set_xlabel

def readSPE(folderpath, prefix, ccd=0, x_min=0, x_max=3300, curvature=False):
    col = 3 if curvature else 4
    x = list()
    y = list()
    if isinstance(prefix, list):
        for p in prefix:
            for name in glob.glob(folderpath+p+'*.h5'):
                file = h5py.File(name, 'r')
                data = file['entry']['data']['data'][:]
                x.extend(data[:,2])
                y.extend(data[:,col])
                file.close()

    elif isinstance(prefix, str):
        for name in glob.glob(folderpath+prefix+'*.h5'):
            file = h5py.File(name, 'r')
            data = file['entry']['data']['data'][:]
            x.extend(data[:,2])
            y.extend(data[:,4])
            file.close()
    
    else:
        print("prefix must be str or list")

    pe = br.PhotonEvents(x, y)
    pe.ylim = [10, 1600]
    
    if ccd==1:
        x_max = 1645
        pe.label='CCD1'
    if ccd==2:
        x_min = 1655
        pe.label='CCD2'
    
    pe.xlim = [x_min, x_max]    
    
    return pe


def read(folderpath, prefix, x_min=0, x_max=3300, n=1, peak=[1000,1000], calib=[1,1], ranges=50):
    x = list()
    y = list()
    if isinstance(prefix, list):
        for p in prefix:
            for name in glob.glob(folderpath+p+'*.h5'):
                file = h5py.File(name, 'r')
                data = file['entry']['data']['data'][:]
                x.extend(data[:,2])
                y.extend(data[:,4])
                file.close()

    elif isinstance(prefix, str):
        for name in glob.glob(folderpath+prefix+'*.h5'):
            file = h5py.File(name, 'r')
            data = file['entry']['data']['data'][:]
            x.extend(data[:,2])
            y.extend(data[:,4])
            file.close()
    
    else:
        print("prefix must be str or list")

    pe = br.PhotonEvents(x, y)
    pe.xlim = [x_min, x_max]
    plt.figure()
    pe.plot()
    plt.show()
    
    pe1= pe.copy()
    pe2 = pe.copy()
    pe1.label='CCD1'
    pe2.label='CCD2'
    pe1.xlim = [x_min, 1645]
    pe2.xlim = [1655, x_max]
    
    ### Photon Events to Zero
    s1 = pe1.calculate_spectrum(nbins = int((max(pe.y)-min(pe.y))*n))
    s2 = pe2.calculate_spectrum(nbins = int((max(pe.y)-min(pe.y))*n))

    s1.xlabel = ''
    s2.xlabel = ''
    
    ### shift reference to zero
    if isinstance(peak, Iterable):
        s1.zero(peak=peak[0], ranges=ranges)
        s2.zero(peak=peak[1], ranges=ranges)
    elif isinstance(peak, int):
        s1.zero(peak=peak, ranges=ranges)
        s2.zero(peak=peak, ranges=ranges)

    ### px to Energy
    if isinstance(calib, Iterable):
        s1.calib = calib[0]
        s2.calib = calib[1]
        
        s1.xlabel = "Energy [eV]" if calib[0] != 1 else "Pixel"
        s2.xlabel = "Energy [eV]" if calib[1] != 1 else "Pixel"
        
    elif isinstance(calib, float):
        s1.calib = calib
        s2.calib = calib
        
        s1.xlabel = "Energy [eV]" if calib != 1 else "Pixel"
        s2.xlabel = "Energy [eV]" if calib != 1 else "Pixel"    
    
    ### shift reference to zero
    elif isinstance(peak, int) or isinstance(peak, Iterable):
        s1.zero(peak=0.0)
        s2.zero(peak=0.0)

    s1.fix_monotonicity()
    s2.fix_monotonicity()
    ss = br.Spectra(s1,s2)
    ss.interp()
    ss.calculate_shift()
    plt.figure()
    ss.plot()
    plt.show()
    
    s = ss.calculate_sum()
    s.label = 'Total Spectrum'
    s.xlabel = s1.xlabel or s2.xlabel or 'Pixel'
    plt.figure()
    s.plot(label=s.label)
    plt.ylabel('Photon counts')
    plt.xlabel(s.xlabel)
    plt.show()
         
    return pe, ss, s


def read2(label, folderpath, prefix, x_min=0, x_max=3300, calib=[1,1], peak=[800,800], ranges=[100,100], n=1, plot=True, xlim=[None,None], zero=True, curvature=False):
    ss = br.Spectra()
    ss1 = br.Spectra()
    ss2 = br.Spectra()
    for i, p in enumerate(glob.glob(folderpath+prefix+'*.h5')):
        file = h5py.File(p, 'r')
        data = file['entry']['data']['data'][:]
        file.close()
        col = 3 if curvature else 4
        pe = br.PhotonEvents(data[:,2], data[:,col])
        pe.xlim = [x_min, x_max]
        
        pe1= pe.copy()
        pe2 = pe.copy()

        pe1.xlim = [x_min, 1650]
        pe2.xlim = [1651, x_max]
        
        ### Photon Events to Zero
        s1 = pe1.calculate_spectrum(nbins = int(1590*n))
        s2 = pe2.calculate_spectrum(nbins = int(1590*n))
        
        s1.label=f'CCD1 {i}'
        s2.label=f'CCD2 {i}'
        
        s1.fix_monotonicity()
        s2.fix_monotonicity()
        
        ss1.append(s1)
        ss2.append(s2)

    ss1.interp()
    ss2.interp()

    ss1.align(ranges=[peak[0]-ranges[0]/2, peak[0]+ranges[0]/2])
    ss2.align(ranges=[peak[1]-ranges[1]/2, peak[1]+ranges[1]/2])
    
    ss1.xlabel = "Energy [eV]" if calib[0] != 1 else "Pixel"
    ss2.xlabel = "Energy [eV]" if calib[1] != 1 else "Pixel"
    
    ss1.ylabel = "Intensity [ph]"
    ss2.ylabel = "Intensity [ph]"

    if plot:
        ss1.sequential_plot()
        ss2.sequential_plot()

    sum_s1 = ss1.calculate_sum()
    sum_s2 = ss2.calculate_sum()
    
    sum_s1.fix_monotonicity()
    sum_s2.fix_monotonicity()

    if zero:
        sum_s1.zero(peak=peak[0], ranges=ranges[0])
        sum_s2.zero(peak=peak[1], ranges=ranges[1])
        
    ### px to Energy
    sum_s1.calib = calib[0]
    sum_s2.calib = calib[1]
    
    ss = br.Spectra(sum_s1, sum_s2)
    ss.fix_monotonicity()
    ss.interp()
    ss.align()
    s = ss.calculate_sum()
    
    s.label = label or 'Total Spectrum'
    s.xlabel = "Energy [eV]" if calib[0] != 1 else "Pixel"
    s.ylabel = "Intensity [ph.eV]"
    if plot:
        plt.figure()
        s.plot(ranges=xlim, label=s.label)
        plt.legend()
        plt.show()
    
    return s, [ss1, ss2]



def curvature(folderpath, prefix, x_min=0, x_max=3300):
    fig, axes = br.subplots(2, 6, sharey='row')
    br.maximize()
    popt = dict()
    for ccd in (1, 2):
        pe = readSPE(folderpath=folderpath, prefix=prefix, ccd=ccd, x_min=x_min, x_max=x_max, curvature=True)
        # binning
        im = pe.binning(nbins=(int(2*1590), 15))

        # Calculate shifts
        im.calculate_roll()

        # fit shifts
        s = br.Spectrum(x=im.x_centers, y=im.calculated_shift)
        popt[ccd], model, R2 = s.polyfit(2)
        print(f'curvature{ccd} = {list(popt[ccd])}')
        
        x = np.linspace(min(pe.x), max(pe.x), 100)
        fit = br.Spectrum(x=x, y=model(x))

        # plot photon events (raw)
        pe.plot(axes[0+(ccd-1)*6], color='black')

        # plot reduced image
        pos = im.plot(axes[2+(ccd-1)*6])

        # plot photon events (raw) with vertical bins
        pe.plot(axes[1+(ccd-1)*6], color='black')
        br.vlines(ax=axes[1+(ccd-1)*6], x=pos.x_edges, color='red', lw=.5)

        # plot horizontal integration of each vertical bin
        cols = im.columns
        cols.switch()
        cols.flip()
        cols.plot(axes[3+(ccd-1)*6])

        # get max and min y (makes plot nicer)
        cols2 = im.columns
        cols2[0].fit_peak()
        ymin = cols2[0].peaks['c'].value - cols2[0].peaks['w'].value*5

        cols2[-1].fit_peak()
        ymax = cols2[-1].peaks['c'].value + cols2[-1].peaks['w'].value*5

        # plot fitting
        offset = cols[0].y[np.argmin(cols[0].x)]
        pe.plot(axes[4+(ccd-1)*6], color='black')
        fit.plot(axes[4+(ccd-1)*6], factor=-1, offset=offset, color='red')

        # set shifts
        pe.plot(axes[5+(ccd-1)*6], color='black')

        pe.set_shift(p=popt[ccd])
        pe.plot(axes[5+(ccd-1)*6], color='red')

        # set y lim
        for i in range(6):
            axes[i+(ccd-1)*6].set_ylim(ymin, ymax)
            axes[i+(ccd-1)*6].set_xlabel(f'x subpixel')

        axes[0+(ccd-1)*6].set_ylabel(f'ccd {ccd}: y subpixel')

    fig.subplots_adjust(top=0.99, bottom=0.05, left=0.05, right=0.99)
    
    return popt

#

def read3(folderpath, prefix, label='', 
          calib1=[1,0], calib2=[1,0], sum_=1, bpp=1,
          peak=[0,0], ranges=[100,100], 
          seq_plot=False, spe_plot=False, spec_plot=True,
          plt_range=[], zero=True, curvature=False, 
          x_min=1, x_max=3300):
    x_all = list()
    y_all = list()
    ss = br.Spectra()
    ss1 = br.Spectra()
    ss2 = br.Spectra()
    E = list()
    TA = list()
    TB = list()
    Ry = list()
    x = list()
    y = list()
    data_ = list()
    files = glob.glob(folderpath+prefix+'*.h5')
    len_ = len(files)
    
    if calib1==[1,0] and calib2==[1,0]:
        zero=False
        peak=[800,800]
        ranges=[200,200]
        

    if sum_ in [0,1,False, 'False','0','1', '']:
        sum_=1

    elif sum_ in [True, 'all']:
        sum_ = len_

    col = 3 if curvature else 4

    for i, p in enumerate(files):
        file = h5py.File(p, 'r')
        data_.append(file['entry']['data']['data'][:])
        NDAttr = file['entry']['instrument']['NDAttributes']
        Ry.append(NDAttr['RIXS_Ry'][:][0])
        E.append(NDAttr['Energy'][:][0])
        TA.append(NDAttr['TempA'][:][0])
        TB.append(NDAttr['TempA'][:][0])
        x_all.extend(data_[i][:,2])
        y_all.extend(data_[i][:,col])
        file.close()
        
    pe_ = br.PhotonEvents()
    if spe_plot:
        pe_.x = x_all
        pe_.y = y_all
        plt.figure()
        pe_.plot()
        plt.show()
    
    E = [np.average(E[i:i+sum_]) for i in range(0, len_, sum_)]
    datas = [data_[i:i+sum_] for i in range(0, len_, sum_)]
    for j, data in enumerate(datas):
        for subdata in data:
            x.extend(subdata[:,2])
            y.extend(subdata[:,col])

        pe = br.PhotonEvents(x, y)
        pe.xlim = [x_min, x_max]
        #pe.ylim = ranges
        
        pe1 = pe.copy()
        pe2 = pe.copy()

        pe1.xlim = [x_min, 1650]
        pe2.xlim = [1651, x_max]
        
        ### Photon Events to Zero:
        s1 = pe1.calculate_spectrum(nbins = int(1590*bpp))
        s2 = pe2.calculate_spectrum(nbins = int(1590*bpp))
        s1.label=f'CCD1 {i}'
        s2.label=f'CCD2 {i}'
        s1.fix_monotonicity()
        s2.fix_monotonicity()
        
        if zero:
            s1.set_shift(value=-(E[j]-calib1[1])/calib1[0])
            s2.set_shift(value=-(E[j]-calib2[1])/calib2[0])

        ss1.append(s1)
        ss2.append(s2)

    ss1.interp()
    ss2.interp()
    
    ss1.align(ranges=[peak[0]-ranges[0]/2, peak[0]+ranges[0]/2])
    ss2.align(ranges=[peak[1]-ranges[1]/2, peak[1]+ranges[1]/2])
    
    ss1.xlabel = "Energy loss[eV]" if calib1[0] != 1 else "Pixel"
    ss2.xlabel = "Energy loss[eV]" if calib2[0] != 1 else "Pixel"
    
    ss1.ylabel = "Intensity [ph]"
    ss2.ylabel = "Intensity [ph]"

    if seq_plot:
        ss1.sequential_plot()
        ss2.sequential_plot()

    sum_s1 = ss1.calculate_sum()
    sum_s2 = ss2.calculate_sum()
    
    sum_s1.fix_monotonicity()
    sum_s2.fix_monotonicity()

    if zero:
        shift1=sum_s1.zero(peak=0, ranges=ranges[0])
        shift2=sum_s2.zero(peak=0, ranges=ranges[1])
        
    ### px to Energy
    sum_s1.calib = -calib1[0]
    sum_s2.calib = -calib2[0]
    
    ss = br.Spectra(sum_s1, sum_s2)
    ss.fix_monotonicity()
    ss.interp()
    ss.align()
    s = ss.calculate_sum()
    
    E = round(np.average(E),2)
    Ry = round(np.average(Ry),2)
    T = round(np.average(TB),2)
    s.label = label or f'E={E}eV T={T}K Ry={Ry}°'
    s.xlabel = "Energy loss[eV]" if calib1[0] != 1 and calib2[0] != 1 else "Pixel"
    s.ylabel = "Photon counts[ph.eV]"
    if spec_plot:
        if not plt_range and zero:
            plt_range=[-50*calib1[0],500*calib1[0]]
        else:
            s.calib = -1
            plt_range=None
 
        plt.figure()
        s.plot(ranges=plt_range, label=s.label)
        plt.xlabel(s.xlabel)
        plt.ylabel(s.ylabel)
        plt.legend()
        plt.show()
    
    return s, [ss1, ss2, pe_]
# %%
=======
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from IPE beamline - Sirius.

Last edited: Carlos Galdino 2024-07-16
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

# %% ========================== rixs metadata ============================= %% #
rixs_attrs = {'ignore': {}, 'raw':{}}

h = rixs_attrs['ignore']
h['modified_date'] = ''
h['scan']          = ''
h['error']         = ''

h = rixs_attrs['raw']
h['Energy']              = 'entry/instrument/NDAttributes/Energy'
h['Energy_SP']           = 'entry/instrument/NDAttributes/Energy_SP'
h['NDArrayEpicsTSSec']   = 'entry/instrument/NDAttributes/NDArrayEpicsTSSec'
h['NDArrayEpicsTSnSec']  = 'entry/instrument/NDAttributes/NDArrayEpicsTSnSec'
h['NDArrayTimeStamp']    = 'entry/instrument/NDAttributes/NDArrayTimeStamp'
h['NDArrayUniqueId']     = 'entry/instrument/NDAttributes/NDArrayUniqueId'
h['PGM_Cff']             = 'entry/instrument/NDAttributes/PGM_Cff'
h['PGM_GR']              = 'entry/instrument/NDAttributes/PGM_GR'
h['PGM_GT']              = 'entry/instrument/NDAttributes/PGM_GT'
h['PGM_MR']              = 'entry/instrument/NDAttributes/PGM_MR'
h['PGM_MT']              = 'entry/instrument/NDAttributes/PGM_MT'
h['RIXSCam_ActualImage'] = 'entry/instrument/NDAttributes/RIXSCam_ActualImage'
h['RIXSCam_NumImages']   = 'entry/instrument/NDAttributes/RIXSCam_NumImages'
h['RIXSCam_exposure']    = 'entry/instrument/NDAttributes/RIXSCam_exposure'
h['RIXS_Ry']             = 'entry/instrument/NDAttributes/RIXS_Ry'
h['RIXS_X']              = 'entry/instrument/NDAttributes/RIXS_X'
h['RIXS_Y']              = 'entry/instrument/NDAttributes/RIXS_Y'
h['RIXS_Z']              = 'entry/instrument/NDAttributes/RIXS_Z'
h['RIXS_GX']             = 'entry/instrument/NDAttributes/RIXS_GX'
h['RIXS_GY']             = 'entry/instrument/NDAttributes/RIXS_GY'
h['RIXS_GZ']             = 'entry/instrument/NDAttributes/RIXS_GZ'
h['RIXS_GRx1']           = 'entry/instrument/NDAttributes/RIXS_GRx1'
h['RIXS_GRz1']           = 'entry/instrument/NDAttributes/RIXS_GRz1'
h['RIXS_DZ']             = 'entry/instrument/NDAttributes/RIXS_DZ'
h['RIXS_DY']             = 'entry/instrument/NDAttributes/RIXS_DY'
h['TempA']               = 'entry/instrument/NDAttributes/TempA'
h['TempB']               = 'entry/instrument/NDAttributes/TempB'
h['Undulator']           = 'entry/instrument/NDAttributes/Undulator'  
                   

# %% =========================== xas metadata ============================= %% #
# xas_attrs = {'ignore': {}, 'raw':{}, 'string':{}, 'round2':{}}                     
# %%

# %% ========================== ascan metadata ============================ %% #
ascan_attrs = {'ignore': {}, 'raw':{}, 'string': {}, 'bool': {}}

h = ascan_attrs['ignore']
h['modified_date'] = ''
h['motors']        = ''
h['error']         = ''
h['snake']         = ''  # mesh
h['motor_x']       = ''  # mesh
h['motor_y']       = ''  # mesh
h['nstep_x']       = ''  # mesh
h['nstep_y']       = ''  # mesh
h['start_x']       = ''  # mesh
h['start_y']       = ''  # mesh
h['stop_x']        = ''  # mesh
h['stop_y']        = ''  # mesh

h = ascan_attrs['string']
h['title']            = 'entry/title'
h['entry_identifier'] = 'entry/entry_identifier'
h['start_time']       = 'entry/start_time'
h['end_time']         = 'entry/end_time'
h['scan_type']        = 'entry/instrument/bluesky/metadata/scan_type'
h['scan']             = 'entry/instrument/bluesky/metadata/scan'

h = ascan_attrs['raw']
h['duration']         = 'entry/duration'
h['num_points']       = 'entry/instrument/bluesky/metadata/num_points'
h['proposal']         = 'entry/instrument/bluesky/metadata/Proposal'

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
        x    = data[:, 2]
        if curv:
            y    = data[:, 4]
        else:
            y    = data[:, 3]
            
        #############################
        # Create PhotonEvent object #
        #############################
        pe1 = br.PhotonEvents(x=x, y=y, xlim=(18, 1650),   ylim=(0, 1608)).crop(18,   1650, None, None)
        pe2 = br.PhotonEvents(x=x, y=y, xlim=(1668, 3300), ylim=(0, 1608)).crop(1668, 3300, None, None)

        #########
        # attrs #
        #########
        metadata = h5.sort_metadata(f=f, attrs_dict=rixs_attrs, verbose=verbose)
        for attr in metadata:
            setattr(pe1, attr, metadata[attr][0])
            setattr(pe2, attr, metadata[attr][0])

        # date
        temp = br.get_modified_date(filepath)
        setattr(pe1, 'modified_date', temp)
        setattr(pe2, 'modified_date', temp)

        # filename
        setattr(pe1, 'filename', filepath.name)
        setattr(pe2, 'filename', filepath.name)

        # scan (this must be romoved when scan number and image number are included as metadata)
        name = str(filepath.name)
        if len(name.split('_')) == 2:
            scan         = name.split('_')[0]
            image_number = name.split('_')[1].split('.')[0]
            setattr(pe1, 'scan', scan)
            setattr(pe2, 'scan', scan)
            setattr(pe1, 'scan', image_number)
            setattr(pe2, 'scan', image_number)
        elif len(name.split('_')) == 3:
            scan         = name.split('_')[1]
            image_number = name.split('_')[2].split('.')[0]
            setattr(pe1, 'scan', scan)
            setattr(pe2, 'scan', scan)
            setattr(pe1, 'scan', image_number)
            setattr(pe2, 'scan', image_number)

        # ccd
        setattr(pe1, 'ccd', 1)
        setattr(pe2, 'ccd', 2)

    return pe1, pe2

def read(fpath, verbose=True, start=0, stop=None, skip=[]):
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
            assert i < start or i > stop, f'skip index ({i}) outside of start ({start}) stop ({start}) image indexes'

        # filter filelist
        filelist = filelist[start:stop + 1]

        # collect data
        x1 = list()
        y1 = list()
        x2 = list()
        y2 = list()
        dummy1 = br.Dummy()
        dummy2 = br.Dummy()
        _pe1 = False
        for j, filepath in enumerate(filelist):
            if j not in skip:
                try:
                    _pe1, _pe2 = _read_rixs(filepath, verbose)
                    x1.extend(_pe1.x)
                    y1.extend(_pe1.y)
                    x2.extend(_pe2.x)
                    y2.extend(_pe2.y)
                    dummy1.append(_pe1)
                    dummy2.append(_pe2)
                except Exception as e:
                    print(f' === ERROR! Image {j} cannot be loaded: {e} ===\n{filepath}')

        # raise error if no image is loaded                    
        assert _pe1 != False, f'no image could be loaded for fpath: {fpath}'      

        # start photon events
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
        pe1.modified_date = dummy1[0].modified_date
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
        pe2.modified_date = dummy2[0].modified_date
        pe2.number_of_images = len(filelist) - len(skip)
        try:
            pe2.scan = int(fpath.name)
        except:
            pass

        return pe1, pe2, dummy1, dummy2

    #########
    # ascan #
    #########
    elif fpath.suffix == '.nxs':
        with h5py.File(Path(fpath), 'r') as f:
            #################
            # get scan type #
            #################
            scan_type = f['entry/instrument/bluesky/metadata/scan_type'][()].decode("utf-8")
            # scan_type = 'mesh'

            ##############
            # get motors #
            ##############
            motors = f['entry/instrument/bluesky/metadata/motors'][()].decode("utf-8").split('\n')[1:-1]
            motors = [m[2:] for m in motors]

            #################
            # get detectors #
            #################
            detectors = f['entry/instrument/bluesky/metadata/detectors'][()].decode("utf-8").split('\n')[:-1]
            detectors = [d[2:] for d in detectors]

            ###################
            # addtional attrs #
            ###################
            metadata_ = {}

            #############
            # read data #
            #############
            x = f['entry/data/' + motors[0]][()]
            if 'RIXS_pd' in detectors:
                PD = br.Spectrum(x=x, y=f['entry/data/RIXS_pd'][()])
            else:
                PD = br.Spectrum()
            if 'RIXS_fy' in detectors:
                TFY = br.Spectrum(x=x, y=f['entry/data/RIXS_fy'][()])
            else:
                TFY = br.Spectrum()
            if 'RIXS_tey in detectors':
                TEY = br.Spectrum(x=x, y=f['entry/data/RIXS_tey'][()])
            else:
                TEY = br.Spectrum()
            if 'RIXS_i0' in detectors:
                I0 = br.Spectrum(x=x, y=f['entry/data/RIXS_i0'][()])
            else:
                I0 = br.Spectrum()        
            _ss  = br.Spectra([TEY, TFY, I0, PD])
            for s in [_ for _ in _ss] + [_ss]:
                    s.EPOCH = f['entry/data/EPOCH'][()]
                    for motor in motors:
                        s.__setattr__('SETPOINT' + '_' + motor, f['entry/data/' + (motor+'_user_setpoint')][()])
        #########
        # ascan #
        #########
        if scan_type == 'ascan':
            assert len(motors) == 1, f'number of motors ({len(motors)}) is not compatible with scan type (ascan)'
            ss = _ss

        ##########
        # a2scan #
        ##########
        elif scan_type == 'a2scan':
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
            ##############
            # mesh attrs #
            ##############
            with h5py.File(Path(fpath), 'r') as f:
                snake = f['entry/instrument/bluesky/metadata/snaking'][()].decode("utf-8").split('\n')[1:-1]
                metadata_['snake']   = [m[2:]=='true' for m in snake]
                metadata_['motor_x'] = f['entry/instrument/bluesky/metadata/motor_x'][()].decode("utf-8")
                metadata_['motor_y'] = f['entry/instrument/bluesky/metadata/motor_y'][()].decode("utf-8")
                metadata_['nstep_x'] = int(f['entry/instrument/bluesky/metadata/nstep_x'][()])
                metadata_['nstep_y'] = int(f['entry/instrument/bluesky/metadata/nstep_y'][()])
                metadata_['start_x'] = f['entry/instrument/bluesky/metadata/start_x'][()]
                metadata_['start_y'] = f['entry/instrument/bluesky/metadata/start_y'][()]
                metadata_['stop_x']  = f['entry/instrument/bluesky/metadata/stop_x'][()]
                metadata_['stop_y']  = f['entry/instrument/bluesky/metadata/stop_y'][()]

            # find `x` motor points (frozen motor --> always motor_y)
            bfinal = np.linspace(metadata_['start_y'], metadata_['stop_y'], metadata_['nstep_y'])
                
            # find `y` motor points
            afinal = np.linspace(metadata_['start_x'], metadata_['stop_x'], metadata_['nstep_x'])

            # reshape data
            ss = br.Dummy([TEY, TFY, I0, PD])
            for j, s in enumerate(_ss):
                # get y
                y = copy.deepcopy(s.y)

                # check if scan finished
                if len(y) < (len(afinal) * len(bfinal)):
                    # y = np.array(list(y) + [y[-1]]*((len(afinal) * len(bfinal))-len(y)))
                    y = np.array(list(y) + [None]*((len(afinal) * len(bfinal))-len(y)))
                    metadata_['error']  = 'mesh interrupted early'

                # reashape (check snake)
                if metadata_['snake'][motors.index(metadata_['motor_x'])] == True:
                    ss[j] = br.Image(data=[row if i%2 == 0 else row[::-1] for i, row in enumerate(y.reshape(metadata_['nstep_y'], metadata_['nstep_x']))])
                else:
                    ss[j] = br.Image(data=[row if i%2 == 0 else row[::] for i, row in enumerate(y.reshape(metadata_['nstep_y'], metadata_['nstep_x']))])
                ss[j].x_centers = afinal
                ss[j].y_centers = bfinal

        #########
        # attrs #
        #########
        with h5py.File(Path(fpath), 'r') as f:
            metadata = h5.sort_metadata(f=f, attrs_dict=ascan_attrs, verbose=verbose)
        metadata.update(metadata_)

        # fix metadata
        metadata['motors'] = motors
        try:
            metadata['modified_date'] = br.get_modified_date(fpath)
            metadata['start_time']    = _str2datetime(metadata['start_time'])
            metadata['start_time']    = _str2datetime(metadata['start_time'])
        except:
            pass

        # set metadata
        for s in [_ for _ in ss] + [ss]:
            for attr in metadata:
                setattr(s, attr, metadata[attr])
        
        return ss

    #######
    # XAS #
    #######
    else:
        ################
        # get metadata #
        ################
        comments = br.load_comments(filepath=fpath, comment_flag='#', stop_flag='#')
        comments = {line.split(':')[0][2:]:line.split(':')[1][1:-1] for line in comments if ':' in line}
        assert 'scan_type' in comments, 'scan type not found. File corrupted'

        #############
        # XAS (old) #
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
            ss  = br.Spectra([TEY, TFY, I0, PD])
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
        
        #######
        # XAS #
        #######
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
            ss  = br.Spectra([TEY, TFY, I0, PD])
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
   
    ############
    # spectrum #
    ############
    s1  = pe1.integrated_rows_vs_y_centers(nrows=sbins)
    s2  = pe2.integrated_rows_vs_y_centers(nrows=sbins)
    ss1 = br.Spectra([s1, s2])
    ss2 = ss1.interp().align()
    s   = ss2.calculate_sum()
    
    #########
    # attrs #
    #########
    s.copy_attrs_from(pe1)
    del s.ccd
    
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
    
    #########
    # calib #
    #########
    if isinstance(calib, Iterable):
        s.calib = calib[0]
        s.shift = -(s.Energy - calib[1])
    elif calib:
        s.calib = calib
    #s.shift = -800*calib

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

# %% =========================== alignment plot =============================== %% #
def alignment(folderpath, scans, atype, sbins=2000, calib=None, norm=False, start=0, stop=None, skip=[]):
    """return rixs spectrum
    """
    
    assert atype in ['r1', 'r2', 'z', 'motor'], "atype must be 'r1', 'r2', 'z' or 'motor'"
    
    basepath = Path(folderpath)
    ss = br.Spectra()
    si = br.Spectra([br.Spectrum()]*2)
    fiti  = br.Spectra([br.Spectrum()]*2)
    popti = [[],[]]
    rixs_z  = list()
    rixs_gz = list()
    rixs_dz = list()
    rixs_dy = list()
    popt = list()
    for scan in scans:
        z  = None
        gz = None
        dz = None
        dy = None
        
        folderpath = basepath/str(scan).zfill(4)
        assert folderpath.exists(), f"Folderpath does not exist, {folderpath}"

        tmp = _process(folderpath, sbins=sbins, calib=calib, norm=norm, start=start, stop=stop, skip=skip)
        s = tmp['s']
        
        
        for i in range(2):
            si[i] = tmp['ss2'][i]
            fiti[i], popt, R2, model = si[i].fit_peak(fixed_m=0, asymmetry=False, limits=[500, 1300])        
            popti[i].append(popt)
        
        if hasattr(s, 'RIXS_Z'):
            z = round(s.RIXS_Z,3)
            rixs_z.append(z)
        if hasattr(s, 'RIXS_GZ'):
            gz = round(s.RIXS_GZ,3)
            rixs_gz.append(gz)
        if hasattr(s, 'RIXS_DZ'):
            dz = round(s.RIXS_DZ,3)
            rixs_dz.append(dz)
        if hasattr(s, 'RIXS_DY'):
            dy = round(s.RIXS_Z,3)
            rixs_dy.append(round(s.RIXS_DY, 3))
            
        if atype == 'z':
            s.label = f'RIXS_Z = {z}'
        elif atype == 'r1':
            s.label = f'RIXS_GZ = {gz}'
        elif atype == 'r2':
            s.label = f'RIXS_DZ = {dz}'
            
        ss.append(s)

    
    if atype == 'z':
        x1 = rixs_z
        x2 = rixs_z
        x1_label = 'RIXS_Z'
        x2_label = 'RIXS_Z'
    elif atype == 'r1':
        x1 = rixs_gz
        x2 = rixs_dz
        x1_label = 'RIXS_GZ'
        x2_label = 'RIXS_DZ'
    elif atype == 'r2':
        x1 = rixs_dz
        x2 = rixs_dy
        x1_label = 'RIXS_DZ'
        x2_label = 'RIXS_DY'
    # elif atype == 'motor':
    #     x1 = rixs_z
    #     x2 = rixs_z
    #     x1_label = 'RIXS_Z'
    #     x2_label = 'RIXS_Z'
    
    # Ordenar x, y e z juntos
    _sorted = sorted(zip(x1, x2, popti[0], popti[1]), key=lambda item: item[0])

    # Separar x, y e z_sorted
    x1_sorted, x2_sorted, popt0_sorted, popt1_sorted = zip(*_sorted)

    # Separar popt_sorted em 4 colunas
    popt0_new = [list(column) for column in zip(*popt0_sorted)]
    popt1_new = [list(column) for column in zip(*popt1_sorted)]

    fig, axs = plt.subplots(3, 2, figsize=(10, 12))
    ax1 = plt.subplot2grid((3, 2), (0, 0), colspan=2)
    ax2 = plt.subplot2grid((3, 2), (1, 0))
    ax3 = plt.subplot2grid((3, 2), (1, 1))
    ax4 = plt.subplot2grid((3, 2), (2, 0))
    ax5 = plt.subplot2grid((3, 2), (2, 1))
    
    fig.suptitle(f'{atype.upper()} - Aligment')
    ss.plot(ax=ax1)
    ax1.legend()
    ax1.set_xlim(min(popt0_new[1])*0.96, max(popt0_new[1])*1.04)
    ax1.set_ylabel('counts/bin')
    if calib is None:
        ax1.set_xlabel('x (pixel)')
    else:
        ax1.set_xlabel('Energy (eV)')

    ax2.set_title('fwhm')
    ax2.plot(x1_sorted, popt0_new[2], 'o-', label='ccd 1')
    ax2.plot(x1_sorted, popt1_new[2], 'o-', label='ccd 2')
    ax2.set_xlabel(x1_label)
    ax22 = ax2.twiny()  # Criando um eixo secundário
    ax22.plot(x2_sorted, popt0_new[2], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
    ax22.plot(x2_sorted, popt1_new[2], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
    ax22.set_xlabel(x2_label)
    ax22.xaxis.set_label_coords(0.5, 0.90)
    ax2.legend()
    if calib is None:
        ax2.set_ylabel('y (pixel)')
    else:
        ax2.set_ylabel('Energy (eV)')
    
    ax3.set_title('center')
    ax3.plot(x1_sorted, popt0_new[1], 'o-', label='ccd 1')
    ax3.plot(x1_sorted, popt1_new[1], 'o-', label='ccd 2')
    ax3.set_xlabel(x1_label)
    ax32 = ax3.twiny()  # Criando um eixo secundário
    ax32.plot(x2_sorted, popt0_new[1], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
    ax32.plot(x2_sorted, popt1_new[1], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
    ax32.set_xlabel(x2_label)
    ax32.xaxis.set_label_coords(0.5, 0.90)
    ax3.set_ylabel('y (pixel)')
    ax3.legend()
    
    ax4.set_title('amplitude')
    ax4.plot(x1_sorted, popt0_new[0], 'o-', label='ccd 1')
    ax4.plot(x1_sorted, popt1_new[0], 'o-', label='ccd 2')
    ax4.set_xlabel(x1_label)
    ax42 = ax4.twiny()  # Criando um eixo secundário
    ax42.plot(x2_sorted, popt0_new[0], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
    ax42.plot(x2_sorted, popt1_new[0], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
    ax42.set_xlabel(x2_label)
    ax42.xaxis.set_label_coords(0.5, 0.90)
    ax4.set_ylabel('counts/bin')
    ax4.legend()
    
    ax5.set_title('offset')
    ax5.plot(x1_sorted, popt0_new[3], 'o-', label='ccd 1')
    ax5.plot(x1_sorted, popt1_new[3], 'o-', label='ccd 2')
    ax5.set_xlabel(x1_label)
    ax52 = ax5.twiny()  # Criando um eixo secundário
    ax52.plot(x2_sorted, popt0_new[3], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
    ax52.plot(x2_sorted, popt1_new[3], 'o-', color='orange', alpha=0)  # Plotando invisível para criar o eixo
    ax52.set_xlabel(x2_label)
    ax52.xaxis.set_label_coords(0.5, 0.90)
    ax5.set_ylabel('y (pixel)')
    ax5.legend()
    
    plt.subplots_adjust(hspace=0.6, wspace=0.4)
    
    return [x1_sorted, x2_sorted], [popt0_new, popt1_new]

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
    fig, axes = br.subplots(1, 6, sharey='row', figsize=figsize)
    # fig.subplots_adjust(top=0.99, bottom=0.05, left=0.05, right=0.99)

    #############
    # read file #
    #############
    pe1, pe2, pes1, pes2 = read(fpath=folderpath)
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
    # print(f'ccd{ccd+1} curvature= {list(popt[ccd])}')
        
    ########
    # plot #
    ########

    # plot photon events (raw)
    pe.plot(axes[0], color='black')

    # plot reduced image
    pos = im.plot(axes[2])

    # plot photon events (raw) with vertical bins
    pe.plot(axes[1], color='black')
    br.axvlines(ax=axes[1], x=pos.x_edges, colors='red', lw=.5)

    # plot horizontal integration of each vertical bin
    cols = im.columns.switch_xy()#.flip_y()
    cols.plot(axes[3])

    # get max and min y (makes plot nicer)
    cols2 = im.columns
    arr100, _popt, err, model = cols2[0].fit_peak()
    ymin = _popt[1] - _popt[2]*12

    arr100, _popt, err, model = cols2[-1].fit_peak()
    ymax = _popt[1] + _popt[2]*12

    # plot fitting
    offset = cols[0].y[np.argmax(cols[0].x)]
    pe.plot(axes[4], color='black')
    fit.plot(axes[4], factor=-1, offset=offset, color='red')

    # set shifts
    pe.plot(axes[5], color='black')

    # fix curvature
    pe = pe.set_vertical_shift_via_polyval(p=popt)
    pe.plot(axes[5], color='red')

    # set y lim
    for i in range(6):
        axes[i].set_ylim(ymin, ymax)
        axes[i].set_xlabel(f'x subpixel')

    axes[0].set_ylabel(f'ccd {ccd}: y subpixel')

    return popt




# %%


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
>>>>>>> main
