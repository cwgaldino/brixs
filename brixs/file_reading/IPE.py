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

# %% ------------------------- Special Imports ---------------------------- %% #
import brixs as br
import h5py

# %% --------------------- matplotlib Configurations ------------------------ %% #
mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['lines.linestyle'] = '-'
mpl.rcParams['lines.markersize'] = 4
mpl.rcParams['lines.marker'] = 'o'
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
        im = pe.binning(nbins=(int(3*1590), 15))

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

def read3(label, folderpath, prefix, x_min=0, x_max=3300, calib1=[1,0], calib2=[1,0],
          peak=[800,800], ranges=[100,100], n=1, 
          plot=True, xlim=[1,3300], zero=True, curvature=False,
          sum_='', plot_spe=False):
    x_all = list()
    y_all = list()
    ss = br.Spectra()
    ss1 = br.Spectra()
    ss2 = br.Spectra()
    E = list()
    TA = list()
    TB = list()
    x = list()
    y = list()
    data_ = list()
    files = glob.glob(folderpath+prefix+'*.h5')
    len_ = len(files)

    if sum_ in [0,1,False, 'False','0','1']:
        sum_=1

    elif sum_ in [True, 'all']:
        sum_ = len_

    col = 3 if curvature else 4

    for i, p in enumerate(files):
        file = h5py.File(p, 'r')
        data_.append(file['entry']['data']['data'][:])
        NDAttr = file['entry']['instrument']['NDAttributes'][:]
        E.append(NDAttr['Energy'])
        TA.append(NDAttr['TempA'])
        TB.append(NDAttr['TempA'])
        x_all.extend(data_[i][:,2])
        y_all.extend(data_[i][:,col])
        file.close()
    
    if plot_spe:
        pe_ = br.PhotonEvents(x_all, y_all)
        plt.figure()
        pe_.plot()
        plt.show()

    datas = [data_[i:i+sum_] for i in range(0, len_, sum_)]
    for data in datas:
        for subdata in data:
            x.extend(subdata[:,2])
            y.extend(subdata[:,col])

        pe = br.PhotonEvents(x, y)
        pe.xlim = [x_min, x_max]
        
        pe1 = pe.copy()
        pe2 = pe.copy()

        pe1.xlim = [x_min, 1650]
        pe2.xlim = [1651, x_max]
        
        ### Photon Events to Zero:
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
    
    ss1.xlabel = "Energy [eV]" if calib1[0] != 1 else "Pixel"
    ss2.xlabel = "Energy [eV]" if calib2[0] != 1 else "Pixel"
    
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
    sum_s1.calib = -calib1[0]
    sum_s2.calib = -calib2[0]
    
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

# %%
