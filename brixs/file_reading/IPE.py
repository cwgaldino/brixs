#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from IPE beamline - Sirius.

Last edited: Felipe Custódio 08-2023
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np
from collections.abc import Iterable
from PIL import Image
import glob

# %% ------------------------- Special Imports ---------------------------- %% #
import brixs as br
import h5py
import matplotlib.pyplot as plt

# %% IPE beamline - SIRIUS - Brazil ============================================
def zero(self:br.Spectrum, peak=None):
    ### Search for the closests peak and shift it to zero
    if peak:
        self.find_peaks()
        c = [self.peaks[i]['c'].value for i in range(self.peaks.number_of_peaks())]
        self.shift= -min(c, key=lambda x: abs(x - peak))
    
def calc_spectrum(self:br.PhotonEvents, n:float):
    return self.calculate_spectrum(nbins = int((max(self.y)-min(self.y))*n))

def set_xlabel(self:br.Spectrum or br.Spectra):
    if self.calib and self.calib != 1:
        self.xlabel = "Energy [eV]"
    else:
        self.xlabel = "Pixel"   

br.PhotonEvents.calc_spectrum = calc_spectrum
br.Spectrum.zero = zero
br.Spectrum.set_xlabel = set_xlabel
br.Spectra.set_xlabel = set_xlabel

def readSPE(folderpath, prefix, ccd=0, x_min=0, x_max=3300):
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
    pe.ylim = [10, 1600]
    
    if ccd==1:
        x_max = 1645
        pe.label='CCD1'
    if ccd==2:
        x_min = 1655
        pe.label='CCD2'
    
    pe.xlim = [x_min, x_max]    
    
    return pe


def read(folderpath, prefix, x_min=0, x_max=3300, n=1, peak=[1000,1000], calib=[1,1]):
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
    pe.ylim = [10, 1600]
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
    
    ### shift reference to zero
    if isinstance(peak, Iterable):
        s1.zero(peak=peak[0])
        s2.zero(peak=peak[1])
    elif isinstance(peak, int):
        s1.zero(peak=peak)
        s2.zero(peak=peak)

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


#

# %%
