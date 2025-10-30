#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""The centroid algorithm returns photon hit positions with sub pixel resolution

This is an approximate description of the algorithm:
    1) enhance image (floor, square, moving average)
    2) get all positions where pixel intensity is above threshold
    3) for positions that are too close (<=n), keep only the one for the brightest pixel
    4) go back to the original image and find brightest pixel around the positions
    suggested by the enhanced image (<=n) 
    5) if two positions are too close (<=_cm_n), assign them as double events
    6) if _cm_bkg, floor the image. See _cm_bkg argument below.
    7) calculate center of mass of position (<=_cm_n)
    8) if the pixel intensity of a position is higher than threshold2, assign
    it as double event (only if threshold2 is not None)

The main function of this implementation is the im.centroid(). See below the arguments of
the centroid method,

    Args:
        nx, ny (int): photon hits candidates that are within n pixels of distance 
            from each other will be considered the same candidate. For better 
            results, set n to be roughly the expected pixel distance between a 
            photon hit and the farthest excited pixel.
        threshold (number): threshold value. Any pixel value above threshold 
            will be considered as a photon hit candidate. If _bkg, _square, and _n, 
            then threshold must be given in terms of pixel intensity of the 
            enhanced image [use im.enhance(n=_n, bkg=_bkg) to get enhanced image]. 
            See enhanced arguments below.
        threshold2 (number, optional): upper limit for threshold. If 
            threshold2 is not None, pixel values
            above threshold2 will be considered double events. If _bkg, _square, and _n,  
            threshold2 must be given in terms of pixel intensity of the 
            enhanced image [use im.enhance(n=_n, bkg=_bkg) to get enhanced image]. 
            Note that, even if threshold2 is None, two photon hit candidates 
            that are closer than _cm_n pixels are still going to be considered 
            double events. Default is None.

        Image enhancement for finding candidates:
        
        _bkg (number, 'auto', or None, optional): the image will be subtracted by _bkg.
            If _bkg='auto', bkg will be defined so the average of the whole 
            image is zero. If None, no offsetting is applied. Default is None. 
            Note that _bkg = 0 implies no flooring (same as _bkg=None).
        _square (bool, optional): If True, the image will be squared. Default 
            is False.
        _nx, _ny (int): Use this to overwrite the size of the moving average window 
            used in the process of enhancing the image, i.e., number of points
            to average. If this is None, the averaging window size will be set 
            to n*2+1. Note that _n=1 implies no moving average. Default is 1.

        Center of mass calculation:

        _cm_bkg (number, 'auto', or None, optional): If _cm_bkg is not None, the image will be 
            subtracted by _cm_bkg before calculating the center of masses.
            If _cm_bkg='auto', _cm_bkg will be defined so the average of the whole 
            image is zero. If _cm_bkg is None, no offsetting of the image is 
            applied. Default is None. Note that _cm_bkg = 0 implies no 
            flooring before calculating center of mass (same as _cm_bkg=None). 
            Note that, center of mass calculation can yield less precise 
            results if image is not floored and _cm_bkg >> n.
        _cm_nx, _cm_ny (int): Use this to overwrite the number of neighbors when
            calculating the center of mass of a photon hit candidate, e.g., if 
            _cm_n=1, only first neighbors. If 'auto', 
            _cm_n will be same as n. Default is 'auto'.     
        _cm_spot_zero_type (str): when calculating the center of mass of a 
            candidate, pixels around the candidate cannot be negative, otherwise
            center of mass calculation can yield to less precise result (the
            result of the center of mass calculation can even be a detector
            position "outside" the range of the detector). Therefore, if a 
            negative pixel intensity is present around a candidate there is 
            two thing we can do 1) set it to zero (which makes sense because
            if it is negative, we assume it is close to zero), or 2) we apply a 
            small intensity offset to the whole spot around the candidate to 
            make all pixels around it positive or zero. Use
            spot_zeroing_type='zero' for 1 and 'offset' for 2. Default is 'zero'
        
        other args:

        MAX_NUMBER_OF_CANDIDATES (int, optional): raises error if number of 
            pixels with intensity higher than the threshold is larger than 
            MAX_NUMBER_OF_CANDIDATES in the enhanced image. Note that this will take into account
            every pixel with intensity above the threshold (without excluding
            pixels around (nx, ny) brightest pixels. Default is 1000.
            Useful for preventing setting up the threshold too low.
    
im.centroid returns a list of photon hits and double events

    Returns:
        pe, pe2
        photon events list, and photon events list for double events
"""

# %% ========================== standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np

# %% =========================== brixs imports =========================== %% #
import brixs as br
import brixs.addons.centroid

# %% ============================== settings ============================= %% #
# matplotlib (optional)
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()

# brixs (optional)
# br.settings.FIGURE_POSITION = (283, 567)
# br.get_window_position()
# %%


# %  ===================================================================== %% #
# %  ======================== ID32 beamline example ====================== %% #
# %% ===================================================================== %% #
# please, refer to brixs/examples/beamlines/id32/centroid.py
# %%

# %  ===================================================================== %% #
# %  ======================== IPE beamline example ====================== %% #
# %% ===================================================================== %% #
# please, refer to brixs/examples/beamlines/ipe/centroid.py
# %%
