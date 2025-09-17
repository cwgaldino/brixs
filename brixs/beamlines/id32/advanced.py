#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced functions for ID32 beamline at ESRF"""

# %% ------------------------- Standard Imports -------------------------- %% #
from pathlib import Path
import numpy as np

# %% ------------------------------ brixs -------------------------------- %% #
import brixs as br
from brixs.beamlines.id32.core import read

# %%

# %% ============================== rixs ================================= %% #
def process(TOP, sample, dataset, scan, nbins=2000, curv=None, calib=None, cosmic=None, centroid=None, include_double_events=True, return_photon_events=False):
    """return rixs spectrum from processed image (centroid and integration mode)

    Warning:
        This function does not do any kind of normalization (exposure time, I0, ...)

    Args:
        TOP (str or path): top/main experiment folderpath where folder such as
            RAW_DATA and PROCESSED_DATA can be found.
        sample (str): sample name.
        dataset (str or None, optional): dataset name. If None, a dictionary
            will be returned with a list of scans for each dataset.
        scan (int): scan number.

        nbins (int, optional): number of bins for converting photon events to 
            spectrum (number of points in the spectrum). Only used if centroid
            is not None. Default is 2000.
        curv (None or list): if not None, curv must be 1D array of polynomial 
            coefficients (including coefficients equal to zero) from 
            highest degree to the constant term 
            [f(x_centers) = curv[n]*x**n + curv[n-1]*x**n-1 + ... + curv[0]]
        calib (number, optional): if not None, the x axis of the final spectrum
            is multiplied by calib. This number must have unit of eV/pixel.

        include_double_events (bool, optional): if True, double events are 
            included for calculating the final spectrum. Only for centroid mode.
            Default is True.
        return_photon_events (bool, optional): if True, returns photon events 
            as a metadata of the spectrum. Only for centroid mode. Default is
            False.
        
        cosmic (dict or None, optional): dict with parameters for cosmic ray removal.
            Valid parameters are:
                n (int): photon hits candidates that are within n pixels of distance 
                    from each other will be considered the same candidate. For better 
                    results, set n to be roughly the expected pixel distance between a 
                    photon hit and the farthest excited pixel.
                threshold (number): threshold value. Any pixel value above threshold 
                    will be considered as a photon hit candidate. If enhance=True, 
                    threshold must be given in terms of pixel intensity of the 
                    enhanced image [use im.enhance(n=_n, bkg=_bkg) to get enhanced image]. 
                threshold2 (number, optional): upper limit for threshold. If 
                    threshold2 is not None, pixel values
                    above threshold2 will be disregarded. If enhance=True, 
                    threshold2 must be given in terms of pixel intensity of the 
                    enhanced image [use im.enhance(n=_n, bkg=_bkg) to get enhanced image]. 
                    Default is None.
                
                Image enhancement args:

                _bkg (number, optional): Use this to overwrite the bkg value for flooring
                    the image. If _bkg is not None, the image will be subtracted by _bkg.
                    If _bkg is None, bkg will be defined so the average of the whole 
                    image is zero. Default is None. Note that _bkg = 0 implies no 
                    flooring.
                _square (bool, optional): If True, the image will be floored (an offset
                    will be applied so avg pixel intensity is zero), squared, 
                    and a moving averaged of size n will be applied. Default is True.
                _n (int): Use this to overwrite the size of the moving average window 
                    used in the process of enhancing the image, i.e., number of points
                    to average. If this is None, the averaging window size will be set 
                    to n*2+1. Note that _n=1 implies no moving average.

                patch args:

                _patch_size (int, optional): Use this to overwrite the patch size (how
                    many pixel around a photon hit will be patched out of the image). 
                    If None, this will be set to be the same size as n. Default is None.
                _patch_value (number, optional): Use this to overwrite the pixel 
                    value of the patches. If None, patches will be dynamically calculated
                    to the average of the surrounding pixels around the patch. Default is
                    None.

                other args:

                MAX_NUMBER_OF_CANDIDATES (int, optional): raises error if number of 
                    photons to be patched out is larger than MAX_NUMBER_OF_CANDIDATES.
                    Useful for preventing too low threshold as patching image is slow.
                    Default is 10. 
            If None, cosmic rays removal is not applied. 

        centroid (dict or None, optional): dict with parameters for centroid.
            Valid parameters are:
                n (int): photon hits candidates that are within n pixels of distance 
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
                
                _bkg (number, optional): Use this to overwrite the bkg value for flooring
                    the image. If _bkg is not None, the image will be subtracted by _bkg.
                    If _bkg is None, bkg will be defined so the average of the whole 
                    image is zero. Default is None. Note that _bkg = 0 implies no 
                    flooring.
                _square (bool, optional): If True, the image will be floored (an offset
                    will be applied so avg pixel intensity is zero), squared, 
                    and a moving averaged of size n will be applied. Default is True.
                _n (int): Use this to overwrite the size of the moving average window 
                    used in the process of enhancing the image, i.e., number of points
                    to average. If this is None, the averaging window size will be set 
                    to n*2+1. Note that _n=1 implies no moving average.

                Center of mass calculation:

                _cm_bkg (bool, optional): If _cm_bkg is not None, the image will be 
                    subtracted by _bkg before calculating the center of masses.
                    If _cm_bkg is None, _cm_bkg will be defined so the average of the whole 
                    image is zero. Default is None. Note that _cm_bkg = 0 implies no 
                    flooring before calculating center of mass. This 
                    is unnecessary if the image is already originally floored. Center of
                    mass calculation can yield less precise results if image is not floored and
                    _cm_bkg >> n.
                _cm_n (int): Use this to overwrite the number of neighbors when
                    calculating the center of mass of a photon hit candidate, e.g., if 
                    _cm_n=1, only first neighbors. _cm_n also defines how close two 
                    candidates need to be to be considered a double event. If None, 
                    _cm_n will be same as n. Default is None.     
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
                    photons to be patched out is larger than MAX_NUMBER_OF_CANDIDATES.
                    Useful for preventing too low threshold as patching image is slow.
                    Default is 10. 

            If None, Integration mode is used.
            
    Return:
        Spectrum
    """
    ims = read(TOP, sample, dataset, scan, processed_rixs=False)

    # remove cosmic rays, centroid, and adding photon events
    if centroid is not None:
        pe  = br.PhotonEvents()
        pe2 = br.PhotonEvents()
        counts = []
        for _im in ims:
            if cosmic is not None:
                _im2, _pec = _im.find_and_patch(**cosmic)
            else:
                _im2 = _im

            # centroid
            _pe, _pe2  = _im2.floor().centroid(**centroid)
            if include_double_events:
                pe += _pe + _pe2
            else:
                pe += _pe
            pe2 += _pe2
            counts.append(len(_pe))

        # curvature correction
        if curv is not None:
            pe  = pe.set_vertical_shift_via_polyval(curv)
            pe2 = pe2.set_vertical_shift_via_polyval(curv)

        # binning
        s  = pe.binning(1, nbins).integrated_rows_vs_y_centers()
        s2 = pe2.binning(1, nbins).integrated_rows_vs_y_centers()

        # calibrate
        if calib is not None:
            s  = s.set_calib(calib)
            s2 = s2.set_calib(calib)
        
        # metadata
        s.number_of_counts = len(pe)
        s.number_of_double = len(pe2)
        s.number_of_counts_per_image = counts
        if return_photon_events:
            s.pe  = pe
            s.pe2 = pe2
            s.s2  = s2
    else:
        if cosmic is not None:
            im, _pec = ims[0].find_and_patch(**cosmic)
            for _im in ims[1:]:
                _im2, _pec = _im.find_and_patch(**cosmic)
                im += _im2
        else:
            im = ims.calculate_sum()

        # curvature correction
        if curv is not None:
            im = im.set_vertical_shift_via_polyval(curv)
        
        # integration
        s = im.integrated_rows_vs_y_centers()

        # calibrate
        if calib is not None:
            s  = s.set_calib(calib)

    # metadata
    s.scan    = scan
    s.sample  = sample
    s.dataset = dataset
    
    try:
        s.command    = ims[0].command
        s.start_time = ims[0].start_time
        s.end_time   = ims[-1].end_time
        s.end_reason = ims[-1].end_reason

        s.pol = ims[0].pol
    except:
        pass

    try:
        s.number_of_images          = len(ims)
        s.aquisition_time_per_image = ims[0].aquisition_time
        s.aquisition_time           = len(ims)*ims[0].aquisition_time
    except:
        pass

    try:
        s.metadata = ims[0].metadata
    except:
        pass

    return s
# %%