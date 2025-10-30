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
        
        cosmic (dict or None, optional): dict with arguments for cosmic ray removal.
            Refer to brixs.addons.centroid.find_and_patch() for documentation on
            possible arguments. If None, cosmic rays removal is not applied. 

        centroid (dict or None, optional): dict with parameters for centroid.
            Refer to brixs.addons.centroid.centroid() for documentation on
            possible arguments. If None, Integration mode is used.
            
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