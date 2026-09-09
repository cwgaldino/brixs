#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced functions for I21 beamline

Todo:
    * implement 'auto' for darkfactor in _process()
"""

# %% ========================== Standard Imports ========================= %% #
from collections.abc import Iterable
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import matplotlib

import matplotlib.lines as lines
import matplotlib.transforms as mtransforms

# %% ========================== Special Imports ========================== %% #
try:
    from scipy.ndimage import median_filter, gaussian_filter
except:
    pass
# %% =============================== brixs =============================== %% #
import brixs as br
from brixs.beamlines.i21.core import get_metadata, readrixs, readxas, readline, settings
import brixs.beamlines.i21.cosmic_rays
# %%

# %% Read
def read(scan=None, start=0, stop=None, verbose=False, folderpath='auto', prefix='auto', filepath='auto'):
    """Returns data from nexus file

    Args:
        scan (int, optional): scan number
        start, stop (int, optional): [Only for rixs scans]
            start and stop index number for images to 
            load. If stop is None, it will get all images. Use this for loading
            scans with large number of images. Start is INCLUSIVE. Stop is 
            EXCLUSIVE. Default is start=0 and stop=None
        verbose (bool, optional): if True, prints the name of metadata is could
            not load. Default is False.
        folderpath, prefix, filepath (str or Path, optional): use this to 
            overwrite i21.settings (FOLDERPATH, PREFIX), or use filepath to 
            directly give the path to a nexus file.

    Returns:
        For rixs:
            dict: im (summed frames), ims (individual frames)
        For xas:
            dict: diff1, fy2, draincurrent, i0
        For line scans:
            dict: diff1, fy2, draincurrent, i0
    """
    try:  # try reading as RIXS scan
        return readrixs(scan=scan, start=start, stop=stop, verbose=verbose, folderpath=folderpath, prefix=prefix, filepath=filepath)
    except KeyError:
        pass

    try:  # try reading as XAS scan
        return readxas(scan=scan, verbose=verbose, folderpath=folderpath, prefix=prefix, filepath=filepath)
    except KeyError:
        pass

    try:  # try reading as line scan
        return readline(scan=scan, verbose=verbose, folderpath=folderpath, prefix=prefix, filepath=filepath)
    except Exception as e:
        raise ValueError('Cannot read scan.\n' + e)

    return

# %% Process
def _process(scan, dark=None, slope=False, calib=False,
            x_start=None, x_stop=None,
            start=None, stop=None,
            i0_norm=True,
            remove_cosmic_rays=True,
            cosmic_rays_contrast_ratio=8,
            cosmic_rays_nrows2bin=2,
            cosmic_rays_median_filter_window=100,
            cosmic_rays_moving_average_window=20,
            cosmic_rays_minimum_intensity=600,
            cosmic_rays_minimum_average=400,
            cosmic_rays_patch_size=4,      
            dark_cosmic_rays_threshold=400,
            dark_cosmic_rays_patch_size=8,  
            dark_smooth=True,
            dark_smooth_median_filter_window=5,
            dark_smooth_sigma=1,
            dark_offset_limits=(1500, 1700), 
            dark_factor_limits=(1100, 1200), 
            folderpath='auto', prefix='auto', filepath='auto'):    

    if slope == 'auto': slope = settings.SLOPE
    if calib == 'auto': calib = settings.CALIB
    d0 = None
    ds0 = None
    ds1 = None
    d2 = None
    dpes1 = None

    cr_im2 = None
    cr_im3 = None
    cr_im6 = None

    im0 = None
    ims0 = None
    ims1 = None
    pes1 = None
    d0 = None
    ds0 = None
    ds1 = None
    dpes1 = None
    d2 = None
    d3 = None
    dark_offsets = None
    dark_factors = None
    ims2 = None
    im2 = None
    im3 = None

    ##########
    # finder #
    ##########
    br.finder.kwargs = vars()
    del br.finder.kwargs['calib']
    s0 = br.finder.search()
    _search_string = br.finder._search_string  # save search_string because the finder for dark will overwrite it
    # if scan was already processed but is still running, re-run the processing
    if s0:
        if bool(s0.metadata['general']['finished']) == False:
            s0 = False

    ##############
    # processing #
    ##############
    if s0 == False:
        ##############
        # dark image #
        ##############
        if dark is not None:
            # check finder for dark image
            br.finder.kwargs = dict(dark=dark, x_start=x_start, x_stop=x_stop,
                                    dark_cosmic_rays_threshold=dark_cosmic_rays_threshold,
                                    dark_cosmic_rays_patch_size=dark_cosmic_rays_patch_size,  
                                    dark_smooth=dark_smooth,
                                    dark_smooth_median_filter_window=dark_smooth_median_filter_window,
                                    dark_smooth_sigma=dark_smooth_sigma)
            br.finder.kwargs['darkimage'] = True
            d3 = br.finder.search()
            if d3 == False:
                # dark image loading and cosmic rays removal
                data = readrixs(scan=dark, folderpath=folderpath, prefix=prefix, filepath=filepath)
                d0  = data['im'].crop(x_start=x_start, x_stop=x_stop)
                ds0 = br.Dummy([_.crop(x_start=x_start, x_stop=x_stop) for _ in data['ims']])
                
                # dark image cosmic rays removal
                ds1   = br.Dummy()
                dpes1 = br.Dummy()
                for _im in ds0:
                    if remove_cosmic_rays:
                        _pe  = _im.detect_outliers_via_threshold(threshold=dark_cosmic_rays_threshold)
                        if len(_pe) > 0:
                            try:
                                _im1 = _im.patch_from_photon_events(_pe, n=dark_cosmic_rays_patch_size)
                            except AssertionError:
                                print(f'ERROR: scan {dark}')
                                _im1 = _im
                        else:
                                _im1 = _im
                        dpes1.append(_pe)  # save cosmics rays in a list for later verification
                    else:
                        _im1 = _im 
                    ds1.append(_im1)

                # average all dark frames
                d2 = ds1[0]
                for _im in ds1[1:]:
                    d2 += _im
                d2 = d2.set_factor(1/len(ds1))

                # smooth
                if dark_smooth:
                    _img_clean = median_filter(d2.data, size=dark_smooth_median_filter_window)
                    if dark_smooth_sigma > 0:
                        _img_smooth = gaussian_filter(_img_clean, sigma=dark_smooth_sigma)
                    else:
                        _img_smooth = _img_clean
                    d3 = br.Image(data=_img_smooth)
                else:
                    d3 = d2.copy()
                
                # finder: save dark image
                br.finder.save(d3)

        ################
        # initial data #
        ################
        data = readrixs(scan=scan, folderpath=folderpath, prefix=prefix, filepath=filepath)

        # check i0
        if i0_norm:
            try:
                i0 = data['im'].metadata['instrument']['i0']
                if i0 is None:
                    i0_norm = False
                    print(f'scan {scan} does not have i0 data. Defaulting to i0_norm=False.')
            except: 
                i0_norm = False
                print(f'scan {scan} does not have i0 data. Defaulting to i0_norm=False.')
            
        # get im and ims
        if start is None and stop is None:
            im0  = data['im'].crop(x_start=x_start, x_stop=x_stop)
            ims0 = br.Dummy([_.crop(x_start=x_start, x_stop=x_stop) for _ in data['ims']])
        else:
            if start is None:
                start = 0
            if stop is None:
                stop = len(data['ims'])
            ims0 = br.Dummy([_.crop(x_start=x_start, x_stop=x_stop) for _ in data['ims'][start:stop]])
            im0 = ims0[0]
            for i, _im in enumerate(ims0[1:]):
                im0 += _im
            im0.metadata = data['im'].metadata
            if i0_norm:
                im0.metadata['instrument']['i0'] = im0.metadata['instrument']['i0'][start:stop]
            
        # cosmic ray removal
        if remove_cosmic_rays:
            ims1 = br.Dummy()
            cr_im2 = br.Dummy()
            cr_im3 = br.Dummy()
            cr_im6 = br.Dummy()
            pes1 = br.Dummy()
            for j, _im in enumerate(ims0):
                _data = _im._detect_outliers_via_row_median(contrast_ratio=cosmic_rays_contrast_ratio, slope=slope, nrows2bin=cosmic_rays_nrows2bin,
                                                        median_filter_window=cosmic_rays_median_filter_window, moving_average_window=cosmic_rays_moving_average_window,
                                                        minimum_intensity=cosmic_rays_minimum_intensity, minimum_average=cosmic_rays_minimum_average)
                _pe = _data['pe1']
                cr_im2.append(_data['im2'])
                cr_im3.append(_data['_im3'])
                cr_im6.append(_data['_im6'])
                if len(_pe) > 0:
                    try:
                        _im1 = _im.patch_from_photon_events(_pe, n=cosmic_rays_patch_size)
                    except AssertionError:
                        print(f'ERROR: scan {scan}:{j}')
                        _im1 = _im
                else:
                    _im1 = _im
                pes1.append(_pe)  # save cosmic rays in a list for later verification
                ims1.append(_im1)
        else:
            ims1 = ims0.copy()
            pes1 = None

        # dark image subtraction
        dark_offsets = []
        dark_factors = []
        if dark is not None:
            ims2 = br.Dummy()
            for _im in ims1:
                imavg0 = _im.crop(y_start=dark_offset_limits[0], y_stop=dark_offset_limits[1]).calculate_average()
                davg0  = d3.crop(y_start=dark_offset_limits[0], y_stop=dark_offset_limits[1]).calculate_average()

                imgavg1 = _im.crop(y_start=dark_factor_limits[0], y_stop=dark_factor_limits[1]).calculate_average()
                davg1   = d3.crop(y_start=dark_factor_limits[0], y_stop=dark_factor_limits[1]).calculate_average()

                dark_factor = (imgavg1 - imavg0)/(davg1 - davg0)
                d3_corrected = d3.set_factor(dark_factor)
                davg2 = d3_corrected.crop(y_start=dark_offset_limits[0], y_stop=dark_offset_limits[1]).calculate_average()

                dark_offset = imavg0-davg2
                dark_offsets.append(dark_offset)
                dark_factors.append(dark_factor)
                ims2.append(_im - d3_corrected.set_offset(dark_offset))

                # imavg = _im.crop(y_start=dark_offset_limits[0], y_stop=dark_offset_limits[1]).calculate_average()
                # davg  = d3.crop(y_start=dark_offset_limits[0], y_stop=dark_offset_limits[1]).calculate_average()
                # dark_offset = imavg-davg
                # dark_offsets.append(dark_offset)
                # ims2.append(_im - d3.set_offset(dark_offset))
        else:
            ims2 = ims1.copy()

        # plot dark image subtraction verification
        if False:
            for i in range(len(ims2)):
                fig, axes = br.subplots(1, 3, figsize=(28, 12), layout='constrained')
                plt.suptitle(f'scan image {i}')
                ims1[i].integrated_rows_vs_y_centers().plot(ax=axes[0])
                d3.set_offset(dark_offsets[i]).integrated_rows_vs_y_centers().plot(ax=axes[0])    
                # d3.integrated_rows_vs_y_centers().plot(ax=axes[0])    
                ims2[i].integrated_rows_vs_y_centers().plot(ax=axes[1])
                ims2[i].plot(ax=axes[2])
        
        # align images
        # <not needed, only single frame scans>
        
        # average all images
        if i0_norm:
            im2 = ims2[0].set_factor(1/im0.metadata['instrument']['i0'][0])
            for i, _im in enumerate(ims2[1:]):
                im2 += _im.set_factor(1/im0.metadata['instrument']['i0'][i])
            im2 = im2.set_factor(1/len(ims2))
        else:
            im2 = ims2[0]
            for i, _im in enumerate(ims2[1:]):
                im2 += _im
            im2 = im2.set_factor(1/len(ims2))
            

        # slope correction
        im3 = im2.set_vertical_shift_via_polyval(p=[-slope, 0])
        
        # calculate spectrum
        s0 = im3.integrated_rows_vs_y_centers()

        # copy metadata
        s0.metadata = im0.metadata

        # finder: save dark image
        br.finder._search_string = _search_string
        br.finder.save(s0)
    
    # energy calibration
    if calib:
        s1 = s0.set_calib(calib)
    else:
        s1 = s0
        
    # text to facilitate verification plots
    text0  = f'dark={dark}, x_start={x_start}, x_stop={x_stop}, threshold={dark_cosmic_rays_threshold}, patch_size={dark_cosmic_rays_patch_size}x{dark_cosmic_rays_patch_size}'
    text1  = f'dark={dark}, smoothening={dark_smooth}, median_filter_window={dark_smooth_median_filter_window}, sigma={dark_smooth_sigma}'
    text2  = f'scan={scan}, x_start={x_start}, x_stop={x_stop}, contrast_ratio={cosmic_rays_contrast_ratio}, slope={slope}, nrows2bin={cosmic_rays_nrows2bin}'
    text2 += '\n' + f'median_filter_window={cosmic_rays_median_filter_window}, moving_average_window={cosmic_rays_moving_average_window}, minimum_intensity={cosmic_rays_minimum_intensity}, minimum_average={cosmic_rays_minimum_average}'
    text2 += '\n' + f'patch_size={cosmic_rays_patch_size}x{cosmic_rays_patch_size}'
    text3  = f'scan={scan}, dark={dark}, dark_offset_limits={dark_offset_limits}, dark_factor_limits={dark_factor_limits}'
    text4  = f'scan={scan}, slope={slope}'

    return {'im0':im0, 'ims0':ims0, 'ims1':ims1, 'pes1':pes1, 
            'cr_im2':cr_im2, 'cr_im3':cr_im3, 'cr_im6':cr_im6, 
            'cr': dict(cosmic_rays_contrast_ratio=cosmic_rays_contrast_ratio, cosmic_rays_nrows2bin=cosmic_rays_nrows2bin, 
                       cosmic_rays_median_filter_window=cosmic_rays_median_filter_window, cosmic_rays_moving_average_window=cosmic_rays_moving_average_window, 
                       cosmic_rays_minimum_intensity=cosmic_rays_minimum_intensity, cosmic_rays_minimum_average=cosmic_rays_minimum_average),
            'd0': d0, 'ds0':ds0, 'ds1':ds1, 'dpes1':dpes1, 
            'd2':d2, 'd3':d3, 'dark_offsets':dark_offsets, 'dark_factors':dark_factors,
            'ims2':ims2, 'im2':im2, 'dark_offset_limits':dark_offset_limits, 'dark_factor_limits':dark_factor_limits,
            'im3':im3, 's0':s0, 's1':s1,
            'text0':text0, 'text1':text1, 'text2':text2, 'text3':text3, 'text4':text4}

def process(scan, dark, slope='auto', calib='auto', **kwargs):
    """Process an I21 RIXS scan and return the final spectrum.

    Processing workflow::

        - Dark frames
            load and crop (d0, ds0)
            └─ optional cosmic-ray removal (ds1, dpes1)
                └─ average (d2)
                    └─ optional smoothing (d3)
        Scan frames
            load and crop (im0, ims0)
            └─ optional cosmic-ray removal (ims1, pes1)
                └─ dark subtraction using d3 (ims2)
                    └─ average images weighted by I0 (im2)
                        └─ curvature correction (im3)
                            └─ calculate spectrum via row integration (s0)
                                └─ energy calibration (s1)

    Notes:
        * Dark-frame offsets are determined automatically by matching the
        average intensity in the detector region defined in dark_offset_limits.
        * Dark-frame step-sizes are determined automatically by matching the
            average intensity in the detector region defined in dark_offset_limits
            and dark_factor_limits.
        * If ``i0_norm=True`` but I0 data are unavailable, normalization is
        automatically disabled.
        * Processed spectra and dark images can be cached through
        ``br.finder`` to avoid repeated processing.

    Args:
        scan (int): Scan number containing the RIXS data.
        dark (int): Scan number containing the dark images.
        slope (float | str, optional): Detector curvature correction slope. 
            If ``'auto'``, ``settings.SLOPE`` is used.
        calib (float | str, optional): Energy calibration in eV/pixel. If ``'auto'``,
            ``settings.CALIB`` is used.
        x_start, x_stop (int, optional): start and stop x-pixel to include in processing. 
            If None, the full detector range is used.
        start, stop (int, optional): First (inclusive) and last (exclusive) frame 
            to include in the processing.
        i0_norm (bool, optional): If True, normalize each image by its I0 value before averaging.

        remove_cosmic_rays (bool, optional): If True, detect and patch cosmic-ray events in individual frames.
        cosmic_rays_contrast_ratio (float, optional): Minimum contrast ratio required for cosmic-ray detection.
        cosmic_rays_nrows2bin (int, optional): Number of detector rows grouped during cosmic-ray detection.
        cosmic_rays_median_filter_window (int, optional): Horizontal median-filter window used to estimate the local
            background.
        cosmic_rays_moving_average_window (int, optional): Moving-average window used for residual normalization during
            cosmic-ray detection.
        cosmic_rays_minimum_intensity (float, optional): Minimum neighbouring-pixel intensity required to validate a
            cosmic-ray candidate.
        cosmic_rays_minimum_average (float, optional): Minimum neighbouring average intensity required to validate a
            cosmic-ray candidate.
        cosmic_rays_patch_size (int, optional): n x n patch size when removing cosmic rays.

        dark_cosmic_rays_threshold (number, optional): Intensity threshold above which a pixel is considered
            an outlier. Only for dark image.
        dark_cosmic_rays_patch_size (int, optional): n x n patch size when removing cosmic rays.
        dark_smooth (bool, optional): If True, smooth dark-frame images before averaging.
        dark_smooth_median_filter_window (int, optional): Median-filter window applied to dark-frame images.
        dark_smooth_sigma (float, optional): Gaussian smoothing sigma applied after median filtering.
            Set to 0 to disable Gaussian smoothing.
        dark_offset_limits (tuple, optional): ymin and ymax row pixel values
            defining the detector region for aligning the lower part of the background.
        dark_factor_limits (tuple, optional): ymin and ymax row pixel values
            defining the detector region for aligning the upper part of the background.

    Returns:
        Dict
    """
    for name in settings.DEFAULT_PROCESSING_PARAMETERS:
        if name not in kwargs:
            kwargs[name] = settings.DEFAULT_PROCESSING_PARAMETERS[name]
    data = _process(scan=scan, dark=dark, slope=slope, calib=calib, **kwargs)
    return data['s1']

# %% line scan
def quick_linescan(scans, figsize=(27, 7), scanned_motor='auto', verbose=True, folderpath='auto', prefix='auto', filepath='auto'):
        
    # load data
    data = {}
    if isinstance(scans, Iterable) == False:
        scans = [scans, ]
    for scan in scans:
        _data = readline(scan=scan, folderpath=folderpath, prefix=prefix, filepath=filepath)
        for det in _data:
            if det not in data:
                data[det] = br.Spectra()
            data[det].append(_data[det])

    # check scanned motor
    commands = []
    key = list(data.keys())[0]
    for s in data[key]:
        commands.append(s.metadata['general']['command'].strip())
    if verbose:
        if br.all_equal([_.split()[1] for _ in commands]) == False:
            print('WARNING: It seems scans do not have the same scanned motor (to stop showing set verbose=False.\n' + str(commands))
    if scanned_motor == 'auto':
        scanned_motor = commands[0].split()[1]

    # colors
    if len(scans) <= 10:
        cmap = plt.get_cmap('tab10')
        colors = [cmap(_) for _ in range(len(scans))]
    else:
        colors = br.get_colors_from_colormap('rainbow', len(scans))
    
    # initialize figure
    fig = br.figure(figsize=figsize, layout='constrained')
    fig.get_layout_engine().set(h_pad=0, w_pad=0)
    nrows = 2
    ncols = 3
    gs = fig.add_gridspec(nrows, ncols, width_ratios=[1, 1, 1], height_ratios=[1, 0.2])
    axes = br.Axes([fig.add_subplot(gs[i, j]) for i in range(nrows-1) for j in range(ncols)], nrows=nrows-1, ncols=ncols)
    ax1 = fig.add_subplot(gs[1, :])

    # plot and set titles
    labels = ('tey', 'fy2', 'diff1')
    titles = ('TEY', 'FY2', 'DIFF1')
    lines = None
    for i, ax in enumerate(axes):
        ax.set_title(titles[i], ha='center', fontsize=22, fontweight='bold', bbox=dict(facecolor='lightgrey', edgecolor='none', pad=0))
        ax.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
        try:
            if lines is None:
                lines = data[labels[i]].plot(ax=ax, labels=scans, color=colors)
            else:
                data[labels[i]].plot(ax=ax, color=colors)
        except:
            pass
        
    # setting up ticks
    for ax in axes:
        ax.tick_params(top=True, right=True, direction='in')
        
    # labels
    axes[0].set_ylabel('Intensity')
    for ax in axes:
        ax.set_xlabel(scanned_motor)

    # legend
    ax1.axis('off')
    ax1.legend(lines, [str(_) for _ in scans], ncols=10, fontsize=8, loc='upper center')
    
    return fig, axes

# %% mesh scan
def mesh(scans, fast_motor='auto', folderpath='auto', prefix='auto', filepath='auto'):
    """Build images from a mesh measurement acquired as repeated line scans.

    This function reconstructs a 2D mesh from a sequence of scans produced
    by repeatedly scanning one motor while stepping another motor between
    scans. 

    Args:
        scans (list[int]): Scan numbers that compose the mesh.
        fast_motor (str, optional):  Motor scanned within each line scan. 
            If ``'auto'``, the motor is inferred from the scan command
            metadata. Valid values are ``'y'``, ``'z'``, and ``'auto'``.
        folderpath, prefix, filepath (str or Path, optional): use this to 
            overwrite i21.settings (FOLDERPATH, PREFIX), or use filepath to 
            directly give the path to a nexus file.

    Returns:
        dict[str, br.Image]:
            Dictionary whose keys are detector names (e.g. ``'diff1'``,
            ``'tey'``, ``'fy2'``) and whose values are reconstructed mesh
            images. Each image contains an additional attribute ``ss`` holding the
            original spectra used to construct the image.

    Notes:
        If ``fast_motor == 'z'``, spectra are stacked as columns and
        the image x-axis corresponds to manipulator ``y`` positions.
        If ``fast_motor == 'y'``, spectra are stacked as rows and
        the image y-axis corresponds to manipulator ``z`` positions.

    Raises:
        ValueError:
            If the scans were not acquired using the same line-scan command.
        ValueError:
            If ``fast_motor`` is not one of ``'y'``, ``'z'``, or ``'auto'``.

    Examples:
        Mesh acquired by scanning z while stepping y (fast_motor=z)::

            for y_ in frange(-1.5, 0.2, 0.05):
                pos y y_
                scan z 1.7 3.5 0.05 diff1 draincurrent fy2

            for y_ in frange(-1.5, 0, 0.1):
                pos y y_
                scan z -0.4 0.52 0.01 diff1 draincurrent fy2
    """
    # load data
    data = {}
    for scan in scans:
        _data = readline(scan=scan, folderpath=folderpath, prefix=prefix, filepath=filepath)
        for det in _data:
            if det not in data:
                data[det] = br.Spectra()
            data[det].append(_data[det])

    # check scanned motor
    commands = []
    key = list(data.keys())[0]
    for s in data[key]:
        commands.append(s.metadata['general']['command'].strip())
    if br.all_equal(commands) == False:
        raise ValueError('Mesh must be created by running the same line scan command.\n' + commands)
    if fast_motor == 'auto':
        fast_motor = commands[0].split()[1]
    else: 
        if fast_motor not in ('y', 'z'):
            raise ValueError('scanned motor must be "y", "z", or "auto".')

    # build image
    final = {}
    for det in data:
        if fast_motor == 'z':
            try:
                final[det] = data[det].interp().stack_spectra_as_columns()
                final[det].x_centers = [s.metadata['manipulator']['y'] for s in data[det]]
                final[det].ss = data[det]
            except TypeError:
                pass
        else:
            try:
                final[det] = data[det].interp().stack_spectra_as_rows()
                final[det].y_centers = [s.metadata['manipulator']['z'] for s in data[det]]
                final[det].ss = data[det]
            except TypeError:
                pass
            # final = {_: final[_].flip() for _ in final}

    return final

def quick_mesh(scans, fast_motor='auto', figsize=(27, 6), vmin=None, vmax=None, folderpath='auto', prefix='auto', filepath='auto'):
    """Plot TEY, FY2, and DIFF1 mesh images side-by-side.

    This is a convenience function for quickly visualizing a mesh scan
    reconstructed with :func:`mesh`. The function creates a figure
    containing TEY, FY2, and DIFF1 detector images, each with its own
    colorbar and a shared interactive crosshair cursor.

    Args:
        scans (list[int]): Scan numbers that compose the mesh.
        fast_motor (str, optional): Motor scanned within each line scan (fast axis of the mesh).
            If ``'auto'``, the motor is inferred from the scan metadata.
        figsize (tuple[float, float], optional): Figure size in cm as 
        ``(width, height)``. Defaults to ``(27, 6)``.
        folderpath, prefix, filepath (str or Path, optional): use this to 
            overwrite i21.settings (FOLDERPATH, PREFIX), or use filepath to 
            directly give the path to a nexus file.

    Returns:
        fig, axes
    """
    # load data
    data = mesh(scans, fast_motor=fast_motor, folderpath=folderpath, prefix=prefix, filepath=filepath)

    # initialize figure
    fig, axes = br.subplots(1, 3, figsize=figsize, sharex=True, sharey=True, layout='constrained')

    # colorbar settings (this function is needed for the callback)
    def lock_offset_text(event_ax):
        event_ax.yaxis.set_offset_position('left')
        event_ax.yaxis.get_offset_text().set_fontsize(6)
        event_ax.yaxis.get_offset_text().set_horizontalalignment('left')

    # plot, set colorbar, and set titles
    labels = ('tey', 'fy2', 'diff1')
    titles = ('TEY', 'FY2', 'DIFF1')
    cbs = []
    for i, ax in enumerate(axes):
        try:
            ax.set_title(titles[i], ha='center', va='top', fontsize=22, fontweight='bold', bbox=dict(facecolor='lightgrey', edgecolor='none', pad=0))
            pos = data[labels[i]].plot(ax=ax, vmin=vmin, vmax=vmax, verbose=False)
            cb = fig.colorbar(pos, ax=ax, fraction=0.046, pad=0.04)
            cb.ax.tick_params(labelsize=6)
            cb.ax.ticklabel_format(style='sci', axis='both', scilimits=(0, 0))
            lock_offset_text(cb.ax)
            cb.ax.callbacks.connect('ylim_changed', lock_offset_text)
            cbs.append(cb)
        except: 
            cbs.append(None)
            pass

    # labels
    for ax in axes:
        ax.set_xlabel('y (mm)')
    axes[0].set_ylabel('z (mm)')

    # fix offset
    # fig.canvas.draw()  # force a draw because Matplotlib computes the offset text lazily
    # for i, ax in enumerate(axes):
    #     try:
    #         cb = cbs[i]
    #         offset = cb.ax.yaxis.get_offset_text()
    #         offset.set_visible(False)
    #         cb.ax.set_title(offset.get_text(), fontsize=6,)
    #     except: 
    #         pass

    return fig, axes
# %%

# %% quick xas
def quick_xas(scans, pre=None, post=None, folderpath='auto', prefix='auto', filepath='auto'):
    """Create an dashboard for rapid XAS data inspection.

    The dashboard displays Diff1, FY2, TEY, and I0 signals for multiple scans,
    showing both raw and I0-normalized spectra. Optional pre-edge offset
    subtraction and post-edge normalization are applied to dedicated panels.

    The command in GDA should be something like ``repeat_xas(530, 560, 4)``

    Args:
        scans (scan number or list of numbers): Scan numbers to load and display.
        pre (tuple[float, float] | None, optional): Energy limits (emin, emax)
            defining the pre-edge region used to calculate and subtract a constant offset.
            If ``None``, no pre-edge correction is applied. Default is None.
        post (tuple[float, float] | None, optional): Energy limits (emin, emax)
            defining the post-edge region used to calculate a normalization factor
            after offset subtraction. The average intensity within this region is 
            scaled to unity. If ``None``, no post-edge normalization is applied.
            Default is None.
        folderpath, prefix, filepath (str or Path, optional): use this to 
            overwrite i21.settings (FOLDERPATH, PREFIX), or use filepath to 
            directly give the path to a nexus file.

    Returns:
        fig, axes
    """
    fig = br.figure(figsize=(28, 18), layout='constrained')
    fig.get_layout_engine().set(h_pad=0, w_pad=0)

    nrows = 4
    ncols = 5
    gs = fig.add_gridspec(nrows, ncols, width_ratios=[1, 1, 0.1, 1, 1], height_ratios=[1, 1, 1, 1], hspace=0)
    axes = br.Axes([fig.add_subplot(gs[i, j]) for i in range(nrows) for j in range(ncols) if j!=2], nrows=nrows, ncols=ncols-1)

    # column header
    title = '           RAW                  RAW                  I0 Norm.           I0 Norm.  '
    fig.suptitle(title, ha='center', va='top', fontsize=22, fontweight='bold', bbox=dict(facecolor='lightgrey', edgecolor='none', pad=0))

    # axes[0].set_title('RAW', ha='center', va='top', pad=4, fontsize=24, fontweight='bold', bbox=dict(facecolor='lightgrey', edgecolor='none', pad=0))
    # axes[1].set_title('RAW', ha='center', va='top', fontsize=24, fontweight='bold', bbox=dict(facecolor='lightgrey', edgecolor='none', pad=1))
    for ax in (axes[1], axes[3]):
        ax.set_title('Pre-edge offset and Post-edge norm.', fontsize=7, ha='center', fontweight='bold', bbox=dict(facecolor='lightgrey', edgecolor='none', pad=0))
    # axes[2].set_title('I0 Norm.', ha='center', va='top', fontsize=24, fontweight='bold',bbox=dict(facecolor='lightgrey', edgecolor='none', pad=1))
    # axes[3].set_title('I0 Norm.', ha='center', va='top', fontsize=24, fontweight='bold',bbox=dict(facecolor='lightgrey', edgecolor='none', pad=1))
    
    # text = plt.suptitle('.', fontsize=28, color='white')
    # fig.text(0.275, 0.999, 'RAW',ha='center',va='top',fontsize=24, fontweight='bold',bbox=dict(facecolor='lightgrey', edgecolor='none', pad=2))
    # fig.text(0.75, 0.999, 'I0 Norm.',ha='center',va='top',fontsize=24,fontweight='bold',bbox=dict(facecolor='lightgrey', edgecolor='none', pad=2),)

    # vertical line
    # vertical line (resizable, but expensive)
    axes2 = fig.add_subplot(gs[0:5, 2])
    transform = mtransforms.blended_transform_factory(axes2.transAxes, fig.transFigure)
    divider_line = fig.add_artist(lines.Line2D([0.5, 0.5], [0, 1], transform=transform, linewidth=3, color='lightgrey'))
    divider_line.set_in_layout(False)
    axes2.set_visible(False)

    # row headers
    row_headers = ('Diff1', 'FY2', 'TEY', 'I0')
    for i, ax in enumerate(axes.cols[0]):
        ax.set_ylabel(row_headers[i], fontsize=24, fontweight='bold',
                bbox=dict(facecolor='lightgrey', edgecolor='none'))
    fig.align_ylabels()

    # x label
    for ax in axes.last_row:
        ax.set_xlabel('Photon Energy (eV)')

    # setting up ticks
    for row in range(2):
        for ax in axes.rows[row]:
            ax.remove_xticklabels()
    axes[8].remove_xticklabels()
    axes[10].remove_xticklabels()
    for ax in axes:
        ax.tick_params(top=True, right=True, direction='in')

    # get data
    data = []
    if isinstance(scans, Iterable) == False:
        scans = [scans, ]
    for scan in scans:
        data.append(read(scan, folderpath=folderpath, prefix=prefix, filepath=filepath))

    # colors
    if len(scans) <= 10:
        cmap = plt.get_cmap('tab10')
        colors = [cmap(_) for _ in range(10)]
    elif len(scans) <= 20:
        cmap = plt.get_cmap('tab20')
        colors = [cmap(_) for _ in range(20)] 
    elif len(scans) <= 40:
        colors = br.get_colors_from_colormap('rainbow', len(scans))
        
    # plot
    lines1 = []
    for i, d in enumerate(data):
        for row, prefix in enumerate(('diff1', 'fy2', 'tey', 'i0')):
            s = d[prefix]

            ax = axes[row * 4 + 0]
            s.plot(ax=ax, color=colors[i])
            if row == 0:
                lines1.append(ax.lines[-1])  # save for the legend later

            if prefix != 'i0':
                ax = axes[row * 4 + 1]
                if pre is None:
                    avg = 0
                else:
                    avg = s.calculate_y_average(limits=(pre[0], pre[1]))
                    if i == 0:
                        ax.axvspan(pre[0], pre[1], color='red', alpha=0.2)
                        br.axhlines(0, ls='--', lw=0.4, color='black', ax=ax)
                if post is None:
                    factor = 1
                else:
                    factor = s.set_offset(-avg).calculate_y_average(limits=(post[0], post[1]))
                    if i == 0:
                        ax.axvspan(post[0], post[1], color='black', alpha=0.2)
                s.set_offset(-avg).set_factor(1/factor).plot(ax=ax, color=colors[i])

            ax = axes[row * 4 + 2]
            s2 = (s/d['i0'])
            s2.plot(ax=ax, color=colors[i])

            if prefix != 'i0':
                ax = axes[row * 4 + 3]
                if pre is None:
                    avg = 0
                else:
                    avg = s2.calculate_y_average(limits=(pre[0], pre[1]))
                    if i == 0:
                        ax.axvspan(pre[0], pre[1], color='red', alpha=0.2)
                        br.axhlines(0, ls='--', lw=0.4, color='black', ax=ax)
                if post is None:
                    factor = 1
                else:
                    factor = s2.set_offset(-avg).calculate_y_average(limits=(post[0], post[1]))
                    if i == 0:
                        ax.axvspan(post[0], post[1], color='black', alpha=0.2)
                s2.set_offset(-avg).set_factor(1/factor).plot(ax=ax, color=colors[i])
                
    # share axes
    br.sharex(axes)

    # textboxes
    axes[13].axis('off')
    text  = ''
    if pre is None:
        pre = ('N/A', 'N/A')
    text += '- Pre-edge (set intensity to 0)\n'
    text += f'- min: {pre[0]}\n- max: {pre[1]}\n\n'
    if post is None:
        post = ('N/A', 'N/A')
    text += '- Post-edge (normalize to 1)\n'
    text += f'- min: {post[0]}\n- max: {post[1]}'
    br.note(text, ax=axes[13], fontsize=8, loc='center', fontweight='bold',
                bbox=dict(facecolor='lightgrey', edgecolor='none'))

    # legend
    axes[15].axis('off')
    if len(scans) <= 8:
        nlcols   = 1
        fontsize = 8
    elif len(scans) <= 16:
        nlcols   = 2
        fontsize = 8  
    elif len(scans) <= 32:
        nlcols   = 3
        fontsize = 6
    else:
        nlcols   = 3
        fontsize = 4
    axes[15].legend(handles=lines1, labels=[str(_) for _ in scans], ncols=nlcols, fontsize=fontsize, loc='center')

    # axes labes ((a), (b), ...)
    br.label_axes(axes[:-3] + [axes[-2], ])

    # Turn off constrained layout because it is expensive
    fig.canvas.draw()
    fig.set_layout_engine('none')

    return fig, axes
# %%

# %% quick curvature correction
def quick_curvature(scan, nbins=3, y_start=550, y_stop=1250, floor_y_limits=None, show=True, **kwargs):
    """Calculate curvature from carbon tape scan.

    Forced:
    i0_norm is set to False
    calib is set to False
    slope is set to False

    If not specified:
    remove_cosmic_rays is False

    If scan has multiple frames, it uses the average of the frames
    """
    # fixing arguments
    kwargs['i0_norm'] = False
    kwargs['calib'] = False
    kwargs['slope'] = False
    if 'x_start' in kwargs: x_start = kwargs['x_start']
    else: x_start = None
    if 'x_stop' in kwargs: x_stop = kwargs['x_stop']
    else: x_stop = None
    if 'remove_cosmic_rays' not in kwargs:
        kwargs['remove_cosmic_rays'] = False

    # process
    previous = br.finder.search_on
    br.finder.search_on = False
    data = _process(scan, **kwargs)
    br.finder.search_on = previous

    im0 = data['im0']
    im1 = data['im2'].crop(y_start=y_start, y_stop=y_stop)
    im2 = im1.binning(nbins)
    if floor_y_limits is not None:
        ss = im2.get_columns().floor(limits=floor_y_limits)
    else:
        ss = im2.get_columns()

    shift = ss.calculate_shift(mode='peak')
    s = br.Spectrum(x=im2.x_centers, y=shift)
    polyfit = s.polyfit(1)
    slope = -polyfit['popt'][0]
    final = im1.set_vertical_shift_via_polyval([-slope, 0]).integrated_rows_vs_y_centers().floor(limits=floor_y_limits)
    fit = final.fit_peak(fixed_m=0)
    amp = fit['popt']['amp']
    fwhm = fit['popt']['fwhm']
    c = fit['popt']['c']
    if show:
        fig, axes = br.subplots(2, 4, figsize=(27, 16), layout='constrained')
        fig.suptitle(f'Curvature correction for scan {scan} (slope = {round(slope, 6)})')

        # axes 0
        im0.plot(ax=axes[0])
        if x_start is None:
            x_start = 0
        if x_stop is None:
            x_stop = 2048
        br.axvlines([x_start, x_stop], color='red', ax=axes[0])
        br.axhlines([y_start, y_stop], color='white', ax=axes[0])

        # axes 1
        im1.plot(ax=axes[1])
        if floor_y_limits is not None:
            if None in floor_y_limits:
                floor_y_limits = list(floor_y_limits)
                if floor_y_limits[0] is None: floor_y_limits[0] = min(im1.y_centers)
                if floor_y_limits[1] is None: floor_y_limits[1] = max(im1.y_centers)
            br.axhlines(floor_y_limits, color='red', ax=axes[1])

        # axes 2
        im2.plot(ax=axes[2])

        # axes 3
        im1.set_vertical_shift_via_polyval([-slope, 0]).plot(ax=axes[3])
        br.axhlines(fit['popt']['c'], color='red', ls='--', ax=axes[3])

        # axes 4
        ss.plot(ax=axes[4])
        if floor_y_limits is not None:
            br.axvlines(floor_y_limits, color='red', ax=axes[4])
        axes[4].set_xlabel('y (pixel)')
        axes[4].set_ylabel('Integrated row intensity')

        # axes 5
        ss.set_shift(shift).plot(ax=axes[5])

        # axes 6
        s.plot(ax=axes[6], marker='o', color='black')
        polyfit['fit'].plot(ax=axes[6], color='red', label=f'm = {round(-slope, 6)}')
        axes[6].legend()
        
        # axes 7
        final.plot(ax=axes[7], color='black')
        fwhm = round(fit['popt']['fwhm'], 2)
        fit['fit'].plot(ax=axes[7], color='red', label=f'FWHM = {fwhm} pixels')
        axes[7].legend()

        # labels
        for ax in axes[:4]:
            ax.set_xlabel('x (pixel)')
            ax.set_ylabel('y (pixel)')
        for ax in axes[4:6] + [axes[-1], ]:
            ax.set_xlabel('y (pixel)')
            ax.set_ylabel('Integrated row intensity')
        axes[6].set_xlabel('x (pixel)')
        axes[6].set_ylabel('Peak center relative to firs column')

        # titles
        fontsize = 9
        axes[0].set_title('1. Full image\n(red: x_start and x_stop)\n(white: y_start and y_stop)', fontsize=fontsize)
        if floor_y_limits is not None:
            axes[1].set_title('2. Croped image\n(red: floor_y_limits)', fontsize=fontsize)
        else:
            axes[1].set_title('2. Croped image', fontsize=fontsize)
        axes[2].set_title('3. Binned image', fontsize=fontsize)
        axes[3].set_title('4. Corrected image', fontsize=fontsize)

        if floor_y_limits is not None:
            axes[4].set_title('5. Columns from binned image\n(red: floor_y_limits)', fontsize=fontsize)
        else:
            axes[4].set_title('5. Columns from binned image', fontsize=fontsize)
        axes[5].set_title('6. Columns from binned image\nafter alignment', fontsize=fontsize)
        axes[6].set_title('7. Shift necessary for alignment', fontsize=fontsize)
        axes[7].set_title('8. Final spectrum', fontsize=fontsize)

    return {'slope': slope, 'amp': amp, 'fwhm':fwhm, 'c':c}

# %% quick cosmic rays removal
def quick_cosmic_rays(scan, slope, frame=0, show=True, **kwargs):
    """

    Forced.
    i0_norm is set to False
    remove_cosmic_rays set to True.
    calib is set to False

    Args:
        same as process()

    Returns:
        photon events object with cosmic rays
    """
    # fixing arguments
    kwargs['remove_cosmic_rays'] = True
    kwargs['i0_norm'] = False
    kwargs['calib'] = False

    # process
    previous = br.finder.search_on
    br.finder.search_on = False
    data = _process(scan, start=frame, stop=frame+1, **kwargs)
    br.finder.search_on = previous

    # im = read(scan, start=frame, stop=frame+1, folderpath=folderpath, prefix=prefix, filepath=filepath)
    # data = im._detect_outliers_via_row_median(contrast_ratio=contrast_ratio, slope=slope, nrows2bin=nrows2bin, median_filter_window=median_filter_window,
    #                                  moving_average_window=moving_average_window, minimum_intensity=minimum_intensity, minimum_average=minimum_average)

    pe1 = data['pes1'][0]
    if show:
        contrast_ratio = data['cr']['cosmic_rays_contrast_ratio']
        median_filter_window = data['cr']['cosmic_rays_contrast_ratio']
        im = data['im0']
        im2 = data['cr_im2'][0]
        _im3 = br.Image(data=data['cr_im3'][0])
        _im6 = br.Image(data=data['cr_im6'][0])
        ss  = im2.get_rows(max_number_of_rows=2048)
        ss3 = _im3.get_rows(max_number_of_rows=2048)
        ss6 = _im6.get_rows(max_number_of_rows=2048)
        
        # initialize figure
        fig, axes = br.subplots(1, 4, figsize=(28, 14), sharex=True, layout='constrained')
        fig.suptitle(f'Image {frame} from Scan {scan}')
        axes[0].sharey(axes[3])

        # axes 0
        im2.plot(ax=axes[0])
        pe1.set_vertical_shift_via_polyval(p=[-slope, 0]).plot(ax=axes[0], s=80, facecolors='none', edgecolors='white')
        hline = axes[0].axhline(0, color='white', lw=2)  # line where mouse was clicked

        # axes 1
        line1 = axes[1].plot(ss[0].x, ss[0].y, color='black', label='Pixel intensity')[0]
        line2 = axes[1].plot(ss3[0].x, ss3[0].y, color='red', label=f'Median (median_filter_window={median_filter_window}')[0]

        # axes 2
        line3 = axes[2].plot(ss6[0].x, ss6[0].y, color='green', label='Contrast')[0]
        br.axhlines(contrast_ratio, ax=axes[2], color='red', label=f'Contrast ratio = {contrast_ratio}')
        
        # axes 3
        im.plot(ax=axes[3])
        pe1.plot(ax=axes[3], s=80, facecolors='none', edgecolors='white')
        
        # titles
        fontsize = 8
        axes[0].set_title(f"Binned and curvature corrected image (im1)\nClick here to inspect rows", fontsize=fontsize)
        axes[1].set_title(f"Pixel intensity profile for row = 0", fontsize=fontsize)
        axes[2].set_title(f"Contrast [(pixel_profile - median)/(moving_average)]", fontsize=fontsize)
        axes[3].set_title(f"Image with final detected outliers (pe1)", fontsize=fontsize)
        
        # leg
        for ax in (axes[1], axes[2]):
            br.leg(ax=ax, fontsize='small')
        
        # click events
        press_event = {"x": None, "y": None}
        def on_press(event):
            if event.inaxes == axes[0]:
                press_event["x"] = event.x
                press_event["y"] = event.y
        
        def on_release(event):
            # ignore clicks outside image
            if event.inaxes != axes[0]:
                return  
            # if no press recorded - ignore
            if press_event["x"] is None:
                return
            # this was a drag - ignore
            dx = abs(event.x - press_event["x"])
            dy = abs(event.y - press_event["y"])
            # threshold in pixels
            if dx > 5 or dy > 5:
                return  
            
            # get clicked row
            y = int(round(event.ydata, 0))
            if y < 0 or y >= (2048-slope*2048-10):
                return

            # update plot axes 0
            hline.set_ydata((y, y))  # clicked line

            # update plot axes 1
            s = ss[int(y/nrows2bin)]
            line1.set_data(s.x, s.y)
            axes[1].set_ylim(np.min(s.y), np.max(s.y))
            axes[1].set_title(f"Pixel intensity profile for row = {y}", fontsize=fontsize)
        
            _s3 = ss3[int(y/nrows2bin)]
            line2.set_data(_s3.x, _s3.y)

            # update plot axes 2
            _s6 = ss6[int(y/nrows2bin)]
            line3.set_data(_s6.x, _s6.y)
            if np.max(_s6.y) < contrast_ratio:
                axes[2].set_ylim(np.min(_s6.y), contrast_ratio+contrast_ratio*0.2)
            else:
                axes[2].set_ylim(np.min(_s6.y), np.max(_s6.y))
        
            _ = fig.canvas.draw_idle()
        
        # Connect event
        _ = fig.canvas.mpl_connect('button_press_event', on_press)
        _ = fig.canvas.mpl_connect('button_release_event', on_release)
    return pe1

# %% quick optimization
def quick_optimization_with_curvature_correction(scans, metadata=None, nbins=4, y_start=300, y_stop=800, floor_y_limits=None, show_curvature=False, show=True, **kwargs):
    """ct as a function of a motor value
    sgm pitch 
    cff

    Args:
        metadata (tuple, optional): (type, name), example, ('instrument', 'sgmpitch')

    """
    data  = [quick_curvature(scan, nbins=nbins, y_start=y_start, y_stop=y_stop, floor_y_limits=floor_y_limits, show=show_curvature, **kwargs) for scan in scans]
    final = {k: [d[k] for d in data] for k in data[0]}

    if metadata is not None:
        values = []
        for scan in scans: 
            values.append(get_metadata(scan)[metadata[0]][metadata[1]])
    else:
        values = np.arange(len(scans))

    if show:
        slope = np.mean(final['slope'])

        fig, axes = br.subplots(1, 2, figsize=(27, 10))
        if metadata is not None:
            fig.suptitle(metadata[0] + '/' + metadata[1])

        ax = axes[0]
        ax.plot(values, final['slope'], color='black', marker='o')
        br.axhlines(slope, color='red', ax=ax)
        ax.set_xlabel('Scan count')
        ax.set_ylabel('Slope')
        
        ax = axes[1]
        ax.plot(values, final['fwhm'], color='black', marker='o')
        ax.set_xlabel('Scan count')
        ax.set_ylabel('FWHM (pixel)')
    return 

def quick_calibration_with_curvature_correction(scan, nbins=4, y_start=300, y_stop=800, floor_y_limits=None, show_curvature=False, show=True, **kwargs):
    """

    """
    settings.METADATA['rixs']['E2'] = ['instrument', 'list_float', '/entry/instrument/energy/value']
    values = get_metadata(scan)['instrument']['E2']

    number_of_frames = len(values)
    data  = [quick_curvature(scan, start=i, stop=i+1, nbins=nbins, y_start=y_start, y_stop=y_stop, floor_y_limits=floor_y_limits, show=show_curvature, **kwargs) for i in range(number_of_frames)]
    final = {k: [d[k] for d in data] for k in data[0]}

    s = br.Spectrum(x=values, y=final['c'])
    fit = s.polyfit(deg=1)
    calib = 1/fit['popt'][0]

    if show:
        slope = np.mean(final['slope'])
        fwhm  = np.mean(final['fwhm'])

        fig, axes = br.subplots(1, 3, figsize=(27, 10))
        fig.suptitle(f'Energy calibration scan {scan}: {round(calib, 6)} ev/pixel')

        ax = axes[0]
        ax.plot(values, final['slope'], color='black', marker='o')
        br.axhlines(slope, color='red', ax=ax)
        ax.set_xlabel('Photon energy (eV)')
        ax.set_ylabel('Slope')
        
        ax = axes[1]
        ax.plot(values, final['fwhm'], color='black', marker='o')
        br.axhlines(fwhm, color='red', ax=ax)
        ax.set_xlabel('Photon energy (eV)')
        ax.set_ylabel('FWHM (pixel)')

        ax = axes[2]
        ax.plot(values, final['c'], color='black', marker='o')
        fit['fit'].plot(ax=ax, color='red')
        ax.set_xlabel('Photon energy (eV)')
        ax.set_ylabel('y pixel (pixel)')
    return {'slope': slope, 'calib': calib, 'fwhm':fwhm*calib}





# %% TENTATIVE: ONLINE PROCESSING
# %% ===================================================================== %% #
# %% ===================================================================== %% #
import time
settings.DARKS = []

def filename2scan(filename):
    """Convert filename to scan number. The filename must have the format 
    '.../i21-1234.nxs' where 1234 is the scan number.

    Args:
        filename (string): name of the hdf5 file.
    
    Returns:
        int: scan number
    """
    return int(filename.split('-')[1].split('.')[0])

def autoprocessing(refresh_time=10, extension='nxs', verbose=True, prefix='auto', folderpath='auto'):
    # watch folder and update sheet
    starttime = time.monotonic()

    if prefix     == 'auto': prefix     = settings.PREFIX
    if folderpath == 'auto': folderpath = settings.FOLDERPATH

    folderpath = Path(folderpath)
    last_processed_scan = 0
    darks = settings.DARKS
    if verbose: print("Watching for new scans...")
    while True:
        # get last scan number from data directory    
        fl = br.filelist(folderpath, string=f'*{extension}')
        scan_numbers_from_directory = np.sort([filename2scan(_.name) for _ in fl])
        if len(scan_numbers_from_directory) == 0:
            last_scan_from_directory = 0
        else:
            last_scan_from_directory = max(scan_numbers_from_directory)

        # check if there are new scans in the directory that are not processed
        # if so process scan
        if last_scan_from_directory > last_processed_scan:
            for scan_index in np.where(scan_numbers_from_directory > last_processed_scan)[0]:
                scan = scan_numbers_from_directory[scan_index]
                try:
                    if verbose: print(f'processing scan {scan}')
                    command = get_metadata(scan)['general']['command']
                    if command.startswith('scan ds'):
                        count_time = int(round(get_metadata(scan)['andor']['count_time'][0]))
                        if count_time not in darks:
                            available_times = list(darks.keys())
                            count_time = available_times[br.index(available_times, count_time, closest=True)]
                        dark = darks[count_time][br.index(darks[count_time], scan, closest=True)]
                        _ = process(scan=scan, dark=dark)
                except Exception as e:
                    if verbose: print(f"Error processing scan {scan}")
                    if verbose: print(e)
                last_processed_scan = scan

        if darks != settings.DARKS:
            darks = settings.DARKS
            last_processed_scan = 0
            print('fffffffasdflsadhfsaodfhlksajdfhsdf\n\n\n')

        # sleep
        time.sleep(float(refresh_time) - ((time.monotonic() - starttime) % float(refresh_time)))
    return

# %% ===================================================================== %% #
# %% ===================================================================== %% #
# %% OLD: Need verification
# Verify
def verify(scan=None, start=0, stop=None, 
           dark=None, darkfactor=1, darkoffset='auto', 
           slope=None, calib=None, norm_exposure=True,
           norm_i0=False, norm_eslit=False, 
           x_start=None, x_stop=None, y_start=None, y_stop=None, 
           verbose=False,
            folderpath='auto', prefix='auto', filepath='auto'):
    """Returns a dict. with data from every processing step and opens a figure with all images

    Note:
        Use keyboard arrow keys to flip detector single images (does not work 
            with inline jupyter plots)


    Note:
        See below the step-by-step processing of the detector image.

                                     d0* ---(crop, norm)---> d1
                                                              |
    ims0 ---(add ims)---> im0 ---(crop, norm)---> im1  ---(dark sub.)---> im2  ---(curv.)---> im3  ---(integration, calib)---> s
    ims0 -----------(crop, norm)----------------> ims1 ---(dark sub.)---> ims2 ---(curv.)---> ims3 ---(integration, calib)---> ss

    *dark image (d0) is assumed to be composed by only one image.

    Args:
        scan (int): scan number.
        folderpath (str or path): folderpath where .nxs files are stored.
        dark (int or Image): dark image is loaded, crop, normalized and
            subtracted from scan. If dark is an br.Image object is
            subtracted (without normalization) directly from scan. Dark image is 
            also cropped with x_start, ..., y_stop and the cropping normalization
            is applied.
        darkfactor (number): dark image intensity will be multiplied by this 
            number.
        darkoffset (number): this value will be added to the dark image 
            intensity. Note that darkfactor is applied before dark offset, so
            one has to account for the multiplicative factor, i.e. 

                dark = dark.set_factor(darkfactor).set_offset(darkoffset)

            if darkoffset='auto', a suitable value will be found as to match
            the average intensity of pixel row integration between y=1500 and 
            y=1800 of the scan and dark image.
        curv (None or list): if not None, curv must be 1D array of polynomial 
            coefficients (including coefficients equal to zero) from 
            highest degree to the constant term 
                
                [f(x_centers) = curv[n]*x**n + curv[n-1]*x**n-1 + ... + curv[0]]
        calib (number, optional): if not None, the x axis of the final spectrum
            is multipled by calib. This number must have unit of eV/pixel.
        norm_i0 (True, optional): if True, image intensity is normalized by
            I0. Default is False.
        norm_exposure (True, optional): if True, image intensity is 
            normalized by image exposure. Default is True.
        norm_eslit (True, optional): if True, image intensity is normalized by
            the size of the exit slit. Default is False.
        x_start, x_stop, y_start, y_stop (int): pixel range in terms of
            x_centers and y_centers. Interval is inclusive. Use None to 
            indicate the edge of the image.
        verbose (bool, optional): if True, a message will print after each
            processing step, also a error message will be printed when 
            metadata cannot be retrieved. Default is False.

    Returns:
        dict {im0, ims0, im1, ims1, d0, d1, im2, ims2, im3, ims3, s, ss}
    """
    ################
    # process data #
    ################
    data = _process(scan=scan, start=start, stop=stop,  dark=dark, 
                    darkfactor=darkfactor, darkoffset=darkoffset, slope=slope, 
                    norm_i0=norm_i0, norm_exposure=norm_exposure,
                    norm_eslit=norm_eslit, calib=calib, 
                    x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop,
                    verbose=verbose,
             folderpath=folderpath, prefix=prefix, filepath=filepath)

    im  = data['im3']
    ims = data['ims3']
    s   = data['s']
    number_of_images = len(ims)

    #######################
    # initial definitions #
    #######################
    ims.__i  = 0

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
    def keyboard(event, ims, axes):
        if event.key == 'right':
            # increase i
            ims.__i = ims.__i + 1
            if ims.__i >= len(ims):
                ims.__i = len(ims) - 1
    
        elif event.key == 'left':# or event.key == 'down':
            # decrease i
            ims.__i = ims.__i - 1
            if ims.__i < 0:
                ims.__i = 0
        else:
            return
            
        # clear axis
        axes[0].cla()
        axes[1].cla()
        axes[2].cla()
        
        # set labels
        # axes[0].set_xlabel('x (pixel)')
        axes[0].set_ylabel('y (pixel)')
        axes[1].set_xlabel('counts/bin')
        
        # change title
        axes[0].set_title('Use left/right keyboard keys to flip through images: ' + str(ims.__i) + '/' + str(number_of_images-1), fontsize='small')

        # plot axes 0
        ims[ims.__i].plot(ax=axes[0])

        # plot axes 1
        ims[ims.__i].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[1])
        
        # plot axes 2
        ims[ims.__i].integrated_columns_vs_x_centers().plot(ax=axes[2])
    
        plt.draw()

    #######################
    # initializing figure #
    #######################
    fig, axes = br.subplots(4, 2, width_ratios=[4, 1], height_ratios=[2, 1, 2, 1], wspace=0.1, hspace=0, figsize=(18, 20))
    plt.subplots_adjust(top=0.97, right=0.99, bottom=0.15)
    
    for i in (4, 5, 6, 7):
        axes[i].ymove(-0.1)
    
    for i in (0, 4):
        axes[i].remove_xticklabels()
    
    for i in (1, 5):
        axes[i].remove_yticklabels()
    
    ##############
    # share axis #
    ##############
    br.sharey([axes[0], axes[1]])
    br.sharey([axes[0], axes[4]])
    br.sharey([axes[0], axes[5]])
    
    br.sharex([axes[0], axes[2]])
    br.sharex([axes[0], axes[4]])
    br.sharex([axes[0], axes[6]])    

    ##################
    # error messages #
    ##################
    # if pe1.RIXSCam_NumImages != len(pes1):
    #     fig.suptitle(f'WARNING: # of images ({data.dims['trainId']}) inside folder is different from # of acquired images ({int(pe1.RIXSCam_NumImages)})', color='red')

    ######################
    # set initial titles #
    ######################
    axes[0].set_title('Use left/right keyboard keys to flip through images: ' + str(0) + '/' + str(number_of_images-1), fontsize='small')
    axes[4].set_title('Summed images', fontsize='small')
    
    ########
    # plot #
    ########
    # plot initial images (axes 0)
    ims[0].plot(ax=axes[0])
    
    # plot initial spectra (axes 1, 2)
    ims[0].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[1])
    ims[0].integrated_columns_vs_x_centers().plot(ax=axes[2])
    
    # plot images summed (axes 4)
    im.plot(ax=axes[4])

    # plot spectra summed (axes 5, 6)
    s.switch_xy().plot(ax=axes[5], label='direct sum')
    _ss = br.Spectra()
    for _im in ims:
        _ss.append(_im.integrated_rows_vs_y_centers())
    _ss.align().calculate_average().switch_xy().plot(ax=axes[5], label='align spectra before sum')
 
    im.integrated_columns_vs_x_centers().plot(ax=axes[6])    

    ##############
    # set labels #
    ##############
    for i in (0, 4):
        axes[i].set_ylabel('y (pixel)', fontsize='x-small')
    
    for i in (1, 5):
        axes[i].set_xlabel('counts/pixel row', fontsize='x-small')
    
    for i in (2, 6):
        axes[i].set_xlabel('x (pixel)', fontsize='x-small')
    
    for i in (2, 6):
        axes[i].set_ylabel('counts/pixel column', fontsize='x-small')
        
    ##########
    # legend #
    ##########
    br.leg(ax=axes[5], fontsize='xx-small')
    
    #################
    # remove unused #
    #################
    axes[3].remove()
    axes[7].remove()

    ######################
    # register callbacks #
    ######################
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ims=ims, axes=axes))
    return data

# verify dark
def verify_dark(scan, start=0, stop=None,
                dark=None, darkfactor=1, darkoffset='auto', 
                norm_i0=False, norm_eslit=False, norm_exposure=True,
                x_start=None, x_stop=None, y_start=None, y_stop=None, 
                verbose=False, folderpath='auto', prefix='auto', filepath='auto'):
    """Opens a figure comparing images with dark image
    
    Note:
        See below the step-by-step processing of the detector image.

                                     d0* ---(crop, norm)---> d1
                                                              |
    ims0 ---(add ims)---> im0 ---(crop, norm)---> im1  ---(dark sub.)---> im2  ---(curv.)---> im3  ---(integration, calib)---> s
    ims0 -----------(crop, norm)----------------> ims1 ---(dark sub.)---> ims2 ---(curv.)---> ims3 ---(integration, calib)---> ss

    *dark image (d0) is assumed to be composed by only one image.

    Args:
        scan (int): scan number.
        folderpath (str or path): folderpath where .nxs files are stored.
        dark (int or Image): dark image is loaded, crop, normalized and
            subtracted from scan. If dark is an br.Image object is
            subtracted (without normalization) directly from scan. Dark image is 
            also cropped with x_start, ..., y_stop and the cropping normalization
            is applied.
        darkfactor (number): dark image intensity will be multiplied by this 
            number.
        darkoffset (number): this value will be added to the dark image 
            intensity. Note that darkfactor is applied before dark offset, so
            one has to account for the multiplicative factor, i.e. 

                dark = dark.set_factor(darkfactor).set_offset(darkoffset)

            if darkoffset='auto', a suitable value will be found as to match
            the average intensity of pixel row integration between y=1500 and 
            y=1800 of the scan and dark image.
        norm_i0 (True, optional): if True, image intensity is normalized by
            I0. Default is False.
        norm_exposure (True, optional): if True, image intensity is 
            normalized by image exposure. Default is True.
        norm_eslit (True, optional): if True, image intensity is normalized by
            the size of the exit slit. Default is False.
        x_start, x_stop, y_start, y_stop (int): pixel range in terms of
            x_centers and y_centers. Interval is inclusive. Use None to 
            indicate the edge of the image.
        verbose (bool, optional): if True, a message will print after each
            processing step, also a error message will be printed when 
            metadata cannot be retrieved. Default is False.

    Returns:
        if darkoffset is 'auto', returns the calculated offset
    """
    ################
    # process data #
    ################
    # data = _process(scan=scan, folderpath=folderpath, dark=dark, darkfactor=1, darkoffset=0, curv=None, norm_i0=norm_i0, norm_exposure=norm_exposure, norm_eslit=norm_eslit, calib=None, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
    data = _process(scan=scan, start=start, stop=stop,  dark=dark, 
                    darkfactor=darkfactor, darkoffset=darkoffset, slope=None, 
                    norm_i0=norm_i0, norm_exposure=norm_exposure,
                    norm_eslit=norm_eslit, calib=None, 
                    x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop,
                    verbose=verbose,
                    folderpath=folderpath, prefix=prefix, filepath=filepath)

    im   = data['im1']
    ims  = data['ims1']
    # s    = data['s']
    number_of_images = len(ims)

    ####################
    # auto dark offset #
    ####################
    if dark is not None:
        if darkoffset == 'auto':
            _temp     = im.integrated_rows_vs_y_centers()
            _bkg      = _temp.crop(1500, 1800).calculate_y_average()
            _bkg_dark = data['d1'].integrated_rows_vs_y_centers().crop(1500, 1800).calculate_y_average()
            darkoffset = (_bkg - _bkg_dark)/len(_temp.x)
            if verbose:
                print(f'calculated dark offset = {darkoffset}')
        data['d1'] = data['d1'].set_factor(darkfactor).set_offset(darkoffset)
    
    #######################
    # initial definitions #
    #######################
    ims.__i  = 0

    ######################
    # change keybindings #
    ######################
    try:
        matplotlib.rcParams['keymap.back'].remove('.')
        matplotlib.rcParams['keymap.forward'].remove(',')
    except ValueError:
        pass

    ###################
    # keyboard events #
    ###################
    def keyboard(event, ims, axes):
        if event.key == 'right':
            # increase i
            ims.__i = ims.__i + 1
            if ims.__i >= len(ims):
                ims.__i = len(ims) - 1
    
        elif event.key == 'left':# or event.key == 'down':
            # decrease i
            ims.__i = ims.__i - 1
            if ims.__i < 0:
                ims.__i = 0
        else:
            return
            
        # clear axis
        axes[0].cla()
        axes[1].cla()
        axes[2].cla()
        
        # set labels
        # axes[0].set_xlabel('x (pixel)')
        axes[0].set_ylabel('y (pixel)')
        axes[1].set_xlabel('counts/bin')
        
        # change title
        axes[0].set_title('Use ,/. keyboard keys to flip through images: ' + str(ims.__i) + '/' + str(number_of_images-1), fontsize='small')

        # plot axes 0
        ims[ims.__i].plot(ax=axes[0])

        # plot axes 1
        ims[ims.__i].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[1])
        data['d1'].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[1])
        # data['d1'].integrated_rows_vs_y_centers().set_factor(darkfactor).set_offset(darkoffset).switch_xy().plot(ax=axes[1])
        
        # plot axes 2
        ims[ims.__i].integrated_columns_vs_x_centers().plot(ax=axes[2])
        data['d1'].integrated_columns_vs_x_centers().plot(ax=axes[2])
        # data['d1'].integrated_columns_vs_x_centers().set_factor(darkfactor).set_offset(darkoffset).plot(ax=axes[2])
    
        plt.draw()

    #######################
    # initializing figure #
    #######################
    fig, axes = br.subplots(4, 2, width_ratios=[4, 1], height_ratios=[2, 1, 2, 1], wspace=0.1, hspace=0, figsize=(18, 20))
    plt.subplots_adjust(top=0.97, right=0.99, bottom=0.15)
    
    for i in (4, 5, 6, 7):
        axes[i].ymove(-0.1)
    
    for i in (0, 4):
        axes[i].remove_xticklabels()

    for i in (1, 5):
        axes[i].remove_yticklabels()
    
    ##############
    # share axis #
    ##############
    br.sharey([axes[0], axes[1]])
    br.sharey([axes[0], axes[4]])
    br.sharey([axes[0], axes[5]])
    
    br.sharex([axes[0], axes[2]])
    br.sharex([axes[0], axes[4]])
    br.sharex([axes[0], axes[6]])    

    ##################
    # error messages #
    ##################
    # if pe1.RIXSCam_NumImages != len(pes1):
    #     fig.suptitle(f'WARNING: # of images ({data.dims['trainId']}) inside folder is different from # of acquired images ({int(pe1.RIXSCam_NumImages)})', color='red')

    ######################
    # set initial titles #
    ######################
    axes[0].set_title('Use left/right keyboard keys to flip through images: ' + str(0) + '/' + str(number_of_images-1), fontsize='small')
    axes[4].set_title('Summed images (note that dark image must align with summed images here - Not individual images up)', fontsize='small')
    
    ########
    # plot #
    ########
    # plot initial images (axes 0)
    ims[0].plot(ax=axes[0])
    
    # plot initial spectra (axes 1, 2)
    ims[0].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[1])
    ims[0].integrated_columns_vs_x_centers().plot(ax=axes[2])

    data['d1'].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[1])
    data['d1'].integrated_columns_vs_x_centers().plot(ax=axes[2])
    # data['d1'].integrated_rows_vs_y_centers().set_factor(darkfactor).set_offset(darkoffset).switch_xy().plot(ax=axes[1])
    # data['d1'].integrated_columns_vs_x_centers().set_factor(darkfactor).set_offset(darkoffset).plot(ax=axes[2])
    
    # plot images summed (axes 4)
    im.plot(ax=axes[4])

    # plot spectra summed (axes 5, 6)
    # s.switch_xy().plot(ax=axes[5], label='direct sum')
    _ss = br.Spectra()
    for _im in ims:
        _ss.append(_im.integrated_rows_vs_y_centers())
    _ss.align().calculate_average().switch_xy().plot(ax=axes[5])
    data['d1'].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[5])
    # data['d1'].integrated_rows_vs_y_centers().set_factor(darkfactor).set_offset(darkoffset).switch_xy().plot(ax=axes[5])
 
    im.integrated_columns_vs_x_centers().plot(ax=axes[6])
    data['d1'].integrated_columns_vs_x_centers().plot(ax=axes[6])
    # data['d1'].integrated_columns_vs_x_centers().set_factor(darkfactor).set_offset(darkoffset).plot(ax=axes[6])
    

    ##############
    # set labels #
    ##############
    for i in (0, 4):
        axes[i].set_ylabel('y (pixel)', fontsize='x-small')
    
    for i in (1, 5):
        axes[i].set_xlabel('counts/pixel row', fontsize='x-small')
    
    for i in (2, 6):
        axes[i].set_xlabel('x (pixel)', fontsize='x-small')
    
    for i in (2, 6):
        axes[i].set_ylabel('counts/pixel column', fontsize='x-small')
        
    ##########
    # legend #
    ##########
    br.leg(ax=axes[5], fontsize='xx-small')
    
    #################
    # remove unused #
    #################
    axes[3].remove()
    axes[7].remove()

    ######################
    # register callbacks #
    ######################
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ims=ims, axes=axes))
    
    return darkoffset

# verify curvature correction
def verify_curv(scan, start=0, stop=None, popt=None, ncols=16, 
                nrows=None, dark=None, darkfactor=1, darkoffset='auto', 
                norm_i0=False, norm_eslit=False, norm_exposure=True,
                x_start=None, x_stop=None, y_start=None, y_stop=None, 
                deg=2,
                folderpath='auto', prefix='auto', filepath='auto' ):
    """Returns a dict. with data from every processing step and opens a figure with all images
    
    Note:
        See below the step-by-step processing of the detector image.

                                     d0* ---(crop, norm)---> d1
                                                              |
    ims0 ---(add ims)---> im0 ---(crop, norm)---> im1  ---(dark sub.)---> im2  ---(curv.)---> im3  ---(integration, calib)---> s
    ims0 -----------(crop, norm)----------------> ims1 ---(dark sub.)---> ims2 ---(curv.)---> ims3 ---(integration, calib)---> ss

    *dark image (d0) is assumed to be composed by only one image.

    Args:
        scan (int): scan number.
        folderpath (str or path): folderpath where .nxs files are stored.
        popt (list or None, optional): Curvature polynomial values. 
            1D array of polynomial 
            coefficients (including coefficients equal to zero) from 
            highest degree to the constant term 
                
                [f(x_centers) = curv[n]*x**n + curv[n-1]*x**n-1 + ... + curv[0]]

            If None, curvature will be calculated by binning the image and 
            aligning pixel columns using cross-correlataion.
        ncols (int or None, optional): number of pixel columns (bins). Default 
            is 16.
        nrows (int or None, optional): number of pixel rows (bins). Default 
            is None (no horizontal binning).
        dark (int or Image): dark image is loaded, crop, normalized and
            subtracted from scan. If dark is an br.Image object is
            subtracted (without normalization) directly from scan. Dark image is 
            also cropped with x_start, ..., y_stop and the cropping normalization
            is applied.
        darkfactor (number): dark image intensity will be multiplied by this 
            number.
        darkoffset (number): this value will be added to the dark image 
            intensity. Note that darkfactor is applied before dark offset, so
            one has to account for the multiplicative factor, i.e. 

                dark = dark.set_factor(darkfactor).set_offset(darkoffset)

            if darkoffset='auto', a suitable value will be found as to match
            the average intensity of pixel row integration between y=1500 and 
            y=1800 of the scan and dark image.
        curv (None or list): if not None, 
        calib (number, optional): if not None, the x axis of the final spectrum
            is multipled by calib. This number must have unit of eV/pixel.
        norm_i0 (True, optional): if True, image intensity is normalized by
            I0. Default is False.
        norm_exposure (True, optional): if True, image intensity is 
            normalized by image exposure. Default is True.
        norm_eslit (True, optional): if True, image intensity is normalized by
            the size of the exit slit. Default is False.
        x_start, x_stop, y_start, y_stop (int): pixel range in terms of
            x_centers and y_centers. Interval is inclusive. Use None to 
            indicate the edge of the image.

    Returns:
        Curvature polynomial values. 1D array of polynomial 
        coefficients (including coefficients equal to zero) from 
        highest degree to the constant term 
                
            [f(x_centers) = curv[n]*x**n + curv[n-1]*x**n-1 + ... + curv[0]]
    """
    ################
    # process data #
    ################
    # data = _process(scan=scan, start=start, stop=stop,  
    #                 dark=dark, darkfactor=darkfactor, darkoffset=darkoffset, 
    #                 norm_i0=norm_i0, norm_exposure=norm_exposure, norm_eslit=norm_eslit, 
    #                 x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, 
    #                 curv=None, calib=None,
    #          folderpath=folderpath, prefix=prefix, filepath=filepath
    #          )
    data = _process(scan=scan, start=start, stop=stop,  dark=dark, 
                    darkfactor=darkfactor, darkoffset=darkoffset, slope=None, 
                    norm_i0=norm_i0, norm_exposure=norm_exposure,
                    norm_eslit=norm_eslit, calib=None, 
                    x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop,
             folderpath=folderpath, prefix=prefix, filepath=filepath)


    im = data['im2']
    
    # if polynomial parameters are not given, calculate curvature
    remove_axes_2 = False
    if popt is None:
        if ncols is None: ncols = im.shape[1]
        if nrows is None: nrows = im.shape[0]

        reduced = im.binning(ncols=ncols, nrows=nrows).crop(y_start=y_start, y_stop=y_stop).floor()

        ss = reduced.columns.floor()
        values = ss.calculate_shift(mode='cc')

        s = br.Spectrum(x=reduced.x_centers, y=values)
        polyfit = s.polyfit(deg=deg)
        fit = polyfit['fit']             
        popt = polyfit['popt']             
    else:
        reduced = im
        fit = br.Spectrum(x=im.x_centers, y=np.polyval(popt, im.x_centers))
        ss  = br.Spectra()
        remove_axes_2 = True


    # apply curvature corredctions
    im2 = im.set_vertical_shift_via_polyval(p=popt)

    #######################
    # initializing figure #
    #######################
    fig, axes = br.subplots(1, 4, width_ratios=[3, 3, 1, 3], sharey=True, wspace=0, figsize=(40, 12))
    plt.subplots_adjust(left=0.05, top=0.945, right=0.9, bottom=0.1)

    for i in (1, 2, 3, 3):
        axes[i].xmove(0.04)
    axes[2].remove_yticklabels()
        

    ##############
    # share axis #
    ##############
    # br.sharey([axes[1], axes[2]])

    ######################
    # set initial titles #
    ######################
    axes[0].set_title('Full image', fontsize='small')
    axes[1].set_title(f'reduced image', fontsize='small')
    axes[3].set_title('corrected image', fontsize='small')

    ########
    # plot #
    ########
    im.plot(ax=axes[0])

    reduced.plot(ax=axes[1])
    ss.switch_xy().plot(ax=axes[2])

    if len(ss) > 0:
        offset = ss[0].get_x_where_y_is_max()
    else:
        offset = im2.integrated_rows_vs_y_centers().get_x_where_y_is_max()
    fit.crop(x_start, x_stop).set_factor(-1).set_offset(offset).plot(ax=axes[1], color='black')
        
    im2.plot(ax=axes[3])

    ##############
    # set labels #
    ##############
    for i in (0, ):
        axes[i].set_ylabel('y (pixel)', fontsize='x-small')

    for i in (2, ):
        axes[i].set_xlabel('counts/pixel row', fontsize='x-small')

    for i in (0, 1, 3):
        axes[i].set_xlabel('x (pixel)', fontsize='x-small')

    # remove axes 2 if curvature was already passed as an argument
    if remove_axes_2:
        axes[2].remove()

    return popt
# %%

# %% OBSOLETE
def _process_OLD1(scan=None, start=0, stop=None, dark=None, darkfactor=1, darkoffset='auto', slope=None,
             norm_exposure=True, norm_i0=False, norm_eslit=False, calib=None, 
             x_start=None, x_stop=None, y_start=None, y_stop=None, verbose=False,
             folderpath='auto', prefix='auto', filepath='auto'):
    """Returns dict. with data from all processing steps from images.

                                      d0* ---(crop, norm)---> d1
                                                              |
    ims0 ---(add ims)---> im0 ---(crop, norm)---> im1  ---(dark sub.)---> im2  ---(curv.)---> im3  ---(integration, calib)---> s
    ims0 -----------(crop, norm)----------------> ims1 ---(dark sub.)---> ims2 ---(curv.)---> ims3 ---(integration, calib)---> ss

    *dark image (d0) is assumed to be composed by only one image.

    Add images: 
        each scan is composed of multiple images. direct sum. No need to align them
    Normalization: 
        image intensity (ims0 and im0) are normalized by number of images and 
        size of image (crop size). Images can also be normalized by I0, 
        exposure time, and exit slit if respective arguments are True.
    Dark image subtraction: 
        Normalized and cropped dark image is subtracted from data image. 
        A multiplicative and additive factors can be applied to the dark image.
    Curvature correction:
        Fix image curvature.
    Integration:
        Horizontal pixel integration
    Energy calibration
        Converts spectrum x axis from pixel to energy loss
    
    Args:
        scan (int or br.Image): scan number or image.
        start, stop (int, optional): [Only for rixs scans]
            start and stop index number for images to 
            load. If stop is None, it will get all images. Use this for loading
            scans with large number of images. Start is INCLUSIVE. Stop is 
            EXCLUSIVE. Default is start=0 and stop=None
        dark (int or Image): dark image is loaded, crop, normalized and
            subtracted from scan. If dark is an br.Image object is
            subtracted (without normalization) directly from scan. Dark image is 
            also cropped with x_start, ..., y_stop and the cropping normalization
            is applied.
        darkfactor (number): dark image intensity will be multiplied by this 
            number.
        darkoffset (number): this value will be added to the dark image 
            intensity. Note that darkfactor is applied before dark offset, so
            one has to account for the multiplicative factor, i.e. 

                dark = dark.set_factor(darkfactor).set_offset(darkoffset)

            if darkoffset='auto', a suitable value will be found as to match
            the average intensity of pixel row integration between y=1500 and 
            y=1800 of the scan and dark image.
        slope (None or list): must be number or auto.
        norm_i0 (True, optional): if True, image intensity is normalized by
            I0. Default is False.
        norm_exposure (True, optional): if True, image intensity is 
            normalized by image exposure. Default is True.
        norm_eslit (True, optional): if True, image intensity is normalized by
            the size of the exit slit. Default is False.
        calib (number, optional): if not None, the x axis of the final spectrum
            is multipled by calib. This number must have unit of eV/pixel.
        x_start, x_stop, y_start, y_stop (int): pixel range in terms of
            x_centers and y_centers. Interval is inclusive. Use None to 
            indicate the edge of the image.
        verbose (bool, optional): if True, a message will print after each
            processing step, also a error message will be printed when 
            metadata cannot be retrieved. Default is False.
        folderpath, prefix, filepath (str or Path, optional): use this to 
            overwrite i21.settings (FOLDERPATH, PREFIX), or use filepath to 
            directly give the path to a nexus file.
        
    Returns:
        dict {im0, ims0, im1, ims1, d0, d1, im2, ims2, im3, ims3, s, ss}
    """
    if slope == 'auto': slope = settings.SLOPE
    if calib == 'auto': calib = settings.CALIB

    # get data
    if verbose: print('read')
    if isinstance(scan, br.Image):
        im0  = scan
        ims0 = scan
    else:
        a = readrixs(scan=scan, start=start, stop=stop, verbose=verbose, folderpath=folderpath, prefix=prefix, filepath=filepath)
        if 'diff1' in a:  raise ValueError('process() is only for RIXS. This is XAS or linescan')
        im0 = a['im']
        ims0 = a['ims']

    # crop and normalization (number of images, image size)
    if verbose: print('crop')
    if x_start is not None or x_stop is not None or y_start is not None or y_stop is not None:
        length_before = im0.shape[1]
        im0 = im0.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
        im0 = im0.set_factor(1/(im0.shape[1]/length_before))
        for i, _im in enumerate(ims0):
            length_before = _im.shape[1]
            ims0[i] = _im.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
            ims0[i] = ims0[i].set_factor(1/(ims0[i].shape[1]/length_before))  # length after crop must be divided by lenght before croping so we can normalize by croped region

    # normalization
    if verbose: print('norm')
    im1  = im0.copy()
    if norm_i0:   
        try:
            im1 = im1.set_factor(1/np.average(im1.metadata['instrument']['i0']))
        except AttributeError:
            print(f'm4c1 does not seem to be recorded for scan {scan}. Cannot normalize by i0. If you want to get rid of this warn, set norm_i0=False')
        except TypeError:
            print(f'm4c1 does not seem to be recorded for scan {scan}. Cannot normalize by i0. If you want to get rid of this warn, set norm_i0=False')
    if norm_exposure:
        try:
            im1 = im1.set_factor(1/np.sum(im1.metadata['andor']['count_time']))
        except AttributeError:
            print(f'exposure_time does not seem to be recorded for scan {scan}. Cannot normalize by exposure_time. If you want to get rid of this warn, set norm_exposure=False')
    if norm_eslit:
        try:
            im1 = im1.set_factor(1/np.sum(im1.exit_slit))
        except AttributeError:
            print(f'exit_slit does not seem to be recorded for scan {scan}. Cannot normalize by exit_slit. If you want to get rid of this warn, set norm_eslit=False')
    ims1 = ims0.copy()
    for i, _im in enumerate(ims0):
        ims1[i] = ims0[i].copy()
        if norm_i0:
            ims1[i] = ims1[i].set_factor(1/np.sum(im1.metadata['instrument']['i0'][i]))
        if norm_exposure:
            ims1[i] = ims1[i].set_factor(1/np.sum(im1.metadata['andor']['count_time'][i]))
        if norm_eslit:
            ims1[i] = ims1[i].set_factor(1/np.sum(ims1[i].exit_slit))

    # get dark, crop, normalize (number of images, image size)
    # optional, normalize by I0, exposure time, exit slit
    if verbose: print('dark')
    if dark is not None:
        if isinstance(dark, br.Image):
            d0  = None
            d1 = dark.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
        else:
            temp = readrixs(scan=dark, folderpath=folderpath, prefix=prefix, filepath=filepath, verbose=verbose)
            d0 = temp['im']
            if x_start is not None or x_stop is not None or y_start is not None or y_stop is not None:
                length_before = d0.shape[1]
                d0 = d0.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
                d0 = d0.set_factor(1/(d0.shape[1]/length_before))  # normalize by size of image
            d1 = d0.copy()
            if norm_i0:
                try:
                    d1 = d1.set_factor(1/np.average(d0.metadata['instrument']['i0']))
                except TypeError:
                    # if i0 haven't been measured for dark image, use the i0 from scan
                    d1 = d1.set_factor(1/np.average(im1.metadata['instrument']['i0']))
                except AttributeError:
                    d1 = d1.set_factor(1/np.average(im1.metadata['instrument']['i0']))
            if norm_exposure:
                d1 = d1.set_factor(1/np.sum(d0.metadata['andor']['count_time']))
            if norm_eslit:
                d1 = d1.set_factor(1/np.sum(d0.exit_slit))
    else:
        d0  = None
        d1  = None

    # apply darkfactor and darkoffset (auto dark offset)
    if dark is not None:
        if darkoffset == 'auto':
            _temp     = im1.integrated_rows_vs_y_centers()
            _bkg      = _temp.crop(1500, 1800).calculate_y_average()
            _bkg_dark = d1.integrated_rows_vs_y_centers().crop(1500, 1800).calculate_y_average()
            darkoffset = (_bkg - _bkg_dark)/len(_temp.x)
        d1 = d1.set_factor(darkfactor).set_offset(darkoffset)
    
    # dark image subtraction
    if dark is not None:
        im2 = im1 - d1

        ims2 = ims1.copy()
        for i, _im in enumerate(ims1):
            ims2[i] = _im - d1
    else:
        im2  = im1.copy()
        ims2 = ims1.copy()

    # curvature correction
    if verbose: print('slope')
    if slope is not None:
        im3 = im2.set_vertical_shift_via_polyval(p=[-slope, 0])

        ims3 = ims2.copy()
        for i, _im in enumerate(ims2):
            ims3[i] = _im.set_vertical_shift_via_polyval(p=[-slope, 0])
    else:
        im3  = im2.copy()
        ims3 = ims2.copy()

    # pixel integration
    if verbose: print('integration')
    s = im3.integrated_rows_vs_y_centers()

    ss = br.Spectra()
    for i, _im in enumerate(ims3):
        ss.append(_im.integrated_rows_vs_y_centers())

    # energy calibration
    if verbose: print('calib')
    if calib is not None:
        if isinstance(calib, Iterable):
            s = s.set_shift(-calib[1])
            s = s.set_calib(calib[0])

            ss.set_shift(-calib[1])
            ss.set_calib(calib[0])
        else:
            s = s.set_calib(calib)
            ss.set_calib(calib)
        
    if verbose: print('done')
    return {'im0':im0, 'ims0':ims0, 
            'im1':im1, 'ims1':ims1,
            'd0': d0, 'd1':d1,
            'im2':im2, 'ims2':ims2,
            'im3':im3, 'ims3':ims3,
            's':s, 'ss':ss}

def _process_OLD2(scan, start=0, stop=None,  
            dark=None, darkfactor=1, darkoffset='auto', 
            slope='auto', calib='auto', norm_exposure=True,
            norm_i0=False, norm_eslit=False,  
            x_start=None, x_stop=None, y_start=None, y_stop=None, 
            verbose=False,
            folderpath='auto', prefix='auto', filepath='auto'):
    """Returns RIXS spectrum from detector image.

                                      d0* ---(crop, norm)---> d1
                                                              |
    ims0 ---(add ims)---> im0 ---(crop, norm)---> im1  ---(dark sub.)---> im2  ---(curv.)---> im3  ---(integration, calib)---> s
    ims0 -----------(crop, norm)----------------> ims1 ---(dark sub.)---> ims2 ---(curv.)---> ims3 ---(integration, calib)---> ss

    *dark image (d0) is assumed to be composed by only one image.

    Note: 
        Image intensity is normalized by number of images and 
        size of image (crop size). Images can also be normalized by I0, 
        exposure time, and exit slit if respective arguments are True.

    Args:
        scan (int): scan number.
        folderpath (str or path): folderpath where .nxs files are stored.
        dark (int or Image): dark image is loaded, crop, normalized and
            subtracted from scan. If dark is an br.Image object is
            subtracted (without normalization) directly from scan. Dark image is 
            also cropped with x_start, ..., y_stop and the cropping normalization
            is applied.
        darkfactor (number): dark image intensity will be multiplied by this 
            number.
        darkoffset (number): this value will be added to the dark image 
            intensity. Note that darkfactor is applied before dark offset, so
            one has to account for the multiplicative factor, i.e. 

                dark = dark.set_factor(darkfactor).set_offset(darkoffset)

            if darkoffset='auto', a suitable value will be found as to match
            the average intensity of pixel row integration between y=1500 and 
            y=1800 of the scan and dark image.
        slope (None or list): must be number or auto.
        calib (number, optional): if not None, the x axis of the final spectrum
            is multipled by calib. This number must have unit of eV/pixel.
        norm_i0 (True, optional): if True, image intensity is normalized by
            I0. Default is False.
        norm_exposure (True, optional): if True, image intensity is 
            normalized by image exposure. Default is True.
        norm_eslit (True, optional): if True, image intensity is normalized by
            the size of the exit slit. Default is False.
        x_start, x_stop, y_start, y_stop (int): pixel range in terms of
            x_centers and y_centers. Interval is inclusive. Use None to 
            indicate the edge of the image.
        verbose (bool, optional): if True, a message will print after each
            processing step, also a error message will be printed when 
            metadata cannot be retrieved. Default is False.
        
    Returns:
        dict {im0, ims0, im1, ims1, d0, d1, im2, ims2, im3, ims3, s, ss}
    """
    return _process_OLD1(scan=scan, start=start, stop=stop,  
                    dark=dark, darkfactor=darkfactor, darkoffset=darkoffset, 
                    slope=slope, calib=calib, norm_exposure=norm_exposure,
                    norm_i0=norm_i0, norm_eslit=norm_eslit, 
                    x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, 
                    verbose=verbose, folderpath=folderpath, prefix=prefix, filepath=filepath)['s']


# %% TENTATIVE: quick xas with variable span
# import matplotlib.lines as lines
# from matplotlib.widgets import MultiCursor
# from matplotlib.widgets import TextBox
def _quick_xas_tentative(scans, pre=(871, 872), post=(915, 918)):
    fig = br.figure(figsize=(28, 18), layout='constrained')
    fig.get_layout_engine().set(h_pad=0, w_pad=0)

    nrows = 4
    ncols = 4
    gs = fig.add_gridspec(nrows, ncols, width_ratios=[1, 1, 1, 1], height_ratios=[1, 1, 1, 1], hspace=0)
    axes = br.Axes([fig.add_subplot(gs[i, j]) for i in range(nrows) for j in range(ncols)], nrows=nrows, ncols=ncols)

    # column headers on ghost axes (preety, but expensive)
    # cols_headers = ('RAW', 'I0 Norm.')
    # axes1 = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 3:5])]
    # for i, ax in enumerate(axes1): 
    #     ax.axis('off')
    #     ax.set_title(cols_headers[i], fontsize=24, fontweight='bold',
    #             bbox=dict(facecolor='lightgrey', edgecolor='none', pad=-1))

    # column header directly on the figure (ugly, but cheap)
    text = plt.suptitle('.', fontsize=28, color='white')
    fig.text(0.275, 0.999, 'RAW',ha='center',va='top',fontsize=24, fontweight='bold',bbox=dict(facecolor='lightgrey', edgecolor='none', pad=2))
    fig.text(0.75, 0.999, 'I0 Norm.',ha='center',va='top',fontsize=24,fontweight='bold',bbox=dict(facecolor='lightgrey', edgecolor='none', pad=2),)

    # vertical line (resizable, but expensive)
    # axes2 = fig.add_subplot(gs[0:5, 2])
    # transform = mtransforms.blended_transform_factory(axes2.transAxes, fig.transFigure)
    # divider_line = fig.add_artist(lines.Line2D([0.5, 0.5], [0, 1], transform=transform, linewidth=3, color='lightgrey'))
    # divider_line.set_in_layout(False)
    # axes2.set_visible(False)

    # vertical line (frozen, but cheap)
    x = axes[1].get_position().x1 + 0.02
    divider_line = lines.Line2D([x, x], [0, 1], transform=fig.transFigure, linewidth=3, color='lightgrey')
    divider_line.set_in_layout(False)
    fig.add_artist(divider_line)

    # row headers
    row_headers = ('Diff1', 'FY2', 'TEY', 'I0')
    for i, ax in enumerate(axes.cols[0]):
        ax.set_ylabel(row_headers[i], fontsize=24, fontweight='bold',
                bbox=dict(facecolor='lightgrey', edgecolor='none'))
    fig.align_ylabels()

    # x label
    for ax in axes.last_row:
        ax.set_xlabel('Photon Energy (eV)')

    # setting up ticks
    for row in range(2):
        for ax in axes.rows[row]:
            ax.remove_xticklabels()
    axes[8].remove_xticklabels()
    axes[10].remove_xticklabels()
    for ax in axes:
        ax.tick_params(top=True, right=True, direction='in')

    # get data
    data = []
    for scan in scans:
        data.append(read(scan))

    # plot
    lines1 = []
    lines2 = []
    for i, d in enumerate(data):
        lines1.append({})
        lines2.append({})
        for row, prefix in enumerate(('diff1', 'fy2', 'tey', 'i0')):
            s = d[prefix]

            ax = axes[row * 4 + 0]
            s.plot(ax=ax)

            if prefix != 'i0':
                ax = axes[row * 4 + 1]
                if pre[0] == pre[1]:
                    avg = 0
                else:
                    avg = s.calculate_y_average(limits=(pre[0], pre[1]))
                lines1[i][prefix] = s.set_offset(-avg).plot(ax=ax)

            ax = axes[row * 4 + 2]
            s2 = (s/d['i0'])
            s2.plot(ax=ax)

            if prefix != 'i0':
                ax = axes[row * 4 + 3]
                if pre[0] == pre[1]:
                    avg = 0
                else:
                    avg = s2.calculate_y_average(limits=(pre[0], pre[1]))
                lines2[i][prefix] = s2.set_offset(-avg).plot(ax=ax)

    # pre-edge span lines
    spans = [0]*6
    zerolines = [0]*6
    if pre[0] != pre[1]:
        for j, ax in enumerate(axes.cols[1][:-1]):
            spans[j] = ax.axvspan(pre[0], pre[1], color='red', alpha=0.2)
            zerolines[j], = br.axhlines(0, ls='--', lw=0.4, color='black', ax=ax)

        for j, ax in enumerate(axes.cols[3][:-1]):
            spans[j+3] = ax.axvspan(pre[0], pre[1], color='red', alpha=0.2)
            zerolines[j+3], = br.axhlines(0, ls='--', lw=0.4, color='black', ax=ax)
    
    # share axes
    br.sharex(axes)



    def update_span(val):
        try:
            # Pull text contents, convert to floats
            xmin = float(editbox[0].text)
            xmax = float(editbox[1].text)

            if xmin == xmax:
                for i, d in enumerate(data):
                    for row, prefix in enumerate(('diff1', 'fy2', 'tey')):
                        s = d[prefix]
                        s2 = s/d['i0']
                        lines1[i][prefix].set_ydata(s.y)
                        lines2[i][prefix].set_ydata(s2.y)
                try:
                    for span in spans:
                        span.remove()
                        span = 0
                except: pass

                for line in zerolines:
                    line.set_visible(False)


                for row, prefix in enumerate(('diff1', 'fy2', 'tey')):
                    for ax in (axes[row * 4 + 1], axes[row * 4 + 3]):
                        ax.relim(visible_only=True)
                        ax.autoscale_view(scalex=False, scaley=True)
            else:
                if xmin >= xmax:
                    _temp = xmin
                    xmin = xmax
                    xmax = _temp
                for i, d in enumerate(data):
                    for row, prefix in enumerate(('diff1', 'fy2', 'tey')):
                        s = d[prefix]
        
                        ax = axes[row * 4 + 1]
                        avg = s.calculate_y_average(limits=(xmin, xmax))
                        lines1[i][prefix].set_ydata(s.set_offset(-avg).y)
                        
                        if i == 0:
                            try:
                                spans[row].remove()
                            except: pass
                            spans[row] = ax.axvspan(xmin, xmax, facecolor='red', alpha=0.2)
                            try:
                                zerolines[row].set_visible(True)
                            except:
                                zerolines[row],  = br.axhlines(0, ls='--', lw=0.4, color='black', ax=ax)
                        ax.relim(visible_only=True)
                        ax.autoscale_view(scalex=False, scaley=True)

                        ax = axes[row * 4 + 3]
                        s2 = s/d['i0']
                        avg = s2.calculate_y_average(limits=(xmin, xmax))
                        lines2[i][prefix].set_ydata(s2.set_offset(-avg).y)
                        ax.relim()
                        ax.autoscale_view(scalex=False, scaley=True)
                        if i == 0:
                            spans[row+3].remove()
                            spans[row+3] = ax.axvspan(xmin, xmax, facecolor='red', alpha=0.2)
                            try:
                                zerolines[row+3].set_visible(True)
                            except:
                                zerolines[row+3], = br.axhlines(0, ls='--', lw=0.4, color='black', ax=ax)
                        
            fig.canvas.draw()
        except ValueError:
            pass 

    
    # textboxes
    axes[13].axis('off')
    axes3 = []
    axes3.append(axes[13].inset_axes([0.4, 0.7, 0.40, 0.12]))
    axes3.append(axes[13].inset_axes([0.4, 0.5, 0.40, 0.12]))
    axes3.append(axes[13].inset_axes([0.4, 0.3, 0.40, 0.12]))
    axes3.append(axes[13].inset_axes([0.4, 0.1, 0.40, 0.12]))
    editbox = []
    editbox.append(TextBox(axes3[0], 'Pre-edge min: ', initial=str(pre[0])))
    editbox.append(TextBox(axes3[1], 'Pre-edge max: ', initial=str(pre[1])))
    editbox.append(TextBox(axes3[2], 'Post-edge min: ', initial=str(post[0])))
    editbox.append(TextBox(axes3[3], 'Post-edge max: ', initial=str(post[1])))
    editbox[0].on_submit(update_span)
    editbox[1].on_submit(update_span)
    editbox[2].on_submit(update_span)
    editbox[3].on_submit(update_span)
    fig.editbox = editbox

    # Turn off constrained layout because it is expensive
    fig.canvas.draw()
    fig.set_layout_engine('none')

    # crosshair
    # for ax in axes:
    #     ax.autoscale(False)
    # fig._crosshair = MultiCursor(None, axes[0:-3], color='red', 
    #                         linestyle='--', linewidth=1,
    #                         horizOn=False, vertOn=True, useblit=True)
    
    return