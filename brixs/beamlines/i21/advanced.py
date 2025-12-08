#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Advanced functions for I21 beamline



Help:
    import brixs.beamlines.i21 as i21

    # read files
    im, ims = i21.read(scan, folderpath)              # For rixs scans
    diff1, TEY, TFY, I0 = i21.read(scan, folderpath)  # For xas and line scans
    
    # return spectrum for rixs scan
    s = i21.process(scan, folderpath, ...)            
    
    # plot quick verification plot (returns dict. with data from all steps)
    steps = i21.verify(scan, folderpath, ...)         
    
    # plot quick verification window for checking dark image subtraction
    darkoffset = i21.verify_dark(scan, folderpath, ...)

    # plot quick verification window for checking curvature correction
    curv = i21.verify_curv(scan, folderpath, ...)

    # returns images from 2D mesh scans
    diff1, TEY, TFY, I0 = i21.mesh(scans, folderpath, motor_scanned='z')

    

Todo:
    * implement 'auto' for darkfactor in _process()
"""

# %% ========================== Standard Imports ========================= %% #
from collections.abc import Iterable
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import matplotlib

# %% =============================== brixs =============================== %% #
import brixs as br
from brixs.beamlines.i21.core import scanlist, _read, _read_xas, _read_line
# %%

# %% Read
def read(scan, folderpath, verbose=False):
    """return rixs image, xas, or linescan data from I21 (nxs file)
    
    Args:
        scan (int): scan number.
        folderpath (str or path): folderpath where .nxs files are stored.
        verbose (bool, optional): if True, print error message when metadata 
            cannot be retrieved. Default is False.
    
    Returns:
        If RIXS: im, ims (summed up final br.Image object and list with 
            individual images)
        If XAS or linescan: diff1, TEY, TFY, I0 (spectra for each detector)
    """
    filepath = Path(folderpath)/ ('i21-' + str(scan) + '.nxs')
    assert filepath.is_file(), f'cannot find scan {scan} inside folder' + '\n' + f'{folderpath}'

    try:  # try reading as RIXS scan
        im, ims = _read(filepath)
        im.scan = scan
        for _im in ims:
            _im.scan = scan
        return im, ims
    except KeyError:
        pass

    try:  # try reading as XAS scan
        diff1, TEY, TFY, I0 = _read_xas(filepath)
        for _s in (diff1, TEY, TFY, I0):
            _s.scan = scan
        return diff1, TEY, TFY, I0
    except KeyError:
        pass

    try:  # try reading as line scan
        diff1, TEY, TFY, I0 = _read_line(filepath)
        for _s in (diff1, TEY, TFY, I0):
            _s.scan = scan
        return diff1, TEY, TFY, I0
    except Exception as e:
        raise ValueError('Cannot read scan.\n' + e)

    return

# %% Process
def _process(scan, folderpath, dark=None, darkfactor=1, darkoffset=0, curv=None,
             norm_i0=False, norm_exposure=True, norm_eslit=False, calib=None, 
             x_start=None, x_stop=None, y_start=None, y_stop=None, verbose=False):
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
        
    Returns:
        dict {im0, ims0, im1, ims1, d0, d1, im2, ims2, im3, ims3, s, ss}
    """
    # get data
    if verbose: print('read')
    if isinstance(scan, br.Image):
        im0  = scan
        ims0 = scan
    else:
        a = read(scan=scan, folderpath=folderpath, verbose=verbose)
        if a[0].type == 'XAS':  raise ValueError('process() is only for RIXS. This is XAS scan')
        if a[0].type == 'line': raise ValueError('process() is only for RIXS. This is line scan')
        im0, ims0 = a

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
            im1 = im1.set_factor(1/np.average(im1.i0))
        except AttributeError:
            print(f'm4c1 does not seem to be recorded for scan {scan}. Cannot normalize by i0. If you want to get rid of this warn, set norm_i0=False')
        except TypeError:
            print(f'm4c1 does not seem to be recorded for scan {scan}. Cannot normalize by i0. If you want to get rid of this warn, set norm_i0=False')
    if norm_exposure:
        try:
            im1 = im1.set_factor(1/np.sum(im1.exposure_time))
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
            ims1[i] = ims1[i].set_factor(1/np.sum(ims1[i].i0[i]))
        if norm_exposure:
            ims1[i] = ims1[i].set_factor(1/np.sum(ims1[i].exposure_time))
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
            d0, _ = read(scan=dark, folderpath=folderpath, verbose=verbose)
            if x_start is not None or x_stop is not None or y_start is not None or y_stop is not None:
                length_before = d0.shape[1]
                d0 = d0.crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
                d0 = d0.set_factor(1/(d0.shape[1]/length_before))  # normalize by size of image
            d1 = d0.copy()
            if norm_i0:
                try:
                    d1 = d1.set_factor(1/np.average(d0.i0))
                except TypeError:
                    # if i0 haven't been measured for dark image, use the i0 from scan
                    d1 = d1.set_factor(1/np.average(im1.i0))
                except AttributeError:
                    d1 = d1.set_factor(1/np.average(im1.i0))
            if norm_exposure:
                d1 = d1.set_factor(1/np.sum(d0.exposure_time))
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
    if verbose: print('curv')
    if curv is not None:
        im3 = im2.set_vertical_shift_via_polyval(p=curv)

        ims3 = ims2.copy()
        for i, _im in enumerate(ims2):
            ims3[i] = _im.set_vertical_shift_via_polyval(p=curv)
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

def process(scan, folderpath, 
            dark=None, darkfactor=1, darkoffset=0, 
            curv=None, calib=None,
            norm_i0=False, norm_exposure=True, norm_eslit=False,  
            x_start=None, x_stop=None, y_start=None, y_stop=None, 
            verbose=False):
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
    return _process(scan=scan, folderpath=folderpath, dark=dark, darkfactor=darkfactor, darkoffset=darkoffset, curv=curv, norm_i0=norm_i0, norm_exposure=norm_exposure, norm_eslit=norm_eslit, calib=calib, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, verbose=verbose)['s']

# %% Verify
def verify(scan, folderpath, 
           dark=None, darkfactor=1, darkoffset=0, 
           curv=None, calib=None, 
           norm_i0=False, norm_exposure=True, norm_eslit=False, 
           x_start=None, x_stop=None, y_start=None, y_stop=None, 
           verbose=False):
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
    data = _process(scan=scan, folderpath=folderpath, dark=dark, 
                    darkfactor=darkfactor, darkoffset=darkoffset, curv=curv, 
                    norm_i0=norm_i0, norm_exposure=norm_exposure, 
                    norm_eslit=norm_eslit, calib=calib, 
                    x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop,
                    verbose=verbose)

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

# %% verify dark
def verify_dark(scan, folderpath, 
                dark=None, darkfactor=1, darkoffset='auto', 
                norm_i0=False, norm_exposure=True, norm_eslit=False, 
                x_start=None, x_stop=None, y_start=None, y_stop=None, 
                verbose=False):
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
    data = _process(scan=scan, folderpath=folderpath, dark=dark, darkfactor=1, darkoffset=0, curv=None, norm_i0=norm_i0, norm_exposure=norm_exposure, norm_eslit=norm_eslit, calib=None, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
    
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

# %% verify curvature correction
def verify_curv(scan, folderpath, popt=None, ncols=16, 
                nrows=None, dark=None, darkfactor=1, darkoffset=0, 
                norm_i0=False, norm_exposure=True, norm_eslit=False, 
                x_start=None, x_stop=None, y_start=None, y_stop=None, 
                deg=2, ):
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
    data = _process(scan=scan, folderpath=folderpath, 
                    dark=dark, darkfactor=darkfactor, darkoffset=darkoffset, 
                    norm_i0=norm_i0, norm_exposure=norm_exposure, norm_eslit=norm_eslit, 
                    x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, 
                    curv=None, calib=None)
    
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

# %% mesh scan
def mesh(scans, folderpath, motor_scanned='z'):
    """Returns images for diff1, TEY, TFY

    Warning:
        Scanned motor must be z, i.e., change y and scan z. Other motor 
        combinations can be implemented in the future.

    Note:
        linecuts can be accessed via 'ss' attribute, i.e., TEY.ss. For
        scanned motor z, each spectrum will be a z scan.

    Args:
        scans (list): ordered list with scan number
        folderpath (str or path): folderpath where .nxs files are stored.
        motor_scanned (str, optional): as for now, only 'z' is possible.
    
    Returns:
        diff1, TEY, TFY, I0
    """
    sss = [br.Spectra(), br.Spectra(), br.Spectra(), br.Spectra()]
    ims = [br.Image(), br.Image(), br.Image(), br.Image()]
    for scan in scans:
        _diff1, _TEY, _TFY, _I0 = read(scan=scan, folderpath=folderpath)
        sss[0].append(_diff1)
        sss[1].append(_TEY)
        sss[2].append(_TFY)
        sss[3].append(_I0)
    for i in range(4):
        try:
            ims[i] = sss[i].interp().stack_spectra_as_columns()
            ims[i].x_centers = [s.sample_y for s in sss[i]]
            ims[i].ss = sss[i]
        except TypeError:
            pass

    return ims
# %%