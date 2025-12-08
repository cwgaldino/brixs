#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from ADRESS beamline - PSI.

Last edited: Carlos Galdino 03-2023
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from datetime import datetime, timedelta
from collections.abc import Iterable
from pathlib import Path
import numpy as np
import copy

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br

# %% ------------------------- Special Imports ---------------------------- %% #
try:
    import h5py
except:
    pass

# %% ------------------------------ XAS ----------------------------------- %% #
def _unpack_attrs_xas(s):
    """Makes ADRESS relevant XAS scan parameters more readable."""

    # polarization
    if s.Polarization == '1':
        s.pol = 'LV'
    else:
        s.pol = 'LH'

    # sample position
    s.SampleX = float(s.Manipulator_X.split()[0])
    s.SampleY = float(s.Manipulator_Y.split()[0])
    s.SampleZ = float(s.Manipulator_Z.split()[0])

    s.T    = round(float(s.Cryostat_temperature[:-1]), 1)
    s.th   = round(float(s.Manipulator_Theta[:-3]), 1)
    s.slit = round(float(s.Exit_slit[:-3]), 1)
    s.label = f'{s.pol}, {s.T} K, th {s.th}, slit {s.slit}'

def _read_xas(filepath):
    """Read xas file from ADRESS beamline.

    Args:
        filepath (str or Path object): filepath.
    
    Returns:
        TEY, TFY, RMU (:py:class:`brixs.Spectrum`)
    
    Last updated: Carlos Galdino 07-2023
    """
    ##################
    # check filepath #
    ##################
    filepath = Path(filepath)
    assert filepath.exists(), 'filepath does not exists'
    assert filepath.is_file(), 'filepath must point to a file'

    #############
    # load data # 
    #############
    data = br.load_data(filepath, force_array=True)
    TEY = br.Spectrum(x=data[:, 0], y=data[:, 1])
    TFY = br.Spectrum(x=data[:, 0], y=data[:, 2])
    RMU = br.Spectrum(x=data[:, 0], y=data[:, 3])

    #################
    # load metadata #
    #################
    header = br.load_comments(filepath)[1:-7]
    nd = {}
    for line in header:
        if line.startswith('# Start-time'):
            split = line[1:].split(':')
            name  = split[0].strip()
            value = ':'.join(split[1:]).split('.')[0].strip()
            nd[name] = value
        else:
            split = line[1:].split(':')
            name  = split[0].strip()
            value = split[1].strip()
            nd[name] = value

    # save attr to object
    for attr in nd:
        attr2 = attr.replace('-', '_')
        setattr(TEY, attr2, nd[attr])
        setattr(TFY, attr2, nd[attr])
        setattr(RMU, attr2, nd[attr])
        
    # additional attrs
    for s in (TEY, TFY, RMU):
        _unpack_attrs_xas(s)
    TEY.mode = 'TEY'
    TFY.mode = 'TFY'
    RMU.mode = 'RMU'

    return TEY, TFY, RMU

def read_xas(*args, **kwargs):
    """Read xas file from ADRESS beamline.

    Usage:
        from brixs.file_reading import ADRESS

        TEY, TFY, RMU = ADRESS.read_xas(filepath=<filepath>)
        TEY, TFY, RMU = ADRESS.read_xas(<filepath>)

        TEY, TFY, RMU = ADRESS.read_xas(folderpath, prefix, scan)
        TEY, TFY, RMU = ADRESS.read_xas(folderpath=<folderpath>, prefix=<prefix>, scan=<scan>)

    Args:
        filepath (str or Path object): filepath to xas file.
        folderpath (str or Path object): folderpath.
        prefix (str): file prefix.
        scan (number): scan number. Example: for 'Cu_0001.xas' the file number is 1.

    Returns:
        TEY, TFY, RMU                 
    """
    ###################################
    # asserting validity of the input #
    ###################################
    error_message = 'Wrong input. Please, use one ' +\
                    'of the examples below to read xas files:\n' +\
                    '\n' +\
                    'TEY, TFY, RMU = ADRESS.read_xas(filepath=<filepath>)\n' +\
                    'TEY, TFY, RMU = ADRESS.read_xas(<filepath>)\n' +\
                    '\n' +\
                    'TEY, TFY, RMU = ADRESS.read_xas(folderpath, prefix, scan)\n' +\
                    'TEY, TFY, RMU = ADRESS.read_xas(folderpath=<folderpath>, prefix=<prefix>, scan=<scan>)\n' +\
                    '\n' +\
                    'filepath and folderpath must be a string or pathlib.Path object' +\
                    'prefix must be a string and scan must be an int'
    if kwargs != {} and args != ():
        raise AttributeError(error_message)
    if kwargs == {} and args == ():
        raise AttributeError(error_message)
    if any([item not in ['filepath', 'folderpath', 'prefix', 'scan'] for item in kwargs.keys()]):
        raise AttributeError(error_message)
    if len(args) > 3 or len(kwargs) > 3:
        raise AttributeError(error_message)
    if len(args) == 2 or len(kwargs) == 2:
        raise AttributeError(error_message)
    if 'folderpath' in kwargs and ('prefix' not in kwargs or 'scan' not in kwargs):
        raise AttributeError(error_message)
    if ('folderpath' in kwargs or 'prefix' in kwargs or 'scan' in kwargs) and 'filepath' in kwargs:
        raise AttributeError(error_message)
    
    ################
    # loading data #
    ################
    # keyword arguments
    if 'filepath' in kwargs:
        filepath = kwargs['filepath']
        return _read_xas(filepath)
    if 'folderpath' in kwargs:
        folderpath = Path(kwargs['folderpath'])
        prefix     = kwargs['prefix']
        scan       = kwargs['scan']

    # positional arguments
    if len(args) == 1:
        filepath = Path(args[0])
        return _read_xas(filepath)
    if len(args) == 3:
        folderpath = Path(args[0])
        prefix     = args[1]
        scan       = args[2]       
    
    # validation
    if folderpath.is_dir() == False:
        raise AttributeError(f'It seems like the folderpath is not a folder\nfolderpath = {folderpath}')
    if prefix is None:
        raise AttributeError('Missing prefix.\nPlease, inform the file prefix.\nExample: for Cu_0001.xas the prefix is `Cu_`')
    if scan is None:
        raise AttributeError('Missing file number (scan number).\nPlease, inform the file number.\nExample: for Cu_0001.xas the file number is 1.')
    if br.numanip.is_integer(scan) == False:
        raise TypeError('scan must be an integer')
    
    filepath = folderpath/(prefix+str(scan).zfill(4)+'.xas')
    return _read_xas(filepath)

# %% ------------------------------ base ---------------------------------- %% #
def _unpack_attrs(s):
    """Makes ADRESS relevant scan parameters more readable."""

    # parameters
    params = dict(pol     = 0,
                 phi      = 0,
                 T        = 0,
                 Tmin     = 0,
                 Tmax     = 0,
                 Tstd     = 0,
                 E        = 0,
                 th       = 0,
                 exposure = 0,
                 split    = 0,
                 slit     = 0)

    # polarization
    if s.PolarMode[0] == 1:
        params['pol'] = 'LV'
    else:
        params['pol'] = 'LH'
    
    # phi
    params['phi'] = round(np.mean(s.SamplePhi), 2)

    # Temperature
    temp = s.SampleTemp
    params['T']    = round(np.mean(temp), 2)
    params['Tmin'] = round(min(temp), 2)
    params['Tmax'] = round(max(temp), 2)
    params['Tstd'] = round(np.std(temp), 2)

    # Photon energy
    params['E']        = round(np.mean(s.PhotonEnergy), 2)

    # theta
    params['th']       = round(np.mean(s.SampleTheta), 2)

    # exposure
    params['exposure'] = round(s.AcquireTime[0], 2)

    # split
    params['split']    = round(s.ExposureSplit[0], 2)

    # slit
    params['slit']     = round(s.ExitSlit[0], 2)


    for name in params:
        setattr(s, name,  params[name])

def _read(filepath, type_='spectrum'):
    """Read one file from ADRESS beamline.

    Args:
        filepath (str or Path object): filepath.
        type_ (str, optional): data to return (case insensitive). Use 
            'spectrum' or 's' for spectrum; 'photon events' or 'pe', for photon
            events list; and 'bad' or 'b' for bad events image.  
        
    
    Returns:
        brixs.Spectrum() or brixs.PhotonEvents()
    
    Last updated: Carlos Galdino 03-2023
    """
    ###################################
    # asserting validity of the input #
    ###################################
    filepath = Path(filepath)
    if filepath.exists == False:
        raise FileNotFoundError(f'file does not exist\nfilepath: {filepath}')

    #############
    # read file #
    #############
    try:
        f = h5py.File(filepath, 'r')
    except NameError:
        raise ModuleNotFoundError('Cannot find h5py. Make sure h5py is installed.') 

    ################
    # get metadata #
    ################
    nd = {key: val[:] for key, val in f['entry']['instrument']['NDAttributes'].items()}
    nd['ModifiedDate'] = br.get_modified_date(filepath)  # modified date
    for key in nd:
        if isinstance(nd[key], np.ndarray):
            nd[key] = list(nd[key])

    ########
    # data # 
    ########
    # spectrum
    if type_.lower() in ['spectrum', 's']:
        # spectrum
        try:
            final = br.Spectrum(f['entry']['analysis']['spectrum'][:])
        except KeyError:
            raise ValueError(f'Cannot find spectrum inside file: {filepath}')
    # photon events
    elif type_.lower() in ['photon events', 'pe']:
        # get array size
        xmax = nd['ArraySizeX'][0]
        ymax = nd['ArraySizeY'][0]
        # get photon events list
        data = f['entry']['analysis']['events'][:]
        final = br.PhotonEvents(x=data[:, 0], y=data[:, 1])
        final.xlim = (0, xmax)
        final.ylim = (0, ymax)      
    # bad events
    elif type_.lower() in ['bad', 'b']:
        # bad events image
        final = br.Image(np.array(f['entry/analysis/bad'][:]))
    else:
        raise ValueError("type_ must be one of 'spectrum', 'photon events', or 'bad'.")

    # save attr to object
    for attr in nd:
        setattr(final, attr, nd[attr])
    _unpack_attrs(final)

    return final
    
def raw(*args, **kwargs):
    """Read raw RIXS files from ADRESS beamline at PSI.

    Usage:
        >>> s   = ADRESS.raw(filepath)
        >>> pe  = ADRESS.raw(filepath, type_='pe')
        >>> bad = ADRESS.raw(filepath, type_='bad')

        >>> data from all 3 ccd's
        >>> ss   = ADRESS.raw(folderpath, prefix, scan)
        >>> ss[0], ss[1], ss[2] will have data from each ccd
        >>> pes  = ADRESS.raw(folderpath, prefix, scan, type_='pe')
        >>> bads = ADRESS.raw(folderpath, prefix, scan, type_='bad')
    
    where pes and bads are lists, while ss is a type br.Spectra of length 3. 
        one can also use keyword arguments instead of positional arguments

    Args:
        filepath (str or Path object): filepath to h5 file (for loading a single file).
        folderpath (str or Path object): folderpath to files.
        prefix (str): file prefix.
        scan (number): scan number. Example: for 'Cu_0001_d1.h5' the file number is 1.
        type_ (str, optional): data to return (not case sensitive). Use 
            'spectrum' or 's' for spectrum; 'photon events' or 'pe', for photon
            events list; 'bad' or 'b' for bad events image; Default is spectrum.

    Returns:
        if ``type_ = 's'``: :py:class:`brixs.Spectrum` or :py:class:`brixs.Spectra`

        if ``type_ = 'pe'``: :py:class:`brixs.PhotonEvents` or a list of :py:class:`brixs.PhotonEvents`.

        if ``type_ = 'bad'``: :py:class:`brixs.Image` or a list of :py:class:`brixs.Image`

    Last updated: Carlos Galdino 07-2023
    """
    ###############
    # select type #
    ###############
    if 'type_' in kwargs: 
        type_ = kwargs.pop('type_')
    else:
        type_ = 's'

    ###################################
    # asserting validity of the input #
    ###################################
    error_message = 'Wrong input. Please, use one ' +\
                    'of the examples below to read RIXS files:\n' +\
                    '\n' +\
                    's   = ADRESS.raw(filepath)\n' +\
                    "pe  = ADRESS.raw(filepath, type_='pe')\n" +\
                    "bad = ADRESS.raw(filepath, type_='bad')\n" +\
                    '\n' +\
                    'ss   = ADRESS.raw(folderpath, prefix, scan)\n' +\
                    "pes  = ADRESS.raw(folderpath, prefix, scan, type_='pe')\n" +\
                    "bads = ADRESS.raw(folderpath, prefix, scan, type_='bad')\n" +\
                    '\n' +\
                    'filepath and folderpath must be a string or pathlib.Path object' +\
                    'prefix must be a string and scan must be an int'
    if kwargs != {} and args != ():
        raise AttributeError(error_message)
    if kwargs == {} and args == ():
        raise AttributeError(error_message)
    if any([item not in ['filepath', 'folderpath', 'prefix', 'scan'] for item in kwargs.keys()]):
        raise AttributeError(error_message)
    if len(args) > 3 or len(kwargs) > 3:
        raise AttributeError(error_message)
    if len(args) == 2 or len(kwargs) == 2:
        raise AttributeError(error_message)
    if 'folderpath' in kwargs and ('prefix' not in kwargs or 'scan' not in kwargs):
        raise AttributeError(error_message)
    if ('folderpath' in kwargs or 'prefix' in kwargs or 'scan' in kwargs) and 'filepath' in kwargs:
        raise AttributeError(error_message)
    
    ################
    # loading data #
    ################
    # keyword arguments
    if 'filepath' in kwargs:
        filepath = kwargs['filepath']
        return _read(filepath, type_=type_)
    if 'folderpath' in kwargs:
        folderpath = Path(kwargs['folderpath'])
        prefix     = kwargs['prefix']
        scan       = kwargs['scan']

    # positional arguments
    if len(args) == 1:
        filepath = Path(args[0])
        return _read(filepath, type_=type_)
    if len(args) == 3:
        folderpath = Path(args[0])
        prefix     = args[1]
        scan       = args[2]       
    
    # validation
    if folderpath.is_dir() == False:
        raise AttributeError(f'It seems like the folderpath is not a folder\nfolderpath = {folderpath}')
    if prefix is None:
        raise AttributeError('Missing prefix.\nPlease, inform the file prefix.\nExample: for Cu_0001.xas the prefix is `Cu_`')
    if scan is None:
        raise AttributeError('Missing file number (scan number).\nPlease, inform the file number.\nExample: for Cu_0001.xas the file number is 1.')
    if br.numanip.is_integer(scan) == False:
        raise TypeError('scan must be an integer')
    
    # folderpath, prefix, scan
    filepaths = [folderpath/(prefix+str(scan).zfill(4)+f'_d{i}.h5') for i in (1, 2, 3)]

    if type_.lower() in ['spectrum', 's']:
        ss  = br.Spectra()
        for i, filepath in enumerate(filepaths):
            _s = _read(filepath)
            _s.scan = scan
            _s.ccd = i
            ss.append(_s)
        ss.scan = scan
        # for attr in ss[0]._get_user_attrs():
        #     setattr(ss, attr, getattr(ss[0], attr))
        # del ss.ccd
        return ss
    elif type_.lower() in ['photon events', 'pe']:
        pes = [0]*3
        for i, filepath in enumerate(filepaths):
            pes[i] = _read(filepath, type_=type_)
            pes[i].scan = scan
            pes[i].ccd = i
        return pes
    elif type_.lower() in ['bad', 'b']:
        bads = [0]*3
        for i, filepath in enumerate(filepaths):
            bads[i] = _read(filepath, type_=type_)
            bads[i].scan = scan
            bads[i].ccd = i
        return bads
    else:
        raise ValueError("type_ must be one of 'spectrum', 'photon events', or 'bad'.")

def read(folderpath, prefix, scan, nbins=None, curvature=None, calib=None, offset=0):
    """Returns RIXS final spectrum from ADRESS beamline at PSI.

    Data from all ccd's are aligned and summed up. If calib is not None, data
        is also calibrated. If nbins is not None, data is rebined. If more than 
        one scan is passed as input, scans are summed up.

    Usage:
        >>> s = ADRESS.read(folderpath, prefix, scan)

    Args:
        folderpath (str or Path object): folderpath to files.
        prefix (str): file prefix.
        scan (number or list): scan number or list with scan numbers.
        nbins (int, optional): if not None, photon events will be rebined to 
            nbins. Default is None. Curvature must be defined.
        curvature (list, optional): list with length 3. Each item should be the
            polynomial coeff. for curvature correction of the corresponding ccd.
            Default is None. Required if nbins is not None.
        calib (number, optional): If not None, calibration factor is
            applied (multiplicative factor). Default is None.
        offset (number, optional): calibration offset. If calib are polynomial
            coeff., offset is added to the const. term. If calib is a number 
            (multiplicative factor), offset is a hard shift on the x-coord. This
            redefine the zero energy position of the detector.
    
    Note:
        If rebining, calib factor must be in terms of subpixels and not bins.

        Use the calibration factor of the first ccd, because the multiplicative 
        term is the same for all ccd's and the cdd's spectra are aligned to the 
        first ccd, so the additive term should be the one for the first ccd.

        Offset is mostly used when you move the detector mid experiment and 
            you don't want to run recalibration again.

    Returns:
        Spectrum
    """
    ##################
    # add many scans #
    ##################
    if isinstance(scan, Iterable):
        ss = br.Spectra()
        for j, s in enumerate(scan):
            if nbins is not None:
                ss1 = br.Spectra()
                if curvature is None:
                    raise ValueError('bining requires curvature')
                pes = raw(folderpath=folderpath, prefix=prefix, scan=s, type_='pe')
                for i, pe in enumerate(pes):
                    pe = pe.set_vertical_shift_via_polyval(curvature[i])
                    ss1.append(pe.calculate_spectrum(nbins=nbins))

                # #####################################
                # # transfer attrs from the first ccd #
                # #####################################
                # for attr in pes[0]._get_user_attrs():
                #     ss1.__setattr__(attr, pes[0].__getattribute__(attr))

            else:
                ss1 = raw(folderpath=folderpath, prefix=prefix, scan=s)

            ##########################
            # align and sum each ccd #
            ##########################
            ss1.align()
            s = ss1.calculate_sum()
            ss.append(s)

            #######################################################
            # transfer attrs from first ccd to the Spectra object #
            #######################################################
            if j == 0:
                for attr in ss1[0].get_attrs():
                    ss.__setattr__(attr, [ss1[0].__getattribute__(attr), ])
            else:
                for attr in ss1[0].get_attrs():
                    ss.__setattr__(attr, ss.__getattribute__(attr) + [ss1[0].__getattribute__(attr), ])
        del ss.ccd

        #######################
        # align and sum scans #
        #######################
        ss.align()
        s = ss.calculate_sum()

        #####################################
        # transfer attrs to summed spectrum #
        #####################################
        for attr in ss.get_attrs():
            s.__setattr__(attr, ss.__getattribute__(attr))
        # print(s.E)
        s.scan = scan

    ############
    # one scan #
    ############
    else:
        if nbins is not None:
            ss = br.Spectra()
            if curvature is None:
                raise ValueError('new bining requires curvature')
            pes = raw(folderpath=folderpath, prefix=prefix, scan=scan, type_='pe')
            for i, pe in enumerate(pes):
                pe = pe.set_vertical_shift_via_polyval(curvature[i])
                ss.append(pe.calculate_spectrum(nbins=nbins))
        else:
            ss = raw(folderpath=folderpath, prefix=prefix, scan=scan)

        ##########################
        # align and sum each ccd #
        ##########################
        s = ss.fix_monotonicity().interp().align().calculate_sum()

        #######################################################
        # transfer attrs from the first ccd to final spectrum #
        #######################################################
        for attr in ss[0].get_attrs():
            s.__setattr__(attr, ss[0].__dict__[attr])
        del s.ccd
        s.scan = scan

    ######################
    # calibrate and zero #
    ######################
    if calib is not None:
        s = s.set_calib(calib)
        # calib = copy.deepcopy(calib)
        # if isinstance(calib, Iterable):
        #     if offset is not None:
        #         calib[-1] += offset
        #     s.set_calib(calib)
        #     if isinstance(scan, Iterable):
        #         s.set_shift(-np.mean(s.E))
        #     else:
        #         s.set_shift(-s.E)
        # else:
        #     s.set_calib(calib)
        #     if offset is not None:
        #         s.set_shift(offset)
    
    #########
    # label #
    #########
    label = (f'scans{s.scan}, '+
             f'{s.T}({s.Tstd})K, '+
             f'{s.E}eV, '+
             f'{s.pol}, '+
             f'phi{s.phi}, '+
             f'th{s.th}, '+
             f'{s.exposure}s, '+
             #f'split={self.split}, '+
             f'slit{s.slit}um'
            )
    s.label = label
    return s

def sequence(folderpath, prefix, scans, nbins=None, curvature=None, calib=None, offset=0):
    """Returns a spectra object with RIXS final spectrum from ADRESS beamline.

    Data from all ccd's are aligned and summed up. If calib is not None, data
        is also calibrated. If nbins is not None, data is rebined. If more than 
        one scan is passed as input, scans are summed up.

    Usage:
        >>> ss = ADRESS.sequence(folderpath, prefix, scans)

        where scans can be 

        >>> scans = (10, 11, 13, 14, 18, 30)

        in this case, the returned object will have 6 spectra (one for each scan).


        >>> scans = ((10, 11), (13, 14), (15, 16))

        In this case, the returned object will have 3 spectra.
        10 and 11 will be added together, so 13 and 14, and finally
        15 and 16.

    Args:
        folderpath (str or Path object): folderpath to files.
        prefix (str): file prefix.
        scan (list): list with scan numbers.
        nbins (int, optional): if not None, photon events will be rebined to 
            nbins. Default is None. Curvature must be defined.
        curvature (list, optional): list with length 3. Each item should be the
            polynomial coeff. for curvature correction of the corresponding ccd.
            Default is None. Required if nbins is not None.
        calib (number or list, optional): If not None, calibration factor is
            applied. Can be a multiplicative factor (number) or a list with
            polynomial coeff. with highest power first. Default is None.
        offset (number, optional): calibration offset. If calib are polynomial
            coeff., offset is added to the const. term. If calib is a number 
            (multiplicative factor), offset is a hard shift on the x-coord. This
            redefine the zero energy position of the detector.
    
    Note:
        If rebining, calib factor must be in terms of subpixels and not bins.

        Use the calibration factor of the first ccd, because the multiplicative 
        term is the same for all ccd's and the cdd's spectra are aligned to the 
        first ccd, so the additive term should be the one for the first ccd.

        Offset is mostly used when you move the detector mid experiment and 
            you don't want to run recalibration again.

    Returns:
        Spectra
    """    
    ss = br.Spectra()
    for i, scan in enumerate(scans):
        ss.append(read(folderpath=folderpath, prefix=prefix, scan=scan, nbins=nbins, curvature=curvature, calib=calib, offset=offset))

    #########
    # attrs #
    #########
    ss.scans    = scans    
    ss.T        = [s.T for s in ss]
    ss.E        = [s.E for s in ss]
    ss.pol      = [s.pol for s in ss]
    ss.phi      = [s.phi for s in ss]
    ss.th       = [s.th for s in ss]
    ss.exposure = [s.exposure for s in ss]
    ss.split    = [s.split for s in ss]
    ss.slit     = [s.slit for s in ss]

    return ss

def calib(folderpath, prefix, start=None, stop=None, scans=None, nbins=None, curvature=None, mode='cc'):
    """Returns the calibration factor (dispersion) of ADRESS.

    if Spectra is given as a function of bin (x=bin, y=Intensity), then
        calib factor comes in units of eV/bin. If Spectra is given in terms of
        detector coordinates (x=um, y=Intensity), then calib = eV/um.

    popt gives the polynomial fitting coeff. of the dispersion curve and can be 
        used to calibrate data absolutely.

    Usage:
        from brixs.file_reading import ADRESS
        popt, sss = ADRESS.calib(folderpath, prefix, start=17, stop=26)

    Args:
        folderpath (str or Path object): folderpath.
        prefix (str): file prefix.
        start (number, optional): first scan number. Default is None.
        stop (number, optional): last scan number. Default is None.
        scans (list, optional): list of scan numbers. Overwrites start and
            stop. Default is None.
        mode (string, optional): method used to calculate the shifts. See 
            br.Spectra.calculate_shift() for mode options.
        nbins (int or tuple): number of bins for calculate spectra.
        curvature (): 1D array of polynomial coefficients 
                (including coefficients equal to zero) from highest degree to 
                the constant term. Polynomial as a function of x_centers (um).

    Returns:
        2 lists:
                    
            polynomial coefficients that fit the dispersion curve. 
            
            Spectra for each ccd.
    """
    # scan numbers
    if scans is None:
        scans = np.arange(start, stop+1)

    # get data
    if nbins is not None:
        sss = [br.Spectra(), br.Spectra(), br.Spectra()]
        for i, scan in enumerate(scans):
            pes = raw(folderpath=folderpath, prefix=prefix, scan=scan, type_='pe')
            for ccd in (0, 1, 2):
                pes[ccd].set_shift(p=curvature[ccd], axis=0)
                s = pes[ccd].calculate_spectrum(nbins=nbins)
                sss[ccd].append(s)
    else:
        sss = [br.Spectra(), br.Spectra(), br.Spectra()]
        for i, scan in enumerate(scans):
            temp = raw(folderpath=folderpath, prefix=prefix, scan=scan)
            for ccd in (0, 1, 2):
                sss[ccd][i].append(temp[ccd])

    # get energies
    energies = [np.mean(s.E) for s in sss[0]]
    
    # save energies to spectra
    for ccd in (0, 1, 2):
        sss[ccd].E = energies

    # calculate calib
    popt = [0, 0, 0]   
    for ccd in (0, 1, 2):
        # calculate shifts
        if mode == 'peaks' or mode == 'peak':
            for ccd2 in (0, 1, 2):
                sss[ccd2].fit_peak()
        sss[ccd].calculate_calib(values=energies, mode=mode, deg=1)

        # get peak of first spectrum
        sss[ccd][0].fit_peak()
        c = sss[ccd][0].peaks[0]['c']
        s = br.Spectrum(x=-sss[ccd].calculated_shift+c, y=energies)
        popt[ccd], model, r2 = s.polyfit(deg=1)        

    return popt, sss

# %% -------------------------- Archiver functions ------------------------ %% #
def str2datetime(string):
    date, time = string.split()
    year, month, day = date.split('-')
    hour, minute, second = time.split(':')
    return datetime(day=int(day), month=int(month), year=int(year), hour=int(hour), minute=int(minute), second=int(second.split('.')[0]))

def read_archiver(filepath):

    # open file
    filepath = Path(filepath)
    with filepath.open("r") as f:
        txt = f.read().split('\n')[:-1]
        
    # create new dict with nearest function
    class MyDict(dict):
        def nearest(self, day, month, year=2023, hour=0, minute=0, second=0):
            pivot = datetime(year, month, day, hour, minute)
            return self['time'].index(min(self['time'], key=lambda x: abs(x - pivot)))
        
        def get_range(self, start=None, stop=None, elapsed=0):
            """elapsed in seconds"""
            assert start is not None or stop is not None, 'start and stop cannot both be None'

            if start is None:
                if type(stop) == str:
                    stop = str2datetime(stop)
                start = stop - timedelta(seconds=elapsed)
            elif stop is None:
                if type(start) == str:
                    start = str2datetime(start)
                stop = start+timedelta(seconds=elapsed)
            else:
                if type(start) == str:
                    start = str2datetime(start)
                if type(stop) == str:
                    stop = str2datetime(stop)
            start = self.nearest(day=start.day, month=start.month, year=start.year, hour=start.hour, minute=start.minute, second=start.second)
            stop  = self.nearest(day=stop.day, month=stop.month, year=stop.year, hour=stop.hour, minute=stop.minute, second=stop.second)
            return {'time': self['time'][start:stop],
                    'RMU': self['RMU'][start:stop],
                    'TFY': self['TFY'][start:stop],
                    'TEY': self['TEY'][start:stop]}

    # sort columns
    data = MyDict(time=[], TEY=[], TFY=[], RMU=[])
    for line in txt:
        if line.startswith('#') or line.startswith('\n') or '#N/A' in line:
            pass
        else:
            temp = line.split('\t')
            if len(temp) == 4:
                data['time'].append(datetime.strptime(temp[0].split('.')[0], '%m/%d/%Y %H:%M:%S'))
                data['TEY'].append(float(temp[1]))
                data['TFY'].append(float(temp[2]))
                data['RMU'].append(float(temp[3]))
        
    return data

# %% -------------------------- high level functions ---------------------- %% #
def _calculate_energy_map(self):
    """Calculates energy map.

    Energy map is saved as an attr:

        ss.energy_map    
    
    Returns:
        None
    """
    ss = self.copy()

    if hasattr(self, 'E') == False:
        try:
            self.E = [s.E for s in self]
        except AttributeError:
            raise AttributeError('spectra does not have ss.E or s.E (energy) attributes.')

    try:
        ss.check_same_x()
    except ValueError:
        ss.interp()

    self.energy_map = ss.calculate_map()
    self.energy_map.x_centers = self.E
br.Spectra.calculate_energy_map = _calculate_energy_map

def _calculate_emission_spectra(self):
    """Calculates emission spectra.

    Emission spectra are saved as an attr:

        ss.emission    
    
    Returns:
        None
    """
    ss = self.copy()

    # x axis from energy loss to emission energy
    if hasattr(self, 'E') == False:
        try:
            for s in ss:
                s.x = s.x + s.E
        except AttributeError:
            raise AttributeError('spectra does not have s.E (energy) attribute.')
    else:
        for i, s in enumerate(ss):
            s.x = s.x + ss.E[i]

    self.emission = ss
br.Spectra.calculate_emission_spectra = _calculate_emission_spectra

def _calculate_emission_map(self):
    """Calculates emission map.

    Emission spectra must be calculated first.

    Emission map are saved as an attr:

        ss.emission_map    
    
    Returns:
        None
    """
    try:
        ss = self.emission.copy()
    except AttributeError:
        raise AttributeError('Must first calculate emission spectra (ss.calculate_emission_spectra())')

    if hasattr(self, 'E') == False:
        try:
            self.E = [s.E for s in self]
        except AttributeError:
            raise AttributeError('spectra does not have s.E (energy) attribute.')

    try:
        ss.check_same_x()
    except ValueError:
        ss.interp()

    self.emission_map = ss.calculate_map(axis=1)
    self.emission_map.y_centers = self.E
    self.emission_map.x_centers = ss.x 
br.Spectra.calculate_emission_map = _calculate_emission_map

def _calculate_th_map(self):
    """Calculates th map.

    Th map is saved as an attr:

        ss.th_map    
    
    Returns:
        None
    """
    ss = self.copy()

    if hasattr(self, 'th') == False:
        try:
            self.th = [s.th for s in self]
        except AttributeError:
            raise AttributeError('spectra does not have ss.th or s.th (theta) attributes.')

    try:
        ss.check_same_x()
    except ValueError:
        ss.interp()

    self.th_map = ss.calculate_map(axis=0)
    self.th_map.x_centers = ss.th 
br.Spectra.calculate_th_map = _calculate_th_map

# ------------------------------- EXPERIMENTAL ---------------------------- %% #
def th2rlu(theta, energy, propagation_vector=(1, 0, 0), lattice_parameters=(3.5066, 4.7485, 7.9341), twotheta=130):
    """[EXPERIMENTAL]"""
    q, qpar, qper = br.momentum_transfer(energy, twotheta-theta, twotheta)
    rlu = br.momentum2rlu(qpar, propagation_vector, *lattice_parameters)
    return rlu

def _calculate_momentum_map(self, energy=None, propagation_vector=None, lattice_parameters=None):
    """Calculates momentum map."""
    ss = self.copy()

    try:
        ss.check_same_x()
    except ValueError:
        ss.interp()

    # th2rlu kwargs
    kwargs = {}
    if propagation_vector is not None:
        kwargs['propagation_vector'] = propagation_vector
    if lattice_parameters is not None:
        kwargs['lattice_parameters'] = lattice_parameters
    if energy is None:
        rlu = [th2rlu(th=s.th, energy=s.energy, **kwargs) for s in ss]
    else:
        kwargs['energy'] = energy
        rlu = [th2rlu(th=s.th, **kwargs) for s in ss]

    self.momentum_map = ss.calculate_map(axis=0)
    self.momentum_map.x_centers = rlu
br.Spectra.calculate_momentum_map = _calculate_momentum_map