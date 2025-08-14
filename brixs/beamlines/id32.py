#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""base functions for ID32 beamline of ESRF"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br

# %% ------------------------- Special Imports ---------------------------- %% #
try:
    import h5py
    from silx.io.specfile import SpecFile
except:
    pass


# %% ============================= support ================================ %% #
def list_of_scans(filepath):
    """"Returns list of available scans"""    
    # list of scans
    try:
        with h5py.File(filepath, 'r') as f:
            # return br.remove_duplicates(list(np.sort([int(float(_)) for _ in list(f.keys()) if 'loopscan' not in f[_]['title'][()].decode()])))
            return np.array(list(np.sort([float(_) for _ in list(f.keys()) if 'loopscan' not in f[_]['title'][()].decode()])))
    except OSError:
        sf = SpecFile(str(filepath))
        # temp = br.remove_duplicates(sf.list())
        temp = np.array([float(_) for _ in sf.keys()])
        sf.close()
        return temp

# %% =============================== XAS ================================== %% #
def _readxas(filepath, scan, order=1, overwrite_moved_motor=False, overwrite_branch=False, verbose=True):
    with h5py.File(filepath, 'r') as f:
        # check scan number
        assert scan in [int(float(_)) for _ in list(f.keys())], f'scan {scan} not found. Available scans: {[int(float(_)) for _ in list(f.keys())]}'
        index = str(scan) + '.' + str(order)
        _orders = [int(_.split('.')[-1]) for _ in list(f.keys()) if _.startswith(str(scan) + '.')]
        assert index in f, f'order {order} not available for scan {scan}. Available orders are {_orders}'
        
        # check if RIXS
        command = f[index]['title'][()].decode()
        if 'loopscan' in command:
            raise ValueError(f'scan {scan} (order={order}) appears to be a RIXS scan')

        # find branch
        if overwrite_branch:
            branch = overwrite_branch
        else:
            if 'it_n' in f[index]["measurement"] or 'ifluo_n' in f[index]["measurement"]:
                branch = 'XMCD'
            elif 'dbig_n' in f[index]["measurement"] or 'sam_n' in f[index]["measurement"]:
                branch = 'RIXS'
            else:
                raise ValueError('branch could not be identified.')
                
        # branch detectors
        if branch == 'XMCD':
            metadata = {'E': 'energy', 'magnetic_field': 'magnet', 'slit': 'exbgap',
                    'th': 'srot', 'sample_y': 'sy', 'sample_z': 'sz'}
            tfy_detector = 'ifluo_n'
            tey_detector = 'it_n'
            i0_detector  = 'i0_n_xmcd'
            moving_motors = {'mtraj': 'energy_enc', 
                             'srot':     'srot', 
                             'sy':       'sy', 
                             'sz':       'sz'}
        elif branch == 'RIXS':
            metadata = {'E': 'energy', 'H': 'H', 'K': 'K', 'L': 'L',
                    'th': 'th', 'tth': 'tth', 'phi': 'phi', 'chi': 'chi', 'slit': 'exagap',
                    'motor_x': 'xsam', 'motor_y': 'ysam', 'motor_z': 'zsam',
                    'sample_x': 'beamx', 'sample_y': 'ysam', 'sample_z': 'beamz', }
            tfy_detector = 'dbig_n'
            tey_detector = 'sam_n'
            i0_detector  = 'mir_rixs'
            moving_motors = {'mtraj': 'energy_enc',
                             'th':     'th', 
                             'xsam':   'beamx', 
                             'ysam':   'beamy', 
                             'zsam':   'zsam'}
        else:
            raise ValueError(f'branch {branch} does not exist. Available branches: [RIXS, XMCD]')
        # print(branch)

        # find moving motor
        x = None
        if overwrite_moved_motor:
            x = f[index]['measurement'][overwrite_moved_motor][()]
        else:
            for _motor in moving_motors:
                # print(command)
                if _motor in command:
                    # print('gg')
                    try:
                        x = f[index]['measurement'][moving_motors[_motor]][()]
                        break
                    except KeyError:
                        pass
        assert x is not None, 'cannot find motor that moved during the scan'

        # spectra
        TEY = br.Spectrum(x=x, y=f[index]['measurement'][tey_detector][()])
        TFY = br.Spectrum(x=x, y=f[index]['measurement'][tfy_detector][()])
        I0  = br.Spectrum(x=x, y=f[index]['measurement'][i0_detector][()])
        TEY.type = 'TEY'
        TFY.type = 'TFY'
        I0.type  = 'I0'

        # metadata
        for _s in (TEY, TFY, I0):
            for attr in metadata:
                try:
                    _s.__setattr__(attr, f[index]['instrument']['positioners'][metadata[attr]][()])
                except KeyError:
                    if verbose: print(f'Cannot load metadata: {attr}')

            _s.end_reason = f[index]['end_reason'][()].decode()
            _s.start_time = f[index]['start_time'][()].decode()
            _s.end_time   = f[index]['end_time'][()].decode()
            _s.command    = command

            # polarization
            polv = f[index]['instrument']['positioners']['hu70ap'][()]
            if polv < -5:
                pol = 'CM'
            elif -5 <= polv < 5:
                pol = 'LH'
            elif 5 <= polv < 30:
                pol = 'CP'
            elif polv >= 31:
                pol = 'LV'  
            _s.pol = pol

            try:
                start_time = _s.start_time.split('T')[-1].split('+')[0].split('.')[0]
                if branch == 'XMCD':
                    _s.label = f'#{scan} {_s.type}, pol={_s.pol}, M={round(_s.magnetic_field, 1)} T, th={round(_s.th, 2)}, y={round(_s.sample_y, 4)}, z={round(_s.sample_z, 4)}, {start_time}'
                elif branch == 'RIXS':
                    _s.label = f'#{scan}, pol={_s.pol}, th={round(_s.th, 2)}, x={round(_s.sample_x, 4)}, y={round(_s.sample_y, 4)}, z={round(_s.sample_z, 4)}, {start_time}'
            except: pass



    return TEY, TFY, I0

# %% ============================== RIXS ================================== %% #
def available_data(filepath, scan):
    """Returns list of available data (columns)"""
    sf = SpecFile(str(filepath))
    assert str(scan)+'.1' in sf.keys(), f'scan {scan} not found. Available scans: {sf.list()}'
    index = sf.index(scan)

    temp = sf.labels(index)
    sf.close()
    
    return temp

def available_motors(filepath, scan):
    """Returns list of available motors (metadata)"""
    sf = SpecFile(str(filepath))
    assert str(scan)+'.1' in sf.keys(), f'scan {scan} not found. Available scans: {sf.list()}'
    index = sf.index(scan)

    temp = sf.motor_names(index)
    sf.close()
    
    return temp

def _read(filepath, scan, order=1, verbose=True):
    """read spec file and return spectrum"""

    # open spec file
    sf = SpecFile(str(filepath))

    # index
    _temp = br.remove_duplicates(sf.list())
    assert scan in _temp, f'scan {scan} not found. Available scans: {_temp}'
    _orders = [int(_.split('.')[-1]) for _ in list(sf.keys()) if _.startswith(str(scan) + '.')]
    assert order in _orders, f'order {order} not available for scan {scan}. Available orders are {_orders}'
    index = sf.index(scan, order)
        
    # spectrum
    # SPC_trend = scan.data_column_by_name("SPC aligned to trend")
    # SPC_shift = scan.data_column_by_name("SPC aligned to shifts")
    # s = br.Spectrum(x=scan.data_column_by_name("SPC aligned to trend"), 
    s = br.Spectrum(x=sf[index].data_column_by_name("Energy (auto)"), 
                    y=sf[index].data_column_by_name("SPC"))
    if 'mir_rixs' in sf.labels(index):
        s.i0 = br.Spectrum(x=sf[index].data_column_by_name("Energy (auto)"), 
                        y=sf[index].data_column_by_name("mir_rixs"))
        metadata = {'E': 'energy', 'H': 'H', 'K': 'K', 'L': 'L',
                'th': 'th', 'tth': 'tth', 'phi': 'phi', 'chi': 'chi', 'slit': 'exagap',
                'motor_x': 'xsam', 'motor_y': 'ysam', 'motor_z': 'zsam',
                'sample_x': 'beamx', 'sample_y': 'ysam', 'sample_z': 'beamz'}
    elif 'mir_xmcd' in sf.labels(index):
        s.i0 = br.Spectrum(x=sf[index].data_column_by_name("Energy (auto)"), 
                        y=sf[index].data_column_by_name("mir_xmcd"))
        metadata = {'E': 'energy', 'magnetic_field': 'magnet', 'slit': 'exbgap', 
                    'th': 'srot', 'sample_y': 'sy', 'sample_z': 'sz'}
    else:
        metadata = {}
        if verbose: print('Cannot load I0')
    s.pixel   = sf[index].data_column_by_name("Pixel")
    s.photons = br.Spectrum(x=sf[index].data_column_by_name("Energy (auto)"), 
                        y=sf[index].data_column_by_name("Photons"))

    s.acquisition_time = sf[index].data_column_by_name("Acquisition time")[0]
    header = [line.split('#C')[-1].strip() for line in sf.scan_header(index) if line.startswith('#C')]
    s.acquisition_time_per_image = float([line for line in header if line.startswith('Acq time')][0].split()[-2])


    # metadata
    for attr in metadata:
        try:
            s.__setattr__(attr, sf[index].motor_position_by_name(name=metadata[attr]))
        except:
            if verbose: print(f'scan {scan} cannot load attr: {attr}')
    
    # other
    s.scan = scan
    s.order = order
    s.command = sf.command(index)
    s.start_time = sf.date(index)

    # polarization
    polv = sf[index].motor_position_by_name(name="hu70ap")
    if polv < -5:
        pol = 'CM'
    elif -5 <= polv < 5:
        pol = 'LH'
    elif 5 <= polv < 30:
        pol = 'CP'
    elif polv >= 31:
        pol = 'LV'  
    s.pol = pol

    try:
        start_time = s.start_time.split('T')[-1].split('+')[0].split('.')[0]
        if 'mir_xmcd' in sf.labels(index):
            s.label = f'#{scan}, E={round(s.E, 2)} eV, pol={s.pol}, slit={s.pol} um, M={round(s.magnetic_field, 1)}, th={round(s.th, 2)}, y={round(s.sample_y, 4)}, z={round(s.sample_z, 4)}, {start_time}'
        elif 'mir_rixs' in sf.labels(index):
            s.label = f'#{scan}, E={round(s.E, 2)} eV, pol={s.pol}, slit={s.pol} um, th={round(s.th, 2)}, tth={int(round(s.tth))}, x={round(s.sample_x, 4)}, y={round(s.sample_y, 4)}, z={round(s.sample_z, 4)}, {start_time}'
    except: pass

    sf.close()
    return s

# %% ============================= FINAL ================================== %% #
def read(filepath, scan, order=1, overwrite_moved_motor=False, overwrite_branch=False, verbose=True):
    """Return RIXS or XAS spectrum

    Args:
        filepath (str or path): filepath to h5 or spec file
        scan (int): scan number
        order (int): not sure what this is yet. Usually, it is 1.
        overwrite_moved_motor (False or str): the function detects the moving
            motor automatically (like, mono, sample_x, ...), in case it fails,
            one can define it manually (motors inside 'measurement' in the h5 file).
        overwrite_branch (False or str): the function detects the branch automatically,
            in case it fails, one can define it manually. Options are: RIXS or XMCD.
        verbose (bool): if True, prints errors.

    Returns:
        for XAS (h5 file)
            TEY, TFY, I0 = read(folderpath, scan, order, overwrite_moved_motor=False, overwrite_branch=False)
        
        for RIXS (spec file)
            s = read(filepath, scan)
    """
    # check if file is h5
    is_h5 = False
    try:
        with h5py.File(filepath, 'r') as f:
            is_h5 = True
    except:
        pass

    # XAS
    if is_h5:
        return _readxas(filepath=filepath, scan=scan, order=order, overwrite_moved_motor=overwrite_moved_motor, overwrite_branch=overwrite_branch, verbose=verbose)
    # RIXS
    else:
        return _read(filepath=filepath, scan=scan, order=order, verbose=verbose)
# %%













#
