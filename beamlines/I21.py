#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""base functions for I21 beamline

Data is saved in hdf and nxs files, but we only need to read nxs files because 
they contain the data as well the metadata



"""

# %% ------------------------- Standard Imports --------------------------- %% #
from pathlib import Path
import numpy as np

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br

# %% ------------------------- Special Imports ---------------------------- %% #
try:
    import h5py
except:
    pass

# %% ================================ read ================================ %% #

# %% ============================== obsolete =============================== %% #
def _h52dict(f):
    """Converts h5 to dict"""
    data = {}
    if type(f) == h5py._hl.files.File:
        data = _h52dict(f['entry'])
    elif type(f) == h5py._hl.group.Group:
        for key in f.keys():
            data[key] = _h52dict(f[key])
    elif type(f) == h5py._hl.dataset.Dataset:
        if f.shape == ():
            return None
        else:
            return f[:]
    else:
        print(f'Cannot recognize type {type(f)}. Skipping {f.name}')
    return data

def _read(filepath):
    """Read files from I21 beamline at Diamond light source.

    Args:
        filepath (str or Path object): filepath.

    Returns:
        :py:class:`brixs.Image`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    filepath = Path(filepath)

    f = h5py.File(Path(filepath), 'r')
    f = _h52dict(f)

    data = np.zeros(f['andor']['data'].shape[1:])
    for temp in f['andor']['data']:
        data += temp
    _ = f['andor'].pop('data')
    _ = f['instrument']['andor'].pop('data')

    im = br.Image(data=data)
    im.nd = f

    return im

def _read1(filepath):
    """Read file from I21 beamline at Diamond light source and return a list.

    Args:
        filepath (str or Path object): filepath.

    Returns:
        list with :py:class:`brixs.Image`

    Last updated: 04/04/2022 by Carlos Galdino
    """
    filepath = Path(filepath)

    f = h5py.File(Path(filepath), 'r')
    f = _h52dict(f)

    data = []
    for temp in f['andor']['data']:
        temp = np.zeros(f['andor']['data'].shape[1:])
        data.append(br.Image(data=temp))

    _ = f['andor'].pop('data')
    _ = f['instrument']['andor'].pop('data')
    for im in data:
        im.nd = f

    return data

# %% HDF5 attributes list ------------------------------------------------------
"""
# script
f['entry']['current_script_name']  # COPYING
# scan
f['entry']['scan_command'][()]  # COPYING
f['entry']['scan_fields'][()]   # not needed
# user data
f['entry']['user01']['facility_user_id'][()]  # COPYING
f['entry']['user01']['name'][()]              # COPYING
# date and time
f['entry']['start_time'][()]  # COPYING
f['entry']['end_time'][()]    # COPYING
f['entry']['duration'][()]    # not needed
# detector
f['entry']['andor'].keys()
f['entry']['andor']['checkbeam'][()]  # COPYING
f['entry']['andor']['data'][()]       # DATA
f['entry']['andor']['ds'][()]         # not needed
# scan2
f['entry']['diamond_scan'].keys() 
f['entry']['diamond_scan']['scan_finished'][()]            # COPYING
f['entry']['diamond_scan']['scan_command'][()]             # repeated
f['entry']['diamond_scan']['point_start_times'][()]        # not needed
f['entry']['diamond_scan']['scan_dead_time'][()]           # not needed
f['entry']['diamond_scan']['scan_dead_time_percent'][()]   # not needed
f['entry']['diamond_scan']['scan_estimated_duration'][()]  # not needed
f['entry']['diamond_scan']['scan_fields'][()]              # repeated
f['entry']['diamond_scan']['scan_rank'][()]                # not needed
f['entry']['diamond_scan']['scan_shape'][()]               # repeated
f['entry']['diamond_scan']['start_time'][()]               # repeated
# whatever 1
f['entry']['m4c1'].keys()  # ['checkbeam', 'ds', 'm4c1']
f['entry']['m4c1']['checkbeam'][()]  # repeated
f['entry']['m4c1']['ds'][()]         # repeated
f['entry']['m4c1']['m4c1'][()]       # not needed???
# whatever 2
f['entry']['m4c1_gain'].keys()  # ['checkbeam', 'ds', 'gain']
f['entry']['m4c1_gain']['checkbeam'][()]  # repeated
f['entry']['m4c1_gain']['ds'][()]         # repeated
f['entry']['m4c1_gain']['gain'][()]       # not needed???
# whatever 3
f['entry']['sample'].keys()  # ['beam']
f['entry']['sample']['beam'].keys()
f['entry']['sample']['beam']['incident_energy'][()]           # COPYING
f['entry']['sample']['beam']['incident_polarization'][()]     # not needed
f['entry']['sample']['beam']['beamExtentScannableName'][()]              # not needed
f['entry']['sample']['beam']['distance'][()]                             # not needed
f['entry']['sample']['beam']['fluxScannableName'][()]                    # not needed
f['entry']['sample']['beam']['incidentBeamDivergenceScannableName'][()]  # not needed
# andor
andor = f['entry']['instrument']['andor']
andor['count_time'][()]    # COPYING    
andor['data'][()]          # repeated
andor['description'][()]   # not needed
andor['local_name'][()]    # not needed
andor['manufacturer'][()]  # not needed
andor['model'][()]         # not needed
# andor settings
andor['andor_settings'].keys()
andor['andor_settings']['accumulation_period'][()]  # not needed
andor['andor_settings']['adc_speed'][()]            # not needed
andor['andor_settings']['binning_x'][()]            # not needed
andor['andor_settings']['binning_y'][()]            # not needed
andor['andor_settings']['cooler_control'][()]       # not needed
andor['andor_settings']['cooler_status'][()]        # not needed
andor['andor_settings']['cooler_temperature'][()]   # not needed
andor['andor_settings']['em_ccd_gain'][()]          # not needed
andor['andor_settings']['image_size_x'][()]         # not needed
andor['andor_settings']['image_size_y'][()]         # not needed
andor['andor_settings']['preamp_gain'][()]          # not needed
andor['andor_settings']['region_size_x'][()]        # not needed
andor['andor_settings']['region_size_y'][()]        # not needed
andor['andor_settings']['region_start_x'][()]       # not needed
andor['andor_settings']['region_start_y'][()]       # not needed
andor['andor_settings']['sensor_size_x'][()]        # not needed
andor['andor_settings']['sensor_size_y'][()]        # not needed
andor['andor_settings']['shutter_ext_TTL'][()]      # not needed
andor['andor_settings']['shutter_mode'][()]         # not needed
andor['andor_settings']['temperature_actual'][()]   # not needed
andor['andor_settings']['vertical_shift_amplitude'][()]  # not needed
andor['andor_settings']['vertical_shift_speed'][()]      # not needed
# poliriser
f['entry']['instrument']['polariser'].keys()
f['entry']['instrument']['polariser']['gamma'][()]  # not needed
f['entry']['instrument']['polariser']['stick'][()]  # not needed
# beamline
f['entry']['instrument']['beamline'][()]  # COPYING
f['entry']['instrument']['name'][()]  # repeated
# ring
f['entry']['instrument']['source'].keys()
f['entry']['instrument']['source']['name'][()]
f['entry']['instrument']['source']['current'][()]
f['entry']['instrument']['source']['energy'][()]
f['entry']['instrument']['source']['probe'][()]
f['entry']['instrument']['source']['type'][()]
# check beam
f['entry']['instrument']['checkbeam'].keys()
f['entry']['instrument']['checkbeam']['checkrc_beamok'][()]          # not needed
f['entry']['instrument']['checkbeam']['checktopup_time_beamok'][()]  # not needed
f['entry']['instrument']['checkbeam']['name'][()]                    # not needed
# whatever
f['entry']['instrument']['ds'].keys()
f['entry']['instrument']['ds']['name'][()]
f['entry']['instrument']['ds']['soft_limit_max'][()]
f['entry']['instrument']['ds']['soft_limit_min'][()]
f['entry']['instrument']['ds']['value'][()]
# fast_shutter
f['entry']['instrument']['fast_shutter'].keys()
f['entry']['instrument']['fast_shutter']['x'][()]
# undulator
f['entry']['instrument']['id'].keys()
f['entry']['instrument']['id']['bottomInner'][()]
f['entry']['instrument']['id']['bottomOuter'][()]
f['entry']['instrument']['id']['enabled'][()]
f['entry']['instrument']['id']['gap'][()]
f['entry']['instrument']['id']['harmonic'][()]
f['entry']['instrument']['id']['mode'][()]
f['entry']['instrument']['id']['polarisation'][()]  # COPYING
f['entry']['instrument']['id']['rowPhase'][()]
f['entry']['instrument']['id']['taper'][()]
f['entry']['instrument']['id']['topInner'][()]
f['entry']['instrument']['id']['topOuter'][()]
f['entry']['instrument']['id']['type'][()]
# mirrors motors
f['entry']['instrument']['m1'].keys()  # 'feedback', 'fine_pitch', 'fpsetpoint', 'height', 'pitch', 'roll', 'x', 'yaw'
f['entry']['instrument']['m2'].keys()  # 'feedback', 'fine_pitch', 'fpsetpoint', 'height', 'pitch', 'roll', 'x', 'yaw'
f['entry']['instrument']['m4'].keys()  # 'femto1', 'femto2', 'longy', 'rx', 'ry', 'rz', 'x', 'y', 'z'
f['entry']['instrument']['m4c1'].keys()  # 'count_time', 'coupling', 'gain', 'local_name', 'm4c1', 'mode'
f['entry']['instrument']['m5'].keys()  # 'hqrx', 'hqry', 'hqrz', 'hqx', 'hqy', 'hqz', 'longy', 'lqrx', 'lqry', 'lqrz', 'lqx', 'lqy', 'lqz', 'tth'
f['entry']['instrument']['pgm'].keys()  # 'b2_shadow', 'cff', 'energy', 'grating_pitch', 'grating_select', 'mirror_pitch', 'mirror_select'
f['entry']['instrument']['sgm'].keys()  # 'grating_select', 'h', 'pitch', 'r1', 'roll', 'wedge_nearside', 'wedge_offside', 'x'
# slit motors (I think s5 is the exit slit)
f['entry']['instrument']['s1'].keys()  # 'x_gap', 'x_pos', 'y_gap', 'y_pos'
f['entry']['instrument']['s2'].keys()  # 'x_gap', 'x_pos', 'y_gap', 'y_pos'
f['entry']['instrument']['s3'].keys()  # 'x_gap', 'x_pos', 'y_gap', 'y_pos'
f['entry']['instrument']['s4'].keys()  # 'lower', 'nearside', 'offside', 'upper', 'x_gap', 'x_pos', 'y_gap', 'y_pos'
f['entry']['instrument']['s5'].keys()  # 'h_gap', 'hdso', 'sut', 'v1_gap', 'v2_gap', 'vdso1', 'vdso2'
f['entry']['instrument']['s6'].keys()  # 'x_gap', 'x_pos', 'y_gap', 'y_pos'
# temperature
f['entry']['instrument']['lakeshore336'].keys()
f['entry']['instrument']['lakeshore336']['cryostat'][()]
f['entry']['instrument']['lakeshore336']['demand'][()]        # COPYING
f['entry']['instrument']['lakeshore336']['heater'][()]
f['entry']['instrument']['lakeshore336']['heater_range'][()]
f['entry']['instrument']['lakeshore336']['sample'][()]        # COPYING
f['entry']['instrument']['lakeshore336']['shield'][()]
# arm
f['entry']['instrument']['spectrometer'].keys()
f['entry']['instrument']['spectrometer']['armtth'][()]     # COPYING
f['entry']['instrument']['spectrometer']['specgamma'][()]
f['entry']['instrument']['spectrometer']['spech'][()]
f['entry']['instrument']['spectrometer']['specl'][()]
# manipulator
f['entry']['instrument']['manipulator'].keys() 
f['entry']['instrument']['manipulator']['chi'][()]  # COPYING
f['entry']['instrument']['manipulator']['phi'][()]  # COPYING
f['entry']['instrument']['manipulator']['th'][()]   # COPYING
f['entry']['instrument']['manipulator']['x'][()]    # COPYING
f['entry']['instrument']['manipulator']['y'][()]    # COPYING
f['entry']['instrument']['manipulator']['z'][()]    # COPYING
f['entry']['instrument']['manipulator']['difftth'][()]
f['entry']['instrument']['manipulator']['draincurrent'][()]
f['entry']['instrument']['manipulator']['sapara'][()]
f['entry']['instrument']['manipulator']['saperp'][()]
# momentum transfer
f['entry']['Q'].keys()
f['entry']['Q']['H'][()]
"""