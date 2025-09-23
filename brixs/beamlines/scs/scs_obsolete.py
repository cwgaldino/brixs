#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for reading files from SCS beamline - XFEL.

Last edited: Carlos Galdino 2024-Sept-03
"""

# %% ------------------------- Standard Imports --------------------------- %% #
from collections.abc import Iterable
from IPython.display import Audio
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import matplotlib
import datetime
import warnings
import copy

# %% ------------------------------ brixs --------------------------------- %% #
import brixs as br
import brixs.addons.centroid

# %% ------------------------- Special Imports ---------------------------- %% #
import toolbox_scs as tb
from toolbox_scs.detectors.hrixs import hRIXS
# %%

# %% ============================= settings =============================== %% #
settings = {'default_proposal': None,
            'has_multiple_energy_max_error': 1.0,
            'has_multiple_delays_max_error': 0.001}
# %%

# %% ========================= useful functions =========================== %% #
def make_sound(freq=200, duration=0.5):
    beep = np.sin(2*np.pi*freq*np.arange(duration*5000*2)/10000)
    return Audio(beep, rate=10000, autoplay=True)

def scanlist(proposal):
    """Return list of scans available for a proposal"""
    # folderpath = Path(folderpath)
    # assert folderpath.exists(), f'fpath does not exist ({folderpath})'
    # return br.parsed_filelist(folderpath, string='*', ref=0, return_type='dict')
    raise NotImplementedError('not implemented')

def verify_params(scan, params=None, proposal=None):
    """plot params with raw signal from a run to check if params is suitable

    params for extracting area of the peaks in signal
    peaks are generated with the x-ray pulses
    params = {'pulseStart': 2176, # start of the first pulse
              'pulseStop':  2195, # end of the first pulse
              'baseStart':  2120, # start of bkg reference
              'baseStop':   2170, # end of bkg reference 
              'period':     None, # 96 for 1.1 MHz (does not matter if bunchPattern is used)
              'npulses':    None} # number of pulses (does not matter if bunchPattern is used)
    
    Args:
        scan (int): scan number.
        params (dict, optional): if None, params will be calculated.
        proposal (int, optional): proposal number. If None, proposal number will be read
            from SCS.settings['default_proposal']. Default is None
    
    Returns
        params used in the plotting
    """
    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    run, ds = _get_dataset_xas(scan, proposal=proposal)

    return tb.check_peak_params(run, 'FastADC2_9raw', params=params)

def has_multiple_energy(ds, max_error=None):
    """Returns True if the photon energy setpoint varies
    
    Returns true if: absolute sum of the photon energy setpoint 
            variation is more then the max_error times the average
            photon energy setpoint
        
    Args:
        ds (xarray.Dataset): dataset
        max_error (number, optional): percentage (of the x step) value 
            of the max error. Default is 0.1 %.
    
    Returns:
        bool                
    """
    if max_error is None:
        max_error = settings['has_multiple_energy_max_error']
        
    if sum(abs(np.diff(ds['nrj_target']))) > ds['nrj_target'].mean()*max_error/100:
        return True
    return False

def _digitize_via_binning(arr, bins, xmin=None, xmax=None):
    """Returns dictionary with bin values (keys) and indexes.
    
    Args:
        arr (list): list with values.
        bins (int, str, list): 
        
            If bins is a list, it is assumed to be the bins edges. 
            
            Bins in inclusive at the bin_edge with smaller value, 
            and excusive at the bin_edge with higher value.
        
            (from np.histogram_bin_edges):
            If bins is an int, it defines the number of equal-width 
            bins in the given range (10, by default). If bins is a sequence, it defines the 
            bin edges, including the rightmost edge, allowing for non-uniform bin widths.

            If bins is a string from the list below, histogram_bin_edges will use the method 
            chosen to calculate the optimal bin width and consequently the number of bins 
            (see the Notes section for more detail on the estimators) from the data that 
            falls within the requested range. While the bin width will be optimal for the 
            actual data in the range, the number of bins will be computed to fill the entire 
            range, including the empty portions. For visualisation, using the ‘auto’ option 
            is suggested. Weighted data is not supported for automated bin size selection.

            ‘auto’
            Minimum bin width between the ‘sturges’ and ‘fd’ estimators. Provides good 
            all-around performance.

            ‘fd’ (Freedman Diaconis Estimator)
            Robust (resilient to outliers) estimator that takes into account data variability 
            and data size.

            ‘doane’
            An improved version of Sturges’ estimator that works better with non-normal datasets.

            ‘scott’
            Less robust estimator that takes into account data variability and data size.

            ‘stone’
            Estimator based on leave-one-out cross-validation estimate of the integrated squared error. 
            Can be regarded as a generalization of Scott’s rule.

            ‘rice’
            Estimator does not take variability into account, only data size. Commonly overestimates 
            number of bins required.

            ‘sturges’
            R’s default method, only accounts for data size. Only optimal for gaussian data and 
            underestimates number of bins for large non-gaussian datasets.

            ‘sqrt’
            Square root (of data size) estimator, used by Excel and other programs for its speed and simplicity.
        xmin, xmax (int, optional): minimum/maximum value for bins (only used if bins are not explicitly 
            given. Default is None. If None, min and max array values are used. 
            
        Returns:
            indexes (dict), bin_edges (list)
    """
    if isinstance(bins, Iterable)==False or isinstance(bins, str):
        bin_edges = np.histogram_bin_edges(arr, bins=20, range=None, weights=None)
    else:
        bin_edges = bins
    temp = np.digitize(arr, bins=bin_edges, right=False)
    
    indexes = {}
    for i in range(len(bin_edges[:-1])):
        # indexes[(bin_edges[i]+bin_edges[i+1])/2] = [_ for _ in np.argwhere(temp==i+1)]
        indexes[(bin_edges[i]+bin_edges[i+1])/2] = [int(_[0]) for _ in list(np.argwhere(temp==i+1))]
        
    return indexes, bin_edges

def _digitize_via_duplicates(arr, max_error=None):
    """Returns dict with mean values (keys) and indexes:

    Args:
        arr (list): list with values.
        max_error (number or None): percentage value of max allowed error between 
            numbers for them to be considered duplicates. if None, numbers are 
            considerate the same only if they are real duplicates.

    Returns:
        indexes (dict)
    """
    indexes = {arr[0]: [0, ]}
    values  = {arr[0]: [arr[0], ]}

    for i, item in enumerate(arr[1:]):
        new = True
        for _item in values:

            average = np.mean(values[_item])
            _min = average - average*max_error/100
            _max = average + average*max_error/100

            # this deals with negative numbers
            if _min > _max:
                _min, _max = (_max, _min)

            if item >= _min and item <= _max:
                indexes[_item].append(i+1)
                values[_item].append(item)
                new = False
                break
        if new:
            indexes[item] = [i+1, ]
            values[item]  = [item, ]  
    return {np.mean(values[key]): indexes[key] for key in values}


# def find_duplicate_delays(scan, max_error=None, proposal=None):

#     #######################
#     # get proposal number #
#     #######################
#     if proposal is None:
#         proposal = settings['default_proposal']
#     assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

#     ####################
#     # get dataset info #
#     ####################
#     run, ds = tb.load(proposal, scan, fields=['PP800_DelayLine'])
#     delay_line = ds['PP800_DelayLine'].values + 0

#     #############
#     # max error #
#     #############
#     if max_error is None:
#         max_error = settings['has_multiple_delay_max_error']

#     return _digitize_via_duplicates(delay_line, max_error=max_error)

# %% =========================== time zero ================================ %% #
# this will force every spectrum object to have the following attrs: 
# t0_ and t0
# is_time_trace
# When t0 changes, s.delay_line is calculated based on s.time_delay if s.time_delay exists
# t0_ and time_delay must be defined for all spectra


# br.Spectrum.is_time_trace = False
# br.Spectrum.t0_ = None
# def getter(self):
#     if hasattr(self, 't0_'):
#         return self.t0_
# def setter(self, value):
#     if hasattr(self, 'time_delay'):
#         # print([float(_) for _ in self.time_delay])
#         # print(tb.delayToPosition([int(_) for _ in self.time_delay], origin=int(self.t0_)))
#         # print(self.t0_)
#         # self.t0_ = 300.0
#         # if self.t0_ is not None:
#         if isinstance(self.time_delay, Iterable):
#             delay_line      = [tb.delayToPosition(_, origin=self.t0_) for _ in self.time_delay]
#             self.time_delay = [tb.positionToDelay(_, origin=value) for _ in delay_line]
#         else:
#             delay_line      = tb.delayToPosition(self.time_delay, origin=self.t0_)
#             self.time_delay = tb.positionToDelay(delay_line, origin=value)
#         if hasattr(self, 'is_time_trace'):
#             if self.is_time_trace:
#                 self.x = self.time_delay
#     self.t0_ = value
#     return
# setattr(br.Spectrum, 't0', property(getter, setter, ))

# def getter(self):
#     if hasattr(self, 'time_delay') and hasattr(self, 't0_'):
#         if isinstance(self.time_delay, Iterable):
#             return [tb.delayToPosition(_, origin=self.t0_) for _ in self.time_delay]
#         else:
#             return tb.delayToPosition(self.time_delay, self.t0_)
# def setter(self, value):
#     raise ValueError('cannot set delay_line value. Change time_delay instead or t0')
# setattr(br.Spectrum, 'delay_line', property(getter, setter, ))

def change_t0(self, t0):
    delayline  = tb.delayToPosition(np.array(self.time_delay), origin=self.t0)
    time_delay = tb.positionToDelay(np.array(delayline), origin=t0)
    if self.is_time_trace:
        self.x = time_delay
    self.time_delay = time_delay
    self.t0 = t0
    return
br.Spectrum.change_t0 = change_t0
# %%

# %% ============================= metadata =============================== %% #
_attrs  = []

_attrs += ['SCS_SA3']                  # ionization chamber after SASE3 monochomator
_attrs += ['SCS_photonFlux']           # ionization chamber before SASE3 monochomator
_attrs += ['nrj']                      # photon energy
_attrs += ['nrj_target']               # energy setpoint (necessary to check if energy is supposed to be changing)
_attrs += ['PP800_DelayLine']          # delay line motor 

_attrs += ['VSLIT', 'HSLIT', 'ESLIT']  # vertical, horizontal, and exit slits
_attrs += ['transmission']             # attenuators (laser)? 
_attrs += ['transmission_setpoint']    # attenuators (laser)? 
_attrs += ['transmission_col2']        # attenuators
_attrs += ['GATT_pressure']            # gas attenuator. After sase3 undulator and slits, before M1
_attrs += ['XRD_STX', 'XRD_STY', 'XRD_STZ']  # manipulator x, y and z positions
_attrs += ['XRD_DRY', 'XRD_SRX', 'XRD_SRY', 'XRD_SRZ']
_attrs += ['XRD_SXT1Y', 'XRD_SXT2Y', 'XRD_SXTX', 'XRD_SXTZ']
# %%

# %% ============================== XAS (Basic) =========================== %% #
def _get_dataset_xas(scan, proposal=None):
    """return run and dataset for a scan (rixs related fields are not loaded)

    Args:
        scan (int): scan number
        proposal (int, optional): proposal number. If None, proposal number will be read
            from SCS.settings['default_proposal']. Default is None

    Returns:
        run, ds
    """
    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    ############
    # get data #
    ############
    # data fields
    fields = _attrs + ['FastADC2_9raw'] 
    
    # load dataset
    run, ds = tb.load(proposal, scan, fields, extract_fadc2=False, extract_fadc=False)
    
    return run, ds

def _dataset2dict_xas(ds2, step=0.1, is_time_trace='auto', raw=False, t0=300, bins_overwrite=None):#, debug=False):
    """return spectrum from dataset

    Args:
        ds2 (dataset)
        step 
        is_time_trace_ (bool, optional): 
        debug (bool): 

    Return:
        spectrum
    """
    ######################################
    # check if delay scan or energy scan # 
    ######################################
    if is_time_trace == 'auto':
        if has_multiple_energy(ds=ds2):
            is_time_trace = False
        else:
            is_time_trace = True
    
    ################
    # get metadata # 
    ################
    metadata = {}

    metadata['raw']  = {attr: ds2[attr].values + 0 for attr in _attrs}
    metadata['stat'] = {attr: {'avg': np.mean(metadata['raw'][attr]),
                               'min': np.min(metadata['raw'][attr]),
                               'max': np.max(metadata['raw'][attr]),
                               'std': np.std(metadata['raw'][attr])} for attr in metadata['raw']}
    if raw == False: del metadata['raw']

    metadata['E']          = round(metadata['stat']['nrj']['avg'], 2)
    metadata['sample_z']   = round(metadata['stat']['XRD_STZ']['avg'], 5)
    metadata['sample_y']   = round(metadata['stat']['XRD_STY']['avg'], 5)
    metadata['sample_x']   = round(metadata['stat']['XRD_STX']['avg'], 5)
    metadata['GATT']       = metadata['stat']['GATT_pressure']['avg']
    metadata['exit_slit']  = round(metadata['stat']['ESLIT']['avg'], 2)

    ##################################################################
    # creating dummy entry for normalized and non-normalized spectra #
    ##################################################################
    ds2['i0_sa3']  = ds2['SCS_SA3']
    ds2['i0_none'] = (('trainId', 'sa3_pId'), np.zeros(ds2['FastADC2_9peaks'].shape) + 1)
    # ds2['i0_photonflux'] = (('trainId', 'sa3_pId'), [[_]*ds2['FastADC2_9peaks'].shape[1] for _ in ds2['SCS_photonFlux'].values])
    
    # params=None, 
    # bunchPattern='sase3'
    ##################################################
    # Calculate delay values as a function of pulses #
    ##################################################
    # delay values
    # ds2['delay_line'] = ds2['PP800_DelayLine']
    # ds2['delay_line'] = tb.delayToPosition(ds2['time_delay'], origin=300) # delay in motor position
    
    ####################
    # creating spectra #
    ####################
    if is_time_trace:
        ds2['time_delay'] = tb.positionToDelay(ds2['PP800_DelayLine'], origin=t0)   # delay in ps
        if step is None:
            step = 0.1
        bins = int((np.amax(ds2['time_delay'])-np.amin(ds2['time_delay']))/step)
        nrjkey = 'time_delay'
    else:
        if step is None:
            step = 0.1
        bins = int((np.amax(ds2['nrj'])-np.amin(ds2['nrj']))/step)
        nrjkey = 'nrj'    

    if bins_overwrite is not None:
        bins = bins_overwrite
    
    # normalized by SCS_SA3
    temp = tb.xas(ds2, bins, Iokey='i0_sa3', Itkey='FastADC2_9peaks', nrjkey=nrjkey, fluorescence=True, plot=False)
    raw  = tb.xas(ds2, bins, Iokey='i0_none', Itkey='FastADC2_9peaks', nrjkey=nrjkey, fluorescence=True, plot=False)
    TFY     = br.Spectrum(x=temp['nrj'], y=temp['muA'])
    TFY.error = br.Spectrum(x=temp['nrj'], y=temp['sterrA'])
    TFY.sigma = br.Spectrum(x=temp['nrj'], y=temp['sigmaA'])
    TFY.i0  = br.Spectrum(x=temp['nrj'], y=temp['muIo'])
    TFY.not_norm = br.Spectrum(x=raw['nrj'],  y=raw['muIo'])

    del ds2['i0_sa3']
    del ds2['i0_none']

    # metadata for final spectrum
    for attr in metadata:
        TFY.__setattr__(attr, metadata[attr])
    TFY.t0 = 300
    if is_time_trace:
        TFY.is_time_trace = True
        TFY.time_delay    = TFY.x
    else:
        TFY.is_time_trace = False
        # TFY.time_delay    = tb.positionToDelay(ds2['PP800_DelayLine'], origin=300)
        # TFY.time_delay    = tb.positionToDelay(ds2['PP800_DelayLine'], origin=300)
        TFY.time_delay    = tb.positionToDelay(metadata['stat']['PP800_DelayLine']['avg'], origin=t0)
        # TFY.delay_line = tb.delayToPosition(TFY.x, origin=300)
    
    return {'ds2':ds2, 'TFY': TFY}

# @br.finder.track
# def _process_xas(scan, onoff=False, params=None, step=0.1, t0=300, proposal=None, is_time_trace='auto', bunchPattern='sase3', raw=False, debug=False):
#     d = _internal_xas(scan=scan, onoff=onoff, params=params, step=step, proposal=proposal, is_time_trace=is_time_trace, bunchPattern=bunchPattern, raw=raw, debug=debug)
#     s = d['TFY']

#     ##########
#     # set t0 #
#     ##########
#     s.t0 = t0

#     return s

# def process_xas(scan, onoff=False, params=None, step=0.1, t0=300, is_time_trace='auto', bunchPattern='sase3', raw=False, proposal=None):
#     #######################
#     # get proposal number #
#     #######################
#     if proposal is None:
#         proposal = settings['default_proposal']
#     assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

#     ###############
#     # get dataset #
#     ###############
#     ds = _get_dataset_xas(scan=scan, proposal=proposal)
#     delay_line = ds['PP800_DelayLine'].values + 0

#     ################################
#     # Does it have multiple delays #
#     ################################
#     if is_tr == 'auto':
#         is_tr = False
#         if br.has_duplicates(delay_line, max_error=settings['has_multiple_delays_max_error']):
#             is_tr = True

#     #####################################
#     # sort multiple delays if necessary #
#     #####################################
#     if is_tr:
#         indexes    = _digitize_via_duplicates(delay_line, max_error=settings['has_multiple_delays_max_error'])
#         delay_line_list = np.sort(list(indexes.keys()))[::-1]

#         ss = br.Spectra()
#         for delay in delay_line_list:
#             _ds = ds.isel({'trainId': indexes[delay]})
#             d = _dataset2dict(_ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)
#             s = d['s']
#             ss.append(s)

#     else:
#         d = _dataset2dict(ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)
#         s = d['s']
#         ss = br.Spectra([s, ])
    
#     return ss.data


    
#     s = _process_xas(scan=scan, onoff=onoff, params=params, step=step, proposal=proposal, is_time_trace=is_time_trace, bunchPattern=bunchPattern, raw=raw, debug=debug)

#     ##########
#     # set t0 #
#     ##########
#     s.t0 = t0

#     return s

def _internal_xas(scan, onoff=False, params=None, step=0.1, t0=300, proposal=None, is_time_trace='auto', bunchPattern='sase3', raw=False, debug=False, bins_overwrite=None):
    """Process XAS and time-trace from 1 run
    
    We are only counting the signal from the detectors if the signal matches in time with 
    a pulse. 
    
    TODO: as for right now, there is no way of integrating the whole signal without 
    discriminating pulses.
    
    TODO: as for right now, there is also no possible way of select which pulses are going 
    to be accounted for.
    
    SCS_SA3        -> ionization chamber after SASE3 monochomator  (function of pulse)
    SCS_photonFlux -> ionization chamber before SASE3 monochomator (function of train)
    
    # params for extracting area of the peaks in signal
    # peaks are generated with the x-ray pulses
    # params = {'pulseStart': 2176, # start of the first pulse
    #           'pulseStop':  2195, # end of the first pulse
    #           'baseStart':  2120, # start of bkg reference
    #           'baseStop':   2170, # end of bkg reference 
    #           'period':     None, # 96 for 1.1 MHz (does not matter if bunchPattern is used)
    #           'npulses':    None} # number of pulses (does not matter if bunchPattern is used)
    
    Note:
        If bunchPattern is used, params['period'] and params['npulses'] are ignored.
    
    Note:
        If bunchPattern is None, and params is given, then params['npulses'] must match the 
        real number. So, it does not matter what number I put in params['npulses']. 
        If it is not right, the function tb.xas() will raise an error
    
    Args:
        scan
        params
        E_step
        proposal
        bunchPattern
    
    Returns:
        TFY1, TFY2, RAW
        TFY1 normalized by SCS_SA3
        TFY2 normalized by SCS_photonFlux
        RAW not normalized
        
    """
    # start time (for debuging only)
    if debug: print('start')
    if debug: start_time = br.start_time()

    
    # load dataset
    run, ds = _get_dataset_xas(scan=scan, proposal=proposal)
    if debug: print('read data: ', br.stop_time_seconds(start_time))
    
    #####################################
    # calculate area of detector signal #
    #####################################
    # only takes into account the signal that matches in time with a x-ray pulse

    # fix bunchPattern argument
    if bunchPattern is None or bunchPattern == 'none':
        bunchPattern = 'None'
    
    # check if params['npulses'] is defined, if not, set it to the right number
    if params is not None:
        if params['npulses'] is None:
            params['npulses'] = ds['SCS_SA3'].shape[1]
    
    # calculation
    _params = copy.deepcopy(params)  # save params (this is necessary because get_digitizer_peaks changes params (facepalm)
    ds2 = tb.get_digitizer_peaks(run, 'FastADC2_9raw', merge_with=ds, integParams=_params, bunchPattern=bunchPattern)
    # ds2['FastADC2_raw'] = ds['FastADC2_raw']
    if debug: print('calc areas: ', br.stop_time_seconds(start_time))
    
    #################################################
    # Calculate delay values as a function of pulses #
    #################################################
    if onoff:
        pumped   = ds2.isel({'sa3_pId': slice(0, None, 2)})
        unpumped = ds2.isel({'sa3_pId': slice(1, None, 2)})
        on   = _dataset2dict_xas(pumped, step=step, is_time_trace=is_time_trace, raw=raw, t0=t0, bins_overwrite=bins_overwrite)
        off  = _dataset2dict_xas(unpumped, step=step, is_time_trace=is_time_trace, raw=raw, t0=t0, bins_overwrite=bins_overwrite)
        ds2 = on['ds2']
        TFY = on['TFY']
        TFY.off = off['TFY']
        
        if debug: print('final: ', br.stop_time_seconds(start_time))
        return {'ds':ds, 'ds2':ds2, 'TFY': TFY}
    
    else:
        d   = _dataset2dict_xas(ds2, step=step, is_time_trace=is_time_trace, raw=raw, t0=t0, bins_overwrite=bins_overwrite)
        ds2 = d['ds2']
        TFY = d['TFY']

        ##################################
        # average signal from all trains #
        ##################################
        # for now, this is only used for plotting
        ds2['FastADC2_9raw_avg'] = ds['FastADC2_9raw'].mean(dim='trainId')

    ############
    # metadata #
    ############
    TFY.scan = scan
            
    if debug: print('final: ', br.stop_time_seconds(start_time))
    return {'ds':ds, 'ds2':ds2, 'TFY': TFY}

# @br.finder.track
# def _process_xas(scan, onoff=False, params=None, step=0.1, proposal=None, is_time_trace='auto', bunchPattern='sase3', raw=False, debug=False):
#     """same as process but without t0"""

#     # # finder
#     # s = br.finder.search(parameters=locals(), folderpath=br.finder.folderpath)
#     # if s is not None:
#     #     return s

#     d = _internal_xas(scan=scan, onoff=onoff, params=params, step=step, proposal=proposal, is_time_trace=is_time_trace, bunchPattern=bunchPattern, raw=raw, debug=debug)
#     s = d['TFY']

#     return s

def process_xas(scan, onoff=False, params=None, step=0.1, t0=300, proposal=None, is_time_trace='auto', bunchPattern='sase3', raw=False, debug=False, bins_overwrite=None):
    # s = _process_xas(scan=scan, onoff=onoff, params=params, step=step, proposal=proposal, is_time_trace=is_time_trace, bunchPattern=bunchPattern, raw=raw, debug=debug)
    
    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    # finder
    _vars = locals()
    # _vars.pop('t0')
    s = br.finder.search(kwargs=_vars, folderpath=br.finder.folderpath)
    if s == False:
        d = _internal_xas(scan=scan, onoff=onoff, params=params, step=step, t0=t0, proposal=proposal, is_time_trace=is_time_trace, bunchPattern=bunchPattern, raw=raw, debug=debug, bins_overwrite=bins_overwrite)
        s = d['TFY']
        s.check_nan()
        if s.has_nan:
            s.remove_nan()
        br.finder.save(obj=s, folderpath=br.finder.folderpath)

    ##########
    # set t0 #
    ##########
    # s.change_t0(t0=t0)

    return s

# %% ============================ XAS (Advanced) ========================== %% #
def sequence_xas(scans, onoff=False, params=None, step=0.1, t0=300, proposal=None, is_time_trace='auto', bunchPattern='sase3', raw=False, debug=False):
    ss = br.Spectra()
    for scan in scans:
        ss.append(process_xas(scan=scan, onoff=onoff, params=params, step=step, t0=t0, proposal=proposal, is_time_trace=is_time_trace, bunchPattern=bunchPattern, raw=raw, debug=debug))
    
    for attr in ss[0].get_attrs():
        try:
            ss.create_attr_from_spectra(attr)
        except:
            pass

    return ss
# %%

# %% ============================= RIXS (Basic) =========================== %% #
def _get_dataset(scan, proposal=None):

    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    ############
    # get data #
    ############
    h  = hRIXS(proposal)
    ds = h.from_run(scan, extra_fields=_attrs + ['hRIXS_exposure', 'hRIXS_delay'])  
    return ds

def _dataset2dict(ds, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, norm=True, include_metadata_for_each_image=False):
    
    ################
    # get metadata # 
    ################
    metadata = {}

    # required metadata
    metadata['number_of_images'] = ds.dims['trainId']  # one image per train
    metadata['exposure']         = ds['hRIXS_exposure'].values + 0
    metadata['total_exposure']   = ds['hRIXS_exposure'].sum().values + 0
    metadata['i0']               = ds['SCS_photonFlux'].values + 0

    # metadata per train
    metadata['raw']  = {attr: ds[attr].values + 0 for attr in _attrs}
    metadata['stat'] = {attr: {'avg': np.mean(metadata['raw'][attr]),
                               'min': np.min(metadata['raw'][attr]),
                               'max': np.max(metadata['raw'][attr]),
                               'std': np.std(metadata['raw'][attr])} for attr in metadata['raw']}
    if include_metadata_for_each_image == False: del metadata['raw']

    # manual assignment
    metadata['E']          = round(metadata['stat']['nrj']['avg'], 2)
    metadata['sample_z']   = round(metadata['stat']['XRD_STZ']['avg'], 5)
    metadata['sample_y']   = round(metadata['stat']['XRD_STY']['avg'], 5)
    metadata['sample_x']   = round(metadata['stat']['XRD_STX']['avg'], 5)
    metadata['GATT']       = metadata['stat']['GATT_pressure']['avg']
    metadata['exit_slit']  = round(metadata['stat']['ESLIT']['avg'], 2)

    ################
    # brixs object #
    ################
    ims = br.Dummy(data=[br.Image(data=_.values).crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop) for _ in ds['hRIXS_det']])    
    
    # metadata for individual images
    if include_metadata_for_each_image:
        for attr in metadata['raw']:
            for i, _im in enumerate(ims):
                ims[i].__setattr__(attr, metadata['raw'][attr][i])
    
    # # scan
    # ims.scan = scan
    # for i, _im in enumerate(ims):
    #     ims[i].scan = scan
    
    # time
    ims.t0 = 300
    ims.time_delay = tb.positionToDelay(metadata['stat']['PP800_DelayLine']['avg'], origin=300)

    ########################
    # curvature correction #
    ########################
    if curv is not None:
        for i, _im in enumerate(ims):
            ims[i] = _im.set_horizontal_shift_via_polyval(p=curv)




    
    ####################
    # aggregate images #
    ####################
    im = br.Image(data=np.zeros(ims[0].shape))
    for i, _im in enumerate(ims):
        if norm:
            factor = 1/metadata['i0'][i]
        else:
            factor = 1
        ims[i]   = _im.set_factor(factor)
        im.data += _im.set_factor(factor).data
    im.x_centers = ims[0].x_centers
    im.y_centers = ims[0].y_centers
    # else:
    #     h     = hRIXS(proposal)
    #     data2 = h.aggregate(ds=ds)
    #     im    = br.Image(data=data2['hRIXS_det'].values).crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
    #     if curv is not None:
    #         im = im.set_horizontal_shift_via_polyval(p=curv)

    #################################
    # metadata for aggregated image #
    #################################
    for attr in metadata:
        im.__setattr__(attr, metadata[attr])
    
    # # scan
    # im.scan = scan 

    # time
    im.t0 = 300
    im.time_delay = ims.time_delay

    ##################
    # final spectrum #
    ##################
    if norm:
        s = im.integrated_columns_vs_x_centers().set_factor(1/metadata['total_exposure'])
    else:
        s = im.integrated_columns_vs_x_centers()

    # metadata for final spectrum
    for attr in metadata:
        s.__setattr__(attr, metadata[attr])

    return {'ds':ds, 'ims':ims, 'im':im, 's':s}







###############################################################################
###############################################################################
###############################################################################
###############################################################################

def _dataset2dict2(ds, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, norm=True, include_metadata_for_each_image=False):
    
    ################
    # get metadata # 
    ################
    metadata = {}

    # required metadata
    metadata['number_of_images'] = ds.dims['trainId']  # one image per train
    metadata['exposure']         = ds['hRIXS_exposure'].values + 0
    metadata['total_exposure']   = ds['hRIXS_exposure'].sum().values + 0
    metadata['i0']               = ds['SCS_photonFlux'].values + 0
    metadata['time_delay2']      = ds['hRIXS_delay'].values + 0

    # metadata per train
    metadata['raw']  = {attr: ds[attr].values + 0 for attr in _attrs}
    metadata['stat'] = {attr: {'avg': np.mean(metadata['raw'][attr]),
                               'min': np.min(metadata['raw'][attr]),
                               'max': np.max(metadata['raw'][attr]),
                               'std': np.std(metadata['raw'][attr])} for attr in metadata['raw']}
    if include_metadata_for_each_image == False: del metadata['raw']

    # manual assignment
    metadata['E']          = round(metadata['stat']['nrj']['avg'], 2)
    metadata['sample_z']   = round(metadata['stat']['XRD_STZ']['avg'], 5)
    metadata['sample_y']   = round(metadata['stat']['XRD_STY']['avg'], 5)
    metadata['sample_x']   = round(metadata['stat']['XRD_STX']['avg'], 5)
    metadata['GATT']       = metadata['stat']['GATT_pressure']['avg']
    metadata['exit_slit']  = round(metadata['stat']['ESLIT']['avg'], 2)

    ################
    # brixs object #
    ################
    ims = br.Dummy(data=[br.Image(data=_.values).crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop) for _ in ds['hRIXS_det']])    
    
    # metadata for individual images
    if include_metadata_for_each_image:
        for attr in metadata['raw']:
            for i, _im in enumerate(ims):
                ims[i].__setattr__(attr, metadata['raw'][attr][i])
    
    # # scan
    # ims.scan = scan
    # for i, _im in enumerate(ims):
    #     ims[i].scan = scan
    
    # time
    ims.t0 = 300
    ims.time_delay = tb.positionToDelay(metadata['stat']['PP800_DelayLine']['avg'], origin=300)

    ########################
    # curvature correction #
    ########################
    if curv is not None:
        for i, _im in enumerate(ims):
            ims[i] = _im.set_horizontal_shift_via_polyval(p=curv)




    
    ####################
    # aggregate images #
    ####################
    im = br.Image(data=np.zeros(ims[0].shape))
    for i, _im in enumerate(ims):
        if norm:
            factor = 1/metadata['i0'][i]
        else:
            factor = 1
        ims[i]   = _im.set_factor(factor)
        im.data += _im.set_factor(factor).data
    im.x_centers = ims[0].x_centers
    im.y_centers = ims[0].y_centers
    # else:
    #     h     = hRIXS(proposal)
    #     data2 = h.aggregate(ds=ds)
    #     im    = br.Image(data=data2['hRIXS_det'].values).crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
    #     if curv is not None:
    #         im = im.set_horizontal_shift_via_polyval(p=curv)

    #################################
    # metadata for aggregated image #
    #################################
    for attr in metadata:
        im.__setattr__(attr, metadata[attr])
    
    # # scan
    # im.scan = scan 

    # time
    im.t0 = 300
    im.time_delay = ims.time_delay

    ##################
    # final spectrum #
    ##################
    if norm:
        s = im.integrated_columns_vs_x_centers().set_factor(1/metadata['total_exposure'])
    else:
        s = im.integrated_columns_vs_x_centers()

    # metadata for final spectrum
    for attr in metadata:
        s.__setattr__(attr, metadata[attr])

    return {'ds':ds, 'ims':ims, 'im':im, 's':s}

###############################################################################
def _process2(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, norm=True, bins=None, t0=300, include_metadata_for_each_image=False, proposal=None):
    """internal function for _process()

    Returns:
        list with spectra
    """    
    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    # finder
    br.finder.kwargs = vars()
    s = br.finder.search()
    if s != False:
        return s
    
    ###############
    # get dataset #
    ###############
    ds = _get_dataset(scan=scan, proposal=proposal)
    # delay_line = ds['PP800_DelayLine'].values + 0

    if bins is None:
        d = _dataset2dict2(ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)
        s = d['s']
        s.t0 = t0
        s.time_delay = tb.positionToDelay(ds['PP800_DelayLine'], origin=t0).values
        s.scan = scan
        
        br.finder.save(s)
        
        return s
    else:
        ss = br.Spectra()
        ss.delays_real_avg = []
        ss.delays_bin_avg  = []
        arr = [round(x, 5) for x in tb.positionToDelay(ds['PP800_DelayLine'], origin=t0).values]
        indexes, bin_edges = _digitize_via_binning(arr, bins, xmin=None, xmax=None)
        for delay in list(indexes.keys())[::2]:
            _ds = ds.isel({'trainId': indexes[delay]})
            d = _dataset2dict2(_ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)
            s = d['s']
            s.t0 = t0
            
            time_delay = tb.positionToDelay(_ds['PP800_DelayLine'], origin=t0).values
            s.time_delay = time_delay
            s.time_delay_avg = np.mean(time_delay)
            s.time_delay_min = np.min(time_delay)#np.mean(indexes[delay])
            s.time_delay_max = np.max(time_delay)#np.mean(indexes[delay])
            s.time_delay_std = np.std(time_delay)#np.mean(indexes[delay])
            s.scan = scan
            
            ss.append(s)
            # ss.delays_real_avg.append(np.mean(indexes[delay]))
            # ss.delays_bin_avg.append(delay)

        br.finder.save(ss)

        return ss

def process2(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, calib=None, norm=True, t0=300, bins=None, include_metadata_for_each_image=False, proposal=None):
    """process rixs images via integration mode

    images -> curv. corr. -> I0 norm. -> sum images to one image (or n images if bins is not None) -> 
    -> integrate columns to get a spectrum -> exposure norm -> calib
    
    Args:
        scan (int or list): scan number. If list of scans, scans will be aligned
            via cross-correlation and summed up.
        x_start, x_stop, y_start, y_stop (int, optional): detector limits to be
            considered in px. If None, the whole image is used.
        curv (list, optional): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term so f(y) = an*y**n + ... + a1*y + a0 gives the horizontal 
                shift (f(y)) for each pixel row (y). If None, no curvature correction
                is applied. Default is None.
        calib (number, optional): if not None, the x axis is multiplied by calib.
            Calib can be a number (multiplicative term) or a tuple with two 
            numbers (linear and constant terms), like calib=[calib, shift].
        norm (bool, optional): if True, images are normalized by I0, while spectra
            is also normalized by total exposure time (number of images and 
            exposure per image. Default is True.
        bins ():
        include_metadata_for_each_image ():
        t0 (number, optional): delay line t0 value. Default is 300.
        proposal (int, optional): proposal number. If None, proposal number will be read
            from SCS.settings['default_proposal']. Default is None

    Returns:
        spectrum
    """
    s = _process2(scan=scan, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, bins=bins, t0=t0, include_metadata_for_each_image=include_metadata_for_each_image, proposal=proposal)
    
    ###############
    # calibration #
    ###############
    if calib is not None:
        s.calib = calib

    return s

def _process3(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, norm=True, bins=None, t0=300, n1=1, n2=4, floor=[850, None, None, None], avg_threshold=900, avg_threshold_max=None, double_threshold=500000, nbins=1400, include_metadata_for_each_image=False, proposal=None):
    """internal function for _process3()
    """    
    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    # finder
    br.finder.kwargs = vars()
    s = br.finder.search()
    if s != False:
        return s
    
    ###############
    # get dataset #
    ###############
    ds = _get_dataset(scan=scan, proposal=proposal)
    # delay_line = ds['PP800_DelayLine'].values + 0

    if bins is None:
        d = _dataset2dict2(ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=None, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)
        
        ims = d['ims']
        pe = br.PhotonEvents()
        pes = br.Dummy() # temporary
        for im in ims:
            _pe, _ = im.floor(*floor).centroid(n1=n1, n2=n2, avg_threshold=avg_threshold, double_threshold=double_threshold, floor=False, avg_threshold_max=avg_threshold_max)
            # print(len(_pe))
            pes.append(_pe)  # temporary
            pe += _pe
        if curv is not None:
            pe = pe.set_horizontal_shift_via_polyval(curv)
        s = pe.integrated_columns_vs_x_centers(ncols=nbins)
        
        s.copy_attrs_from(d['s'])
        s.t0 = t0
        s.time_delay = tb.positionToDelay(ds['PP800_DelayLine'], origin=t0).values
        s.scan = scan

        # temporary
        s.ims = ims
        s.pes = pes
        
        br.finder.save(s)
        
        return s
    else:
        ss = br.Spectra()
        ss.delays_real_avg = []
        ss.delays_bin_avg  = []
        arr = [round(x, 5) for x in tb.positionToDelay(ds['PP800_DelayLine'], origin=t0).values]
        # arr = [round(x, 5) for x in ds['hRIXS_delay'].values]
        # print(arr)
        indexes, bin_edges = _digitize_via_binning(arr, bins, xmin=None, xmax=None)
        
        for delay in list(indexes.keys())[::2]:
            _ds = ds.isel({'trainId': indexes[delay]})
            d = _dataset2dict2(_ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=None, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)

            # s = d['s']
            ims = d['ims']
            pe = br.PhotonEvents()
            pes = br.Dummy() # temporary
            for im in ims:
                _pe, _ = im.floor(*floor).centroid(n1=n1, n2=n2, avg_threshold=avg_threshold, double_threshold=double_threshold, floor=False, avg_threshold_max=avg_threshold_max)
                # print(len(_pe))
                pes.append(_pe)  # temporary
                pe += _pe
            if curv is not None:
                pe = pe.set_horizontal_shift_via_polyval(curv)
            s = pe.integrated_columns_vs_x_centers(ncols=nbins)

            s.copy_attrs_from(d['s'])
            s.t0 = t0
            time_delay = tb.positionToDelay(_ds['PP800_DelayLine'], origin=t0).values
            s.time_delay = time_delay
            s.time_delay_avg = np.mean(time_delay)
            s.time_delay_min = np.min(time_delay)#np.mean(indexes[delay])
            s.time_delay_max = np.max(time_delay)#np.mean(indexes[delay])
            s.time_delay_std = np.std(time_delay)#np.mean(indexes[delay])
            s.scan = scan

            # temporary
            s.ims = ims
            s.pes = pes
            
            ss.append(s)
            # ss.delays_real_avg.append(np.mean(indexes[delay]))
            # ss.delays_bin_avg.append(delay)

        br.finder.save(ss)

        return ss

def process3(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, calib=None, norm=True, t0=300, bins=None, n1=1, n2=4, floor=[850, None, None, None], avg_threshold=900, avg_threshold_max=None, double_threshold=500000, nbins=1400,  include_metadata_for_each_image=False, proposal=None):
    """same as process2 but with centroid
    """
    s = _process3(scan=scan, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, bins=bins, t0=t0, n1=n1, n2=n2, floor=floor, avg_threshold=avg_threshold, avg_threshold_max=avg_threshold_max, double_threshold=double_threshold, nbins=nbins, include_metadata_for_each_image=include_metadata_for_each_image, proposal=proposal)
    
    ###############
    # calibration #
    ###############
    if calib is not None:
        s.calib = calib

    return s




















def _process4(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, norm=True, bins=None, t0=300, n1=1, n2=4, floor=[850, None, None, None], avg_threshold=900, avg_threshold_max=None, double_threshold=500000, nbins=1400, include_metadata_for_each_image=False, proposal=None, _dummy=0):
    """internal function for _process3()
    """    
    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    # finder
    br.finder.kwargs = vars()
    s = br.finder.search()
    if s != False:
        return s
    
    ###############
    # get dataset #
    ###############
    ds = _get_dataset(scan=scan, proposal=proposal)
    # delay_line = ds['PP800_DelayLine'].values + 0

    if bins is None:
        d = _dataset2dict2(ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=None, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)
        
        ims = d['ims']
        pe = br.PhotonEvents()
        pes = br.Dummy() # temporary
        for im in ims:
            _pe, _ = im.floor(*floor).centroid(n1=n1, n2=n2, avg_threshold=avg_threshold, double_threshold=double_threshold, floor=False, avg_threshold_max=avg_threshold_max)
            # print(len(_pe))
            pes.append(_pe)  # temporary
            pe += _pe
        if curv is not None:
            pe = pe.set_horizontal_shift_via_polyval(curv)
        s = pe.integrated_columns_vs_x_centers(ncols=nbins)
        
        s.copy_attrs_from(d['s'])
        s.t0 = t0
        s.time_delay = tb.positionToDelay(ds['PP800_DelayLine'], origin=t0).values
        s.scan = scan

        # temporary
        s.ims = ims
        s.pes = pes
        
        br.finder.save(s)
        
        return s
    else:
        ss = br.Spectra()
        ss.delays_real_avg = []
        ss.delays_bin_avg  = []
        # arr = [round(x, 5) for x in tb.positionToDelay(ds['PP800_DelayLine'], origin=t0).values]
        arr = [round(x, 5) for x in ds['hRIXS_delay'].values]
        # print(arr)
        indexes, bin_edges = _digitize_via_binning(arr, bins, xmin=None, xmax=None)
        
        for delay in list(indexes.keys())[::2]:
            _ds = ds.isel({'trainId': indexes[delay]})
            d = _dataset2dict2(_ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=None, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)

            # s = d['s']
            ims = d['ims']
            pe = br.PhotonEvents()
            pes = br.Dummy() # temporary
            for im in ims:
                _pe, _ = im.floor(*floor).centroid(n1=n1, n2=n2, avg_threshold=avg_threshold, double_threshold=double_threshold, floor=False, avg_threshold_max=avg_threshold_max)
                # print(len(_pe))
                pes.append(_pe)  # temporary
                pe += _pe
            if curv is not None:
                pe = pe.set_horizontal_shift_via_polyval(curv)
            s = pe.integrated_columns_vs_x_centers(ncols=nbins)

            s.copy_attrs_from(d['s'])
            s.t0 = t0
            time_delay = _ds['hRIXS_delay'].values
            s.time_delay = time_delay
            s.time_delay_avg = np.mean(time_delay)
            s.time_delay_min = np.min(time_delay)#np.mean(indexes[delay])
            s.time_delay_max = np.max(time_delay)#np.mean(indexes[delay])
            s.time_delay_std = np.std(time_delay)#np.mean(indexes[delay])
            s.scan = scan

            # temporary
            s.ims = ims
            s.pes = pes
            
            ss.append(s)
            # ss.delays_real_avg.append(np.mean(indexes[delay]))
            # ss.delays_bin_avg.append(delay)

        br.finder.save(ss)

        return ss

def process4(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, calib=None, norm=True, t0=300, bins=None, n1=1, n2=4, floor=[850, None, None, None], avg_threshold=900, avg_threshold_max=None, double_threshold=500000, nbins=1400,  include_metadata_for_each_image=False, proposal=None):
    """same as process2 but with centroid
    """
    s = _process4(scan=scan, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, bins=bins, t0=t0, n1=n1, n2=n2, floor=floor, avg_threshold=avg_threshold, avg_threshold_max=avg_threshold_max, double_threshold=double_threshold, nbins=nbins, include_metadata_for_each_image=include_metadata_for_each_image, proposal=proposal, _dummy=0)
    
    ###############
    # calibration #
    ###############
    if calib is not None:
        s.calib = calib

    return s
    













# @br.finder.track2
def _process(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, norm=True, is_tr='auto', include_metadata_for_each_image=False, proposal=None):
    """internal function for _process()

    Returns:
        list with spectra
    """    
    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    ###############
    # get dataset #
    ###############
    ds = _get_dataset(scan=scan, proposal=proposal)
    delay_line = ds['PP800_DelayLine'].values + 0

    ################################
    # Does it have multiple delays #
    ################################
    if is_tr == 'auto':
        is_tr = False
        if br.has_duplicates(delay_line, max_error=settings['has_multiple_delays_max_error']):
            is_tr = True

    #####################################
    # sort multiple delays if necessary #
    #####################################
    if is_tr:
        indexes    = _digitize_via_duplicates(delay_line, max_error=settings['has_multiple_delays_max_error'])
        delay_line_list = np.sort(list(indexes.keys()))[::-1]

        ss = br.Spectra()
        for delay in delay_line_list:
            _ds = ds.isel({'trainId': indexes[delay]})
            d = _dataset2dict(_ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)
            s = d['s']
            ss.append(s)

    else:
        d = _dataset2dict(ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image)
        s = d['s']
        ss = br.Spectra([s, ])
    
    return ss.data

def process(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, calib=None, norm=True, t0=300, is_tr='auto', include_metadata_for_each_image=False, proposal=None):
    """process rixs images via integration mode

    images -> curv. corr. -> I0 norm. -> sum images to one image -> 
    -> integrate columns to get a spectrum -> exposure norm -> calib
    
    Args:
        scan (int or list): scan number. If list of scans, scans will be aligned
            via cross-correlation and summed up.
        x_start, x_stop, y_start, y_stop (int, optional): detector limits to be
            considered in px. If None, the whole image is used.
        curv (list, optional): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term so f(y) = an*y**n + ... + a1*y + a0 gives the horizontal 
                shift (f(y)) for each pixel row (y). If None, no curvature correction
                is applied. Default is None.
        calib (number, optional): if not None, the x axis is multiplied by calib.
            Calib can be a number (multiplicative term) or a tuple with two 
            numbers (linear and constant terms), like calib=[calib, shift].
        norm (bool, optional): if True, images are normalized by I0, while spectra
            is also normalized by total exposure time (number of images and 
            exposure per image. Default is True.
        t0 (number, optional): delay line t0 value. Default is 300.
        proposal (int, optional): proposal number. If None, proposal number will be read
            from SCS.settings['default_proposal']. Default is None

    Returns:
        spectrum
    """
    if scan is Iterable:
        _ss = br.Spectra()
        for _scan in scan:
            _temp = _process(scan=_scan, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, is_tr=is_tr, include_metadata_for_each_image=include_metadata_for_each_image, proposal=proposal)
            if len(_temp) > 1:
                raise ValueError(f'scan {_scan} has multiple delays. Cannot add this scan to other scans {scan}')
            else:
                _ss.append(_temp[0])

        for attr in _ss[0].get_attrs():
            try:
                _ss.create_attr_from_spectra(attr)
            except:
                pass
        s = _ss.align().interp().calculate_average()
    else:
        _temp = _process(scan=scan, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, is_tr=is_tr, include_metadata_for_each_image=include_metadata_for_each_image, proposal=proposal)
        if len(_temp) > 1:
            s = br.Spectra(data=_temp)
        else:
            s = _temp[0]

    ###############
    # calibration #
    ###############
    if calib is not None:
        s.calib = calib

    ############
    # metadata #
    ############
    if isinstance(s, br.Spectra):
        # print('multiple delays')
        for _s in s:
            # scan
            _s.scan = scan 
            # time
            _s.t0 = t0
        for attr in s[0].get_attrs():
            try:
                s.create_attr_from_spectra(attr)
            except:
                pass
    else:
        s.t0 = t0
    s.scan = scan 

    #################
    # normalization #
    #################
    # if norm:
    #     s = s.set_factor(1/s.total_exposure)

    return s

# %% ============================ RIXS (Advanced) ========================= %% #
def sequence(scans, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, calib=None, norm=True, t0=300, include_metadata_for_each_image=False, proposal=None):
    """return a list with rixs spectra

    Example:

    >>> # ss1 will have 3 spectra
    >>> ss1 = sequence([100, 101, 102])
    >>>
    >>> # ss2 will have 3 spectra. The middle one will a sum of 2 spectra
    >>> ss2 = sequence([100, [101, 102], 103])

    Args:
        scans (list): list of rixs scan number. Replace a scan number for a list to 
            sum spectra inside list.
        x_start, x_stop, y_start, y_stop (int, optional): detector limits to be
            considered in px. If None, the whole image is used.
        curv (list, optional): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term so f(y) = an*y**n + ... + a1*y + a0 gives the horizontal 
                shift (f(y)) for each pixel row (y). If None, no curvature correction
                is applied. Default is None.
        calib (number, optional): if not None, the x axis is multiplied by calib.
            Calib can be a number (multiplicative term) or a tuple with two 
            numbers (linear and constant terms), like calib=[calib, shift].
        norm (bool, optional): if True, images are normalized by I0, while spectra
            is also normalized by total exposure time (number of images and 
            exposure per image. Default is True.
        t0 (number, optional): delay line t0 value. Default is 300.
        proposal (int, optional): proposal number. If None, proposal number will be read
            from SCS.settings['default_proposal']. Default is None

    Returns:
        Spectra
    """
    ss = br.Spectra()
    for scan in scans:
        ss.append(process(scan=scan, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, calib=calib, norm=norm, is_tr=False, t0=t0, include_metadata_for_each_image=include_metadata_for_each_image, proposal=proposal))
    
    for attr in ss[0].get_attrs():
        try:
            ss.create_attr_from_spectra(attr)
        except:
            pass

    return ss

def verify(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, norm=True, include_metadata_for_each_image=False, proposal=None, **kwargs):
    """plot separate images for a run
    
    Args:
        scan (int): scan number
        x_start, x_stop, y_start, y_stop (int, optional): detector limits to be
            considered in px. If None, the whole image is used.
        curv (list, optional): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term so f(y) = an*y**n + ... + a1*y + a0 gives the horizontal 
                shift (f(y)) for each pixel row (y). If None, no curvature correction
                is applied. Default is None.
        norm (bool, optional): if True, images are normalized by I0, while spectra
            is also normalized by total exposure time (number of images and 
            exposure per image. Default is True.
        proposal (int, optional): proposal number. If None, proposal number will be read
            from SCS.settings['default_proposal']. Default is None

    Returns:
        dictionary {'ds':ds, 'ims':ims, 'im':im, 's':s}
    """    
    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    ###############
    # get dataset #
    ###############
    ds = _get_dataset(scan=scan, proposal=proposal)
    
    ################
    # process data #
    ################
    d   = _dataset2dict(ds=ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image, proposal=proposal)
    ims = d['ims']
    im  = d['im']
    s   = d['s']
    number_of_images = s.number_of_images
    
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
        ims[ims.__i].plot(ax=axes[0], **kwargs)

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
    
    ax2 = axes[6].twinx()
    ax2.tick_params(axis='y', labelcolor='red')
    for i in (4, 5, 6, 7):
        axes[i].ymove(-0.1)
    # ax2.ymove(-0.1)
    
    for i in (0, 4):
        axes[i].remove_xticklabels()
    
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
    # plot initial photon events (axes 0)
    ims[0].plot(ax=axes[0], **kwargs)
    
    # plot initial spectra (axes 1, 2)
    ims[0].integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[1])
    ims[0].integrated_columns_vs_x_centers().plot(ax=axes[2])
    
    # plot photon events summed (axes 4)
    im.plot(ax=axes[4])

    # plot spectra summed (axes 5 and 6)
    im.integrated_rows_vs_y_centers().switch_xy().plot(ax=axes[5])
    curveA = im.integrated_columns_vs_x_centers().plot(ax=axes[6], label='direct sum')
    
    _ss = br.Spectra()
    for _im in ims:
        _ss.append(_im.integrated_columns_vs_x_centers())
    curveB = _ss.align().calculate_sum().plot(ax=ax2, color='red', label='align spectra before sum (no norm.)')

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
    br.leg(ax=axes[6], handles=[curveA, curveB])
    
    #################
    # remove unused #
    #################
    axes[3].remove()
    axes[7].remove()

    ######################
    # register callbacks #
    ######################
    fig.canvas.mpl_connect('key_press_event', lambda event: keyboard(event, ims=ims, axes=axes))
    return d

def verify_curv(scan, popt=None, ncols=None, nrows=None, deg=2, x_start=None, x_stop=None, y_start=None, y_stop=None, proposal=None):
    """plot step-by-step curvature correction process
    
    Args:
        scan (int): scan number
        popt (list): If None, curvature will be calculated via cross-correlation.
            If list, curvature correction will be applied where popt is assumed
             to be the 1D array of polynomial coefficients (including 
            coefficients equal to zero) from highest degree to the constant 
            term so f(y) = an*y**n + ... + a1*y + a0 gives the horizontal 
            shift (f(y)) for each pixel row (y). Default is None.
        ncols, nrows (int, optional): number of row and columns to bin (reduce)
            the image. If None, the max number of pixels will used. It is 
            recommended nrows to be a number around 10.
        deg (int, optional): Degree of the curvature fitting polynomial. 
                Default is 2.
        x_start, x_stop, y_start, y_stop (int, optional): detector limits to be
            considered in px. If None, the whole image is used.
        proposal (int, optional): proposal number. If None, proposal number will be read
            from SCS.settings['default_proposal']. Default is None

    Returns:
        curvature correction optimized parameters (list)
    """    
    #######################
    # get proposal number #
    #######################
    if proposal is None:
        proposal = settings['default_proposal']
    assert proposal is not None, 'proposal cannot be None. Either provide a proposal number of set default proposal number at SCS.settings["default_proposal"]'

    ###############
    # get dataset #
    ###############
    ds = _get_dataset(scan=scan, proposal=proposal)
    
    ################
    # process data #
    ################
    d   = _dataset2dict(ds=ds, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=None, norm=False, include_metadata_for_each_image=False, proposal=proposal)
    ims = d['ims']
    im  = d['im']
    s   = d['s']

    # if polynomial parameters are not given, calculate curvature
    if popt is None:
        if ncols is None: ncols = im.shape[1]
        if nrows is None: nrows = im.shape[0]

        reduced = im.binning(ncols=ncols, nrows=nrows).crop(y_start=y_start, y_stop=y_stop).floor()

        ss = reduced.rows.floor()
        values = ss.calculate_shift(mode='cc')

        s = br.Spectrum(x=reduced.y_centers, y=values)
        fit, popt, R2, model = s.polyfit(deg=deg)             
    else:
        reduced = im
        fit = br.Spectrum(x=im.y_centers, y=np.polyval(popt, im.y_centers))
        ss  = br.Spectra()

    # apply curvature corredctions
    im2 = im.set_horizontal_shift_via_polyval(p=popt)

    #######################
    # initializing figure #
    #######################
    fig, axes = br.subplots(4, 1, height_ratios=[2, 2, 1, 2], hspace=0, figsize=(12, 20))
    plt.subplots_adjust(top=0.97, right=0.99, bottom=0.25)

    for i in (1, 2, 3, 3):
        axes[i].ymove(-0.08)

    for i in (1, ):
        axes[i].remove_xticklabels()

    ##############
    # share axis #
    ##############
    br.sharex([axes[1], axes[2]])

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
    ss.plot(ax=axes[2])

    if len(ss) > 0:
        shift = ss[0].get_x_where_y_is_max()
    else:
        shift = im2.integrated_columns_vs_x_centers().get_x_where_y_is_max()
    fit.crop(y_start, y_stop).set_factor(-1).switch_xy().plot(shift=shift, ax=axes[1], color='black')
        

    im2.plot(ax=axes[3])

    ##############
    # set labels #
    ##############
    for i in (0, 1, 3):
        axes[i].set_ylabel('y (pixel)', fontsize='x-small')

    for i in (2, ):
        axes[i].set_ylabel('counts/pixel row', fontsize='x-small')

    for i in (0, 1, 2, 3):
        axes[i].set_xlabel('x (pixel)', fontsize='x-small')

    return popt

def verify_calib(scans, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, norm=True, proposal=None):
    """plot step-by-step calibration value calculation process

    Args:
        scan (int): scan number
        x_start, x_stop, y_start, y_stop (int, optional): detector limits to be
            considered in px. If None, the whole image is used.
        curv (list, optional): 1D array of polynomial coefficients (including 
                coefficients equal to zero) from highest degree to the constant 
                term so f(y) = an*y**n + ... + a1*y + a0 gives the horizontal 
                shift (f(y)) for each pixel row (y). If None, no curvature correction
                is applied. Default is None.
        norm (bool, optional): if True, images are normalized by I0, while spectra
            is also normalized by total exposure time (number of images and 
            exposure per image. Default is True.
        proposal (int, optional): proposal number. If None, proposal number will be read
            from SCS.settings['default_proposal']. Default is None

    Returns:
        calibratio value
    """    
    ss = sequence(scans=scans, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, calib=None, norm=norm, include_metadata_for_each_image=False, proposal=proposal)
    # ss = br.Spectra()
    # for scan in scans:
    #     ss.append(process(scan=scan, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, proposal=proposal))
     
    ss = ss.floor()
    s = ss.calculate_calib(values=[_.E for _ in ss], mode='cc', deg=1)
    
    #######################
    # initializing figure #
    #######################
    fig, axes = br.subplots(2, 1, height_ratios=[2, 2], hspace=0, sharex=True, figsize=(10, 10))
    plt.subplots_adjust(top=0.97, right=0.99, left=0.2, bottom=0.15)
    
    ########
    # plot #
    ########
    # plot initial photon events (axes 0)
    ss.plot(ax=axes[0])
    
    # plot initial spectra (axes 1, 2)
    _s, popt, err, f = ss[0].fit_peak()
    s.plot(ax=axes[1], shift=popt[1], marker='o', color='black')
    s.fit.plot(ax=axes[1], shift=popt[1], color='red')

    ##############
    # set labels #
    ##############
    for i in (0, ):
        axes[i].set_ylabel('Intensity', fontsize='x-small')
    
    for i in (1, ):
        axes[i].set_xlabel('x (pixel)', fontsize='x-small')
    
    for i in (1, ):
        axes[i].set_ylabel('Photon energy (eV)', fontsize='x-small')
        
    return s.popt
# %%

# def pa():
#     s = _dict2spectrum(scan=scan, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, include_metadata_for_each_image=include_metadata_for_each_image, proposal=proposal)

#     ################
#     # get metadata # 
#     ################
#     metadata = {}
#     metadata['number_of_images'] = ds.dims['trainId']  # one image per train
#     metadata['exposure']         = ds['hRIXS_exposure'].values + 0
#     metadata['total_exposure']   = ds['hRIXS_exposure'].sum().values + 0
#     metadata['i0']               = ds['SCS_photonFlux'].values + 0
#     # metadata['t0']               = None
#     # metadata['time_delay']       = None
#     # i0_sa3           = ds['SCS_SA3'].sum(dim='sa3_pId').values + 0
#     # i0_sa3           = ds['SCS_SA3'].mean(dim='sa3_pId').values + 0

#     metadata['raw']  = {attr: ds[attr].values + 0 for attr in _attrs}
#     metadata['stat'] = {attr: {'avg': np.mean(metadata['raw'][attr]),
#                                'min': np.min(metadata['raw'][attr]),
#                                'max': np.max(metadata['raw'][attr]),
#                                'std': np.std(metadata['raw'][attr])} for attr in metadata['raw']}
#     if raw == False: del metadata['raw']

#     metadata['E']          = round(metadata['stat']['nrj']['avg'], 2)
#     # metadata['delay_line'] = metadata['stat']['PP800_DelayLine']['avg']
#     # metadata['time_delay'] = tb.metadata['stat']['PP800_DelayLine']['avg']
#     metadata['sample_z']   = round(metadata['stat']['XRD_STZ']['avg'], 5)
#     metadata['sample_y']   = round(metadata['stat']['XRD_STY']['avg'], 5)
#     metadata['sample_x']   = round(metadata['stat']['XRD_STX']['avg'], 5)
#     metadata['GATT']       = metadata['stat']['GATT_pressure']['avg']
#     metadata['exit_slit']  = round(metadata['stat']['ESLIT']['avg'], 2)

#     ################
#     # brixs object #
#     ################
#     ims = br.Dummy(data=[br.Image(data=_.values).crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop) for _ in ds['hRIXS_det']])    
    
#     # metadata for individual images
#     if 'raw' in metadata:
#         for attr in metadata['raw']:
#             for i, _im in enumerate(ims):
#                 ims[i].__setattr__(attr, metadata['raw'][attr][i])
    
#     # scan
#     ims.scan = scan
#     for i, _im in enumerate(ims):
#         ims[i].scan = scan
    
#     # time
#     ims.t0_ = 300
#     ims.time_delay = tb.positionToDelay(metadata['stat']['PP800_DelayLine']['avg'], origin=300)

#     ########################
#     # curvature correction #
#     ########################
#     if curv is not None:
#         for i, _im in enumerate(ims):
#             ims[i] = _im.set_horizontal_shift_via_polyval(p=curv)
    
#     ####################
#     # aggregate images #
#     ####################
#     if norm:
#         im = br.Image(data=np.zeros(ims[0].shape))
#         for i, _im in enumerate(ims):
#             ims[i]   = _im.set_factor(1/metadata['i0'][i])
#             im.data += _im.set_factor(1/metadata['i0'][i]).data
#         im.x_centers = ims[0].x_centers
#         im.y_centers = ims[0].y_centers
#     else:
#         data2 = h.aggregate(ds=ds)
#         im    = br.Image(data=data2['hRIXS_det'].values).crop(x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop)
#         if curv is not None:
#             im = im.set_horizontal_shift_via_polyval(p=curv)

#     #################################
#     # metadata for aggregated image #
#     #################################
#     for attr in metadata:
#         im.__setattr__(attr, metadata[attr])
    
#     # scan
#     im.scan = scan 

#     # time
#     im.t0_ = 300
#     im.time_delay = ims.time_delay

#     ##################
#     # final spectrum #
#     ##################
#     if norm:
#         s = im.integrated_columns_vs_x_centers().set_factor(1/metadata['total_exposure'])
#     else:
#         s = im.integrated_columns_vs_x_centers()

#     # metadata for final spectrum
#     for attr in metadata:
#         s.__setattr__(attr, metadata[attr])

#     return {'ds':ds, 'ims':ims, 'im':im, 's':s}

# @br.finder.track
# def _process(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, 
#              curv=None, norm=True, proposal=None, raw=False):
#     """same as internal but with finder"""
#     d = _internal(scan=scan, proposal=proposal, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, raw=raw)
#     return d['s']

# def process(scan, x_start=None, x_stop=None, y_start=None, y_stop=None, curv=None, 
#             calib=None, norm=True, t0=300, proposal=None, raw=False):
#     """process rixs images via integration mode

#     images -> curv. corr. -> I0 norm. -> sum images to one image -> 
#     -> integrate columns to get a spectrum -> exposure norm -> calib
    
#     Args:
#         scan (int or list): scan number. If list of scans, scans will be aligned
#             via cross-correlation and summed up.
#         x_start, x_stop, y_start, y_stop (int, optional): detector limits to be
#             considered in px. If None, the whole image is used.
#         curv (list, optional): 1D array of polynomial coefficients (including 
#                 coefficients equal to zero) from highest degree to the constant 
#                 term so f(y) = an*y**n + ... + a1*y + a0 gives the horizontal 
#                 shift (f(y)) for each pixel row (y). If None, no curvature correction
#                 is applied. Default is None.
#         calib (number, optional): if not None, the x axis is multiplied by calib.
#             Calib can be a number (multiplicative term) or a tuple with two 
#             numbers (linear and constant terms), like calib=[calib, shift].
#         norm (bool, optional): if True, images are normalized by I0, while spectra
#             is also normalized by total exposure time (number of images and 
#             exposure per image. Default is True.
#         t0 (number, optional): delay line t0 value. Default is 300.
#         proposal (int, optional): proposal number. If None, proposal number will be read
#             from SCS.settings['default_proposal']. Default is None

#     Returns:
#         spectrum
#     """
#     if scan is Iterable:
#         _ss = br.Spectra()
#         for _scan in scan:
#             s = _process(scan=_scan, proposal=proposal, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm, raw=raw)
#             _ss.append(s)

#         for attr in _ss[0].get_attrs():
#             try:
#                 _ss.create_attr_from_spectra(attr)
#             except:
#                 pass

#         s = _ss.align().interp().calculate_sum()


#     else:
#         s = _process(scan=scan, proposal=proposal, x_start=x_start, x_stop=x_stop, y_start=y_start, y_stop=y_stop, curv=curv, norm=norm)
    
#     ###############
#     # calibration #
#     ###############
#     if calib is not None:
#         if isinstance(calib, br.Iterable):
#             s.shift = -calib[1]
#             s.calib = calib[0]
#         else:
#             s.calib = calib

#     ##########
#     # set t0 #
#     ##########
#     s.t0 = t0

#     return s



# %% ============================= Centroid =============================== %% #
def _centroid(self, n1, n2, avg_threshold, double_threshold, floor=False, avg_threshold_max=None, include_doubles=True, MAX_PHOTONS=100000):
        """Returns a PhotonEvents object.
        
        This function is adapted from the ESRF RIXS toolbox and from XFEL SCS toolbox.
        
        Warning:
            for the algorithm to work, the background average must be zero. One can acomplish that by
            using the function im.floor(x_start, x_stop, y_start, y_stop), where start and stop args
            can be used to select a region of the image which is only background. One can also use the
            argument floor=True in this function. This will `floor` the whole image (the average of the
            whole image is zero). This should be ok for most images with low number of photons.
        
        Note:
            image x_centers and y_centers must be monotonic
            
        Args:
            n1 (int): number of pixels to average to detect a photon hit candidate, e.g., a photon hit candidate 
            is selected if the average of the intensities within a n1-by-n1 square exceeds avg_threshold.
            n2 (int): For a photon hit candidates the 'center-of-mass' of the photon hit is calculated 
                within a (n1+n2)-by-(n1+n2) square. n2 must be an even number to ensure that the 
                `center of mass` of a photon hit is calculated in a square where the pixel with the 
                photon hit is the central pixel.
            avg_threshold (number): any pixel with intensity higher than avg_threshold in the n1-by-n1 
                averaged image will be selected as a photon-hit-candidate position.
            double_threshold (number): any a photon-hit-candidate position where the sum of the surounding
                n1+n2 -by- n1+n2 square is higher than double_threshold is considered a double hit
            floor (bool, optional): if True, an intensity offset is added to the image such as the average
                intensity of the whole image is zero. Default is False.
            avg_threshold_max (number or None, optional): any pixel with intensity higher than avg_threshold_max
                in the n1xn1 averaged image will removed as a photon-hit-candidate position. If None,
                this functionaly is disabled and only the lower limit avg_threshold is used. Default is None.
            include_doubles (bool, optional): if False, double photon events will be removed from the
                final result
            MAX_PHOTONS (number or False): if number, this function will raise an error if the number of 
                detected photons exceed MAX_PHONTOS
        
        Note:
            The averaging n1 prevents that noise is counted as a candidate. By averaging, we are
            requiring that, not only the `central` pixel is iluminated, but also the surrounding 
            pixels. n1=2 seems like a good averaging for most cases.
        
        Note:
            n1+n2 must be roughly the same number of pixel which a photon can excite. For example,
            if n1+n2 = 4, we are expecting that a photon will excite at most a 4x4 array of pixels.
        
        Note:
            if n1+n2 is too large, than the algorithm might mistakenly assign double photon hits as
            it will think that multiple photons will be falling whitin the n1+n2 squares.
            
        Note:
            I think a good metric for seting up the avg_threshold is to choose the highest threshold possible
            such as one has no photon at `negative` energy loss. One can do this by slowly increasing 
            avg_threshold, and/or by visually inspecting the image using
            
            >>> br.figure()
            >>> im.moving_average(n1).plot()
            
        Note:
            For determining a resonable double_threshold, I think the best option is to get a pixel which
            for sure is a single hit, then draw a n1+n2 square around it, and sum the intensities inside 
            this square and multiply by 1.5. For example, given the following matrix with a single photon
            hit in the center (note that the center is slightly shifted to the left and top because the 
            lenght of the square is a even number):
            
                [[ 4,  6, -4,  8],
                 [-9, 79, 52, 10],
                 [-2, 19, 17, -2],
                 [-6, -9,  2, -8]]
            
            we see that the photon hit is on 79, then for n1+n2=4 we have that the sum is 157 (see also example
            below). A suitable double hit treshold would be 1.5 * 157 = 235.
            
            
        Example:
        
            This example goes over the step-by-step algorithm for detecting a candidate and assigning a photon hit.
        
            Given the follwoing matrix:

                [[14, -7,  4,  6, -4,  8],
                 [14,  4, -9, 79, 52, 10],
                 [ 1, 17, -2, 19, 17, -2],
                 [-1, -4, -6, -9,  2, -8],
                 [-4,  3,  2,  7, -9, -3]]

            Its approx. rounded averaged counterpart for n1=2 is:

                [[ 6, -4, 18, 33, 16,  1],
                 [ 9,  1, 20, 42, 19,  0],
                 [ 3,  1,  0,  6,  2, -5],
                 [-2, -1, -2, -5, -7, -6],
                 [ 0,  4,  4, -2, -6,  0]]
                 
            (note that in this example we are using a aprox rounded average to facilitate the visualization
            of the matrix, but the script does an exact calculation).
                 
            For avg_treshold = 25, we have two spots which are candidates: (x, y) = (3, 0) and (3, 1).
            
            Between these two spots, (3, 0) will be disregarded, because the intensity of (3, 1) is 
            brighest.
            
            Going back to the original matrix, the position (3, 1) yields intensity 79. We then get
            pixels surrounding it (n1 x n1) like:
            
            [[79, 52],
             [19, 17]]
             
            (these are the pixels that yielded 42 in the averaged matrix).
            
            We then add pixel rows and cols to the left/right/top/bottom so to make a square of size
            n1+n2. For n2 = 2 we have
            
            [[ 4,  6, -4,  8],
             [-9, 79, 52, 10],
             [-2, 19, 17, -2],
             [-6, -9,  2, -8]]
            
            this is what we call a `spot`. Now we only have to determine if this spot is a double or
            single hit. 
            
            the sum of the spot is 157. Let's say that double_threshold = 255, therefore, this spot is
            a single hit.
            
            Here is an example with a double hit
            
            Original:
            [[-6, -1,  8, -4,  6,  0],
             [ 4, -9, 14, 23,  0,  7],
             [ 4,  4, 99, 41, -9, -9],
             [37, 81, 50, 12, 13, -9],
             [-3, 46, 16, -6,  5, -9]]
            
            aprox. rounded averaged n1=2
            [[-4,  2, 10,  6,  3,  1],
             [ 0, 27, 45, 14, -4, -1],
             [31, 59, 51, 14, -6, -3],
             [40, 48, 18,  6, -1, -1],
             [12, 17,  3, -1,  0, -2]]
             
            From the averaged matrix we have many candidates (intensity > avg_treshold = 25). Taking the
            brighest one, we have x, y = (1, 2), which has intensity 59
            
            the spot for this pixel is (given n2=2)
            
            [[ 4, -9, 14, 23],
             [ 4,  4, 99, 41],
             [37, 81, 50, 12],
             [-3, 46, 16, -6]]
             
            For this example, we know that a photon hit threshold should be around 200. Therefore, a 
            suitable double hit treshold would be 1.5 * 200 = 300.
             
            This sum is 413, so it will be regarded as a double photon hit.
            
            The position of a photon hit is then calculated by the 2D weighted sum of the intensities 
            (the mean x and y values within a spot).
            
            
            
            
                
                
        Obsolete:
            Note:
                ESRF algorithm: n1=1, n2=2, mode='energy', energy=<photon-energy>
                XFEL algorithm: n1=2, n2=2, mode='std' [note that XFEL cannot be 100% 
                reproduced with this code because in the XFEL algorithm the brightest pixel
                within a `spot` is detected on the averaged (moving average) image instead
                of on the original image. Mostly, the original XFEL code leads to some candidates 
                at the edges to be left out]
            
            Args:
                mode (str): mode for calculating threshold
                    `manual`: 'treshold' must be passed as kwargs

                    `std` [XFEL]: 'treshold' must be given in terms of `number of 
                    standard deviations`, i.e., 'treshold'=3.5 (default) will assume
                    that pixels with intensity higher than intensity_avg+3.5*intensity_std
                    is a candidate for photon hit.

                    `energy` [ESRF]: 'threshold' is defined as (0.2) * (photon_energy/3.6/1.06)
                    where 3.6 and 1.06 are a factor defined by the ADU of the detector. 
                    Photon energy must be passed with the kwarg 'energy'.
                double_factor (number, optional): the treshold factor in which to indentify double photon events, e.g.,
                    if the summed intensity of a spot [(n1+n2)-by-(n1+n2) square] is bigger than 
                    treshold*double_factor, then this spot is identified as a double photon event.
        
        Returns:
            PhotonEvents, double PhotonEvents in terms of x_centers and y_centers
        """
        ###################
        # check n1 and n2 #
        ###################
        assert n1 > 0,  f'n1 must be an positive integer'
        assert n2 >= 0,  f'n2 must be an positive even integer'
        assert n2%2 == 0, f'n2 must be an even number to ensure that the `center of mass` of a photon hit is calculated in a square where the pixel with the photon hit is the central one'
        
        #################################
        # check x_centers and y_centers #
        #################################
        if self.x_monotonicity is None:
            self.check_x_monotonicity()
        if self.y_monotonicity is None:
            self.check_y_monotonicity()
        
        ###############
        # floor image #
        ###############
        if floor:
            image = self.floor()
        else:
            image = self.copy()
        
        ##################
        # moving average #
        ##################
        if n1 > 1:
            image = image.moving_average(n1)

        ##########################
        # check double threshold #
        ##########################
        assert double_threshold  > avg_threshold, 'avg_threshold must be smaller than double_threshold'
        
        # if mode == 'manual':
        #     assert 'treshold' in kwargs, f'treshold missing. mode=`manual` requires kwarg `treshold`'
        #     low_treshold  = kwargs['treshold']
        #     if 'high_treshold' in kwargs:
        #         high_treshold = kwargs['high_treshold']
        #     else:
        #         high_treshold = image.max()
                
#         elif mode == 'std':
#             if 'treshold' not in kwargs:
#                 kwargs['treshold'] = 3.5
#             # assert 'treshold' in kwargs, f'treshold missing. mode=`std` requires kwarg `treshold`'
#             low_treshold = image.calculate_average() + kwargs['treshold'] * image.calculate_sigma()
#             if 'high_treshold' in kwargs:
#                 high_treshold = kwargs['high_treshold']
#             else:
#                 high_treshold = image.max()
            
#         elif mode == 'energy':
#             assert 'energy' in kwargs, f'Photon energy missing. mode=`energy` requires kwarg `energy`'
                    
#             # Multiplication factor * ADU/photon
#             photons       = energy/3.6/1.06
#             low_treshold  = 0.2 * photons
#             high_treshold = 1.0 * photons
            
#             # SpotLOW    = 0.2 * photons
#             # SpotHIGH   = 1.5 * photons
#         else:
#             raise ValueError(f'mode={mode} not valid. Valid modes: `manual`, `std`, `energy`')

        ###################
        # find candidates #
        ###################
        # remove the edges of image (because one cannot calculated the `center of mass` at the edges)
        # select pixels with intensity between low_treshold and high_treshold
        if avg_threshold_max is not None:
            assert avg_threshold_max > avg_threshold, 'avg_threshold must be smaller than avg_threshold_max'
            cp = np.argwhere((image.data[n2//2:-n2//2, n2//2:-n2//2] > avg_threshold)*(image.data[n2//2 : -n2//2, n2//2 : -n2//2] < avg_threshold_max))
        else:
            cp = np.argwhere((image.data[n2//2:-n2//2, n2//2:-n2//2] > avg_threshold))
            
        #############################################################
        # shift photon position because we removed the edges before # 
        #############################################################
        # when we were looking for candidates, we excluded the edges
        # now, to have the candidates position in terms of positions 
        # in the original image (image) we have to add back the edges
        # that we excluded
        cp += np.array((n2//2, n2//2))
        
        #################################
        # max allowed number of photons #
        #################################
        if MAX_PHOTONS:
            if len(cp) > MAX_PHOTONS:
                raise RuntimeError(f'number of detected photons is too high (> {MAX_PHOTONS}). Maybe try and change MAX_PHOTONS or Threshold.')
    
        ########################
        # centroid algorithm 1 #
        ########################
        # runs slower than algorithm 2 (I guess)
        # has the problem that the brightest pixel must have n1 x n1 average above threshold for this pixel to be counted
#         res  = []
#         dres = []
#         for i, (y, x) in enumerate(cp):

#             # isolate the spot where photon hit is in the center
#             spot = self.data[y-n2//2:y+n1+n2//2, x-n2//2:x+n1+n2//2]

#             # check if the central spot is the brightest one
#             if (spot > self.data[y, x]).sum() == 0:

#                 # calculate x center of mass
#                 mx = np.average(np.arange(x-n2//2, x+n1+n2//2), weights=spot.sum(axis=0))

#                 # calculate y center of mass
#                 my = np.average(np.arange(y-n2//2, y+n1+n2//2), weights=spot.sum(axis=1))

#                 # check if spot is a double event
#                 # if (spot.sum() >= SpotLOW) and (spot.sum() < SpotHIGH):         
#                 if (spot.sum() <= double_threshold):               
#                     res.append((mx, my))
#                 else:
#                     if include_doubles:
#                         res.append((mx, my))
#                         res.append((mx, my))
#                     dres.append((mx, my))
          
        ########################
        # centroid algorithm 2 #
        ########################
        # an improved version of algorithm 1
        flag = []
        res  = []
        dres = []
        for i, (y, x) in enumerate(cp):       
            
            # isolate the spot where photon hit is in the center
            spot = self.data[y-n2//2:y+n1+n2//2, x-n2//2:x+n1+n2//2]
            
            # get brightest pixel
            _y, _x = np.unravel_index(np.argmax(spot), np.array(spot).shape)
            bx = x + (_x - n2//2)
            by = y + (_y - n2//2)
                             
            if (bx, by) not in flag:
                # flag that this point has been accounted for
                flag.append((bx, by))
                
                # calculate x center of mass
                mx = np.average(np.arange(x-n2//2, x+n1+n2//2), weights=spot.sum(axis=0))
                
                # calculate y center of mass
                my = np.average(np.arange(y-n2//2, y+n1+n2//2), weights=spot.sum(axis=1))

                # check if spot is a double event
                if (spot.sum() <= double_threshold):               
                    res.append((mx, my))
                else:
                    if include_doubles:
                        res.append((mx, my))
                        res.append((mx, my))
                    dres.append((mx, my))

        ############################
        # convert pixel to centers #
        ############################
        data = np.array([[self.x_centers[int(x)] + (self.x_centers[int(x)+1]-self.x_centers[int(x)])*(x-int(x)), self.y_centers[int(y)] + (self.y_centers[int(y)+1]-self.y_centers[int(y)])*(y-int(y))] for x, y in res])        
        if len(data) == 0:
            pe = br.PhotonEvents()
        else:
            pe = br.PhotonEvents(x=data[:, 0], y=data[:, 1], xlim=(min(self.x_centers), max(self.x_centers)), ylim=(min(self.y_centers), max(self.y_centers)))
        pe.copy_attrs_from(self)
        
        data = np.array([[self.x_centers[int(x)] + (self.x_centers[int(x)+1]-self.x_centers[int(x)])*(x-int(x)), self.y_centers[int(y)] + (self.y_centers[int(y)+1]-self.y_centers[int(y)])*(y-int(y))] for x, y in dres])        
        if len(data) == 0:
            bad = br.PhotonEvents()
        else:
            bad = br.PhotonEvents(x=data[:, 0], y=data[:, 1])
        bad.copy_attrs_from(self)

        return pe, bad
br.Image._centroid = _centroid            
# %%

# %% ============================= Obsolete =============================== %% #
# _main  = []
# _main += ['FastADC2_9raw']            # Detector (photon diode? )
# _main += ['SCS_SA3']                  # ionization chamber after SASE3 monochomator
# _main += ['SCS_photonFlux']           # ionization chamber before SASE3 monochomator
# _main += ['nrj']                      # photon energy
# _main += ['nrj_target']               # energy setpoint (necessary to check if energy is supposed to be changing)
# _main += ['PP800_DelayLine']          # delay line motor 
    
# _metadata  = []
# _metadata += ['VSLIT', 'HSLIT', 'ESLIT']  # vertical, horizontal, and exit slits
# _metadata += ['transmission']             # attenuators (laser)? 
# _metadata += ['transmission_setpoint']    # attenuators (laser)? 
# _metadata += ['transmission_col2']        # attenuators
# _metadata += ['GATT_pressure']            # gas attenuator. After sase3 undulator and slits, before M1
# _metadata += ['XRD_STX', 'XRD_STY', 'XRD_STZ']  # manipulator x, y and z positions
# _metadata += ['XRD_DRY', 'XRD_SRX', 'XRD_SRY', 'XRD_SRZ']
# _metadata += ['XRD_SXT1Y', 'XRD_SXT2Y', 'XRD_SXTX', 'XRD_SXTZ']

# _metadata_rixs = []
# _metadata_rixs += ['hRIXS_exposure', 'nrj', 'nrj_target',  'SCS_SA3', 'SCS_photonFlux', 'PP800_DelayLine']

# # temporary
# # _metadata += ['laser', 'npulses_laser', 'maindump', 'PP800_SynchDelayLine', 'tpi', 'bunchPatternTable']  # apparently useless

# # metadata simple names will also create attrs with simplified names
# _metadata_simple_names = {'ESLIT':           'exit_slit'
#                          ,'GATT_pressure':   'GATT'
#                          ,'XRD_STX':         'sample_x'
#                          ,'XRD_STY':         'sample_y'
#                          ,'XRD_STZ':         'sample_z'
#                          ,'PP800_DelayLine': 'delay_line'
#                          ,'nrj':             'E'}

# _extra  = []
# # _extra += ['PP800_PhaseShifter', 'PP800_SynchDelayLine', 'PP800_HalfWP', 'PP800_FocusLens']  # other delay fields
# # _extra += ['sase3', 'sase2', 'sase1']
# # _extra += ['bunchpattern', 'bunchPatternTable']
# # _extra += ['npulses_sase3', 'npulses_sase1']
# # _extra += ['mono_order', 'M2BEND', 'tpi']
# # _extra += ['UND', 'UND2', 'UND3', 'MAG_CHICANE_DELAY']
# # _extra += ['AFS_DelayLine', 'AFS_FocusLens']
# # _extra += ['LIN_FocusLens', 'FFT_FocusLens']
# # _extra += ['FastADC0raw', 'FastADC3raw', 'FastADC5raw', 'FastADC9raw'] # 
