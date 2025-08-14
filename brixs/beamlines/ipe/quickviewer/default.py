# %% ============================= USER PROPOSAL ========================== %% #
proposal = '20232839'
# %%

# %% ==========================  STANDARDS IMPORTS ======================== %% #
from collections.abc import Iterable
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import sys
import os

# %%

# %% =============================== BASE ================================= %% #
if os.name == 'nt':  # Windows
    BASE = Path('Z:/')
elif os.name == 'posix':  # Linux
    BASE = Path('/ibira')
# %%
    
# %% =========================== FOLDERPATHS ============================== %% #
IBIRA = BASE/'lnls/beamlines/ipe'

IPE   = IBIRA/'apps/IPE'

USER  = IBIRA/'proposals'/str(proposal)
PROC  = USER/'proc'
DATA  = PROC/'DATA'

RIXS  = USER/'data/RIXS/SPE'
TIF   = USER/'data/RIXS/TIF'
XAS   = DATA/'XAS'
ASCAN = DATA/'ASCAN'
MESH  = DATA/'MESH'

BRIXS  = PROC/'brixs'

OUT    = PROC/'out'   # output
TMP    = PROC/'tmp'   # temporary files (can be deleted)
STORE  = PROC/'store' # storage (not output file, cannot be deleted)
# %%

# %% =========================== BRIXS IMPORT ============================= %% #
sys.path.insert(0, str(BRIXS))
import brixs as br

# from brixs.beamlines import IPE as _ipe
import brixs.beamlines.IPE2 as ipe
import brixs.addons.fitting
# %%

# %% ========================== FINDER SETTINGS =========================== %% #
br.finder.folderpath = TMP
# %%

# %% ======================== SUPPORT FUNCTIONS =========================== %% #
def xas(scan, verbose=True):
    """return TEY, TFY, I0, PD for XAS

    Args:
        scan (int): xas scan number
        verbose (bool, optional): turn verbose on. Default is True.

    Returns:
        TEY, TFY, I0, PD
    """
    return ipe.read(fpath=XAS/(str(scan).zfill(4) + '.dat'), verbose=verbose)

def ascan(scan, verbose=True):
    """return TEY, TFY, I0, PD for ASCAN

    Args:
        scan (int): ascan scan number
        verbose (bool, optional): turn verbose on. Default is True.

    Returns:
        TEY, TFY, I0, PD
    """
    return ipe.read(fpath=ASCAN/(str(scan).zfill(4) + '.nxs'), verbose=verbose)

def mesh(scan, verbose=True):
    """return TEY, TFY, I0, PD for MESH

    Args:
        scan (int): ascan scan number
        verbose (bool, optional): turn verbose on. Default is True.

    Returns:
        TEY, TFY, I0, PD
    """
    return ipe.read(fpath=MESH/(str(scan).zfill(4) + '.nxs'), verbose=verbose)

def rixs(scan, sbins, calib=None, norm=True, start=0, stop=None, skip=[]):
    """return rixs spectrum
    
    PhotonEvents for each Image are summed up. The summed up PhotonEvents is 
    turned into one spectrum. The spectrum for each CCD is then aligned and summed. 
    
    Args:
        scan (int): rixs scan number
        sbins (int): number of bins for converting photon events to spectrum (number 
            of points in the spectrum).
        calib (number, optional): if not None, the x axis is multipled by calib
        norm (bool, optional): if True, spectrum is divided by the exposure time,
            number of images, and number of bins (sbins).
        start, stop, skip (list, optional): For RIXS only. Start and stop are 
            the indexes for the first and last image to sum (inclusive). 
            Default start is 0 and the default for stop is the None, which
            will get up to the last image available. skip should be a list with
            image number indexes to not read (skip). Default is an empty list [].
    
    Returns:
        Spectrum
    """
    return ipe.process(RIXS/str(scan).zfill(4), sbins=sbins, calib=calib, norm=norm, start=start, stop=stop, skip=skip)

def sequence(scans, sbins, calib=True, norm=True):
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
        sbins (int): number of bins for converting photon events to spectrum (number 
            of points in the spectrum).
        calib (number, optional): if not None, the x axis is multipled by calib
        norm (bool, optional): if True, spectrum is divided by the exposure time,
            number of images, and number of bins (sbins).
        
    Returns:
        Spectra
    """
    ss = br.Spectra()
    for scan in scans:
        if isinstance(scan, Iterable):
            _ss = br.Spectra()
            for _scan in scan:
                _ss.append(read(scan, sbins=sbins, calib=calib, norm=norm))
            # fix attrs
            for attr in _ss[0].get_attrs():
                try:
                    _ss.create_attr_from_spectra(attr)        

                    if attr.endswith('_min'):
                        _ss.__setattr__(attr, min(_ss.__getattr__(attr)))
                    elif attr.endswith('_max') or attr.endswith('_sigma'):
                        _ss.__setattr__(attr, max(_ss.__getattr__(attr)))
                    else:
                        _ss.__setattr__(attr, np.mean(_ss.__getattr__(attr)))
                except:
                    pass
            ss.append(_ss.align().interp().calculate_sum())

        else:
            ss.append(read(scan, sbins=sbins, calib=calib, norm=norm))    

    #########
    # attrs #
    #########
    for attr in ss[0].get_attrs():
        try:
            ss.create_attr_from_spectra(attr)
        except:
            pass

    return ss

def verify(scan, sbins, calib=None, norm=True, **kwargs):
    """open a figure with step-by-step rixs data processing

    Args:
        scan (int): rixs scan number
        sbins (int): number of bins for converting photon events to spectrum (
            number of points in the spectrum).
        calib (number, optional): if not None, the x axis is multipled by calib
        norm (bool, optional): if True, spectrum is divided by the exposure time,
            number of images, and number of bins (sbins).
        **kwargs are passed to the scatter plot that plots photon events.
        
    Note:
        Use the argument s=10, to increase the marker size of photon evets plots.
    
    Returns:
        dict {'pe1':pe1, 'pe2':pe2, 'pes1':pes1, 'pes2':pes2, 'ss1':ss1, 'ss2':ss2, 's':s}
    """
    return ipe.verify(RIXS/str(scan).zfill(4), sbins=sbins, calib=calib, norm=norm, **kwargs)
# %%

# %% ================================ PRINT =============================== %% #
dtext  = '====== IPE beamline default imports, variables, and functions =========' + '\n'
dtext += '\n'
dtext += '1) Standard imports' + '\n'
dtext += 'import matplotlib.pyplot as plt' + '\n'
dtext += 'from pathlib import Path' + '\n'
dtext += 'import numpy as np' + '\n'
dtext += 'import sys' + '\n'
dtext += 'import os' + '\n'
dtext += '\n'
dtext += '2) Brixs imports' + '\n'
dtext += 'import brixs as br' + '\n'
dtext += 'from brixs.beamlines import IPE as _ipe' + '\n'
dtext += 'import brixs.beamlines.IPE2 as ipe' + '\n'
dtext += 'import brixs.addons.fitting' + '\n'
dtext += '\n'
dtext += '3) Proposal number' + '\n'
dtext += f'proposal = {proposal}' + '\n'
dtext += '\n'
dtext += '4) Top paths' + '\n'
dtext += f'BASE  = {BASE}' + '\n'
dtext += f'IBIRA = {IBIRA}' + '\n'
dtext += f'IPE   = {IPE}' + '\n'
dtext += f'USER  = {USER}' + '\n'
dtext += f'PROC  = {PROC}' + '\n'
dtext += f'DATA  = {DATA}' + '\n'
dtext += '\n'
dtext += '4) Raw data paths' + '\n'
dtext += f'RIXS  = {RIXS}' + '\n'
dtext += f'TIF   = {TIF}' + '\n'
dtext += f'XAS   = {XAS}' + '\n'
dtext += f'ASCAN = {ASCAN}' + '\n'
dtext += f'MESH  = {MESH}' + '\n'
dtext += '\n'
dtext += '5) brixs paths' + '\n'
dtext += f'BRIXS = {BRIXS}' + '\n'
dtext += 'br.__file = ' + str((br.__file__)) + '\n'
dtext += '\n'
dtext += '6) analysis paths' + '\n'
dtext += f'OUT     = {OUT}' + '\n'
dtext += f'TMP     = {TMP}' + '\n'
dtext += f'STORE   = {STORE}' + '\n'
dtext += '\n'
dtext += '7) functions' + '\n'
dtext += f'tey, tfy, i0, pd = xas(scan)'   + '\n'
dtext += f'tey, tfy, i0, pd = ascan(scan)' + '\n'
dtext += f'tey, tfy, i0, pd = mesh(scan)'  + '\n'
dtext += f's    = rixs(scan, sbins, calib=None, norm=True, start=0, stop=None, skip=[])' + '\n'
dtext += f'ss   = sequence(scans, sbins, calib=None, norm=True)' + '\n'
dtext += f'temp = verify(scan, sbins, calib=None, norm=True, **kwargs)' + '\n'
dtext += '\n'



