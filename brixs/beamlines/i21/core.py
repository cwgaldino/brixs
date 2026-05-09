# %% ========================== Standard Imports ========================= %% #
from pathlib import Path
import numpy as np

# %% =============================== brixs =============================== %% #
import brixs as br

# %% ========================== Special Imports ========================== %% #
try:
    import h5py
except:
    pass

# %% ============================= settings ============================= %% #
from .settings import _settings
settings = _settings()
from .support import _str2datetime

# %% ============================== scanlist ============================= %% #
def scanlist(folderpath='auto'):
    """Returns list with all available scan numbers inside folderpath
    
    Args:
        folderpath (str or Path, optional): folderpath where .nxs files can be 
            found. If 'auto' it uses the FOLDERPATH defined in i21.settings.

    Returns:
        list
    """
    if folderpath == 'auto': folderpath = settings.FOLDERPATH

    # return [int(_.name.split('.')[0].split('-')[1]) for _ in br.parsed_filelist(folderpath, string='*.nxs', ref=1)]
    return br.parsed_filelist(folderpath, string='*.nxs', ref=1, return_type='dict')
# %%

# %% ============================== readrixs ============================= %% #
def readrixs(scan=None, start=0, stop=None, verbose=False, folderpath='auto', prefix='auto', filepath='auto'):
    """Returns data from nexus file containing rixs data

    Args:
        scan (int, optional): scan number
        start, stop (int, optional): start and stop index number for images to 
            load. If stop is None, it will get all images. Use this for loading
            scans with large number of images. Start is INCLUSIVE. Stop is 
            EXCLUSIVE. Default is start=0 and stop=None
        verbose (bool, optional): if True, prints the name of metadata is could
            not load. Default is False.
        folderpath, prefix, filepath (str or Path, optional): use this to 
            overwrite i21.settings (FOLDERPATH, PREFIX), or use filepath to 
            directly give the path to a nexus file.

    Returns:
        dict: im (summed frames), ims (individual frames)
    """
    # assert scan is int
    # assert start is int
    # assert stop is int or None
    # assert verbose is bool
    # assert filepath is exists and points to a file

    # get filepath
    if prefix     == 'auto': prefix     = settings.PREFIX
    if folderpath == 'auto': folderpath = settings.FOLDERPATH
    if filepath   == 'auto': filepath   = folderpath/ (prefix + str(scan) + '.nxs')       

    # open nexus file
    with h5py.File(Path(filepath), 'r', libver='latest', swmr=True) as f:
        # check if data comes from andor or maranaX
        detector = None
        if "/entry/andor" in f:
            detector = 'andor'
        elif "/entry/maranax" in f:
            detector = 'maranax'
        else:
            raise KeyError('Cannot find which rixs detector was used (Andor, MaranaX)')

        # check number of images
        expected_number_of_images    = f['entry/scan_shape'][()]
        number_of_images_in_the_file = len(f['entry'][detector]['data'])
        isfinished = bool(f['/entry/diamond_scan/scan_finished'][()])
        # if isfinished:
        #     expected_number_of_images = number_of_images_in_the_file
                          
        if stop is None: stop = number_of_images_in_the_file
        assert stop >= start, 'stop must be higher than start'
        if stop-start > settings.MAX_IMAGES_TO_LOAD_AT_ONCE:
            _text = 'This scan has too many images and it might crash the kernel\n' +\
                    'Please, load slices of this scan by using `start` and `stop` arguments\n' +\
                    'Or, use `settings.MAX_IMAGES_TO_LOAD_AT_ONCE` to increase the maximum allowed number of images to be loaded at the same time'
            raise ValueError(_text)

        # checks if loaded images represents the full scan
        # full = False
        # if number_of_images_in_the_file == expected_number_of_images:
        #     full = True

        # get data
        ims = br.Dummy()
        for data in f['entry'][detector]['data'][start:stop]:
            ims.append(br.Image(data=data))

        # sum images
        data = np.zeros(ims[0].shape)
        for _im in ims:
            data += _im.data
        im = br.Image(data)

        # get attrs
        metadata = {}
        if True:
            for name in settings.METADATA['rixs']:
                key     = settings.METADATA['rixs'][name][0]
                proc    = settings.METADATA['rixs'][name][1]
                address = settings.METADATA['rixs'][name][2]                
                if key not in metadata:
                    metadata[key] = {}
                try: 
                    if proc == 'string':         _value = f[address][()].decode("utf-8")
                    if proc == 'int':            _value = int(f[address][()])
                    if proc == 'number':         _value = f[address][()]
                    if proc.startswith('round'): _value = round(f[address][()], int(proc.split('round')[1]))
                    if proc == 'datetime':       _value = _str2datetime(f[address][()].decode("utf-8"))
                    if proc == 'bool':           _value = f[address][()][0] == 1
                except Exception as e:
                    if verbose: print(e)
                    _value = None
                metadata[key][name] = _value
        im.metadata   = metadata
        im.isfinished = isfinished
        im.nexpected  = expected_number_of_images
        im.scan = scan
        
        return dict(im=im, ims=ims)
# %%

# %% ============================== readxas ============================== %% #
def readxas(scan=None, verbose=False, folderpath='auto', prefix='auto', filepath='auto'):
    """Returns data from nexus file containing rixs data

    Args:
        scan (int, optional): scan number
        verbose (bool, optional): if True, prints the name of metadata is could
            not load. Default is False.
        folderpath, prefix, filepath (str or Path, optional): use this to 
            overwrite i21.settings (FOLDERPATH, PREFIX), or use filepath to 
            directly give the path to a nexus file.

    Returns:
        dict: diff1, fy2, draincurrent, i0
    """
    # assert scan is int
    # assert verbose is bool
    # assert filepath is exists and points to a file

    # get filepath
    if prefix     == 'auto': prefix     = settings.PREFIX
    if folderpath == 'auto': folderpath = settings.FOLDERPATH
    if filepath   == 'auto': filepath   = folderpath/ (prefix + str(scan) + '.nxs')       

    # open nexus file
    with h5py.File(Path(filepath), 'r', libver='latest', swmr=True) as f:
        
        if 'entry/diff1_c' not in f and 'entry/fy2_c' not in f and 'entry/draincurrent_c' not in f:
            raise KeyError('This scan does not seem to be a XAS scan')
        
        # data
        data = {}
        try: 
            data['diff1'] = br.Spectrum(x=f[f'entry/diff1_c/energy'][()], y=f['entry/diff1_c/data'][()])
        except Exception as e:
            if verbose: print(e)
        try: 
            data['fy2'] = br.Spectrum(x=f['entry/fy2_c/energy'][()], y=f['entry/fy2_c/data'][()])
        except Exception as e:
            if verbose: print(e) 
        try:           
            data['tey'] = br.Spectrum(x=f['entry/draincurrent_c/energy'][()], y=f['entry/draincurrent_c/data'][()])
        except Exception as e:
            if verbose: print(e)
        try: 
            data['i0'] = br.Spectrum(x=f['entry/m4c1_c/energy'][()], y=f['entry/m4c1_c/data'][()])
        except Exception as e:
            if verbose: print(e)
        
        # get attrs
        metadata = {}
        if True:
            for name in settings.METADATA['xas']:
                key     = settings.METADATA['xas'][name][0]
                proc    = settings.METADATA['xas'][name][1]
                address = settings.METADATA['xas'][name][2]                
                if key not in metadata:
                    metadata[key] = {}
                try: 
                    if proc == 'string':         _value = f[address][()].decode("utf-8")
                    if proc == 'int':            _value = int(f[address][()])
                    if proc == 'number':         _value = f[address][()]
                    if proc.startswith('round'): _value = round(f[address][()], int(proc.split('round')[1]))
                    if proc == 'datetime':       _value = _str2datetime(f[address][()].decode("utf-8"))
                    if proc == 'bool':           _value = f[address][()][0] == 1
                except Exception as e:
                    if verbose: print(e)
                    _value = None
                metadata[key][name] = _value
        for key in data:
            data[key].metadata = metadata
            data[key].scan = scan

    return data
# %%

# %% ============================ readlinescan =========================== %% #
def readline(scan=None, verbose=False, folderpath='auto', prefix='auto', filepath='auto'):
    """Returns data from nexus file containing rixs data

    Args:
        scan (int, optional): scan number
        verbose (bool, optional): if True, prints the name of metadata is could
            not load. Default is False.
        folderpath, prefix, filepath (str or Path, optional): use this to 
            overwrite i21.settings (FOLDERPATH, PREFIX), or use filepath to 
            directly give the path to a nexus file.

    Returns:
        dict: diff1, fy2, draincurrent, i0
    """
    # assert scan is int
    # assert verbose is bool
    # assert filepath is exists and points to a file

    # get filepath
    if prefix     == 'auto': prefix     = settings.PREFIX
    if folderpath == 'auto': folderpath = settings.FOLDERPATH
    if filepath   == 'auto': filepath   = folderpath/ (prefix + str(scan) + '.nxs')       

    # open nexus file
    with h5py.File(Path(filepath), 'r', libver='latest', swmr=True) as f:
        dets = ['diff1', 'fy2', 'draincurrent', 'm4c1']
        det_name = {'diff1': 'diff1', 'fy2': 'fy2', 'draincurrent': 'tey', 'm4c1': 'i0'}

        data = {}
        for det in dets:
            try:
                if f'entry/{det}_i' in f:
                    name = list(f[f'entry/{det}_i'].keys())
                    name.remove('{det}_i')
                    name = name[0]
                    data[det_name[det]] = br.Spectrum(x=f[f'entry/{det}_i/{name}'][()], y=f[f'entry/{det}_i/{det}_i'][()])
            except:
                pass

        # for scan y -2.9 3.1 0.1 diff1 fy2 draincurrent
        if f'entry/diff1' in f:
            name = [item for item in list(f['entry/diff1'].keys()) if item not in dets][0]
            for det in dets:
                try:
                    data[det_name[det]] = br.Spectrum(x=f[f'entry/diff1/{name}'][()], y=f[f'entry/diff1/{det}'][()])
                except:
                    pass
        
        # get attrs
        metadata = {}
        if True:
            for name in settings.METADATA['xas']:
                key     = settings.METADATA['xas'][name][0]
                proc    = settings.METADATA['xas'][name][1]
                address = settings.METADATA['xas'][name][2]                
                if key not in metadata:
                    metadata[key] = {}
                try: 
                    if proc == 'string':         _value = f[address][()].decode("utf-8")
                    if proc == 'int':            _value = int(f[address][()])
                    if proc == 'number':         _value = f[address][()]
                    if proc.startswith('round'): _value = round(f[address][()], int(proc.split('round')[1]))
                    if proc == 'datetime':       _value = _str2datetime(f[address][()].decode("utf-8"))
                    if proc == 'bool':           _value = f[address][()][0] == 1
                except Exception as e:
                    if verbose: print(e)
                    _value = None
                metadata[key][name] = _value
        for key in data:
            data[key].metadata = metadata
            data[key].scan = scan

    return data
# %%