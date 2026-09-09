from pathlib import Path
import numpy as np
import platform
import fnmatch
import time
import re
import os

# %% =============================== brixs =============================== %% #
import brixs as br
from brixs.beamlines.i21.core import get_metadata, settings
from brixs.beamlines.i21.advanced import process
settings.DARKS = []

# %% ===================================================================== %% #
# %% ========================= operating system ========================== %% #
# %% ===================================================================== %% #
def operating_system():
    """Return string with name of operating system (windows, linux, or mac)."""
    system = platform.system().lower()
    is_windows = system == 'windows'
    is_linux = system == 'linux'
    is_mac = system == 'darwin'
    if is_windows:
        return 'windows'
    elif is_linux:
        return 'linux'
    elif is_mac:
        return 'mac'
    else:
        raise ValueError('OS not recognized')
is_windows = operating_system() == 'windows'
is_linux   = operating_system() == 'linux'
is_mac     = operating_system() == 'mac'

# %% ===================================================================== %% #
# %% ============================== support ============================== %% #
# %% ===================================================================== %% #  
def filelist(dirpath='.', string='*', case_sensitive=True):
    """Returns a list with all the files containing `string` in its name.

    Note:
        List is sorted by the filename.

    Args:
        dirpath (str or pathlib.Path, optional): list with full file directory
            paths.
        string (str, optional): pattern. only filenames with this string will be considered.
            Use '*' for matching anything. Default is '*'.
        case_sensitive (bool or None, optional): Default is True. If None, 
            case_sensitive will matches paths using platform-specific casing rules.
            Typically, case-sensitive on POSIX, and case-insensitive on Windows.

    Return:
        list
    """
    dirpath = Path(dirpath)

    if '*' not in string:
        string = '*' + string + '*'

    # on linux and mac, glob is naturally case sensitive
    if (is_linux or is_mac) and case_sensitive == False:
        rule = re.compile(fnmatch.translate(string), re.IGNORECASE)
        temp = [dirpath/name for name in os.listdir(dirpath) if rule.match(name)]
    # on windows, glob is naturally case INsensitive
    elif is_windows and case_sensitive:
        temp = list(dirpath.glob(pattern=string))
        match = re.compile(fnmatch.translate(str(dirpath/string))).match
        temp  = [path for path in temp if match(str(path))]
    else:
        temp = list(dirpath.glob(pattern=string))

    # sort and return
    temp2 = [filepath.name for filepath in temp]
    return [x for _,x in sorted(zip(temp2, temp))]

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
    if verbose: print("Watching for new scans...")
    while True:
        # get last scan number from data directory    
        fl = filelist(folderpath, string=f'*{extension}')
        scan_numbers_from_directory = np.sort([filename2scan(_.name) for _ in fl])
        if len(scan_numbers_from_directory) == 0:
            last_scan_from_directory = 0
        else:
            last_scan_from_directory = max(scan_numbers_from_directory)

        # check if there are new scans in the directory that are not processed
        # if so process scan
        if last_scan_from_directory > last_processed_scan:
            # if verbose: print("New scan detected, updating google sheet...")
            for scan_index in np.where(scan_numbers_from_directory > last_processed_scan)[0]:
                scan = scan_numbers_from_directory[scan_index]
                dark = settings.DARKS[br.index(settings.DARKS, scan, closest=True)]
                try:
                    command = get_metadata(scan)['general']['command']
                    if command.startswith('scan ds'):
                        if verbose: print(f'processing scan {scan}')
                        _ = process(scan=scan, dark=dark)
                except Exception as e:
                    if verbose: print(f"Error processing scan {scan}")
                    if verbose: print(e)
                last_processed_scan = scan

        # sleep
        time.sleep(float(refresh_time) - ((time.monotonic() - starttime) % float(refresh_time)))
    return