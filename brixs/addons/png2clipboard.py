#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Function for copying png image to clipboard"""

# %% ------------------------- Standard Imports --------------------------- %% #
from io import BytesIO
import subprocess

# %% -------------------------- operating system ------------------------ %% #
import platform
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

# %% ========================= png2clipboard ============================== %% #
if is_windows:
    import win32clipboard
    from PIL import Image as _image
def png2clipboard(filepath):
    """Copy png pictures to clipboard.

    On linux (tested on ubuntu) it uses ``xsel`` package (``sudo apt-get install -y xsel``)
    On windows it uses the python package win32clipboard and PIL
    On mac it will raise a error (not implemented yet).

    """
    if is_windows:
        # assert win32clipboardok, 'png2clipboard() cannot copy figure\nError: python package `win32clipboard` not found\nmaybe install it via ``pip install win32clipboard``' 
        # assert PILok, 'png2clipboard() cannot copy figure\nError: python package `PIL` not found\nmaybe install it via ``pip install PIL``' 
        win32clipboard.OpenClipboard()
        win32clipboard.EmptyClipboard()
        image = _image.open(filepath)
        output = BytesIO()
        image.convert("RGB").save(output, "DIB")
        data = output.getvalue()
        output.close()
        win32clipboard.SetClipboardData(win32clipboard.CF_DIB, data)
        win32clipboard.CloseClipboard()
        # raise NotImplementedError('This function is not implemented on Windows yet.')
    elif is_linux:
        p = subprocess.Popen([f'xclip -selection clipboard -t image/png -i {filepath}'], shell=True)  # ctrl+V
    elif is_mac:
        raise NotImplementedError('This function is not implemented on Mac yet.')
    return
