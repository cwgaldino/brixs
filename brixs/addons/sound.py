#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""functions to produce sound"""

# %% ------------------------- Standard Imports --------------------------- %% #
import os

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
# %%

# %% ============================== sound ================================= %% #
if is_windows:
    import winsound

def make_sound(duration=1, freq=440):
    """Make a sound.

    On linux (tested on ubuntu) it uses sox (``sudo apt install sox``).
    On windows it uses the python package winsound (``pip install winsound``).
    On mac it will raise a error (not implemented yet).

    Args:
        duration (int, optional): duration in seconds.
        freq (int, optional): frequency of sound in Hertz
    """
    if is_windows:
        # assert winsoundok, 'make_sound() cannot generate sound\nError: python package `winsound` not found\nmaybe install it via ``pip install winsound``' 
        duration = duration*1000  # milliseconds
        if type(duration) == int:
            pass
        else:
            duration = int(duration)
        if duration <= 0:
            raise ValueError(f'Duration cannot be negative nor zero.\nduration={duration}')
        winsound.Beep(freq, duration)
    elif is_linux:
        os.system('play -nq -t alsa synth {} sine {}'.format(duration, freq))
    elif is_mac:
        raise NotImplementedError('make_sound() not supported on mac yet.')

def say(message):
    """Text-to-speech function.

    On Linux ``speech-dispatcher`` package is necessary
    (``sudo apt-get install -y speech-dispatcher``).

    On windows it uses ``wsay`` package
    (for more information see https://github.com/p-groarke/wsay).

    On mac it will raise a error (not implemented yet).

    Args:
        message (str): message to be said out loud.

    Return:
        None
    """
    if is_windows:
        try:
            os.system('wsay "' + str(message) + '"')
        except:
            pass
    elif is_linux:
        try:
            os.system('spd-say "' + str(message) + '"')
        except:
            pass
    elif is_mac:
        raise NotImplemented('This function is not implemented on mac yet')
    return
# %%
