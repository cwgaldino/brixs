#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Everyday use functions for user interaction."""

import os
import sys
import warnings

try:
    import winsound
except ModuleNotFoundError:
    pass


def makeSound(duration=1, freq=440):
    """Make a sound.

    On linux (tested on ubuntu) it uses sox (``sudo apt install sox``).
    On windows it uses the python package winsound (``pip install winsound``).

    Args:
        duration (int, optional): duration in seconds.
        freq (int, optional): frequence of sound in Hertz
    """
    uname = sys.platform.lower()
    if os.name == 'nt':
        uname = 'win'
    if uname.startswith('linux'):
        uname = 'linux'

    duration = duration*1000  # milliseconds
    freq = 440  # Hz

    if uname == 'win':
        try: winsound.Beep(freq, duration)
        except: warnings.warn('Cannot generate sound.')
    else:
        try: os.system('play -nq -t alsa synth {} sine {}'.format(duration/1000, freq))
        except: warnings.warn('Cannot generate sound.')


def sayOutLoud(message):
    """Generate speech like sounds.

    On Linux (tested only on ubuntu 18.04) speech-dispatcher package is necessary
    (``sudo apt-get install -y speech-dispatcher``).

    Args:
        message (str): message to be said out loud.
    """
    os.system('spd-say "' + str(message) + '"')


def query_yes_no(question, default="yes"):
    """Ask a yes/no question and return answer.

    Note:
        It accepts many variations of yes and no as answer, like, "y", "YES", "N", ...

    Args:
        question (str): string that is presented to the user.
        default ('yes', 'no' or None): default answer if the user just hits
            <Enter>. If None, an answer is required of the user.

    Returns:
        True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True, "Y": True, "YES": True, "YE": True,
             "no": False, "n": False, "No":True, "NO":True, "N":True}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt + '\n')
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")
