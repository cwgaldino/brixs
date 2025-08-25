#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> query and copy2clipboard"""

# %% ------------------------- Standard Imports --------------------------- %% #
import subprocess
import platform
import time
import sys

# %% -------------------------- operating system ------------------------ %% #
def operating_system():
    """Return string with name of operating system (windows, linux, or mac)."""
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: filemanip

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

# %% ============================== query ================================= %% #
def query(question, default="yes"):
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
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: filemanip

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
                             "('y' or 'n').\n")
    return
# %%

# %% =============================== time ================================= %% #
def start_time():
    """returns tuple with (perf counter, process time) in ns

    perf counter will run like a 'stopwatch' (almost like real time measurement).
    
    process time will return the time spent by the computer for the 
    current process. This way, process time does not include:
     sleep-time, OS-scheduling artifacts, time spent on 'external resources', 
     'system call' ... In general, if and only if, your function is entirely 
    CPU bound, does not access 
    external resources, nor cause any system call whatsoever, 
    process_time is the more accurate and better choice.
    
    Usage:

        >>> start_time = br.start_time()
        >>> run......
        >>> print(br.stop_time(start_time))
    
    Returns:
        tuple (perf counter, process time) in ns
    """
    return (time.perf_counter_ns(), time.process_time_ns())

def stop_time(start_time):
    """returns (perf counter, process time) subtracted by the start time in ns

    perf counter will run like a 'stopwatch' (almost like real time measurement).
    
    process time will return the time spent by the computer for the 
    current process. This way, process time does not include:
     sleep-time, OS-scheduling artifacts, time spent on 'external resources', 
     'system call' ... In general, if and only if, your function is entirely 
    CPU bound, does not access 
    external resources, nor cause any system call whatsoever, 
    process_time is the more accurate and better choice.
    
    Usage:

        >>> start_time = br.start_time()
        >>> run......
        >>> print(br.stop_time(start_time))
    
    Returns:
        tuple (elapsed perf counter, elapsed process time) in ns
    """
    return (time.perf_counter_ns() - start_time[0], time.process_time_ns() - start_time[1])

def stop_time_seconds(start_time):
    """returns (perf counter, process time) subtracted by the start time in seconds

    perf counter will run like a 'stopwatch' (almost like real time measurement).
    
    process time will return the time spent by the computer for the 
    current process. This way, process time does not include:
     sleep-time, OS-scheduling artifacts, time spent on 'external resources', 
     'system call' ... In general, if and only if, your function is entirely 
    CPU bound, does not access 
    external resources, nor cause any system call whatsoever, 
    process_time is the more accurate and better choice.
    
    Usage:

        >>> start_time = br.start_time()
        >>> run......
        >>> print(br.stop_time_seconds(start_time))
    
    Returns:
        tuple (elapsed perf counter, elapsed process time) in seconds
    """
    return (round((time.perf_counter_ns() - start_time[0])/1e9, 2), round((time.process_time_ns() - start_time[1])/1e9, 2))

# %% =========================== clipboard ================================ %% #
def copy2clipboard(txt):
    """Copy text to clipboard.

    on linux it uses ``xclip`` package (``sudo apt install xclip``).
    """
    # backpack developers note!!!!
    # if this function changes, it needs to be copied to these files: figmanip

    if is_windows:
        # cmd='echo ' + txt.strip() + ' | clip'
        cmd=f'echo|set /p={txt.strip()}| clip'
        # cmd='echo ' + txt.strip() + '| Set-Clipboard -Value {$_.Trim()}'
        subprocess.check_call(cmd, shell=True)
    elif is_linux:
        # p = subprocess.Popen(['xsel','-bi'], stdin=subprocess.PIPE)
        # p.communicate(input=bytes(txt.strip()).encode())
        # p.communicate(input=bytes(txt.strip(), encoding='utf-8'))
        # p.communicate(input=txt.strip())
        with subprocess.Popen(['xclip','-selection', 'clipboard'], stdin=subprocess.PIPE) as pipe:
            pipe.communicate(input=txt.strip().encode('utf-8'))
    elif is_mac:
        cmd='echo '+ txt.strip() + ' | pbcopy'
        subprocess.check_call(cmd, shell=True)
    return

def svg2clipboard(filepath):
    """Copy svg pictures to clipboard.

    Warning:
        Only implemented on Linux.

    On linux it uses ``xsel`` package (``sudo apt-get install -y xsel``)."""
    if is_windows:
        raise NotImplementedError('This function is not implemented on Windows yet.')
    elif is_linux:
        p = subprocess.Popen([f'xclip -selection clipboard -t image/svg+xml -i {filepath}'], shell=True)  # ctrl+V
        p = subprocess.Popen([f'echo {filepath}  |xclip -in -selection primary -target text/uri-list'], shell=True)  # ctrl+V
    elif is_mac:
        raise NotImplementedError('This function is not implemented on Mac yet.')
    return
