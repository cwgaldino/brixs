#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Useful functions for everyday use ---> use/machine interaction"""

# %% ------------------------- Standard Imports --------------------------- %% #
from io import BytesIO
import subprocess
import platform
import sys
import os

# %% -------------------------- operating system ------------------------ %% #
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

# %% ============================== email ================================= %% #
smtplibok = False
getpassok = False
try:
    import smtplib
    smtplibok = True
except ModuleNotFoundError:
    pass
try:
    from getpass import getpass
    getpassok = True
except ModuleNotFoundError:
    pass
def send_email(recipient, subject, body, sender, password=None):
    """Send email.

    If using a gmail sender email turn ``Allow less secure apps`` to ON.
    Be aware that this makes it easier for others to gain access to your account.

    Warning:
        Be aware that writing a password on a python script is a bad practice.
        It is recommended creating a development email just for usage with this
        function and never share sensitive information via this account.

    Args:
        recipient (str): email of the recipient.
        subject (str): email subject.
        body (str): email body.
        sender (str): sender email.
        password (str, optional): password of the sender email. If ``None``, it will
            securely ask the user to input the password.

    Returns:
        None
    """
    assert smtplibok, 'send_email() cannot send email\nError: python package `smtplib` not found\nmaybe install it via ``pip install smtplib``' 
    if password is None:
        assert getpassok, 'send_email() cannot ask for password securely\nError: python package `getpass` not found\nmaybe install it via ``pip install getpass``' 
        password = getpass()
    server = smtplib.SMTP_SSL('smtp.gmail.com', 465)
    server.ehlo()
    server.login(sender, password)

    # Enter the headers of the email
    headers = '\r\n'.join(['from: ' + sender
                           , 'subject: ' + subject
                           , 'to: ' + recipient
                           , 'mime-version: 1.0'
                           , 'content-type: text/html'
                           ])
    # Enter the text of the body of the email
    body_of_email = body

    # Tie the headers and body together into the email's content
    content = headers + '\r\n\r\n' + body_of_email

    server.sendmail(sender, recipient, content)
    server.close()
    return
# %%

# %% ============================== sound ================================= %% #
winsoundok = False
if is_windows:
    try:
        import winsound
        winsoundok = True
    except ModuleNotFoundError:
        pass
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
        assert winsoundok, 'make_sound() cannot generate sound\nError: python package `winsound` not found\nmaybe install it via ``pip install winsound``' 
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

# %% ============================= objects ================================ %% #
def get_attributes(object):
    """returns a dictionary with attributes (name and value) of an object."""
    return {name: attr for name, attr in object.__dict__.items()
            if not name.startswith("__")
            and not callable(attr)
            and not type(attr) is staticmethod}

# %% =========================== clipboard ================================ %% #
def copy2clipboard(txt):
    """Copy text to clipboard.

    on linux it uses ``xclip`` package (``sudo apt install xclip``).
    """
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

PILok = False
win32clipboardok = False
if is_windows:
    try:
        import win32clipboard
        win32clipboardok = True
    except ModuleNotFoundError:
        pass
    try:
        from PIL import Image as _image
        PILok = True
    except:
        pass
def png2clipboard(filepath):
    """Copy png pictures to clipboard.

    Warning:
        Only implemented on Linux.

    On linux it uses ``xsel`` package (``sudo apt-get install -y xsel``)."""
    if is_windows:
        assert win32clipboardok, 'png2clipboard() cannot copy figure\nError: python package `win32clipboard` not found\nmaybe install it via ``pip install win32clipboard``' 
        assert PILok, 'png2clipboard() cannot copy figure\nError: python package `PIL` not found\nmaybe install it via ``pip install PIL``' 
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

# %% ============================= Experimental =========================== %% #
IPythonok = False
try:
    from IPython import get_ipython
    IPythonok = True
except ModuleNotFoundError:
    pass
def is_notebook():
    """[EXPERIMENTAL] Return True if running from a ipython interactive terminal or jupyter notebook."""
    assert IPythonok, 'is_notebook() cannot check for notebook\nError: python package `IPython` not found\nmaybe install it via ``pip install IPython``' 
    try:
        shell = get_ipython().__class__.__name__
        if shell == 'ZMQInteractiveShell':
            return True   # Jupyter notebook or qtconsole
        elif shell == 'TerminalInteractiveShell':
            return True  # Terminal running IPython
        else:
            return False  # Other type (?)
    except NameError:
        return False      # Probably standard Python interpreter

