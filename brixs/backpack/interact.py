#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Functions for user interaction."""

import os
import sys
import warnings
import subprocess
import platform

try:
    import winsound
except ModuleNotFoundError:
    pass

try:
    import smtplib
except ModuleNotFoundError:
    pass

try:
    from getpass import getpass
except ModuleNotFoundError:
    pass


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

def make_sound(duration=1, freq=440):
    """Make a sound.

    On linux (tested on ubuntu) it uses sox (``sudo apt install sox``).
    On windows it uses the python package winsound (``pip install winsound``).

    Args:
        duration (int, optional): duration in seconds.
        freq (int, optional): frequence of sound in Hertz
    """

    duration = duration*1000  # milliseconds
    freq = 440  # Hz

    if is_windows:
        try: winsound.Beep(freq, duration)
        except: warnings.warn('Cannot generate sound.')
    elif is_linux:
        try: os.system('play -nq -t alsa synth {} sine {}'.format(duration/1000, freq))
        except: warnings.warn('Cannot generate sound.')

def say(message):
    """Text-to-speech function.

    On Linux ``speech-dispatcher`` package is necessary
    (``sudo apt-get install -y speech-dispatcher``).

    On windows it uses ``wsay`` package
    (for more information see https://github.com/p-groarke/wsay).

    Args:
        message (str): message to be said out loud.
    """
    if is_windows:
        os.system('wsay "' + str(message) + '"')

    elif is_linux:
        os.system('spd-say "' + str(message) + '"')

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

def copy2clipboard(txt):
    """Copy text to clipboard.

    on linux it uses ``xsel`` package (``sudo apt-get install -y xsel``)."""
    if is_windows:
        cmd='echo ' + txt.strip() + ' | clip'
        subprocess.check_call(cmd, shell=True)
    elif is_linux:
        p = subprocess.Popen(['xsel','-bi'], stdin=subprocess.PIPE)
        p.communicate(input=bytes(txt.strip()).encode())
    elif is_mac:
        cmd='echo '+ txt.strip() + ' | pbcopy'
        subprocess.check_call(cmd, shell=True)
    return

def png2clipboard(filepath):
    """Copy png pictures to clipboard.

    Warning:
        Only implemented on Linux.

    On linux it uses ``xsel`` package (``sudo apt-get install -y xsel``)."""
    if is_windows:
        raise NotImplementedError('This function is not implemented on Windows yet.')
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
        p = Popen([f'xclip -selection clipboard -t image/svg+xml -i {filepath}'], shell=True)  # ctrl+V
        p = Popen([f'echo {filepath}  |xclip -in -selection primary -target text/uri-list'], shell=True)  # ctrl+V
    elif is_mac:
        raise NotImplementedError('This function is not implemented on Mac yet.')
    return

def send_email(recipient, subject, body, sender, password=None):
    """Send email.

    If using a gmail sender email turn ``Allow less secure apps to ON``.
    Be aware that this makes it easier for others to gain access to your account.

    Warning:
        Be aware that writting a password on a python script is a bad practice.
        It is recomended creating a development email just for usage with this
        function and never share sensitive information via this account.

    Args:
        recipient (str): email of the recipient.
        subject (str): email subject.
        body (str): email body.
        sender (str): sender email.
        password (str): password of the sender email. If ``None``, it will
            securely ask the user to input the password.

    """
    if password is None:
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
