#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Support functions for manipulating Libreoffice Calc files.

- Install libreoffice. Follow the instructions.
In this case we are in Ubuntu 18.04 and libreoffice version 7.0.4

tar -xvzf file
cd folder
cd DEBS
sudo dpkg -i *.deb


sudo apt-get install python3-uno  >> import uno
import sys
sys.path.append('/opt/libreoffice7.0/program')
import uno

local = uno.getComponentContext()
resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)
context = resolver.resolve("uno:socket,host=localhost,port=8100;urp;StarOffice.ComponentContext")

desktop = context.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", context)
document = desktop.getCurrentComponent()

cursor = document.Text.createTextCursor()
document.Text.insertString(cursor, "This text is being added to openoffice using python and uno package.", 0)

TO DO:
    - method to delete rows and cols
    - method in sheet: get sheet name
    - conditional formatting is not wroking (it seems)

Used methods from unotools.component.calc.Calc:
    # 'set_rows_str', OK---
    # 'set_columns_str', OK---
    # 'set_rows_value', OK---
    # 'set_columns_value', OK---
    # 'set_rows_formula', OK---
    # 'set_columns_formula', OK---
    # 'get_cell_by_position', OK---
    # 'get_cell_range_by_position', OK---

Not used:
    # 'set_rows_cell_data',
    # 'set_columns_cell_data',
    # 'set_rows',
    # 'set_columns',
    # 'get_cell_range_by_name',
    # 'charts',
    # 'get_charts_count',
    # 'add_charts_new_by_name',
    # 'get_chart_by_index',
    # 'get_chart_by_name',

"""

# # standard imports
# import numpy as np
# from pathlib import Path
# import os
# import inspect
# import psutil
# import signal
# import subprocess
# import time
# import warnings
# import re
# import copy
# from collections.abc import Iterable
#
# # import uno
# import sys
# # from unotools import Socket, connect
# # from unotools.component.calc import Calc
# # from unotools.unohelper import convert_path_to_url
#
# from backpack.intermanip import query_yes_no
#
# # %%
# import uno
# local = uno.getComponentContext()
# resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)
# context = resolver.resolve("uno:socket,host=localhost,port=8100;urp;StarOffice.ComponentContext")
# desktop = context.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", context)
# document = desktop.getCurrentComponent()
# document.getURL()
# sheets = document.getSheets()
# sheet = sheets.getByName('Sheet1')
# # # %% WORKS =============
# #
# # import subprocess
# # cmd = '/opt/libreoffice7.0/program/python'
# #
# # proc = subprocess.Popen([cmd, "-c", 'print(5+5)'], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# # stdout, stderr = proc.communicate()
# # print(stdout.decode('utf-8'))
# #
# #
# #
# # # %% WORKS =============
# # import subprocess
# # cmd = '/opt/libreoffice7.0/program/python'
# #
# # proc = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# # stdout, stderr = proc.communicate(input=b'print(5+3)')
# # print(stdout.decode('utf-8'))
#
# # %% WORKS =============
# def start_python(python_exe):
#     return subprocess.Popen([python_exe, '-q', '-i'],
#                             stdin=subprocess.PIPE,
#                             stdout=subprocess.PIPE,
#                             stderr=subprocess.PIPE)
#
# def read(process):
#     return process.stdout.readline().decode("utf-8").strip()
#
# def write(process, message):
#     process.stdin.write(f"{message.strip()}\n".encode("utf-8"))
#     process.stdin.flush()
#
#
# def terminate(process):
#     process.stdin.close()
#     process.terminate()
#     process.wait(timeout=0.2)
#
# # %%
# import subprocess
# cmd = '/opt/libreoffice7.0/program/python'
# output_buffer = sys.stdout
# # output_buffer = open('file', 'wb')
# s=subprocess.Popen([cmd, '-q', '-i'],
#                         stdin=subprocess.PIPE,
#                         stdout=output_buffer,
#                         stderr=subprocess.PIPE)
# a = write(s, '2')
# a
# write(s, 'print(5+1)')
# write(s, '2+')
# write(s, "print('aaaafff')")
#
#
# print("Input: ", end="", file=output_buffer, flush=True)
# read(2)
#
# cmd = '/opt/libreoffice7.0/program/python'
# p = start_python(cmd)
# p.stdin.open()
# write(p, "a = [1,2,3]")
# write(p, "print(a)")
# read(p)
#
#
# while True:
#     line = p.stdout.readline()
#     if not line:
#         break
#     print(line)
#
#
# p.stdout.read().decode("utf-8").strip()
# b = read(p)
# b
# write(p, "print(5+)")
# p.stdout.readline()
# write(p, "2")
# message = "print(a)"
#
# out, err = p.communicate(f"{message.strip()}\n".encode("utf-8"))
# for line in p.stdout:
#     print(line)
#
#
#
#
# # %% subprocess vs pexpect =====================================================
# import subprocess
# import pexpect
# import time
#
# # open libreoffice
# libreoffice_folder = '/opt/libreoffice7.0/program'
# port = 8100
# p0 = subprocess.Popen([f"{libreoffice_folder}/soffice --nodefault --norestore --nologo --accept='socket,host=localhost,port={port};urp;'"], shell=True, close_fds=True)
#
# # connect via subprocess
# cmd = f'{libreoffice_folder}/python'
# p = start_python(cmd)
# write(p, 'import sys')
# read(p)
# write(p, f'sys.stderr=open("/home/galdino/.libremanip/errorlog.text", "w")')
# write(p, '2+222+')
# f = open("/home/galdino/.libremanip/errorlog.text", "a")
# f.write('='*20)
# f.close()
# write(p, 'import uno')
# write(p, 'local = uno.getComponentContext()')
# write(p, 'resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)')
# write(p, f'context = resolver.resolve("uno:socket,host=localhost,port={port};urp;StarOffice.ComponentContext")')
# write(p, 'desktop = context.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", context)')
# write(p, 'document = desktop.getCurrentComponent()')
# start_time = time.time()
# write(p, 'cursor = document.Text.createTextCursor()')
# write(p, 'document.Text.insertString(cursor, "This text is being added to openoffice using python and uno package.", 0)')
# print("--- %s seconds ---" % (time.time() - start_time))
#
# # connect via subprocess - connect every time (~25 times slower than previous method)
# start_time = time.time()
# p2 = subprocess.Popen([f"{cmd}", '-q', '/home/galdino/github/py-backpack/backpack/test.py'])#, shell=True, close_fds=True)
# print("--- %s seconds ---" % (time.time() - start_time))
#
# # connect via pexpect
# cmd = f'{libreoffice_folder}/python'
# c = pexpect.spawn(cmd)
# c.write('import uno\n')
# c.interact(escape_character='$')
# c.write('local = uno.getComponentContext()\n')
# c.write('resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)\n')
# c.write(f'context = resolver.resolve("uno:socket,host=localhost,port={port};urp;StarOffice.ComponentContext")\n')
# c.write('desktop = context.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", context)\n')
# c.write('document = desktop.getCurrentComponent()\n')
# c.write('cursor = document.Text.createTextCursor()\n')
# c.write('document.Text.insertString(cursor, "This text is being added to openoffice using python and uno package.", 0)\n')
#
# # %% ===========================================================================
#
# import uno
# local = uno.getComponentContext()
# resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)
# context = resolver.resolve("uno:socket,host=localhost,port=8100;urp;StarOffice.ComponentContext")
#
# desktop = context.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", context)
# document = desktop.getCurrentComponent()
#
# cursor = document.Text.createTextCursor()
# document.Text.insertString(cursor, "This text is being added to openoffice using python and uno package.", 0)
#
# def redirect_to_file(text):
#     original = sys.stdout
#     sys.stdout = open('/home/galdino/Desktop/redirect.txt', 'w')
#     print('This is your redirected text:')
#     print(text)
#     sys.stdout = original
# redirect_to_file(2+)
# sys.stderr = open('/home/galdino/Desktop/errorlog.txt', 'w')
# with open('/home/galdino/Desktop/errorlog.txt', 'w')
#
#
# with redirect_stdout(f):
#     print('foobar')
#     print(12)
# print('Got stdout: "{0}"'.format(f.getvalue()))
#
# cmd = f'{libreoffice_folder}/python'
# c = pexpect.spawn(cmd)
# c.sendline('print(5+1)')
# c.expect('\n')
# c.before
# c.after
# c.interact(escape_character='$')
# c.read()
# # %%
# import pexpect
#
# cmd = '/opt/libreoffice7.0/program/python'
# c = pexpect.spawn(cmd)
#
# f = open('file', 'wb')
# c.logfile = f
#
# c.write('a = [1,2,3, 4]\n')
# c.write('print(a)\n')
# c.readline()
#
#
#
# c.before
# c.after
# # %%
# from __future__ import absolute_import
# from __future__ import print_function
# from __future__ import unicode_literals
# import pexpect
# import sys
# cmd = '/opt/libreoffice7.0/program/python'
# # c = pexpect.spawnu(cmd)#, encoding='')
# c = pexpect.spawn(cmd, encoding='utf-8')
#
# f = open('file', 'wb')
# c.logfile = sys.stdout
# # f_send = open('send', 'wb')
# # c.logfile_send = f_send
# # f_read = open('read', 'wb')
# # c.logfile_read = f_read
#
#
# c.write('a = [1,2,3, 4]\n')
# c.write('print(a)\n')
#
# c.interact(escape_character='$')
#
#
# c.write('a = [1,2,3]\n'.encode('utf-8'))
#
#
# c.expect('>>>')
# c.expect('\n')
#
# c.write('a = [1,2,3]\n')
# c.before
# c.after
# c.write('print(a)\n')
# c.read()
# '\x1d'.encode('utf8')
# c.interact(escape_character='$')
#
#
# import unicodedata
# unicodedata.normalize('NFKD', U+001B).encode('ascii', 'ignore')
#
# # %%
# process = start(cmd)
# write(process, "print(5+1)")
# print(read(process))
# terminate(process)
#
#
# cmd = '/opt/libreoffice7.0/program/python'
# p = start_python(cmd)
# write(p, 'local = uno.getComponentContext()')
# write(p, 'resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)')
# read(p)
# write(p, "print(5+1)")
# read(p)
# write(p, "print(5+)")
# p.stdout.read()
# p.stdout.readline()
# a = p.stderr.read()
# write(p, 'context = resolver.resolve("uno:socket,host=localhost,port=8100;urp;StarOffice.ComponentContext")')
# write(p, '')
# write(p, '')
# write(p, '')
# write(p, '')
#
# # %%
# import sys
# sys.path.append('/opt/libreoffice7.0/program')
# import uno
#
# local = uno.getComponentContext()
# resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)
# context = resolver.resolve("uno:socket,host=localhost,port=8100;urp;StarOffice.ComponentContext")
#
# desktop = context.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", context)
# document = desktop.getCurrentComponent()
#
# cursor = document.Text.createTextCursor()
# document.Text.insertString(cursor, "This text is being added to openoffice using python and uno package.", 0)
# document.get_sheets_count()


# %%
# standard imports
import numpy as np
from pathlib import Path
import os
import inspect
import psutil
import signal
import subprocess
import time
import warnings
import re
import copy
from collections.abc import Iterable
import sys

from backpack.intermanip import query_yes_no

def start_python(python_exe, output_filepath):
    with open(output_filepath, 'w') as f:
        return subprocess.Popen([python_exe, '-q', '-i'],
                                stdin=subprocess.PIPE,
                                # stdout=subprocess.PIPE,
                                stdout=f,
                                stderr=subprocess.PIPE)

def read(file):
    return process.stdout.readline().decode("utf-8").strip()

def write(process, message):
    process.stdin.write(f"{message.strip()}\n".encode("utf-8"))
    process.stdin.flush()

def terminate(process):
    process.stdin.close()
    process.terminate()
    process.wait(timeout=0.2)

def check_error_file(old_error_log, error_filepath):
    f=open(error_filepath, 'r')
    error_log = f.read()
    f.close()

    if old_error_log == error_log:
        return (False,0, 0 )
    else:
        return (True, error_log, error_log.replace(old_error_log, '', 1))

# %%
class soffice():

    def __init__(self, port=8100, norestore=False, libreoffice_folder=None, config_folder=None, default_waiting_time=0.01, default_timeout=10):
        self.pid_previous = self._libreoffice_pid_list()
        self.default_waiting_time = default_waiting_time
        self.default_timeout = default_timeout
        self.port = port

        # config folder
        if config_folder is None:
            self.config_folder = Path('~/.libremanip').expanduser()
        else:
            self.config_folder = Path(config_folder).expanduser()
        if self.config_folder.is_dir():
            pass
        else:
            self.config_folder.mkdir()
            print(f'Configuration folder not found.')
            print(f'New folder created: {self.config_folder}')

        # set libreoffice folder
        if libreoffice_folder is not None:
            libreoffice_folder = Path(libreoffice_folder)
        else:
            libreoffice_folder = Path('/opt/libreoffice7.0/program/')

        # initialize libreoffice
        if norestore:
            self.process = subprocess.Popen([f"{libreoffice_folder/'soffice'} --nodefault --norestore --nologo --accept='socket,host=localhost,port={port};urp;'"], shell=True, close_fds=True)
        else:
            self.process = subprocess.Popen([f"{libreoffice_folder/'soffice'} --nodefault --nologo --accept='socket,host=localhost,port={port};urp;'"], shell=True, close_fds=True)
        time.sleep(0.1)


        # initialize python from libreoffice
        self.python_exe = f'{libreoffice_folder}/python'
        self.output_filepath = self.config_folder/"output.text"
        self.python = start_python(self.python_exe, self.output_filepath)
        print(f'Output log at: {self.output_filepath}')
        self.output_file = open(self.output_filepath, 'r')

        # clean error file
        self.error_filepath = str(self.config_folder/"errorlog.text")
        # with open(self.error_filepath, 'w') as f:
        #     f.truncate(0)

        # redirect stderr
        write(self.python, 'import sys')
        write(self.python, f'sys.stderr=open("{self.error_filepath}", "w")')
        print(f'Error log at: {self.error_filepath}')
        time.sleep(0.1)
        self.error_file = open(self.error_filepath, 'r')

        # imports
        dummy = self.send('import uno; print(uno)', require_answer=True)
        dummy = self.send('from com.sun.star.beans import PropertyValue; print(PropertyValue)', require_answer=True)  # used for saving files

        # initialize communication
        dummy = self.send('local = uno.getComponentContext(); print(local)', require_answer=True)
        dummy = self.send('resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local); print(resolver)', require_answer=True)
        dummy = self.send(f'context = resolver.resolve("uno:socket,host=localhost,port={port};urp;StarOffice.ComponentContext"); print(context)', require_answer=True)
        dummy = self.send('desktop = context.ServiceManager.createInstanceWithContext("com.sun.star.frame.Desktop", context); print(desktop)', require_answer=True)


        time.sleep(1)
        self.pid = self.process.pid
        self.pid_children = self._get_children_pid(self.pid)
        self.apps = []

    def send(self, message, require_answer=False, timeout=10):
        write(self.python, message=message)
        time.sleep(self.default_waiting_time)
        out = self.output_file.read()
        err = self.error_file.read()
        if require_answer:
            start_time = time.time()
            # time.sleep(self.default_waiting_time)
            while out == '' and err == '':
                out = self.output_file.read()
                err = self.error_file.read()
                time.sleep(0.1)
                if time.time() > start_time+ timeout:
                    raise TimeoutError(f'timeout: {timeout}s, message: {message}')
                    break
        if out != '':
            return out.strip()
        elif err != '':
            raise Exception(err)
        else:
            return None

    def receive(self):
        out = self.output_file.read()
        err = self.error_file.read()
        if out != '':
            return out.strip()
        elif err != '':
            raise Exception(err)
        else:
            return None

    def set_default_waiting_time(self, waiting_time):
        self.default_waiting_time = waiting_time

    def get_default_waiting_time(self):
        return self.default_waiting_time


    def openCalc(self, filepath=None, force_new_file=False):
        children_old = self._get_children_pid(self.pid)

        calcObject = calc(filepath, libreoffice=self, force_new_file=force_new_file)
        calcObject.pid = [child for child in self._get_children_pid(self.pid) if child not in children_old]

        self.pid_children += calcObject.pid
        # self.pid_children = [item for sublist in self.pid_children for item in sublist]
        self.apps.append(calcObject)

        return calcObject


    def _get_children_pid(self, pid):
        output = subprocess.check_output(["bash", "-c", f"pstree -p -n  {pid}"])
        pattern = re.compile(r"\((\d+)\)")
        pid_children = pattern.findall(output.decode('utf-8'))
        return [int(pid) for pid in pid_children]


    def _get_pid_by_name(self, string):
        """Get a list of all the PIDs of all the running process whose name contains
        string.

        Args:
            string (str): string.
        """

        listOfProcessObjects = []

        # Iterate over the all the running process
        for proc in psutil.process_iter():
           try:
               pinfo = proc.as_dict(attrs=['pid', 'name', 'create_time'])
               # Check if process name contains the given name string.
               if string.lower() in pinfo['name'].lower() :
                   listOfProcessObjects.append(pinfo)
           except (psutil.NoSuchProcess, psutil.AccessDenied , psutil.ZombieProcess):
               pass

        return [item['pid'] for item in listOfProcessObjects];


    def _libreoffice_pid_list(self):
        """Return a list of processes associated with libreoffice.

        Note:
            This function try to match the name of a process with names tipically
            related with libreoffice processes ('soffice.bin' or 'oosplash').
            Therefore, it might return processes that are not related to
            libreoffice if their name mathces with words: 'soffice'
            and 'oosplash'.

        Returns:
            list.
        """
        process_list = []
        for proc in self._get_pid_by_name('soffice'):
            process_list.append(proc)
        for proc in self._get_pid_by_name('oosplash'):
            process_list.append(proc)
        return [i for sublist in [self._get_children_pid(item) for item in process_list] for i in sublist]


    def kill_libreoffice_processes(self):
        """Kill libreoffice processes.

        Note:
            It will close ALL processes that are related to libreoffice (processes
            that have 'soffice.bin' or 'oosplash' in their name)."""

        pid_list = self._libreoffice_pid_list()
        for proc in process_list:
            os.kill(proc['pid'], signal.SIGKILL)


    def terminate(self, ask=True):

        if ask:
            if len(self.apps)>0:
                print('You appear to have opened libreoffice apps.')
                for app in self.apps:
                    if app.filepath is None:
                        msg = f'{app.object} is not associated with a file. Do you wish to discard this work?'
                    else:
                        msg = f'Progress not saved on {app.filepath} will be discarted. Do you wish to continue?'
                    if not query_yes_no(msg, 'no'):
                        return
            msg = 'All instances of libreoffice (even if not initialized via python) will be discarted and progress not saved will be lost. Do you wish to continue?'
            if not query_yes_no(msg, 'no'):
                return

        for pid in self._libreoffice_pid_list():
            try:
                os.kill(int(pid), 9)
                print(f'{int(pid)} killed.')
            except ProcessLookupError:
                print(f'{int(pid)} does not exist.')
        print('Done!')


# %%

class calc():


    def __init__(self, filepath=None, libreoffice=None, force_new_file=False):
        self.pid = []
        # self.object = None
        self.libreoffice = libreoffice
        self.default_waiting_time = self.libreoffice.default_waiting_time
        self.python = self.libreoffice.python

        if force_new_file:
            self.send("URL = 'private:factory/scalc'")
            print(f'Opening a new file.')
            self.send('document = desktop.loadComponentFromURL(URL, "_default", 0, ())')
            title = self.send("print(document.getTitle())", require_answer=True)
            # title = self.receive()
            print(f'Connected with opened file: {title}')
            self.filepath = None
        else:
            if filepath is None: # if a filepath is not defined
                print('Filepath was not given.')
                self.send("document =  desktop.getCurrentComponent()")
                # is there something open???
                isOpen = self.send("print('False') if document is None else print('True')", require_answer=True)
                if str2bool(isOpen):  # if something is open what it is? is it a calc instance
                    print('However, something is open.')
                    type = self.send("print(document.getImplementationName())", require_answer=True)
                    time.sleep(0.1)
                    if type == 'ScModelObj': # yes, it is a calc instance
                        # print('And it is a calc instance')
                        hasFilepath = self.send("print('False') if document.getURL() is None else print('True')", require_answer=True)
                        if str2bool(hasFilepath): # check if file has filepath
                            url = self.send("print(uno.fileUrlToSystemPath(document.getURL()))", require_answer=True)
                            # url = self.receive()
                            print(f'Connected with opened file: {url}')
                            self.filepath = Path(str(url))
                        else:
                            title = self.send("print(document.getTitle())", require_answer=True)
                            # title = self.receive()
                            print(f'Connected with opened file: {title}')
                            self.filepath = None
                    else: # if it is not calc, check if there is a calc opened
                        print('Searching for opened Calc instances...')
                        self.send("documents = desktop.getComponents().createEnumeration()")
                        open_new = True
                        hasMoreElements = str2bool(self.send("print(documents.hasMoreElements())", require_answer=True))
                        # hasMoreElements = str2bool(self.receive())
                        while hasMoreElements:
                            self.send("document = documents.nextElement()")
                            type = self.send("print(document.getImplementationName())", require_answer=True)
                            if type == 'ScModelObj': # yes, it is a calc instance
                                hesFilepath = self.send("print('False') if document.getURL() is None else print('True')", require_answer=True)
                                if str2bool(hesFilepath): # check if file has filepath
                                    title = self.send("print(document.getTitle())", require_answer=True)
                                    # title = self.receive()
                                    print(f'Connected with opened file: {title}')
                                    self.filepath = None
                                else:
                                    url = self.send("print(uno.fileUrlToSystemPath(document.getURL()))", require_answer=True)
                                    # url = self.receive()
                                    print(f'Connected with opened file: {url}')
                                    self.filepath = Path(str(url)[1:])
                                open_new = False
                                break
                            hasMoreElements = str2bool(self.send("print(documents.hasMoreElements())"))
                            # hasMoreElements = str2bool(self.receive())
                        if open_new: # if there is none, open a new one.
                            print(f'No opened Calc instances were found.')
                            self.send("URL = 'private:factory/scalc'")
                            print(f'Opening a new file...')
                            self.send('document = desktop.loadComponentFromURL(URL, "_default", 0, ())')
                            title = self.send("print(document.getTitle())", require_answer=True)
                            # title = self.receive()
                            print(f'Connected with opened file: {title}')
                            self.filepath = None
                else:   # if nothing is open, open a new calc instance
                    self.send("URL = 'private:factory/scalc'")
                    print(f'Nothing is opened.')
                    print(f'Opening a new file.')
                    self.send('document = desktop.loadComponentFromURL(URL, "_default", 0, ())')
                    title = self.send("print(document.getTitle())", require_answer=True)
                    # title = self.receive()
                    print(f'Connected with opened file: {title}')
                    self.filepath = None
            else: # if filepath was given
                print('Filepath was given.')
                # open file if it is not opened
                self.filepath = Path(filepath).absolute()
                self.send(f"URL = uno.systemPathToFileUrl('{self.filepath}')")
                self.send('document = desktop.loadComponentFromURL(URL, "_default", 0, ())')
                url = self.send("print(uno.fileUrlToSystemPath(document.getURL()))", require_answer=True)
                # url = self.receive()
                print(f'Connected with file: {url}')

    def send(self, message, require_answer=False):
        return self.libreoffice.send(message=message, require_answer=require_answer)

    def receive(self):
        return self.libreoffice.receive()

    def set_default_waiting_time(self, waiting_time):
        self.libreoffice.default_waiting_time = waiting_time

    def get_default_waiting_time(self):
        return self.libreoffice.default_waiting_time


    def save(self, filepath=None):
        """Save ods file.

        Note:
            If filepath have no suffix, it adds '.ods' at the end of the filename.

        Args:
            filepath (string or pathlib.Path, optional): filepath to save file.
        """
        if filepath is None and self.filepath is None:
            temporary_path = Path.cwd()/'Untitled.ods'
            if not query_yes_no(f'Filepath not defined. Wish to save at {temporary_path}?'):
                filepath = temporary_path
                return
        elif filepath is not None:
            self.filepath = Path(filepath)


        self.filepath = Path(self.filepath)

        # fix extension
        if self.filepath.suffix != '.ods':
            self.filepath = Path(self.filepath).with_suffix('.ods')

        # save
        self.send(f'URL = uno.systemPathToFileUrl("{self.filepath}")')
        # u = self.send('print(URL)', True)
        # print(f'fffff:{u}')
        self.send(f"properties = ( PropertyValue('FilterName', 0, 'writer8', 0), )")
        # u = self.send('print(properties)', True)
        # print(f'fffff:{u}')
        # time.sleep(1)
        self.send(f"document.storeAsURL(URL, properties)")
        print(f'Saved at: {self.filepath}')

    def close(self):
        """Close window."""
        self.send(f'document.close(True)')

    def terminate(self, ask=True):
        self.libreoffice.terminate(ask)

    def get_sheets_count(self):
        return int(self.send("print(len(document.getSheets()))", require_answer=True))

    def get_sheets_name(self):
        """Returns the sheets names in a tuple."""
        n = self.send("print(document.getSheets().ElementNames)", require_answer=True)
        return eval(n)



    def insert_sheets(self, name, position=None):
        """position starts from 1. If position = 1, the sheet will be the first one.
        """
        if position is None:
            position = self.get_sheets_count()+1

        if type(name) == str:
            name = [name]

        existing_names = [name2 for name2 in name if name2 in self.get_sheets_name()]
        if len(existing_names) != 0:
            raise SheetNameExistError(existing_names)

        self.object.insert_multisheets_new_by_name(name, position-1)


    def remove_sheets(self, name):
        self.remove_sheets_by_name(name)


    def remove_sheets_by_name(self, name):

        if type(name) == str:
            name = [name]

        not_existing_names = [name2 for name2 in name if name2 not in self.get_sheets_name()]
        if len(not_existing_names) != 0:
            raise SheetNameDoNotExistError(not_existing_names)

        for n in name:
            if len(self.get_sheets_name()) == 1:
                raise SheetRemoveError(n)
            self.object.remove_sheets_by_name(n)


    def remove_sheets_by_position(self, position):
        names = self.get_sheets_name()

        if position > len(names) or position < 1:
            raise IndexError('Position outside range.')

        if len(names) == 1:
            raise SheetRemoveError(names[position-1])

        self.object.remove_sheets_by_name(names[position-1])


    def get_sheet_by_name(self, name):
        return sheet(name, self)

    def get_sheets(self, name=None):
        if name is None:
            name = self.get_sheets_name()
        else:
            if type(name) == str:
                name = [name]

            if type(name) == int:
                name = [name]

        sheet_objects = []
        for n in name:
            if type(n) == str:
                sheet_objects.append(self.get_sheet_by_name(n))
            if type(n) == int:
                sheet_objects.append(self.get_sheets_by_position(n))
        if len(sheet_objects) == 1:
            return sheet_objects[0]
        else:
            return sheet_objects

    def get_sheets_by_position(self, position):
        names = self.get_sheets_name()

        if type(position) == int:
            position = [position]

        outside_range = []
        for p in position:
            if p > len(names) or p < 1:
                outside_range.append(p)
        if len(outside_range) > 0:
            raise IndexError(f'Positions {outside_range} outside range.')

        if len(position) == 1:
            return sheet(names[position[0]-1], self)
        else:
            return [sheet(names[p-1], self) for p in position]



class SheetNameExistError(Exception):

    # Constructor or Initializer
    def __init__(self, sheet_names):
        if type(sheet_names) == str:
                sheet_names = [sheet_names]
        self.names = sheet_names

    # __str__ is to print() the value
    def __str__(self):
        msg = ''
        for name in self.names:
            msg += f"'{name}' already exists\n"
        return(msg)


class SheetNameDoNotExistError(Exception):

    # Constructor or Initializer
    def __init__(self, sheet_names):
        if type(sheet_names) == str:
                sheet_names = [sheet_names]
        self.names = sheet_names

    # __str__ is to print() the value
    def __str__(self):
        msg = ''
        for name in self.names:
            msg += f"'{name}' does not exists\n"
        return(msg)


class SheetRemoveError(Exception):

    # Constructor or Initializer
    def __init__(self, sheet_names):
        if type(sheet_names) == str:
                sheet_names = [sheet_names]
        self.names = sheet_names

    # __str__ is to print() the value
    def __str__(self):
        msg = ''
        for name in self.names:
            msg += f"'{name}' cannot be removed because it is the only existing sheet.\n"
        return(msg)

# %%
class sheet():

    def __init__(self, name, calc):
        # check if sheet name exists
        #
        #
        #

        self.name = name
        self.calc = calc
        self.libreoffice = self.calc.libreoffice
        self.python = self.libreoffice.python

    def send(self, message, require_answer=False):
        return self.libreoffice.send(message=message, require_answer=require_answer)

    def receive(self):
        return self.libreoffice.receive()

    def _get_sheet(self):
        dummy = self.send('sheets = document.getSheets(); print(sheets)', require_answer=True)
        dummy = self.send(f'sheet = sheets.getByName("{self.name}"); print(sheet)', require_answer=True)

    def set_default_waiting_time(self, waiting_time):
        self.libreoffice.default_waiting_time = waiting_time

    def get_default_waiting_time(self):
        return self.libreoffice.default_waiting_time

    def get_name(self):
        self._get_sheet()
        # time.sleep(0.2)
        return self.send('print(sheet.getName())', require_answer=True)

    def set_name(self, name):
        self._get_sheet()
        self.send(f"sheet.setName('{name}')")


    def get_last_row(self):
        # higher_edited_row = len(self.object.getRowDescriptions()) + self.object.queryVisibleCells().Count -1
        # higher_edited_col = len(self.object.getColDescriptions()) + self.object.queryVisibleCells().Count -1
        #
        # row = higher_edited_row
        # if not any(v for v in self.get_row_values(row=row, col_stop=higher_edited_col)):
        #
        #
        # higher_edited_row
        self._get_sheet()
        idx = self.send("print(sheet.getRowDescriptions()[-1].split(' ')[-1])", True)
        #print(idx)
        idx = int(idx)
        # idx = int(self.send("print(idx)", True))
        # idx = int(self.receive())

        visible = False
        while visible == False:
            visible = self.send(f"print(sheet.getRows().getByIndex({idx}).IsVisible)", True)
            visible = str2bool(visible)
            idx += 1
        return idx-1

    def get_last_col(self):
        self._get_sheet()
        idx = _letter2num(self.send("print(sheet.getColumnDescriptions()[-1].split(' ')[-1])", True))
        # idx = _letter2num(self.send("print(idx)", True))
        # idx = _letter2num(self.receive())

        visible = False
        while visible == False:
            visible = self.send(f"print(sheet.getColumns().getByIndex({idx}).IsVisible)", True)
            visible = str2bool(visible)
            idx += 1
        return idx-1


    def set_col_width(self, width, col=None):
        if col is None:
            col = np.arange(0, self.get_last_col())
        else:
            col = _check_col_value(col)
        self._get_sheet()
        self.send("colsObject = sheet.getColumns()")
        for c in col:
            self.send(f"colsObject[{c}].setPropertyValue('Width', {width})")

    def get_col_width(self, col=None):
        if col is None:
            col = [0, ]
        else:
            col = _check_col_value(col)
        self._get_sheet()
        self.send("colsObject = sheet.getColumns()")

        width = []
        for c in col:
            width.append(int(self.send(f"print(colsObject[{c}].Width)", True)))
            # width.append(int(self.receive()))
        if len(col) == 1:
            return width[0]
        else:
            return width


    def set_row_height(self, height, row=None):
        if row is None:
            row = np.arange(0, self.get_last_row())
        else:
            row = _check_row_value(row)
        self._get_sheet()
        self.send("rowsObject = sheet.getRows()")
        for r in row:
            self.send(f"rowsObject[{r}].setPropertyValue('Height', {height})")

    def get_row_height(self, row=None):
        if row is None:
            row = [0, ]
        else:
            row = _check_row_value(row)
        self._get_sheet()
        self.send("rowsObject = sheet.getRows()")

        height = []
        for r in row:
            height.append(int(self.send(f"print(rowsObject[{r}].Height)", True)))
        if len(row) == 1:
            return height[0]
        else:
            return height


    def set_cell_value(self, value, row, col, format='formula'):
        """
        """
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        self._get_sheet()

        if format == 'formula':
            self.send(f"sheet.getCellByPosition({col}, {row}).setFormula('{value}')")
        elif format == 'string':
            self.send(f"sheet.getCellByPosition({col}, {row}).setString('{value}')")
        elif format == 'number':
            self.send(f"sheet.getCellByPosition({col}, {row}).setValue({value})")
        else:
            raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")

    def get_cell_value(self, row, col, format='string'):
        """
        """
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        self._get_sheet()

        if format == 'formula':
            value = self.send(f"sheet.getCellByPosition({col}, {row}).getFormula()", True)[1:-1]
            # value = self.receive()[1:-1]
        elif format == 'string':
            value = self.send(f"sheet.getCellByPosition({col}, {row}).getString()", True)[1:-1]
            # value = self.receive()[1:-1]
        elif format == 'number':
            value = float(self.send(f"sheet.getCellByPosition({col}, {row}).getValue()", True))
            # value = float(self.receive())
        else:
            raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")

        return value


    def set_cells_value(self, data, row_start=1, col_start=1, format='formula'):
        """
        if data is 1d array or list, data is placed in a row.

        formula set formulas, but also set numbers fine. Dates and time not so much because it changes the formating (if setting date and time iwth formula you might wanna format the
        cell like date or time using copy_cells to copy formatting).

        string (data) works fine with date, time and number, but formulas are set as string. Therefore, formulas do not work.

        value (data_number) works fine for numbers ONLY.

        """
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]

        try:
            row_count, col_count = np.shape(data)
            col_stop = col_start+col_count-1
            row_stop = row_start+row_count-1
        except ValueError:
            col_count = len(data)
            row_count = 0
            data = [data,]
            col_stop = col_start+col_count-1
            row_stop = row_start

        self._get_sheet()
        self.send(f"sheet_data = sheet.getCellRangeByPosition({col_start}, {row_start}, {col_stop}, {row_stop})")
        self.send(f"data = {data}")

        if format == 'formula':
            self.send("sheet_data.setFormulaArray(data)")
        elif format == 'string':
            self.send("sheet_data.setDataArray(data)")
        elif format == 'number':
            self.send("sheet_data.setDataArray(data)")
        else:
            raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")

    def get_cells_value(self, row_start=1, col_start=1, row_stop=None, col_stop=None, format='string'):
        """
        """
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]

        if row_stop is None:
            row_stop = self.get_last_row()
        if col_stop is None:
            col_stop = self.get_last_col()

        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        if col_stop < col_start:
            raise ValueError('col_start cannot be bigger than col_stop')
        if row_stop < row_start:
            raise ValueError('row_start cannot be bigger than row_stop')

        self._get_sheet()
        self.send(f"sheet_data = sheet.getCellRangeByPosition({col_start}, {row_start}, {col_stop}, {row_stop})")

        if format == 'formula':
            self.send(f"sheet_data = list(sheet_data.getFormulaArray())")
            data = self.send(f"print(sheet_data)", True)
            sheet_data = eval(data)

        elif format == 'string':
            self.send(f"sheet_data = list(sheet_data.getDataArray())")
            data = self.send(f"print(sheet_data)", True)
            sheet_data = eval(data)

        elif format == 'number':
            self.send(f"sheet_data = list(sheet_data.getDataArray())")
            data = self.send(f"print(sheet_data)", True)
            sheet_data = eval(data)

        else:
            raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")

        # transform in list
        for row_number, row_data in enumerate(sheet_data):
            sheet_data[row_number] = list(row_data)

            # if one column or one row data, transform in vector
            if col_start == col_stop:
                sheet_data[row_number] = row_data[0]
        if row_start == row_stop:
            sheet_data = sheet_data[0]

        return sheet_data


    def get_row_values(self, row, col_start=1, col_stop=None, format='string'):
        return self.get_cells_value(row_start=row, col_start=col_start, row_stop=row, col_stop=col_stop, format=format)

    def set_row_values(self, data, row, col_start=1, format='formula'):
        """data must be a 1D list."""
        try:
            row_count, col_count = np.shape(data)
            warnings.warn('Data must be a 1D list.')
        except ValueError:
            self.set_cells_value(data, row_start=row, col_start=col_start, format=format)


    def get_col_values(self, col, row_start=1, row_stop=None, format='string'):
        return self.get_cells_value(row_start=row_start, col_start=col, row_stop=row_stop, col_stop=col, format=format)

    def set_col_values(self, data, col, row_start=1, format='formula'):
        """data must be a 1D list."""
        try:
            row_count, col_count = np.shape(data)
            warnings.warn('Data must be a 1D list.')
        except ValueError:
            self.set_cells_value(transpose(data), row_start=row_start, col_start=col, format=format)






    def list_properties(self):
        return self.object._show_attributes()

    def list_cell_properties(self, filter=None, row=1, col=1):
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        cellObject = self.object.get_cell_by_position(col, row)
        if filter is None:
            # return cellObject._show_attributes()
            return cellObject.as_raw().__dir__()
        else:
            p = cellObject.as_raw().__dir__()
            matching = [s for s in p if filter.lower() in s.lower()]
            return matching

    def get_cell_property(self, property, row, col):
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        if type(property) == str:
            property = [property]

        object_cell = self.object.get_cell_by_position(col, row)
        attr = getattr(object_cell, property[0])

        for idx in range(1, len(property)):
            attr = getattr(attr, property[idx])

        # if type(attr) != int and type(attr) != float and type(attr) != bool:
        #     return attr, attr.value.__dir__()
        # else:
        #     return attr, None
        try:
            return attr, attr.value.__dir__()
        except AttributeError:
            return attr, None

    def set_cell_property(self, property, value, row, col):
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        if type(property) == str:
            property = [property]

        self.object.get_cell_by_position(col, row).setPropertyValue(property[0], value)

    def get_cells_properties(self, property, row_start, col_start, row_stop, col_stop):
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]
        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        if type(property) == str:
            property = [property]

        object_cell = self.object.get_cell_range_by_position(col_start, row_start, col_stop, row_stop)
        attr = getattr(object_cell, property[0])

        for idx in range(1, len(property)):
            attr = getattr(attr, property[idx])

        try:
            return attr, attr.value.__dir__()
        except AttributeError:
            return attr, None

    def set_cells_properties(self, property, value, row_start, col_start, row_stop, col_stop):
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]
        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        if type(property) == str:
            property = [property]

        self.object.get_cell_range_by_position(col_start, row_start, col_stop, row_stop).setPropertyValue(property[0], value)

    def get_cell_object(self, row=1, col=1):
        row = _check_row_value(row)[0]
        col = _check_col_value(col)[0]

        return self.object.get_cell_by_position(col, row)

    def get_cell_formating(self, row, col, extra=None):
        """font : bold, font name, font size, italic, color
        text: vertical justify,horizontal justify
        cell: border, color
        conditional formating: conditional formating
        """

        p = ['FormatID']
        p += ['CharWeight', 'CharFontName', 'CharHeight', 'CharPosture', 'CharColor']
        p += ['VertJustify', 'HoriJustify']
        p += ['CellBackColor', 'TableBorder', 'TableBorder2']
        p += ['ConditionalFormat']

        if extra is not None:
            p += extra

        p_list = []
        for property in p:
            obj, _ = self.get_cell_property(property, row=row, col=col)
            p_list.append(obj)

        return p_list

    def set_cell_formating(self, obj_list, row, col, extra=None):
        """font : bold, font name, font size, italic, color
        text: vertical justify,horizontal justify
        cell: border, color
        conditional formating: conditional formating
        """
        p = ['FormatID']
        p += ['CharWeight', 'CharFontName', 'CharHeight', 'CharPosture', 'CharColor']
        p += ['VertJustify', 'HoriJustify']
        p += ['CellBackColor', 'TableBorder', 'TableBorder2']
        p += ['ConditionalFormat']

        if extra is not None:
            p += extra

        for obj, property in zip(obj_list, p):
            self.set_cell_property(property, value=obj, row=row, col=col)

    def get_cells_formatting(self, row_start=1, col_start=1, row_stop=None, col_stop=None, extra=None):
        if row_stop is None:
            row_stop = self.get_last_row()
        if col_stop is None:
            col_stop = self.get_last_col()

        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        cell_formatting_list = []
        for idx, row in enumerate(range(row_start, row_stop+2)):
            cell_formatting_list.append([])
            for col in range(col_start, col_stop+2):
                p_list = self.get_cell_formating(row=row, col=col, extra=extra)
                cell_formatting_list[idx].append(p_list)
        return cell_formatting_list

    def set_cells_formatting(self, cell_formatting_list, row_start=1, col_start=1, extra=None):

        for idx, row in enumerate(range(row_start, row_start+len(cell_formatting_list))):
            for idx2, col in enumerate(range(col_start, col_start+len(cell_formatting_list[idx]))):
                self.set_cell_formating(cell_formatting_list[idx][idx2], row=row, col=col, extra=extra)

    def get_merged(self, ):
        merged_ranges = []
        #################
        start=[]
        stop = []
        #################

        ucf = self.object.getUniqueCellFormatRanges()
        for ranges in ucf:
            rgtest = ranges.getByIndex(0)
            if rgtest.getIsMerged():
                for rg in ranges:
                    oCursor = rg.getSpreadsheet().createCursorByRange(rg)
                    oCursor.collapseToMergedArea()
                    ######################
                    row_start = int(oCursor.getRowDescriptions()[0].split(' ')[-1])-1
                    row_end = int(oCursor.getRowDescriptions()[-1].split(' ')[-1])-1
                    for row in range(row_end-row_start+1):
                        col_start = backpack.libremanip._letter2num(oCursor.getColumnDescriptions()[0].split(' ')[-1])-1
                        col_end = backpack.libremanip._letter2num(oCursor.getColumnDescriptions()[-1].split(' ')[-1])-1
                        for col in range(col_end-col_start+1):
                            if oCursor.getCellByPosition(col, row).IsMerged:
                                # print(row+row_start, col+col_start)
                                start.append([row+row_start, col+col_start])
                            else:
                                stop.append([row+row_start, col+col_start])
                                # print(f'not: {row+row_start}, {col+col_start}')
                    ######################
                    addr = oCursor.getRangeAddress()

                    col_start = addr.StartColumn+1
                    row_start = addr.StartRow+1
                    col_stop = addr.EndColumn+1
                    row_stop = addr.EndRow+1
                    merged_ranges.append([row_start, col_start, row_stop, col_stop])
        return merged_ranges

    def merge(self, row_start, col_start, row_stop, col_stop):
        row_start = _check_row_value(row_start)[0]
        col_start = _check_col_value(col_start)[0]
        row_stop = _check_row_value(row_stop)[0]
        col_stop = _check_col_value(col_stop)[0]

        sheet_data = self.object.get_cell_range_by_position(col_start, row_start, col_stop, row_stop)
        sheet_data.merge(True)

# libreoffice = soffice()
# calcObject = libreoffice.openCalc(current_file=True)
# calcObject.get_sheets_count()
# calcObject.get_sheets_name()

# s = calcObject.get_sheet_by_name('Sheet1')
# s.get_name()
# s.set_name('dddd')
# s.set_name('Sheet1')
#
# s.get_last_row()
# s.get_last_col()
# s.get_cell_value(1, 1)
# s.get_cell_value(1, 1, format='number')
# s.get_cell_value(1, 1, format='formula')
#
# s.get_cell_value(3, 2)
# s.get_cell_value(3, 2, format='number')
# s.get_cell_value(3, 2, format='formula')
#
# s.set_cell_value(342, 1, 4, format='string')
# s.set_cell_value(34, 2, 4, format='number')
# s.set_cell_value('=D1+D2', 3, 4, format='formula')
#
# a = s.get_cells_value(row_start=1, col_start=1, row_stop=3, col_stop=3, format='string')
# s.get_col_values(col=3, format='string')
# s.get_row_values(row=3, col_start=1, col_stop=None, format='string')
#
# a = s.get_cells_value(row_start=1, col_start=1, row_stop=3, col_stop=1, format='string')
# a = transpose(a)
# s.set_cells_value(a, row_start=1, col_start=1, format='formula')
#
# a = s.get_cells_value(row_start=1, col_start=1, row_stop=3, col_stop=1, format='string')
# s.set_col_values(a, row_start=1, col='G', format='formula')
#
# s.get_row_height()
# s.set_col_width(2000)



def str2bool(string):
    if string == 'False':
        return False
    else:
        return True

def transpose(l):
    """Transpose lists."""
    try:
        row_count, col_count = np.shape(l)
        return [list(x) for x in list(zip(*l))]
    except ValueError:
        return [[x] for x in l]

def transpose2(a):
    r = len(a)
    c = len(a[0])
    if c>1:
        b = []
        for i in range(c):
            b.append([0]*r)
        for i in range(r):
            for j in range(c):
              b[j][i] = a[i][j]
    else:
        b = [[] for x in range(r)]
        for i in range(r):
            b[i].append(a[i])
    return b



def _iterable(obj):
    return isinstance(obj, Iterable)

def _letter2num(string):

    string = string.lower()
    alphabet = 'abcdefghijklmnopqrstuvwxyz'

    col = 0
    for idx, s in enumerate(string):
        col += alphabet.index(s)+(idx)*26+1

    return(col)

def _check_row_value(row):
    if _iterable(row) == False:
        row2 = [row]
    else:
        row2 = copy.deepcopy(row)

    row2 = [r-1 for r in row2]

    outside_range = []
    for r in row2:
        if r < 0:
            outside_range.append(r)
    if len(outside_range) > 0:
        raise IndexError(f'Row number {outside_range} outside range.')

    return row2

def _check_col_value(col):
    if _iterable(col)==False or type(col)==str:
        col2 = [col]
    else:
        col2 = copy.deepcopy(col)

    for idx, c in enumerate(col2):
        if type(c) == str:
            col2[idx] = _letter2num(c)-1
        else:
            col2[idx] = c-1
            # raise TypeError(f'Col value {c} must be an int or str.')

    outside_range = []
    for c in col2:
        if type(c) == int:
            if c < 1:
                outside_range.append(c)
    if len(outside_range) > 1:
        raise IndexError(f'Col number {outside_range} outside range.')

    return col2

def get_cell_value_from_sheets(sheetObject_list, row, col, format='string'):
    """
    """
    values = []
    for sheetObject in sheetObject_list:
        values.append(sheetObject.get_cell_value(row=row, col=col, format=format))
    return values
