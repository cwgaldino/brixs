#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tmux
tmux ls
tmux attach -t 0
tmux attach-session -t 0
tmux rename-session -t 0 new_name
tmux new -s session_name
tmux kill-session -t name

panes
ctrl+b %                create panes vertical
ctrl+b <arrow-keys>     move to panes
ctrl+b "               create pane horizontal

windows:
ctrl+b c               new window
ctrl+b 0               move to window 0
ctrl+b 1               move to window 1
ctrl+b ,               rename window
exit                   close window

sessions:


ctrl+b d detach


windows:
...

"""
# %%
import libtmux
import numpy as np
from pathlib import Path
import time
import psutil
import os
import signal
import re
import subprocess
import warnings
import sys

def start_soffice(port=8100, folder='/opt/libreoffice7.0', nodefault=True, norestore=False, nologo=False, tmux_config=True, timeout=10):

    # set libreoffice folder
    folder = Path(folder)/'program'

    # tmux server
    if tmux_config:
        tmux_server = libtmux.Server(config_file=str(Path(os.path.dirname(__file__))/'tmux.conf'))
    else:
        tmux_server = libtmux.Server()

    # new session
    if tmux_server.has_session('libreoffice-wrapper'):
        tmux_session = tmux_server.find_where({'session_name':'libreoffice-wrapper'})
    else:
        tmux_session = tmux_server.new_session('libreoffice-wrapper')

    # new pane
    if tmux_session.find_where({'window_name':'soffice'}) is None:
        tmux_pane = tmux_session.new_window(attach=False, window_name='soffice').panes[0]
    else:
        tmux_pane = tmux_session.find_where({'window_name':'soffice'}).panes[0]
    time.sleep(1)
    tmux_pane.capture_pane = lambda: [x.rstrip() for x in tmux_pane.cmd('capture-pane', '-p', '-J').stdout]

    # initialize libreoffice
    print('Initialing soffice (it may take a few moments.)')
    tmux_pane.send_keys(f'cd {folder}', enter=True, suppress_history=True)
    time.sleep(0.1)

    options = ''
    if nodefault:
        options += ' --nodefault'
    if norestore:
        options += ' --norestore'
    if nologo:
        options += ' --nologo'

    tmux_pane.send_keys(f"./soffice{options} --accept='socket,host=localhost,port={port};urp;' &", enter=True, suppress_history=True)
    time.sleep(0.1)

    pid = 0
    start_time = time.time()
    while time.time() < start_time + timeout:
        try:
            pid = int(tmux_pane.capture_pane()[-2].split()[-1])
            if pid != 0:
                print(f'soffice started\nProcess pid {pid}')
                if not _has_pid(pid):
                    warnings.warn(f'\nsoffice seems to be already running. \nCommunication stabilished. \nNote that killing this process (pid={pid}) will not close soffice', stacklevel=2)
                return pid
        except ValueError:
            pass
    raise TimeoutError(f'soffice taking more than {timeout_soffice} seconds to load. Maybe try again.')

def _has_pid(pid):
    for proc in psutil.process_iter():
        if proc.as_dict(attrs=['pid'])['pid'] == pid: return True
    return False

def _get_children(pid):
    output = subprocess.check_output(["bash", "-c", f"pstree -p -n  {pid}"])
    pattern = re.compile(r"\((\d+)\)")
    pid_children = pattern.findall(output.decode('utf-8'))
    return [int(pid) for pid in pid_children]

def kill(*args, recursive=True):
    """Kill processes."""
    try:
        if recursive:
            pid2 = []
            for pid in args:
                for p in _get_children(pid):
                    pid2.append(p)
        else:
            pid2 = args
    except TypeError:
        if recursive:
            pid2 = []
            for pid in args[0]:
                for p in _get_children(pid):
                    pid2.append(p)
        else:
            pid2 = args

    for pid in pid2:
        try:
            os.kill(pid, signal.SIGKILL)
        except ProcessLookupError:
            pass

def kill_tmux(session_name='libreoffice-wrapper'):
    tmux_server = libtmux.Server()
    if tmux_server.has_session('libreoffice-wrapper'):
        tmux_server.kill_server()

def str2bool(string):
    if string == 'False':
        return False
    else:
        return True

def _name2ImplementationName(name):
    if name == 'scalc' or name == 'calc':
        return 'ScModelObj'
    elif name == 'swriter' or name == 'writer':
        return 'SwXTextDocument'
    elif name == 'simpress' or name == 'impress' or name == 'draw':
        return 'SdXImpressDocument'
    elif name == 'smath' or name == 'math':
        return 'com.sun.star.comp.Math.FormulaDocument'
    elif name == 'base' or name == 'sbase' or name == 'obase':
        return 'com.sun.star.comp.dba.ODatabaseDocument'
    else:
        raise ValueError('type must be calc, writer, impress, draw, math, or base')

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
                             "('y' or 'n').\n")

def _letter2num(string):

    string = string.lower()
    alphabet = 'abcdefghijklmnopqrstuvwxyz'

    n = 0
    for idx, s in enumerate(string):
        n += alphabet.index(s)+(idx)*26

    return n

def cell2num(string):
    if not isinstance(string, str):
        raise TypeError('value is not a cell.')
    if ':' in string:
        raise ValueError('cell seems to be a range.')
    try:
        temp = re.compile("([a-zA-Z]+)([0-9]+)")
        res = temp.match(string).groups()
        return _letter2num(res[0]), int(res[1])-1
    except AttributeError:
        raise ValueError('value is not a cell.')

def range2num(string):
    if ':' not in string:
        raise ValueError('range seems to be a cell.')

    return flatten((cell2num(string.split(':')[0]), cell2num(string.split(':')[1])))

def _check_row_value(row):
    if type(row) == str:
        row = int(row)-1

    if row < 0:
        raise ValueError('row cannot be negative')

    return row

def _check_column_value(column):
    if type(column) == str:
        column = _letter2num(column)

    if column < 0:
        raise ValueError('column cannot be negative')

    return column

def flatten(x):
    """Returns the flattened list or tuple."""
    if len(x) == 0:
        return x
    if isinstance(x[0], list) or isinstance(x[0], tuple):
        return flatten(x[0]) + flatten(x[1:])
    return x[:1] + flatten(x[1:])

def transpose(l):
    """Transpose lists."""
    try:
        row_count, col_count = np.shape(l)
        return [list(x) for x in list(zip(*l))]
    except ValueError:
        return [[x] for x in l]

def chunk(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

def partitionate(value, max_string=4000, i=1):
    # find i
    while len(str(value)) > max_string:
        i += 1

        if i >= len(value):
            raise ValueError(f'It seemns like length of rows are too big. Try setting shorter rows at a time and/or setting shorter data (less characters)')

        g = chunk(value, int(len(value)/i))
        while True:
            try:
                if len(str(next(g))) > max_string:
                    return partitionate(value, max_string=max_string, i=i)
            except StopIteration:
                break
        return chunk(value, int(len(value)/i))
    return chunk(value, int(len(value)/i))

def _parse_args(args, kwargs, req_args=[], opt_args=[]):
    """
    opt_args must be given in kwargs

    Range
    cell_start, cell_stop
    column_start, row_start, column_stop, row_stop

    cell
    column, row
    """
    req = {k:False for k in req_args}
    opt = {k:None for k in opt_args}

    for i, k in enumerate(opt_args):
        if k in kwargs:
            opt[k] = kwargs[k]

    lookout = 0
    del_list = []
    for i, k in enumerate(req_args):
        try:
            req[k] = kwargs[k]
            del_list.append(k)
        except KeyError:
            lookout += 1
    for k in del_list:
        del req_args[req_args.index(k)]

    for i in range(lookout):
        req[req_args[::-1][i]] = args[::-1][i]
    if lookout != 0:
        args = args[:-lookout]

    isRange = False
    if len(args) == 0:
        isRange = -1
    elif len(args) == 1:
        try:
            column_start, row_start = cell2num(args[0])
        except ValueError:
            column_start, row_start, column_stop, row_stop = range2num(args[0])
            isRange = True
    elif len(args) == 2:
        try:
            column_start, row_start = cell2num(args[0])
            column_stop, row_stop   = cell2num(args[1])
            isRange = True
        except TypeError:
            column_start = args[0]
            row_start    = args[1]
        except ValueError:
            column_start = args[0]
            row_start    = args[1]
    elif len(args) == 3:
        raise ValueError('cannot parse arguments.')
    elif len(args) == 4:
        column_start = args[0]
        row_start    = args[1]
        column_stop  = args[2]
        row_stop     = args[3]
        isRange = True

    if 'column' in kwargs:
        column_start = kwargs['column']
        isRange = False
    if 'row' in kwargs:
        row_start = kwargs['row']
        isRange = False
    if 'column_start' in kwargs:
        column_start = kwargs['column_start']
        isRange = False
    if 'row_start' in kwargs:
        row_start = kwargs['row_start']
        isRange = False
    if 'column_stop' in kwargs:
        column_stop = kwargs['column_stop']
        isRange = True
    if 'row_stop' in kwargs:
        row_stop = kwargs['row_stop']
        isRange = True
    if 'cell' in kwargs:
        column_start, row_start = cell2num(kwargs['cell'])
        isRange = False
    if 'cell_start' in kwargs:
        column_start, row_start = cell2num(kwargs['cell_start'])
        isRange = False
    if 'cell_stop' in kwargs:
        column_stop, row_stop = cell2num(kwargs['cell_stop'])
        isRange = True
    if 'Range' in kwargs:
        column_start, row_start, column_stop, row_stop = range2num(kwargs['Range'])
        isRange = True

    if isRange == -1:
        return -1, (), req, opt
    else:
        column_start = _check_column_value(column_start)
        row_start = _check_row_value(row_start)
        if isRange:
            column_stop = _check_column_value(column_stop)
            row_stop = _check_row_value(row_stop)

            if column_stop < column_start:
                raise ValueError(f'column_start ({column_start}) cannot be bigger than column_stop ({column_stop})')
            if row_stop < row_start:
                raise ValueError(f'row_start ({row_start}) cannot be bigger than row_stop ({row_stop})')

            return True, (column_start, row_start, column_stop, row_stop), req, opt
        else:
            return False, (column_start, row_start), req, opt

operator_map = {0:('is equal to', 'equal', '=', '==', 'EQUAL'),
                1:('is less than', 'less', '<', 'LESS'),
                2:('is greater than', 'greater', '>', 'GREATER'),
                3:('is less than or equal to', 'less-equal', 'less equal', '<=', 'LESS_EQUAL'),
                4:('is greater than or equal to', 'greater-equal', 'greater equal', '>=', 'GREATER_EQUAL'),
                5:('is not equal to', 'not-equal', 'not equal', '<>', '!=', 'NOT_EQUAL'),
                6:('is between', 'between', 'BETWEEN'),
                7:('is not between', 'not-between', 'not between', 'NOT_BETWEEN'),
                8:('is duplicated', 'duplicated', 'DUPLICATED', ''),
                9:('is not duplicated', 'not-duplicated', 'not duplicated', 'NOT_DUPLICATED', ''),
                10:('is in top N elements', 'is in top n elements', 'top', 'top n', 'top-n', 'TOP', ''),
                11:('is in bottom N elements', 'is in bottom n elements', 'bottom', 'bottom n', 'bottom-n', 'BOTTOM', ''),
                12:('is in top N percent', 'is in top n percent', 'top percent', 'top-percent', 'top n percent', 'top-n-percent', 'TOP_PERCENT', ''),
                13:('is in bottom N percent', 'is in bottom n percent', 'bottom percent', 'bottom-percent', 'bottom n percent', 'bottom-n-percent', 'BOTTOM_PERCENT', ''),
                14:('is above average', 'above', 'above average', 'above-average', 'ABOVE', ''),
                15:('is below average', 'below', 'below average', 'below-average', 'BELOW', ''),
                16:('is above or equal average', 'above equal', 'above-equal', 'above equal average', 'above-equal-average', 'ABOVE-EQUAL', ''),
                17:('is below or equal average', 'below equal', 'below-equal', 'below equal average', 'below-equal-average', 'BELOW-EQUAL', ''),
                18:('is error', 'error', 'ERROR', ''),
                19:('is not error', 'not error', 'not-error', 'NOT-ERROR', ''),
                20:('begins with', 'begins', 'BEGINS', ''),
                21:('ends with', 'ends', 'ENDS', ''),
                22:('contains', 'CONTAINS', ''),
                23:('does not contain', 'not contain', 'not-contain', 'NOT_CONTAIN', ''),
                24:('formula', 'Formula', 'FORMULA'),
               }


class soffice():

    def __init__(self, tmux_config=True, port=8100, folder='/opt/libreoffice7.0'):
        """
        server -> session (pycalc) -> window (soffice-python) -> pane (0)
        """
        timeout_python = 10
        timeout_soffice = 60

        # set libreoffice folder
        self.folder = Path(folder)/'program'

        # tmux server
        if tmux_config:
            self.tmux_server = libtmux.Server(config_file=str(Path(os.path.dirname(__file__))/'tmux.conf'))
        else:
            self.tmux_server = libtmux.Server()

        # new session
        if self.tmux_server.has_session('libreoffice-wrapper'):
            self.tmux_session = self.tmux_server.find_where({'session_name':'libreoffice-wrapper'})
        else:
            self.tmux_session = self.tmux_server.new_session('libreoffice-wrapper')

        # new pane
        if self.tmux_session.find_where({'window_name':'python'}) is None:
            self.tmux_pane = self.tmux_session.new_window(attach=False, window_name='python').panes[0]
        else:
            self.tmux_pane = self.tmux_session.find_where({'window_name':'python'}).panes[0]
        time.sleep(1)

        self.tmux_pane.capture_pane = lambda starting_line=0: [x.rstrip() for x in self.tmux_pane.cmd('capture-pane', '-p', '-J', f'-S {starting_line}').stdout]

        # initialize python
        print('Initialing python (it may take a few moments.)')
        self.tmux_pane.send_keys(f'cd {self.folder}', enter=True, suppress_history=True)
        time.sleep(0.1)
        self.tmux_pane.send_keys('./python', enter=True, suppress_history=True)
        start_time = time.time()
        flag = True
        while time.time() < start_time + timeout_python:
            if self.tmux_pane.capture_pane()[-1] == '>>>':
                print('done')
                flag = False
                break
        if flag:
            raise TimeoutError(f'Python taking more than {timeout_python} seconds to load. Maybe try again.')

        # python imports
        print('Importing modules (it may take a few moments.)')
        self.tmux_pane.send_keys(f"""try: exec(\"import traceback\")""" + '\n'\
                          """except Exception as e: print(e)""" + '\n'
                          , enter=True, suppress_history=False)
        start_time = time.time()
        while time.time() < start_time + 10:
            if self.check_running():
                pass
            else:
                break
        self.write(f"""import uno\n"""+\
                   f"""from com.sun.star.beans import PropertyValue\n"""
                   )
        time.sleep(0.1)
        print('done')

        # initialize communication
        print('Initializing comunications (it may take a few moments.)')
        self.write(f"""local = uno.getComponentContext()\n""" +\
                   f"""resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)\n""" +\
                   f"""context = resolver.resolve("uno:socket,host=localhost,port={port};urp;StarOffice.ComponentContext")\n""" +\
                   f"""desktop = context.ServiceManager.createInstanceWithContext('com.sun.star.frame.Desktop', context)\n"""
                   )
        time.sleep(0.1)

        # self.write("local = uno.getComponentContext()")
        # time.sleep(0.1)
        # self.write("""resolver = local.ServiceManager.createInstanceWithContext("com.sun.star.bridge.UnoUrlResolver", local)""")
        # time.sleep(0.1)
        # self.write(f"""context = resolver.resolve("uno:socket,host=localhost,port={port};urp;StarOffice.ComponentContext")""")
        # time.sleep(0.1)
        # self.write("desktop = context.ServiceManager.createInstanceWithContext('com.sun.star.frame.Desktop', context)")
        # time.sleep(0.1)
        try:
            self._get_tag()
        except soffice_python_error:
            self._set_tag(0)
        time.sleep(0.1)
        print('done')

    def _get_tag(self):
        return int(self.write("print(tag)")[0])

    def _set_tag(self, tag):
        self.write(f"tag = {tag}")

    def kill(self, force_kill=False):

        if self.tmux_server.has_session('libreoffice-wrapper') and force_kill == False:
            tmux_session = self.tmux_server.find_where({'session_name':'libreoffice-wrapper'})
            if tmux_session.find_where({'window_name':'soffice'}) is None:
                self.tmux_server.kill_server()
            else:
                self.write("documents = desktop.getComponents().createEnumeration()")
                hasMoreElements = str2bool(self.write("print(documents.hasMoreElements())")[0])
                if hasMoreElements:
                    print('Killing this process will close (without saving) the following opened files')
                    while hasMoreElements:
                        self.write("document = documents.nextElement()")
                        print(self.write(f"print(document.getImplementationName())")[0] + ': ' + self.write(f"document.getTitle()")[0])
                        hasMoreElements = str2bool(self.write("print(documents.hasMoreElements())")[0])
                if query_yes_no('Continue?', 'no'):
                    self.tmux_server.kill_server()
        else:
            self.tmux_server.kill_server()

    def check_running(self):
        # if self.tmux_pane.capture_pane()[-1] in ['>>>', '... >>>', '... ... >>>']:
        if self.tmux_pane.capture_pane()[-1].endswith('>>>'):
            return False
        else:
            return True

    def write(self, string, timeout=10, max_output=1000):
        # check quotes

        keys2send = f"try: exec(\"\"\"{string}\"\"\")" + '\n'\
                     """except: print('calc-ERROR'); print(traceback.format_exc())"""
        self.tmux_pane.send_keys(keys2send, enter=True, suppress_history=False)
        self.tmux_pane.enter()
        # time.sleep(0.05)
        time.sleep(0.001)

        start_time = time.time()
        while time.time() < start_time + timeout:
            if self.check_running():
                pass
            else:

                # output = self.tmux_pane.capture_pane()
                # i = output[::-1].index('>>> ' + keys2send.split('\n')[0])

                output = self.tmux_pane.capture_pane()
                starting_line = 0
                while '>>> ' + keys2send.split('\n')[0] not in output and starting_line>-max_output:
                    starting_line -= 10
                    output = self.tmux_pane.capture_pane(starting_line)
                if starting_line<-max_output:
                    raise valueError(f'output seems to be bigger than max_output {max_output}.')

                i = output[::-1].index('>>> ' + keys2send.split('\n')[0])
                if 'calc-ERROR' in output[-i:]:
                    f = output[::-1].index('>>>')
                    # try:
                    #     i = output[::-1].index('.calc-ERROR')
                    # except ValueError:
                        # i = [i for i, s in enumerate(output[-i::-1]) if 'calc-ERROR' in s][0]
                    raise soffice_python_error('\n'.join(output[-i:-1+f]))
                else:
                    return output[-i+2:-1]

                # i = output[::-1].index('>>> ' + keys2send.split('\n')[0])
                # if output[-i+2] == 'calc-ERROR':
                #     f = output[::-1].index('>>>')
                #     raise soffice_python_error('\n'.join(output[-i+3:-1+f]))
                # elif output[-i+2].startswith('...'):
                #     if output[-i+3] == 'calc-ERROR':
                #         f = output[::-1].index('>>>')
                #         raise soffice_python_error('\n'.join(output[-i+4:-1+f]))
                #     else:
                #         return output[-i+3:-1]
                # else:
                #     return output[-i+2:-1]

        raise TimeoutError(f'timeout: {timeout}s. Process is still running.')

    def read(self, timeout=10, max_output=1000):
        start_time = time.time()
        while time.time() < start_time + timeout:
            if self.check_running():
                pass
            else:
                output = self.tmux_pane.capture_pane()
                starting_line = 0
                while '>>>' not in output[:-1] and starting_line>-max_output:
                    starting_line -= 10
                    output = self.tmux_pane.capture_pane(starting_line)
                if starting_line<-max_output:
                    raise valueError(f'output seems to be bigger than max_output {max_output}.')

                output = output[[i for i, x in enumerate(output[:-1]) if x.startswith('>>>')][-1]+3:-1]
                if len(output) > 0:
                    if output[0] == 'calc-ERROR':
                        raise soffice_python_error('\n'.join(output))
                    else:
                        return output
                else:
                    return output
        raise TimeoutError(f'timeout: {timeout}s. Process is still running.')

    def Calc(self, filepath=None, new_file=False):

        type = 'scalc'
        if new_file:
            tag = self._new_file(type=type)
        elif filepath is None: # if a filepath is not defined
            print('Filepath was not given.')
            if self._has_open():
                print('However, something is open.')
                tag = self._connect_with_open(type=type)
                if not tag:
                    print(f'Opening a new file...')
                    tag = self._new_file(type=type)
            else:   # if nothing is open, open a new calc instance
                print(f'Nothing is opened.')
                tag = self._new_file(type=type)
        else:
            tag = self._open_file(filepath)

        return Calc(tag, soffice=self)

    def Writer(self, filepath=None, new_file=False):
        type = 'swriter'
        if new_file:
            tag = self._new_file(type=type)
        elif filepath is None: # if a filepath is not defined
            print('Filepath was not given.')
            if self._has_open():
                print('However, something is open.')
                tag = self._connect_with_open(type=type)
                if not tag:
                    print(f'Opening a new file...')
                    tag = self._new_file(type=type)
            else:   # if nothing is open, open a new calc instance
                print(f'Nothing is opened.')
                tag = self._new_file(type=type)
        else:
            tag = self._open_file(filepath)

        return Writer(tag, soffice=self)

    def Impress(self, filepath=None, new_file=False):
        type = 'simpress'
        if new_file:
            tag = self._new_file(type=type)
        elif filepath is None: # if a filepath is not defined
            print('Filepath was not given.')
            if self._has_open():
                print('However, something is open.')
                tag = self._connect_with_open(type=type)
                if not tag:
                    print(f'Opening a new file...')
                    tag = self._new_file(type=type)
            else:   # if nothing is open, open a new calc instance
                print(f'Nothing is opened.')
                tag = self._new_file(type=type)
        else:
            tag = self._open_file(filepath)

        return Impress(tag, soffice=self)

    def Draw(self, filepath=None, new_file=False):
        type = 'sdraw'
        if new_file:
            tag = self._new_file(type=type)
        elif filepath is None: # if a filepath is not defined
            print('Filepath was not given.')
            if self._has_open():
                print('However, something is open.')
                tag = self._connect_with_open(type=type)
                if not tag:
                    print(f'Opening a new file...')
                    tag = self._new_file(type=type)
            else:   # if nothing is open, open a new calc instance
                print(f'Nothing is opened.')
                tag = self._new_file(type=type)
        else:
            tag = self._open_file(filepath)

        return Draw(tag, soffice=self)

    def Math(self, filepath=None, new_file=False):
        type = 'smath'
        if new_file:
            tag = self._new_file(type=type)
        elif filepath is None: # if a filepath is not defined
            print('Filepath was not given.')
            if self._has_open():
                print('However, something is open.')
                tag = self._connect_with_open(type=type)
                if not tag:
                    print(f'Opening a new file...')
                    tag = self._new_file(type=type)
            else:   # if nothing is open, open a new calc instance
                print(f'Nothing is opened.')
                tag = self._new_file(type=type)
        else:
            tag = self._open_file(filepath)

        return Math(tag, soffice=self)

    def Base(self, filepath=None, new_file=False):
        raise NotImplementedError('Base not implemented yet')

    def save(self, tag, type, extension, filepath=None):
        """Save  file.

        Args:
            filepath (string or pathlib.Path, optional): filepath to save file.
        """
        if filepath is None and self.get_filepath(tag)=='':  # ok
            temporary_path = Path.cwd()/self.get_title(tag)
            temporary_path = temporary_path.with_suffix(extension)
            if query_yes_no(f'Filepath not defined. Wish to save at {temporary_path}?'):
                filepath = temporary_path
            else:
                print('File not saved')
                return
        elif filepath is None and self.get_filepath(tag)!='':
            filepath = Path(self.get_filepath(tag)).with_suffix(extension)
        else:
            filepath = Path(filepath).with_suffix(extension)

        # save
        self.write(f'URL = uno.systemPathToFileUrl("{filepath}")')
        self.write(f"properties = (PropertyValue('FilterName', 0, '{type}', 0),)")
        self.write(f"document_{tag}.storeAsURL(URL, properties)")
        time.sleep(0.1)
        print(f'Saved at: {filepath}')
        return self.get_filepath(tag)

    def close(self, tag):
        """Close window."""
        self.write(f'document_{tag}.close(True)')
        self.write(f'del document_{tag}')

    def get_filepath(self, tag):
        hasFilepath = str2bool(self.write(f"print(document_{tag}.hasLocation())")[0])
        if hasFilepath: # check if file has filepath
            return self.write(f"print(uno.fileUrlToSystemPath(document_{tag}.getURL()))")[0]
        else: return ''

    def get_title(self, tag):
        return self.write(f"print(document_{tag}.getTitle())")[0]

    def _new_file(self, type):
        """type = calc or writer
        https://wiki.documentfoundation.org/Macros/Basic/Documents
        """

        if type == 'scalc' or type == 'calc':
            self.write("URL = 'private:factory/scalc'")
        elif type == 'swriter' or type == 'writer':
            self.write("URL = 'private:factory/swriter'")
        elif type == 'simpress' or type == 'impress':
            self.write("URL = 'private:factory/SdXImpressDocument'")
        elif type == 'sdraw' or type == 'draw':
            self.write("URL = 'private:factory/sdraw'")
        elif type == 'smath' or type == 'math':
            self.write("URL = 'private:factory/smath'")
        elif type == 'base' or type == 'obase' or type == 'sbase':
            raise NotImplementedError('base not implemented yet')
        else:
            raise ValueError('type must be calc, writer, impress, draw, math, or base')
        tag = self._get_tag() + 1
        self._set_tag(tag)
        self.write(f"document_{tag} = desktop.loadComponentFromURL(URL, '_default', 0, ())")
        # time.sleep(0.1)
        title = self.get_title(tag=tag)
        time.sleep(0.1)
        start_time = time.time()
        while title is [] and time.time() < start_time + 10:
            time.sleep(0.1)
            title = self.get_title(tag=tag)

        print(f'Connected with opened file: {title}')
        return tag

    def _open_file(self, filepath):
        tag = self._get_tag() + 1
        self._set_tag(tag)
        print('Filepath was given.')
        filepath = Path(filepath).absolute()
        self.write(f"URL = uno.systemPathToFileUrl('{filepath}')")
        self.write(f'document_{tag} = desktop.loadComponentFromURL(URL, "_default", 0, ())')
        url = self.write(f"print(uno.fileUrlToSystemPath(document_{tag}.getURL()))")[0]
        print(f'Connected with file: {url}')
        return tag

    def _has_open(self):
        # is there something open???
        self.write("document =  desktop.getCurrentComponent()")
        return str2bool(self.write("print('False') if document is None else print('True')")[0])

    def _connect_with_open(self, type):

        self.write(f"document = desktop.getCurrentComponent()")
        implementationName = self.write(f"print(document.getImplementationName())")
        if _name2ImplementationName(type) in implementationName:
            return self._connect_with_current()
        else:
            print(f'current document is not type {type}.')
            print(f'Searching for any opened {type} instances...')
            self.write("documents = desktop.getComponents().createEnumeration()")


            open_new = True
            hasMoreElements = str2bool(self.write("print(documents.hasMoreElements())")[0])
            while hasMoreElements:
                self.write("document = documents.nextElement()")
                implementationName = self.write("print(document.getImplementationName())")
                if _name2ImplementationName(type) in implementationName:
                    tag = self._get_tag() + 1
                    self._set_tag(tag)
                    self.write(f"document_{tag} = document")

                    if self.get_filepath(tag) == '':
                        title = self.get_title(tag)
                        print(f'Connected with opened file: {title}')
                    else:
                        print(f'Connected with opened file: {self.get_filepath(tag)}')
                    return tag
                hasMoreElements = str2bool(self.write("print(documents.hasMoreElements())")[0])
            print(f'No opened {type} instances were found.')
            return False

    def _connect_with_current(self):
        tag = self._get_tag() + 1
        self._set_tag(tag)
        self.write(f"document_{tag} = desktop.getCurrentComponent()")

        print(f'Connecting with current file.')

        time.sleep(0.1)
        if self.get_filepath(tag) != '':
            print(f'Connected with opened file: {self.get_filepath(tag)}')
        else:
            title = self.get_title(tag)
            time.sleep(0.1)
            start_time = time.time()
            while title is [] and time.time() < start_time + 10:
                time.sleep(0.1)
                title = self.get_title(tag)
            print(f'Connected with opened file: {title}')
        return tag


class soffice_python_error(Exception):

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class Writer():

    def __init__(self, tag, soffice):

        self.soffice = soffice
        self.tag = tag
        self.sheet_tags = []
        self.filepath = self.get_filepath()

    def write(self, string, timeout=10):
        return self.soffice.write(string=string, timeout=timeout)

    def read(self, timeout=10):
        return self.soffice.read(timeout=timeout)

    def get_filepath(self):
        return self.soffice.get_filepath(self.tag)

    def get_title(self):
        return self.soffice.get_title(self.tag)

    def save(self, filepath=None, type='odt'):
        if type in ['odt', '.odt', 'swriter', 'writer', 'Writer']:
            type = 'writer8'
            extension = '.odt'
        else:
            raise ValueError('type not recognised')

        if filepath is None and self.filepath!='' and self.get_filepath()!='':
            if self.filepath != self.get_filepath():
                print('File was last saved in a different filepath than the one stored here.')
                print('Last saved path:' + self.get_filepath())
                print('Stored path:' + self.filepath)
                if query_yes_no(f'Wish to save at {self.get_filepath()}?'):
                    # self.soffice.Calc(filepath()
                    filepath = self.get_filepath()
                else:
                    if query_yes_no(f'Wish to save at {self.filepath}?'):
                        filepath = self.filepath
                    else:
                        print('File not saved')
                        return

        self.filepath = self.soffice.save(tag=self.tag, type=type, extension=extension, filepath=filepath)

    def close(self):
        """Close window."""
        for tag in self.sheet_tags:
            try:
                self.write(f'del document_{tag}')
            except soffice_python_error:
                pass
        self.soffice.close(self.tag)


class Impress():

    def __init__(self, tag, soffice):

        self.soffice = soffice
        self.tag = tag
        self.sheet_tags = []
        self.filepath = self.get_filepath()

    def write(self, string, timeout=10):
        return self.soffice.write(string=string, timeout=timeout)

    def read(self, timeout=10):
        return self.soffice.read(timeout=timeout)

    def get_filepath(self):
        return self.soffice.get_filepath(self.tag)

    def get_title(self):
        return self.soffice.get_title(self.tag)

    def save(self, filepath=None, type='odp'):
        if type in ['odp', '.odt', 'simpress', 'impress', 'Impress']:
            type = 'impress8'
            extension = '.odp'
        else:
            raise ValueError('type not recognised')

        if filepath is None and self.filepath!='' and self.get_filepath()!='':
            if self.filepath != self.get_filepath():
                print('File was last saved in a different filepath than the one stored here.')
                print('Last saved path:' + self.get_filepath())
                print('Stored path:' + self.filepath)
                if query_yes_no(f'Wish to save at {self.get_filepath()}?'):
                    # self.soffice.Calc(filepath()
                    filepath = self.get_filepath()
                else:
                    if query_yes_no(f'Wish to save at {self.filepath}?'):
                        filepath = self.filepath
                    else:
                        print('File not saved')
                        return

        self.filepath = self.soffice.save(tag=self.tag, type=type, extension=extension, filepath=filepath)

    def close(self):
        """Close window."""
        for tag in self.sheet_tags:
            try:
                self.write(f'del document_{tag}')
            except soffice_python_error:
                pass
        self.soffice.close(self.tag)


class Draw():

    def __init__(self, tag, soffice):

        self.soffice = soffice
        self.tag = tag
        self.sheet_tags = []
        self.filepath = self.get_filepath()

    def write(self, string, timeout=10):
        return self.soffice.write(string=string, timeout=timeout)

    def read(self, timeout=10):
        return self.soffice.read(timeout=timeout)

    def get_filepath(self):
        return self.soffice.get_filepath(self.tag)

    def get_title(self):
        return self.soffice.get_title(self.tag)

    def save(self, filepath=None, type='odf'):
        if type in ['odf', '.odf', 'simpress', 'impress', 'Impress']:
            type = 'impress8'
            extension = '.odf'
        else:
            raise ValueError('type not recognised')

        if filepath is None and self.filepath!='' and self.get_filepath()!='':
            if self.filepath != self.get_filepath():
                print('File was last saved in a different filepath than the one stored here.')
                print('Last saved path:' + self.get_filepath())
                print('Stored path:' + self.filepath)
                if query_yes_no(f'Wish to save at {self.get_filepath()}?'):
                    # self.soffice.Calc(filepath()
                    filepath = self.get_filepath()
                else:
                    if query_yes_no(f'Wish to save at {self.filepath}?'):
                        filepath = self.filepath
                    else:
                        print('File not saved')
                        return

        self.filepath = self.soffice.save(tag=self.tag, type=type, extension=extension, filepath=filepath)

    def close(self):
        """Close window."""
        for tag in self.sheet_tags:
            try:
                self.write(f'del document_{tag}')
            except soffice_python_error:
                pass
        self.soffice.close(self.tag)


class Math():

    def __init__(self, tag, soffice):

        self.soffice = soffice
        self.tag = tag
        self.sheet_tags = []
        self.filepath = self.get_filepath()

    def write(self, string, timeout=10):
        return self.soffice.write(string=string, timeout=timeout)

    def read(self, timeout=10):
        return self.soffice.read(timeout=timeout)

    def get_filepath(self):
        return self.soffice.get_filepath(self.tag)

    def get_title(self):
        return self.soffice.get_title(self.tag)

    def save(self, filepath=None, type='odf'):
        if type in ['odf', '.odf', 'smath', 'math', 'Math']:
            type = 'math8'
            extension = '.odf'
        else:
            raise ValueError('type not recognised')

        if filepath is None and self.filepath!='' and self.get_filepath()!='':
            if self.filepath != self.get_filepath():
                print('File was last saved in a different filepath than the one stored here.')
                print('Last saved path:' + self.get_filepath())
                print('Stored path:' + self.filepath)
                if query_yes_no(f'Wish to save at {self.get_filepath()}?'):
                    # self.soffice.Calc(filepath()
                    filepath = self.get_filepath()
                else:
                    if query_yes_no(f'Wish to save at {self.filepath}?'):
                        filepath = self.filepath
                    else:
                        print('File not saved')
                        return

        self.filepath = self.soffice.save(tag=self.tag, type=type, extension=extension, filepath=filepath)

    def close(self):
        """Close window."""
        for tag in self.sheet_tags:
            try:
                self.write(f'del document_{tag}')
            except soffice_python_error:
                pass
        self.soffice.close(self.tag)


class Base():
    def __init__(self):
        raise NotImplementedError('Base not implemented yet')


class Calc():

    def __init__(self, tag, soffice):

        self.soffice = soffice
        self.tag = tag
        self.sheet_tags = []
        self.filepath = self.get_filepath()

        self.write(f"""def get_cell_property_recursively(tag, column, row, name, attrs=[]):\n""" +\
                   f"""    f = dict()\n""" +\
                   f"""    sheet = eval(f"sheet_{{tag}}")\n""" +\
                   f"""    i = [x.Name for x in sheet.getCellByPosition(column, row).getPropertySetInfo().Properties].index(name)\n""" +\
                   f"""    if attrs == []:\n""" +\
                   f"""        try:\n""" +\
                   f"""            if type(sheet.getCellByPosition(column, row).getPropertyValue(name).value) is int:\n""" +\
                   f"""                return sheet.getCellByPosition(column, row).getPropertyValue(name).value\n""" +\
                   f"""            elif type(sheet.getCellByPosition(column, row).getPropertyValue(name).value) is str:\n""" +\
                   f"""                return sheet.getCellByPosition(column, row).getPropertyValue(name).value\n""" +\
                   f"""            else:\n""" +\
                   f"""                keys = sheet.getCellByPosition(column, row).getPropertyValue(name).value.__dir__()\n""" +\
                   f"""        except AttributeError:\n""" +\
                   f"""            return sheet.getCellByPosition(column, row).getPropertyValue(name)\n""" +\
                   f"""    else:\n""" +\
                   f"""        try:\n""" +\
                   f"""            t = ''\n""" +\
                   f"""            for attr in attrs:\n""" +\
                   f"""                t += f"__getattr__('{{attr}}')."\n""" +\
                   f"""            if type(eval(f"sheet.getCellByPosition(column, row).getPropertyValue(name).{{t}}value")) is str:\n""" +\
                   f"""                return eval(f"sheet.getCellByPosition(column, row).getPropertyValue(name).{{t}}value")\n""" +\
                   f"""            else:\n""" +\
                   f"""                keys = eval(f"sheet.getCellByPosition(column, row).getPropertyValue(name).{{t}}value.__dir__()")\n""" +\
                   f"""        except AttributeError:\n""" +\
                   # f"""            print(t)\n""" +\
                   f"""            return eval(f"sheet.getCellByPosition(column, row).getPropertyValue(name).{{t[:-1]}}")\n""" +\
                   f"""    for key in keys:\n""" +\
                   f"""        f[key] = get_cell_property_recursively(tag, column, row, name, attrs+[key])\n""" +\
                   f"""    return f""")


        self.write(f"""def get_cells_property_recursively(tag, column_start, row_start, column_stop, row_stop, name, attrs=[]):\n""" +\
                   f"""    f = dict()\n""" +\
                   f"""    sheet = eval(f"sheet_{{tag}}")\n""" +\
                   f"""    i = [x.Name for x in sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertySetInfo().Properties].index(name)\n""" +\
                   f"""    if attrs == []:\n""" +\
                   f"""        try:\n""" +\
                   f"""            if type(sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name).value) is int:\n""" +\
                   f"""                return sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name).value\n""" +\
                   f"""            elif type(sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name).value) is str:\n""" +\
                   f"""                return sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name).value\n""" +\
                   f"""            else:\n""" +\
                   f"""                keys = sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name).value.__dir__()\n""" +\
                   f"""        except AttributeError:\n""" +\
                   f"""            return sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name)\n""" +\
                   f"""    else:\n""" +\
                   f"""        try:\n""" +\
                   f"""            t = ''\n""" +\
                   f"""            for attr in attrs:\n""" +\
                   f"""                t += f"__getattr__('{{attr}}')."\n""" +\
                   f"""            if type(eval(f"sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name).{{t}}value")) is str:\n""" +\
                   f"""                return eval(f"sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name).{{t}}value")\n""" +\
                   f"""            else:\n""" +\
                   f"""                keys = eval(f"sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name).{{t}}value.__dir__()")\n""" +\
                   f"""        except AttributeError:\n""" +\
                   # f"""            print(t)\n""" +\
                   f"""            return eval(f"sheet.getCellRangeByPosition(column_start, row_start, column_stop, row_stop).getPropertyValue(name).{{t[:-1]}}")\n""" +\
                   f"""    for key in keys:\n""" +\
                   f"""        f[key] = get_cells_property_recursively(tag, column_start, row_start, column_stop, row_stop, name, attrs+[key])\n""" +\
                   f"""    return f""")

        self.write(f"""def _letter2num(string):\n"""+\
                    f"""    string = string.lower()\n"""+\
                    f"""    alphabet = 'abcdefghijklmnopqrstuvwxyz'\n"""+\
                    f"""    n = 0\n"""+\
                    f"""    for idx, s in enumerate(string):\n"""+\
                    f"""        n += alphabet.index(s)+(idx)*26\n"""+\
                    f"""    return n""")

        self.write(f"""from com.sun.star.table import TableBorder, BorderLine2, ShadowFormat\n""" +\
                  f"""from com.sun.star.awt import Size, Point\n""" +\
                  f"""from com.sun.star.lang import Locale\n""" +\
                  f"""from com.sun.star.util import CellProtection\n"""+\
                  f"""from com.sun.star.sheet.ConditionOperator import NONE""")

        for k in operator_map:
            if operator_map[k][-1] != '':
                self.write(f"from com.sun.star.sheet.ConditionOperator import {operator_map[k][-1]}")

    def write(self, string, timeout=10):
        return self.soffice.write(string=string, timeout=timeout)

    def read(self, timeout=10):
        return self.soffice.read(timeout=timeout)

    def get_filepath(self):
        return self.soffice.get_filepath(self.tag)

    def get_title(self):
        return self.soffice.get_title(self.tag)

    def save(self, filepath=None, type='ods'):
        if type in ['ods', '.ods', 'calc8', 'calc', 'Calc']:
            type = 'calc8'
            extension = '.ods'
        elif type in ['excel', 'xlsx', '.xlsx']:
            type = 'Calc MS Excel 2007 XML'
            extension = '.xlsx'
        elif type in ['xls', '.xls']:
            type = 'MS Excel 97'
            extension = '.xls'
        elif type in ['csv', '.csv', 'text', 'txt', '.txt']:
            type = 'Text - txt - csv (StarCalc)'
            extension = '.csv'
        else:
            raise ValueError('type not recognised')

        if filepath is None and self.filepath!='' and self.get_filepath()!='':
            if self.filepath != self.get_filepath():
                print('File was last saved in a different filepath than the one stored here.')
                print('Last saved path:' + self.get_filepath())
                print('Stored path:' + self.filepath)
                if query_yes_no(f'Wish to save at {self.get_filepath()}?'):
                    # self.soffice.Calc(filepath()
                    filepath = self.get_filepath()
                else:
                    if query_yes_no(f'Wish to save at {self.filepath}?'):
                        filepath = self.filepath
                    else:
                        print('File not saved')
                        return

        self.filepath = self.soffice.save(tag=self.tag, type=type, extension=extension, filepath=filepath)

    def close(self):
        """Close window."""
        for tag in self.sheet_tags:
            try:
                self.write(f'del document_{tag}')
            except soffice_python_error:
                pass
        self.soffice.close(self.tag)




    def get_sheets_count(self):
        return int(self.write(f"print(document_{self.tag}.Sheets.getCount())")[0])

    def get_sheets_name(self):
        """Returns the sheets names in a tuple."""
        return eval(self.write(f"print(document_{self.tag}.Sheets.getElementNames())")[0])

    def get_sheet_position(self, name):
        try:
            return self.get_sheets_name().index(name)
        except ValueError:
            raise SheetNameDoNotExistError(f'{name} does not exists')

    def get_sheet_name_by_position(self, position):
        names = self.get_sheets_name()
        if position > len(names) or position < 0:
            raise IndexError('Position outside range.')
        return names[position]

    def get_sheet_by_position(self, position):
        name = self.get_sheet_name_by_position(position)
        return self.get_sheet(name)

    def get_sheet(self, name):
        if name in self.get_sheets_name():
            tag = self.soffice._get_tag() + 1
            self.write(f"sheet_{tag} = document_{self.tag}.Sheets.getByName('{name}')")
            self.sheet_tags.append(tag)
            return Sheet(tag, self)
        else:
            raise SheetNameDoNotExistError(f'{name} does not exists')

    def insert_sheet(self, name, position=None):
        """name can be a string or a list

        position starts from 0. If position = 0, the sheet will be the first one.
        """
        if position is None:
            position = self.get_sheets_count()+1

        if name in self.get_sheets_name():
            raise SheetNameExistError(f'{name} already exists')

        self.write(f"document_{self.tag}.Sheets.insertNewByName('{name}', {position})")

    def remove_sheet(self, name):
        if name in self.get_sheets_name():
            if self.get_sheets_count() > 1:
                self.write(f"document_{self.tag}.Sheets.removeByName('{name}')")
            else:
                raise SheetRemoveError(f"{name} cannot be removed because it is the only existing sheet")
        else:
            raise SheetNameDoNotExistError(f'{name} does not exists')

    def remove_sheets_by_position(self, position):
        self.remove_sheet(self.get_sheet_name_by_position(position))

    def move_sheet(self, name, position):
        if position < 0:
            raise ValueError('position cannot be negative')
        if name in self.get_sheets_name():
            self.write(f"document_{self.tag}.Sheets.moveByName('{name}', {position})")
        else:
            raise SheetNameDoNotExistError(f'{name} does not exists')

    def copy_sheet(self, name, new_name, position):
        if position < 0:
            raise ValueError('position cannot be negative')
        if name in self.get_sheets_name():
            self.write(f"document_{self.tag}.Sheets.copyByName('{name}', '{new_name}', {position})")
        else:
            raise SheetNameDoNotExistError(f'{name} does not exists')

    


    def get_styles(self):
        return self.write(f"""print([x.Name for x in document_{self.tag}.getStyleFamilies()['CellStyles']])""")[0]

    def remove_style(self, name):
        if name in self.get_styles():
            self.write(f"""document_{self.tag}.getStyleFamilies()['CellStyles'].removeByName('{name}')""")
        else:
            raise valueError(f'{name} is not a style')

    def new_style(self, name, properties):
        """
        name = 'MyStyle'
        overwrite = True
        properties = {'CellBackColor':-1, 'TopBorder':{'Color':16776960}}
        """

        keys_string = 'keys = ('
        values_string = 'values = ('
        for name0 in properties:
            keys_string += "'" + name0 +"', "
            sheet = self.get_sheet_by_position(0)
            func = sheet._get_property_function(0, 0, name0)
            if func is not None:
                if isinstance(properties[name0], dict):
                    d = sheet.get_cell_property(0, 0, name0)
                    for k in d:
                        if k in properties[name0]:
                            d[k] = properties[name0][k]
                        else:
                            raise ValueError(f'{k} is not present in {name0}.')
                    f_string = f"{func}(" + ','.join([f"{k}={v}" for k, v in d.items()]) + ")"
                    values_string += f_string +', '

                else:
                    raise ValueError(f'property for {name0} must e a dict')
            else:
                values_string += str(properties[name0]) +', '

        if name in self.get_styles():
            raise ValueError(f'{name} already exists')

        self.write(f"""new_style = document_{self.tag}.createInstance('com.sun.star.style.CellStyle')\n"""+\
                f"""document_{self.tag}.getStyleFamilies()['CellStyles'].insertByName('{name}', new_style)\n"""+\
                f"""{keys_string[:-2]})\n"""+\
                f"""{values_string[:-2]})\n"""+
                f"""new_style.setPropertyValues(keys, values)""")


class SheetNameExistError(Exception):

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class SheetNameDoNotExistError(Exception):

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class SheetRemoveError(Exception):

    def __init__(self, message):
        self.message = message
        super().__init__(self.message)


class Sheet():

    def __init__(self, tag, Calc):
        self.Calc = Calc
        self.tag = tag

    def write(self, string, timeout=10):
        return self.Calc.soffice.write(string=string, timeout=timeout)

    def read(self, timeout=10):
        return self.Calc.soffice.read(timeout=timeout)

    def get_name(self):
        return self.write(f'print(sheet_{self.tag}.getName())')[0]

    def set_name(self, name):
        self.write(f"sheet_{self.tag}.setName('{name}')")

    def isVisible(self):
        return str2bool(self.write(f"print(sheet_{self.tag}.IsVisible)")[0])

    def move(self, position):
        self.Calc.move_sheet(self.get_name(), position)

    def remove(self):
        self.write(f'del sheet_{self.tag}')
        self.Calc.remove_sheet(self.get_name())

    def get_last_row(self):
        """starts from 0"""
        l = self.write(f"""row_name = sheet_{self.tag}.getRowDescriptions()[-1]\n"""+\
                           f"""idx = int(row_name.split()[-1])\n"""+\
                           f"""visible = False\n"""+\
                           f"""while visible == False:\n"""+\
                           f"""    visible = sheet_{self.tag}.getRows()[idx].IsVisible\n"""+\
                           f"""    idx += 1\n"""+\
                           f"""print(idx-2)""")[-1]
        return int(l)

    def get_last_column(self):
        """starts from 0"""
        l = self.write(f"""col_name = sheet_{self.tag}.getColumnDescriptions()[-1]\n"""+\
                           f"""idx = int(_letter2num(col_name.split()[-1]))\n"""+\
                           f"""visible = False\n"""+\
                           f"""while visible == False:\n"""+\
                           f"""    visible = sheet_{self.tag}.getColumns()[idx+1].IsVisible\n"""+\
                           f"""    idx += 1\n"""+\
                           f"""print(idx-1)""")[-1]
        return int(l)

    def get_row_length(self, row):
        row = _check_row_value(row)
        l = self.write(f"""col_name = sheet_{self.tag}.getColumnDescriptions()[-1]\n"""+\
                           f"""idx = int(_letter2num(col_name.split()[-1]))\n"""+\
                           f"""visible = False\n"""+\
                           f"""while visible == False:\n"""+\
                           f"""    visible = sheet_{self.tag}.getColumns()[idx+1].IsVisible\n"""+\
                           f"""    idx += 1\n"""+\
                           f"""lc = idx-1\n"""+\
                           f"""if len(sheet_{self.tag}.getCellRangeByPosition(0, {row}, lc, {row}).queryEmptyCells()) == 0:\n"""+\
                           f"""    row_length = lc+1\n"""+\
                           f"""else:\n"""+\
                           f"""    startColumn = sheet_{self.tag}.getCellRangeByPosition(0, {row}, lc, {row}).queryEmptyCells()[-1].RangeAddress.StartColumn\n"""+\
                           f"""    endColumn = sheet_{self.tag}.getCellRangeByPosition(0, {row}, lc, {row}).queryEmptyCells()[-1].RangeAddress.EndColumn\n"""+\
                           f"""    if endColumn == lc:\n"""+\
                           f"""        row_length = sheet_{self.tag}.getCellRangeByPosition(0, {row}, lc, {row}).queryEmptyCells()[-1].RangeAddress.StartColumn\n"""+\
                           f"""    else:\n"""+\
                           f"""        row_length = lc + 1\n"""+\
                           f"""print(row_length)""")[-1]
        return int(l)

    def get_column_length(self, column):
        column = _check_column_value(column)
        l = self.write(f"""row_name = sheet_{self.tag}.getRowDescriptions()[-1]\n"""+\
                           f"""idx = int(row_name.split()[-1])\n"""+\
                           f"""visible = False\n"""+\
                           f"""while visible == False:\n"""+\
                           f"""    visible = sheet_{self.tag}.getRows()[idx].IsVisible\n"""+\
                           f"""    idx += 1\n"""+\
                           f"""lr = idx-1\n"""+\
                           f"""if len(sheet_{self.tag}.getCellRangeByPosition({column}, 0, {column}, lr).queryEmptyCells()) == 0:\n"""+\
                           f"""    column_length = lr+1\n"""+\
                           f"""else:\n"""+\
                           f"""    startRow = sheet_{self.tag}.getCellRangeByPosition({column}, 0, {column}, lr).queryEmptyCells()[-1].RangeAddress.StartRow\n"""+\
                           f"""    endRow = sheet_{self.tag}.getCellRangeByPosition({column}, 0, {column}, lr).queryEmptyCells()[-1].RangeAddress.EndRow\n"""+\
                           f"""    if endRow == lr:\n"""+\
                           f"""        column_length = sheet_{self.tag}.getCellRangeByPosition({column}, 0, {column}, lr).queryEmptyCells()[-1].RangeAddress.StartRow\n"""+\
                           f"""    else:\n"""+\
                           f"""        column_length = lr + 1\n"""+\
                           f"""print(column_length)""")[-1]
        return int(l)

    def set_column_width(self, column, width):
        """column can be single value or list
        width can be single value or list
        """
        if isinstance(column, list) or isinstance(column, tuple):
            column = [_check_column_value(c) for c in column]
            if isinstance(width, list) or isinstance(width, tuple):
                if len(column) != len(width):
                    raise ValueError("Width and column must be the same length.")
                for w in width:
                    if w < 0:
                        raise ValueError('width cannot be negative')
            else:
                if width < 0:
                    raise ValueError('width cannot be negative')
                width = [width]*len(column)
            self.write(f"""for r, w in zip({column}, {width}):\n"""+\
                       f"""    sheet_{self.tag}.getColumns()[r].setPropertyValue('Width', int(w))""")
        else:
            column = _check_column_value(column)
            if width < 0:
                raise ValueError('width cannot be negative')
            self.write(f"sheet_{self.tag}.getColumns()[{column}].setPropertyValue('Width', {int(width)})")

    def get_column_width(self, column):
        column = _check_column_value(column)
        return int(self.write(f"print(sheet_{self.tag}.getColumns()[{column}].Width)")[0])

    def set_row_height(self, row, height):
        """row can be single value or list
        height can be single value or list
        """
        if isinstance(row, list) or isinstance(row, tuple):
            row = [_check_row_value(r) for r in row]
            if isinstance(height, list) or isinstance(height, tuple):
                if len(row) != len(height):
                    raise ValueError("Height and row must be the same lenght.")
                for h in height:
                    if h < 0:
                        raise ValueError('height cannot be negative')
            else:
                if height < 0:
                    raise ValueError('height cannot be negative')
                height = [height]*len(row)
            self.write(f"""for r, h in zip({row}, {height}):\n"""+\
                       f"""    sheet_{self.tag}.getRows()[r].setPropertyValue('Height', h)""")
        else:
            row = _check_row_value(row)
            if height < 0:
                raise ValueError('height cannot be negative')
            self.write(f"sheet_{self.tag}.getRows()[{row}].setPropertyValue('Height', {height})")

    def get_row_height(self, row):
        row = _check_row_value(row)
        return int(self.write(f"print(sheet_{self.tag}.getRows()[{row}].Height)")[0])

    def set_value(self, *args, **kwargs):
        """
        value has no maximum number of rows.
        value limits the number of columns.
        The maximum number of columns is determined by the length of the values in each cell (number of characters).
        """
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=['value'], opt_args=['format'])
        value = req['value']

        if opt['format'] is None:
            format = 'formula'
        else:
            format = opt['format']

        if isRange == -1:
            raise ValueError('missing range/cell')
        elif isRange == False:  # check if value is a range
            if type(value) is not str:
                try:
                    if len(value) > 0:
                        if type(value[0]) is not str:
                            if len(value) > 0:
                                isRange = True
                                c = list(c)
                                c.append(c[0] + len(value[0])-1)
                                c.append(c[1] + len(value)-1)
                except TypeError:
                    pass

        if isRange == True:
            if c[2] - c[0] != len(value[0])-1 or c[3] - c[1] != len(value)-1:
                raise ValueError('cell range does not match data/value size.')

            if len(value) > 1:
                if sum((len(value[0]) - len(v) for v in value)) != 0:
                    raise ValueError('all columns/rows must have the same length.')


            # deal with big columns
            value = [list(x) for x in value]
            g = partitionate(value)
            self.write(f"""value = {next(g)}""")
            while True:
                try:
                    self.write(f"""value.extend({next(g)})""")
                except StopIteration:
                    break

            if format == 'formula':
                self.write(f"sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).setFormulaArray(value)")
            elif format == 'string':
                self.write(f"sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).setDataArray(value)")
            elif format == 'number':
                self.write(f"sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).setDataArray(value)")
            else:
                raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")
        elif isRange == False:
            if format == 'formula':
                self.write(f"sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).setFormula('{value}')")
            elif format == 'string':
                self.write(f"sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).setString('{value}')")
            elif format == 'number':
                self.write(f"sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).setValue({value})")
            else:
                raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")

    def get_value(self, *args, **kwargs):
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=[], opt_args=['format'])
        if opt['format'] is None:
            format = 'string'
        else:
            format = opt['format']

        if isRange == -1:
            raise ValueError('missing range/cell')
        elif isRange == True:
            if format == 'formula':
                return list(list(x) for x in eval(self.write(f"print(sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).getFormulaArray())")[0]))
            elif format == 'string':
                return list(list(x) for x in eval(self.write(f"print(sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).getDataArray())")[0]))
            elif format == 'number':
                return list(list(x) for x in eval(self.write(f"print(sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).getDataArray())")[0]))
            else:
                raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")
        elif isRange == False:
            if format == 'formula':
                value =  self.write(f"print(sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).getFormula())")
                if len(value) == 1: return value[0]
                else: return ''
            elif format == 'string':
                value = self.write(f"print(sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).getString())")
                if len(value) == 1: return value[0]
                else: return ''
            elif format == 'number':
                return float(self.write(f"print(sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).getValue())")[0])
            else:
                raise ValueError(f"{format} is not a valid format (valid formats: 'formula', 'string', 'number').")

    def clear(self, *args, **kwargs):
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=[], opt_args=[])

        if isRange == -1:
            if query_yes_no('No cell/range was given. Wish to clear all data?', 'no'):
                c = (0, 0, self.get_last_column(), self.get_last_row())
                isRange = True
            else:
                return
        if isRange == True:
            value = [['']*(c[2]-c[0]+1)]*(c[3]-c[1]+1)
            self.set_value(c[0], c[1], c[2], c[3], value)
        elif isRange == False:
            self.set_value(c[0], c[1], '')

    def set_column(self, column, value, row_start=0, format='formula'):
        if type(column) is str:
            try:
                column, row_start = cell2num(column)
            except ValueError:
                column = _check_column_value(column)

        try:
            if len(value[0]) > 0 and type(value[0]) is not str:
                pass
            else:
                value = transpose(value)
        except TypeError:
            value = transpose(value)

        row_start = _check_row_value(row_start)

        # clear
        l = self.get_column_length(column)
        row_stop = l if row_start < l else row_start
        self.clear(column_start=column, row_start=row_start, column_stop=column, row_stop=row_stop)

        return self.set_value(value=value, column_start=column, column_stop=column, row_start=row_start, row_stop=row_start + len(value)- 1, format=format)

    def get_column(self, column, row_start=0, format='string'):
        """
        """
        if type(column) is str:
            try:
                column, row_start = cell2num(column)
            except ValueError:
                column = _check_column_value(column)

        row_stop  = self.get_column_length(column)-1
        return transpose(self.get_value(column_start=column, column_stop=column, row_start=row_start, row_stop=row_stop, format=format))[0]

    def set_row(self, row, value, column_start=0, format='formula'):
        if type(row) is str:
            try:
                column_start, row = cell2num(row)
            except ValueError:
                row = _check_row_value(row)

        try:
            if len(value[0]) > 0 and type(value[0]) is not str:
                pass
            else:
                value = [value]
        except TypeError:
           value = [value]


        column_start = _check_column_value(column_start)
        l = self.get_row_length(row)
        column_stop = l if column_start < l else column_start
        self.clear(column_start=column_start, row_start=row, column_stop=column_stop, row_stop=row)
        return self.set_value(value=value, column_start=column_start, column_stop=column_start + len(value[0])-1, row_start=row, row_stop=row, format=format)

    def get_row(self, row, column_start=0, format='string'):
        """
        """
        if type(row) is str:
            try:
                column_start, row = cell2num(row)
            except ValueError:
                row = _check_row_value(row)

        column_stop  = self.get_row_length(row)-1
        return self.get_value(column_start=column_start, column_stop=column_stop, row_start=row, row_stop=row, format=format)[0]

    def clear_column(self, column, row_start=0):
        if type(column) is str:
            try:
                column, row_start = cell2num(column)
            except ValueError:
                column = _check_column_value(column)

        row_start = _check_row_value(row_start)

        # clear
        l = self.get_column_length(column)
        row_stop = l if row_start < l else row_start
        self.clear(column_start=column, row_start=row_start, column_stop=column, row_stop=row_stop)

    def clear_row(self, row, column_start=0):
        if type(row) is str:
            try:
                column_start, row = cell2num(row)
            except ValueError:
                row = _check_row_value(row)

        column_start = _check_column_value(column_start)
        l = self.get_row_length(row)
        column_stop = l if column_start < l else column_start
        self.clear(column_start=column_start, row_start=row, column_stop=column_stop, row_stop=row)

    def merge(self, *args, **kwargs):
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=[], opt_args=[])
        if isRange == True:
            self.write(f"sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).merge(True)")
        elif isRange == False or isRange == -1:
            raise ValueError('no range to merge')

    def unmerge(self, *args, **kwargs):
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=[], opt_args=[])
        if isRange == True:
            self.write(f"sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).merge(False)")
        elif isRange == False or isRange == -1:
            raise ValueError('no range to unmerge')

    def cell_properties(self, *args, **kwargs):
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=[], opt_args=[])

        if isRange == -1:
            return eval(self.write(f"print([x.Name for x in sheet_{self.tag}.getCellByPosition(0, 0).getPropertySetInfo().Properties])")[0])
        else:
            return eval(self.write(f"print([x.Name for x in sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).getPropertySetInfo().Properties])")[0])

def _letter2num(letter):
    """Returns position of letter in the alphabet starting from 0 (A, B, ..., Z, AA, AB, ..., ZZ, AAA, ...)"""
    letter   = letter.lower()
    alphabet = 'abcdefghijklmnopqrstuvwxyz'
    n = 0
    for idx, s in enumerate(letter):
        n += alphabet.index(s)+(idx)*26
    return n

import re
def _str2num(string):
    """Returns col, row number based on string like 'A2', 'AJ20', ..."""

    # check input
    if isinstance(string, str) == False:
        msgbox(f'_cell2num error: input must be type str\n\ninput = {string}\ninput type={type(string)}')
        return
    if ':' in string:
        start = _str2num(string.split(':')[0])
        stop  = _str2num(string.split(':')[1])
        return start[0], start[1], stop[0], stop[1]
        msgbox(f'_cell2num error: input must be a cell position, not range\n\ninput = {string}')
        return
    
    # cell2num
    try:
        temp = re.compile("([a-zA-Z]+)([0-9]+)")
        res = temp.match(string).groups()
        return _str2num(res[0]), int(res[1])-1
    except AttributeError:
        msgbox(f'_cell2num error: cannot covert input to number\n\ninput = {string}')
        return



def _check_row_value(row):
    if type(row) == str:
        row = int(row)-1

    if row < 0:
        raise ValueError('row cannot be negative')

    return row

def _check_column_value(column):
    if type(column) == str:
        column = _letter2num(column)

    if column < 0:
        raise ValueError('column cannot be negative')
    return column

def _parse_args(args, kwargs, req_args=[], opt_args=[]):
    """
    opt_args must be given in kwargs

    Range
    cell_start, cell_stop
    column_start, row_start, column_stop, row_stop

    cell
    column, row
    """
    req = {k:False for k in req_args}
    opt = {k:None for k in opt_args}

    for i, k in enumerate(opt_args):
        if k in kwargs:
            opt[k] = kwargs[k]

    lookout = 0
    del_list = []
    for i, k in enumerate(req_args):
        try:
            req[k] = kwargs[k]
            del_list.append(k)
        except KeyError:
            lookout += 1
    for k in del_list:
        del req_args[req_args.index(k)]

    for i in range(lookout):
        req[req_args[::-1][i]] = args[::-1][i]
    if lookout != 0:
        args = args[:-lookout]

    isRange = False
    if len(args) == 0:
        isRange = -1
    elif len(args) == 1:
        try:
            column_start, row_start = cell2num(args[0])
        except ValueError:
            column_start, row_start, column_stop, row_stop = range2num(args[0])
            isRange = True
    elif len(args) == 2:
        try:
            column_start, row_start = cell2num(args[0])
            column_stop, row_stop   = cell2num(args[1])
            isRange = True
        except TypeError:
            column_start = args[0]
            row_start    = args[1]
        except ValueError:
            column_start = args[0]
            row_start    = args[1]
    elif len(args) == 3:
        raise ValueError('cannot parse arguments.')
    elif len(args) == 4:
        column_start = args[0]
        row_start    = args[1]
        column_stop  = args[2]
        row_stop     = args[3]
        isRange = True

    if 'column' in kwargs:
        column_start = kwargs['column']
        isRange      = False
    if 'row' in kwargs:
        row_start = kwargs['row']
        isRange   = False
    if 'column_start' in kwargs:
        column_start = kwargs['column_start']
        isRange      = False
    if 'row_start' in kwargs:
        row_start = kwargs['row_start']
        isRange   = False
    if 'column_stop' in kwargs:
        column_stop = kwargs['column_stop']
        isRange     = True
    if 'row_stop' in kwargs:
        row_stop = kwargs['row_stop']
        isRange  = True
    if 'cell' in kwargs:
        column_start, row_start = cell2num(kwargs['cell'])
        isRange = False
    if 'cell_start' in kwargs:
        column_start, row_start = cell2num(kwargs['cell_start'])
        isRange = False
    if 'cell_stop' in kwargs:
        column_stop, row_stop = cell2num(kwargs['cell_stop'])
        isRange = True
    if 'Range' in kwargs:
        column_start, row_start, column_stop, row_stop = range2num(kwargs['Range'])
        isRange = True

    if isRange == -1:
        return -1, (), req, opt
    else:
        column_start = _check_column_value(column_start)
        row_start    = _check_row_value(row_start)
        if isRange:
            column_stop = _check_column_value(column_stop)
            row_stop    = _check_row_value(row_stop)

            if column_stop < column_start:
                raise ValueError(f'column_start ({column_start}) cannot be bigger than column_stop ({column_stop})')
            if row_stop < row_start:
                raise ValueError(f'row_start ({row_start}) cannot be bigger than row_stop ({row_stop})')

            return True, (column_start, row_start, column_stop, row_stop), req, opt
        else:
            return False, (column_start, row_start), req, opt




    def _get_property_function(self, column, row, name):

        string = self.write(f"print(sheet_{self.tag}.getCellByPosition({column}, {row}).getPropertyValue('{name}'))")

        if len(string) > 0:
            if '{' in string[0]:
                return string[0].split('{')[0][1:-1].split('.')[-1]
            else: return None
        else:
            return None

    def set_property(self, *args, **kwargs):
        """
        name = 'property.subproperty'
        value = simple value

        name = 'property'
        value = dict
        """
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=['name', 'value'], opt_args=[])
        name = req['name']
        value = req['value']

        if isRange == -1:
            raise ValueError('missing range/cell')

        if name.count('.') == 0:
            if isinstance(value, list) or isinstance(value, np.ndarray) or isinstance(value, tuple):
                raise ValueError(f'Property ({name}) requires a value that is not a list, tuple, or array.')
            elif isinstance(value, dict):

                d = self.get_property(c[0], c[1], name)
                for k, v in value.items():
                    if k in d:
                        d[k] = v
                    else:
                        raise ValueError(f'{k} is not a valid option of {name}.')

                func = self._get_property_function(c[0], c[1], name)
                f_string = f"{func}(" + ','.join([f"{k}={v}" for k, v in d.items()]) + ")"

                if isRange == False:
                    self.write(f"sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).setPropertyValue('{name}', {f_string})")
                else:
                    self.write(f"sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).setPropertyValue('{name}', {f_string})")
            else:
                if isRange == False:
                    self.write(f"sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).setPropertyValue('{name}', {value})")
                else:
                    self.write(f"sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).setPropertyValue('{name}', {value})")


        elif name.count('.') == 1:
            name0 = name.split('.')[0]
            name1 = name.split('.')[1]

            func = self._get_property_function(c[0], c[1], name0)

            if func is None:
                raise ValueError(f'Property name given ({name}) has a dot (.), which mean it is a nested property, but actual property is not nested.')
            else:
                d = self.get_property(c[0], c[1], name0)
                if name1 in d:
                    d[name1] = value
                else:
                    raise ValueError(f'{name1} is not a valid option of {name0}.')

                f_string = f"{func}(" + ','.join([f"{k}={v}" for k, v in d.items()]) + ")"
                if isRange == False:
                    self.write(f"sheet_{self.tag}.getCellByPosition({c[0]}, {c[1]}).setPropertyValue('{name0}', {f_string})")
                else:
                    self.write(f"sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).setPropertyValue('{name0}', {f_string})")
        else:
            raise ValueError(f'Only single nested properties can be edited. Note that complex properties can always be edited through other simple nested properties.')

    def get_property(self, *args, **kwargs):
        """
        cannot get values from here:

        UserDefinedAttributes
        NumberingRules
        Validation
        ValidationLocal
        ValidationXML
        ConditionalFormat
        ConditionalFormatLocal
        ConditionalFormatXML
        """
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=['name'], opt_args=[])
        name = req['name']

        if isRange == -1:
            raise ValueError('missing range/cell')
        elif isRange == True:
            t = self.write(f"print(type(get_cells_property_recursively({self.tag}, {c[0]}, {c[1]}, {c[2]}, {c[3]}, '{name}')))")[0]
            try:
                output = self.write(f"print(get_cells_property_recursively({self.tag}, {c[0]}, {c[1]}, {c[2]}, {c[3]}, '{name}'))")[0]
            except IndexError:
                return None
        elif isRange == False:
            t = self.write(f"print(type(get_cell_property_recursively({self.tag}, {c[0]}, {c[1]}, '{name}')))")[0]
            try:
                output = self.write(f"print(get_cell_property_recursively({self.tag}, {c[0]}, {c[1]}, '{name}'))")[0]
            except IndexError:
                return None

        if  t == "<class 'int'>" :
            return int(output)
        elif t == "<class 'bool'>" :
            return str2bool(output)
        elif t == "<class 'dict'>":
            return eval(output)
        else:
            return output

    def remove_conditional_format(self, *args, **kwargs):
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=[], opt_args=['range_index', 'index'])

        if opt['index'] is not None:
            if isRange == -1:
                d = self.get_conditional_formats()
                for r in d[range_index]['Ranges']:
                    self.write(f"""cf = sheet_{self.tag}.getCellRangeByPosition({r[0]}, {r[1]}, {r[2]}, {r[3]}).ConditionalFormat\n"""+\
                                f"""cf.removeByIndex({index})\n"""+\
                                f"""sheet_{self.tag}.getCellRangeByPosition({r[0]}, {r[1]}, {r[2]}, {r[3]}).ConditionalFormat = cf""")
            elif isRange == True:
                self.write(f"""cf = sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).ConditionalFormat\n"""+\
                            f"""cf.removeByIndex({index})\n"""+\
                            f"""sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).ConditionalFormat = cf""")
            elif isRange == False:
                self.write(f"""cf = sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[0]}, {c[1]}).ConditionalFormat\n"""+\
                            f"""cf.removeByIndex({index})\n"""+\
                            f"""sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[0]}, {c[1]}).ConditionalFormat = cf""")
            else:
                raise valueError('range not defined.')
        else:
            if isRange == -1:
                if opt['range_index'] is None:
                    raise ValueError('range (or range_index) not defined.')
                else:
                    range_index = opt['range_index']
                d = self.get_conditional_formats()
                for r in d[range_index]['Ranges']:
                    self.write(f"""cf = sheet_{self.tag}.getCellRangeByPosition({r[0]}, {r[1]}, {r[2]}, {r[3]}).ConditionalFormat\n"""+\
                                f"""cf.clear()\n"""+\
                                f"""sheet_{self.tag}.getCellRangeByPosition({r[0]}, {r[1]}, {r[2]}, {r[3]}).ConditionalFormat = cf""")
            elif isRange == True:
                self.write(f"""cf = sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).ConditionalFormat\n"""+\
                            f"""cf.clear()\n"""+\
                            f"""sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).ConditionalFormat = cf""")
            elif isRange == False:
                self.write(f"""cf = sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[0]}, {c[1]}).ConditionalFormat\n"""+\
                            f"""cf.clear()\n"""+\
                            f"""sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[0]}, {c[1]}).ConditionalFormat = cf""")
            else:
                raise valueError('range not defined.')

    def get_conditional_formats(self, *args, **kwargs):
        output = self.write(f"""d = dict()\n"""+\
                    f"""for i in range(sheet_{self.tag}.ConditionalFormats.getLength()):\n"""+\
                    f"""    r = [(r.StartColumn, r.StartRow, r.EndColumn, r.EndRow) for r in sheet_{self.tag}.ConditionalFormats.ConditionalFormats[i].Range.getRangeAddresses()]\n"""+\
                    f"""    id = sheet_{self.tag}.ConditionalFormats.ConditionalFormats[i].ID-1\n"""+\
                    f"""    d[id] = {{'ConditionalFormats':{{}}, 'Ranges':r}}\n"""+\
                    f"""    for j in range(sheet_{self.tag}.ConditionalFormats.ConditionalFormats[i].getCount()):\n"""+\
                    f"""        d[id]['ConditionalFormats'][j] = dict()\n"""+\
                    f"""        cf = sheet_{self.tag}.ConditionalFormats.ConditionalFormats[i].getByIndex(j)\n"""+\
                    f"""        d[id]['ConditionalFormats'][j]['Operator'] = cf.getPropertyValue('Operator')\n"""+\
                    f"""        d[id]['ConditionalFormats'][j]['Formula1'] = cf.getPropertyValue('Formula1')\n"""+\
                    f"""        d[id]['ConditionalFormats'][j]['Formula2'] = cf.getPropertyValue('Formula2')\n"""+\
                    f"""        d[id]['ConditionalFormats'][j]['StyleName'] = cf.getPropertyValue('StyleName')\n"""+\
                    # f"""        d[id]['ConditionalFormats'][j]['SourcePosition'] = cf.getPropertyValue('SourcePosition')\n"""+\
                    f"""print(d)""")[-1]
        d = eval(output)
        for i in d:
            for j in d[i]['ConditionalFormats']:
                d[i]['ConditionalFormats'][j]['Operator'] = operator_map[d[i]['ConditionalFormats'][j]['Operator']][0]

        if len(args) == 0 and len(kwargs)==0:
            return d
        else:
            isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=[], opt_args=[])
            if isRange == True:
                column_start, row_start, column_stop, row_stop = c
                for i in d:
                    for r in d[i]['Ranges']:
                        if (column_start, row_start, column_stop, row_stop) == r:
                            d_filtered[i] = d[i]
                return d_filtered

            elif isRange == False:
                column, row = c
                d_filtered = {}
                for i in d:
                    for r in d[i]['Ranges']:
                        if column >= r[0] and column <= r[2]:
                            if row >= r[1] and row <= r[3]:
                                d_filtered[i] = d[i]
                return d_filtered

    def new_conditional_format(self, *args, **kwargs):
        isRange, c, req, opt = _parse_args(args=args, kwargs=kwargs, req_args=['Operator', 'Formula1', 'StyleName'], opt_args=['index_range', 'Formula2'])

        # adjust conditinal parameters
        if type(req['Operator']) is str:
            Operator = [operator_map[k][-1] for k in operator_map if req['Operator'] in operator_map[k]][0]
            if Operator == '':
                # req['Operator'] = 'NONE'
                Opertor = req['Operator']
                raise NotImplementedError(f'conditional format of type "{Operator}" not implemented yet.')
            else:
                req['Operator'] = Operator
        req['Formula1'] = str(req['Formula1'])
        if opt['Formula2'] is not None:
             Formula2 = str(opt['Formula2'])
        else:
            Formula2 = '0'


        Opertor = req['Operator']
        Formula1 = req['Formula1']
        StyleName = req['StyleName']
        if isRange == -1:
            if 'index_range' not in opt:
                raise valueError('range/cell not defined.')
            else:
                d = self.get_conditional_formats()
                for c in d[range_index]['Ranges']:
                    self.write(f"""args = dict(Operator = {Operator}, Formula1 = '{Formula1}', Formula2 = '{Formula2}', StyleName = '{StyleName}')\n"""+\
                                f"""pv = [PropertyValue(Name=n, Value=v) for n, v in args.items()]\n"""+\
                                f"""cf = sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).ConditionalFormat\n"""+\
                                f"""cf.addNew(pv)\n"""+\
                                f"""sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).ConditionalFormat = cf""")
        elif isRange == True:
            self.write(f"""args = dict(Operator = {Operator}, Formula1 = '{Formula1}', Formula2 = '{Formula2}', StyleName = '{StyleName}')\n"""+\
                        f"""pv = [PropertyValue(Name=n, Value=v) for n, v in args.items()]\n"""+\
                        f"""cf = sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).ConditionalFormat\n"""+\
                        f"""cf.addNew(pv)\n"""+\
                        f"""sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[2]}, {c[3]}).ConditionalFormat = cf""")
        elif isRange == False:
            self.write(f"""args = dict(Operator = {Operator}, Formula1 = '{Formula1}', Formula2 = '{Formula2}', StyleName = '{StyleName}')\n"""+\
                        f"""pv = [PropertyValue(Name=n, Value=v) for n, v in args.items()]\n"""+\
                        f"""cf = sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[0]}, {c[1]}).ConditionalFormat\n"""+\
                        f"""cf.addNew(pv)\n"""+\
                        f"""sheet_{self.tag}.getCellRangeByPosition({c[0]}, {c[1]}, {c[0]}, {c[1]}).ConditionalFormat = cf""")
        else:
            raise valueError('range/cell not defined.')
