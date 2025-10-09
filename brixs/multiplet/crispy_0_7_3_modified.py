#! /usr/bin/env python3
# -*- coding: utf-8 -*-
r"""Multiplet calculation of XPS, XES, XAS, and RIXS via crispy and Quanty.

Notes:
    This module requires a Quanty executable (Quanty.exe) that can be downloaded
    at https://www.quanty.org/index.html

Warning:
    This module was only tested on windows and for a limited number of calculation parameters.
    It is good practice to compare the results from this module with results 
    obtained with the latest version of Crispy.

Notes:
    If you use this module on publications. Please, reference Crispy 0.7.3 as this
    module is a modified version of Crispy.

Usage:
    See brixs/examples/multiplet folder
    
TO DO:
    ( ) Check if broaden works with non-monotonic data.
    ( ) implement new shapes for q.plot_geometry() 
    ( ) q.plot.geometry() brakes depending on kin, LHin, and LVin
    ( ) Implement linear RIXS for 1s2p and 1s3p
    ( ) Implement circular RIXS for 1s2p and 1s3p
"""

from collections.abc import Iterable, MutableMapping
from pathlib import Path
import collections
import subprocess
import tempfile
import json
import copy
import os

import numpy as np

import brixs as br


# %% ============================== settings ============================= %% #
class _settings():

    def __init__(self):
        self._QUANTY_FILEPATH      = ''
        self._WSL_NTHREADS         = 1
        self._WSL_QUANTY_FILEPATH  = ''
        self._WSL_TEMP_FILEPATH    = ''
        self._TEMPLATES_FOLDERPATH = ''
        self._PARAMETERS_FILEPATH  = ''
        self._PARAMETERS           = None
        self._ELEMENTS             = None

    @property
    def QUANTY_FILEPATH(self):
        return self._QUANTY_FILEPATH
    @QUANTY_FILEPATH.setter
    def QUANTY_FILEPATH(self, value):
        value = Path(value)
        assert value.exists(), 'Cannot find filepath'
        assert value.is_file(), 'filepath does not point to a file'
        self._QUANTY_FILEPATH = value
    @QUANTY_FILEPATH.deleter
    def QUANTY_FILEPATH(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def WSL_NTHREADS(self):
        return self._WSL_NTHREADS
    @WSL_NTHREADS.setter
    def WSL_NTHREADS(self, value):
        value = int(round(value))
        vmax = os.cpu_count()
        assert value<=vmax and value>1, f'Invalid (WSL_NTHREADS={value}). It must be a integer higher or equal to 1 and less or equal to the max number of available threads on the machine (vmax={vmax})'
        self._WSL_NTHREADS = value
    @WSL_NTHREADS.deleter
    def WSL_NTHREADS(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def WSL_QUANTY_FILEPATH(self):
        return self._WSL_QUANTY_FILEPATH
    @WSL_QUANTY_FILEPATH.setter
    def WSL_QUANTY_FILEPATH(self, value):
        value = Path(value)
        assert value.exists(), 'Cannot find filepath'
        assert value.is_file(), 'filepath does not point to a file'
        self._WSL_QUANTY_FILEPATH = value
    @WSL_QUANTY_FILEPATH.deleter
    def WSL_QUANTY_FILEPATH(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def WSL_TEMP_FILEPATH(self):
        return self._WSL_TEMP_FILEPATH
    @WSL_TEMP_FILEPATH.setter
    def WSL_TEMP_FILEPATH(self, value):
        value = Path(value)
        assert value.parent.exists(), 'Cannot find folderpath'
        # assert value.is_file(), 'filepath does not point to a file'
        self._WSL_TEMP_FILEPATH = value
    @WSL_TEMP_FILEPATH.deleter
    def WSL_TEMP_FILEPATH(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def TEMPLATES_FOLDERPATH(self):
        return self._TEMPLATES_FOLDERPATH
    @TEMPLATES_FOLDERPATH.setter
    def TEMPLATES_FOLDERPATH(self, value):
        value = Path(value)
        assert value.exists(), 'Cannot find folderpath'
        assert value.is_dir(), 'folderpath does not point to a folder'
        self._TEMPLATES_FOLDERPATH = value
    @TEMPLATES_FOLDERPATH.deleter
    def TEMPLATES_FOLDERPATH(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def PARAMETERS_FILEPATH(self):
        return self._PARAMETERS_FILEPATH
    @PARAMETERS_FILEPATH.setter
    def PARAMETERS_FILEPATH(self, value):
        value = Path(value)
        assert value.exists(), 'Cannot find filepath'
        assert value.is_file(), 'filepath does not point to a file'
        self._PARAMETERS_FILEPATH = value
    @PARAMETERS_FILEPATH.deleter
    def PARAMETERS_FILEPATH(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def PARAMETERS(self):
        return self._PARAMETERS
    @PARAMETERS.setter
    def PARAMETERS(self, value):
        raise AttributeError('Cannot modify object.')
    @PARAMETERS.deleter
    def PARAMETERS(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def ELEMENTS(self):
        return self._ELEMENTS
    @ELEMENTS.setter
    def ELEMENTS(self, value):
        raise AttributeError('Cannot modify object.')
    @ELEMENTS.deleter
    def ELEMENTS(self):
        raise AttributeError('Cannot delete object.')
    

    def __str__(self):
        text  = ''
        text += f'QUANTY_FILEPATH:      {self.QUANTY_FILEPATH}\n' +\
                f'WSL_NTHREADS:         {self.WSL_NTHREADS}\n' +\
                f'WSL_QUANTY_FILEPATH:  {self.WSL_QUANTY_FILEPATH}\n' +\
                f'WSL_TEMP_FILEPATH:    {self.WSL_TEMP_FILEPATH}\n' +\
                f'TEMPLATES_FOLDERPATH: {self.TEMPLATES_FOLDERPATH}\n' +\
                f'PARAMETERS_FILEPATH:  {self.PARAMETERS_FILEPATH}\n' +\
                f'ELEMENTS:             {self.ELEMENTS}\n' +\
                f'PARAMETERS:           returns full parameter list for each valid element'
                # f'default: {self.FIGURE_FORCE_NEW_WINDOW}'
        return text

    def __repr__(self):
        return self.__str__()

    def help(self):
        print(self._help)
    
    def pretty_print(self):
        print(self.__str__())
settings = _settings()
# %%

# %% ========================== Operating system ========================= %% #
import platform
system = platform.system().lower()
is_windows = system == 'windows'
is_linux   = system == 'linux'
is_mac     = system == 'darwin'
# %%

# %% =============================== odict =============================== %% #
class odict(collections.OrderedDict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value
# %%

# %% ========================== support classes ========================== %% #
class _hamiltonianState(MutableMapping):
    def __init__(self, initial, parent):
        self.store = dict(initial)
        self.parent = parent

        # change from int to bool
        for name in self.store:
            if self.store[name]:
                self.store[name] = True
            else:
                self.store[name] = False

    def __str__(self):
        return str(self.store).replace(', ', '\n ')

    def __repr__(self):
        # return str(self.store)
        return str(self.store).replace(', ', '\n ')

    def __getitem__(self, name):
        return self.store[name]

    def __setitem__(self, name, value):
        if name not in self:
            raise ValueError(f'New terms cannot be created.\nInvalid term: {name}\nValid terms are: {list(self.keys())}')
        assert isinstance(value, bool), f'hamiltonianState[{name}] must be a bool (True or False), not type: {type(value)}'

        if self.store[name] == False and value == True:  # from False to True
            
            # check if LMCT and MLCT are True at the same time
            if 'LMCT' in name:
                # find MLCT
                for _opposite in self:
                    if 'MLCT' in _opposite:
                        break
                
                # check if opposite is also True
                if self.store[_opposite]:
                    raise ValueError('Ligands Hybridization LMCT and MLCT cannot be both True.\nIf you want to set LMCT to True, you have to set MLCT to False first')
            elif 'MLCT' in name:
                # find MLCT
                for _opposite in self:
                    if 'LMCT' in _opposite:
                        break

                # check if opposite is also True
                if self.store[_opposite]:
                    raise ValueError('Ligands Hybridization MLCT and LMCT cannot be both True.\nIf you want to set MLCT to True, you have to set LMCT to False first')
            
            # set value
            self.store[name] = True

            # Determine the maximum number of allowed configurations.
            if 'LMCT' in name:
                if 'd' in self.parent.block:
                    self.parent._nConfigurationsMax = 10 - self.parent.nElectrons + 1
                elif 'f' in self.parent.block:
                    self.parent._nConfigurationsMax = 14 - self.parent.nElectrons + 1
                self.parent.nConfigurations = 2
                print(f'LMCT changed from False to True: nConfigurationsMax changed from 1 to {self.parent.nConfigurationsMax}')
                print(f'LMCT changed from False to True: nConfigurations changed from 1 to 2')
            elif 'MLCT' in name:
                self.parent._nConfigurationsMax = self.parent.nElectrons + 1
                self.parent.nConfigurations = 2
                print(f'MLCT changed from False to True: nConfigurationsMax changed from 1 to {self.parent.nConfigurationsMax}')
                print(f'MLCT changed from False to True: nConfigurations changed from 1 to 2')

        elif self.store[name] == True and value == False:  # from True to False
            # set value
            self.store[name] = False
            
            # Determine the maximum number of allowed configurations.
            if 'LMCT' in name:
                print(f'LMCT changed from True to False: nConfigurationsMax changed from {self.parent.nConfigurationsMax} to 1')
                print(f'LMCT changed from True to False: nConfigurations changed from {self.parent.nConfigurations} to 1')
                self.parent._nConfigurationsMax = 1
                self.parent.nConfigurations    = 1
            elif 'MLCT' in name:
                self.parent._nConfigurationsMax = 1
                self.parent.nConfigurations    = 1
                print(f'MLCT changed from True to False: nConfigurationsMax changed from {self.parent.nConfigurationsMax} to 1')
                print(f'MLCT changed from True to False: nConfigurations changed from {self.parent.nConfigurations} to 1')
        else:
            pass

    def __delitem__(self, key):
        raise AttributeError('Items cannot be deleted')

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len({name:self.store[name] for name in self.store if self.store[name] is not None})

    def export(self):
        temp = {}
        for item in self.store:
            if self.store[item]:
                temp[item] = 2
            else:
                temp[item] = 0
        return temp

class _hamiltonianData(MutableMapping):
    def __init__(self, initial):
        self.store = dict(copy.deepcopy(initial))
        for key in self.store:
            if type(self.store[key]) != float and type(self.store[key]) != list:
                self.store[key] = _hamiltonianData(self.store[key])
        self.initial = dict(copy.deepcopy(initial))

        self.sync = False

    def __str__(self):
        # toprint = str(self.store).replace('\':', '\':\n'+' '*10).replace(', \'', ', \n'+' '*0+'\'')
        # toprint = toprint.replace('(\'Final', '\n'+' '*11+'\'Final').replace('), (\'', ',\n'+' '*20+'\'')
        # toprint = toprint.replace('odict([(', '').replace('Hamiltonian\', \'', 'Hamiltonian\',\n'+' '*20+'\'')
        # toprint = toprint.replace(')]))])', '')
        # toprint = toprint.replace(')]))', '')
        # return toprint

        i = 8
        toprint = str(self.store).replace('\': {\'Initial', '\':\n'+' '*i+'{\'Initial')  # initial
        toprint = toprint.replace('}, \'Intermediate', '\':\n'+' '*i+'\'Intermediate')  # Intermediate
        toprint = toprint.replace('}, \'Final', '},\n'+' '*(i+1)+'\'Final')  # final
        toprint = toprint.replace('}}, \'', '\n \'')  # terms
        toprint = toprint.replace(', \'', '\n'+' '*i*2 +'\'')  # values
        toprint = toprint.replace('\': {\'', '\':\n'+' '*i*2 +'\'')  # values from the first line
        return toprint

    def __repr__(self):
        return str(self.store)

    def __getitem__(self, name):
        return self.store[name]

    def __setitem__(self, name, value):
        if name not in self:
            raise ValueError(f'New terms cannot be created.\nInvalid term: {name}\nValid terms are: {list(self.keys())}')
        if type(self.store[name]) != type(value):
            if (type(self.store[name]) == float and type(value) == int) or (type(self.store[name]) == int and type(value) == float):
                pass
            else:
                raise ValueError(f'New value does not have a valid type.\nOld value: {self.store[name]}\nOld type: {type(self.store[name])}\nNew value: {value}\nNew type: {type(value)}')
        self.store[name] = value

    def __delitem__(self, key):
        raise AttributeError('Itens cannot be deleted')

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len({name:self.store[name] for name in self.store if self.store[name] is not None})

    def reset(self):
        self.store = copy.deepcopy(self.initial)
# %%

# %% ========================= support functions ========================= %% #
def quanty(filepath, wsl=False):
    """Run Quanty.

    Args:
        filepath (string or pathlib.Path): path to file.

    Returns:
        Calculation output (stdout).
    """
    if wsl:
        assert is_windows, 'Failed to recognized windows machine. One must be on a windows machine to run quanty with wsl'
        assert settings.WSL_QUANTY_FILEPATH != '', 'settings.WSL_QUANTY_FILEPATH not set'
        assert settings.WSL_TEMP_FILEPATH != '', 'settings.WSL_TEMP_FILEPATH not set'

        quanty_exe = Path(settings.WSL_QUANTY_FILEPATH).absolute()
        filepath   = Path(filepath).absolute()
        temporary  = Path(settings.WSL_TEMP_FILEPATH).absolute()

        # converting filepath from windows to wsl
        drive = quanty_exe.home().drive[0].lower()  # get drive letter
        quanty_exe = '/'.join(['/mnt', drive] + list(quanty_exe.parts[1:]))
        drive = filepath.home().drive[0].lower()  # get drive letter
        filepath = '/'.join(['/mnt', drive] + list(filepath.absolute().parts[1:]))
        drive = temporary.home().drive[0].lower()  # get drive letter
        temporary = '/'.join(['/mnt', drive] + list(temporary.absolute().parts[1:]))
        

        # threads
        NTHREADS = int(settings.WSL_NTHREADS)

        # run
        # print(quanty_exe)
        # print(filepath)
        quanty = subprocess.Popen(r'C:\Windows\System32\wsl.exe -e sh -c "export OMP_NUM_THREADS=' + f'{NTHREADS};' + f' {quanty_exe} {filepath}' + r'>' + f'{temporary}' + r'"', shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         # quanty = subprocess.Popen(['wsl', f'export OMP_NUM_THREADS={NTHREADS}', f"{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         # quanty = subprocess.Popen([f'wsl ~ -e sh -c ' + f'export OMP_NUM_THREADS={NTHREADS}; {quanty_exe} {filepath}'], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         # quanty = subprocess.Popen([f'wsl ~ -e sh -c ' + f'export OMP_NUM_THREADS={NTHREADS}; {quanty_exe}'], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         # quanty = subprocess.Popen([f'wsl ~ -e sh -c "' + f'export OMP_NUM_THREADS={NTHREADS}' + '"'], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         import subprocess
#         # quanty = subprocess.Popen([f'C:\Windows\System32\wsl.exe ~'], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         quanty = subprocess.Popen([r'C:\Windows\System32\wsl.exe'], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         quanty = subprocess.Popen([r'C:\Windows\System32\wsl.exe'], close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         quanty = subprocess.Popen([r'C:\Windows\System32\wsl.exe'], close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         quanty = subprocess.Popen(r'C:\Windows\System32\wsl.exe', shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         quanty = subprocess.Popen(r'C:\Windows\System32\wsl.exe -e sh -c "export OMP_NUM_THREADS=2; ' + r'/mnt/c/Users/galdin_c/github/quanty/quanty_lin/Quanty' + '"', shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# r'/mnt/c/Users/galdin_c/github/quanty/quanty_lin/Quanty'
# C:\Windows\System32\wsl.exe -e sh -c "export OMP_NUM_THREADS=2; /mnt/c/Users/galdin_c/github/quanty/quanty_lin/Quanty /mnt/c/Users/galdin_c/AppData/Local/Temp/tmpfh8gkcgi"
# r'/mnt/c/Users/galdin_c/AppData/Local/Temp/tmpfh8gkcgi
# quanty = subprocess.Popen(r'C:\Windows\System32\wsl.exe -e sh -c "export OMP_NUM_THREADS=' + f'{16};' + f' /mnt/c/Users/galdin_c/github/quanty/quanty_lin/Quanty /mnt/c/Users/galdin_c/AppData/Local/Temp/tmpfh8gkcgi' + r' > /mnt/c/Users/galdin_c/Downloads/temp.txt"', shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# export OMP_NUM_THREADS=2
# /mnt/c/Users/galdin_c/github/quanty/quanty_lin/Quanty /mnt/c/Users/galdin_c/AppData/Local/Temp/tmpfh8gkcgi
    else:
        assert settings.QUANTY_FILEPATH != '', 'settings.QUANTY_FILEPATH not set'

        quanty_exe = Path(settings.QUANTY_FILEPATH).absolute()
        filepath = Path(filepath).absolute()

        if is_windows:
            quanty = subprocess.Popen([f"{quanty_exe}", f"{filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # if quanty_exe.is_absolute():
            #     quanty = subprocess.Popen([f"{quanty_exe}", f"{filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # else:
            #     quanty = subprocess.Popen([f"./{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif is_linux:
            quanty = subprocess.Popen([f"{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # if quanty_exe.is_absolute():
            #     quanty = subprocess.Popen([f"{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # else:
            #     quanty = subprocess.Popen([f"./{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        elif is_mac:
            quanty = subprocess.Popen([f"./{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if wsl:
        output = quanty.stdout.read().decode("utf-8")
        with Path(settings.WSL_TEMP_FILEPATH).open() as f:
            output = f.read()
        # br.rm(Path(settings.WSL_TEMP_FILEPATH))
    else:
        output = quanty.stdout.read().decode("utf-8")
    error  = quanty.stderr.read().decode("utf-8")
    if error != '':
        raise RuntimeError(f"Error while reading file: {filepath}\nERROR: {error}")
    if 'Error while loading the script' in output:
        error = output[output.find('Error while loading the script')+len('Error while loading the script:')+1:]
        raise ValueError(f'Error while loading file: {filepath}. \n {error}')
    return output

def load_calculation(filepath):
    """Load calculation.

    Args:
        filepath (string or pathlib.Path): filepath.

    Returns:
        None
    """
    filepath = Path(filepath)

    with open(str(filepath), 'r') as file:
        par = json.load(file)

    hamiltonian = par.pop('hamiltonian')
    hamiltonianData  = hamiltonian.pop('hamiltonianData')
    hamiltonianState = hamiltonian.pop('hamiltonianState')

    initial = par.pop('initial')
    q = Calculation(element    = initial.pop('element'),
                    charge     = initial.pop('charge'),
                    symmetry   = initial.pop('symmetry'),
                    experiment = initial.pop('experiment'),
                    edge       = initial.pop('edge'),)
    geometry = par.pop('geometry')
    q._tth      = geometry.pop('tth')
    q._R        = geometry.pop('R')
    q._kin      = geometry.pop('kin')
    q._LHin     = geometry.pop('LHin')
    q._LVin     = geometry.pop('LVin')
    q._kout     = geometry.pop('kout')
    q._LHout    = geometry.pop('LHout')
    q._LVout    = geometry.pop('LVout')
    q._kin_cf   = geometry.pop('kin_cf')
    q._LVin_cf  = geometry.pop('LVin_cf')
    q._LHin_cf  = geometry.pop('LHin_cf')
    q._kout_cf  = geometry.pop('kout_cf')
    q._LVout_cf = geometry.pop('LVout_cf')
    q._LHout_cf = geometry.pop('LHout_cf')
    q._cf       = geometry.pop('cf')

    experiment = par.pop('experiment')
    q._temperature              = experiment.pop('temperature')
    q._magneticField            = experiment.pop('magneticField')
    q._magneticFieldOrientation = experiment.pop('magneticFieldOrientation')
                    
    calculation = par.pop('calculation')
    q._templateName   = calculation.pop('templateName')
    q._templatePath   = calculation.pop('templatePath')
    q._configurations = calculation.pop('configurations')
    q._nConfigurationsMax = calculation.pop('nConfigurationsMax')
    q._block          = calculation.pop('block')
    q._nElectrons     = calculation.pop('nElectrons')
    q._polarization   = calculation.pop('polarization')
    q._nPsis          = calculation.pop('nPsis')
    q._nPsisMax       = calculation.pop('nPsisMax')
    q._nPsisAuto      = calculation.pop('nPsisAuto')
    q._xLabel         = calculation.pop('xLabel')
    q._yLabel         = calculation.pop('yLabel')
    q._xMin           = calculation.pop('xMin')
    q._xMax           = calculation.pop('xMax')
    q._xNPoints       = calculation.pop('xNPoints')
    q._resonance      = calculation.pop('resonance')
    q._gamma1         = calculation.pop('gamma1')
    q._gamma2         = calculation.pop('gamma2')
    q._verbosity      = calculation.pop('verbosity')
    q._denseBorder    = calculation.pop('denseBorder')

    for k1 in hamiltonianState:
        q.hamiltonianState[k1] = hamiltonianState[k1]
    for k1 in hamiltonianData:
        for k2 in hamiltonianData[k1]:
            for k3 in hamiltonianData[k1][k2]:
                q.hamiltonianData[k1][k2][k3] = hamiltonianData[k1][k2][k3]
    return q

def remove_greek_letters_from_hamiltonianData(hamiltonianData):
    """Remove greek letters from hamiltonian data"""
    a = hamiltonianData
    new = {}
    for k1 in a:
        new[k1] = {}
        for k2 in a[k1]:
            new[k1][k2] = {}
            for k3 in a[k1][k2]:
                _k3 = k3.replace('ζ', 'zeta')
                _k3 = _k3.replace('Δ', 'Delta')
                _k3 = _k3.replace('σ', 'sigma')
                _k3 = _k3.replace('τ', 'tau')
                _k3 = _k3.replace('μ', 'mu')
                _k3 = _k3.replace('ν', 'nu')
                new[k1][k2][_k3] = a[k1][k2][k3]
    return new

def put_back_greek_letters_from_hamiltonianData(hamiltonianData):
    """Remove greek letters from hamiltonian data"""
    a = hamiltonianData
    new = {}
    for k1 in a:
        new[k1] = {}
        for k2 in a[k1]:
            new[k1][k2] = {}
            for k3 in a[k1][k2]:
                _k3 = k3.replace('zeta', 'ζ')
                _k3 = _k3.replace('Delta', 'Δ')
                _k3 = _k3.replace('sigma', 'σ')
                _k3 = _k3.replace('tau', 'τ')
                _k3 = _k3.replace('mu', 'μ')
                _k3 = _k3.replace('nu', 'ν')
                new[k1][k2][_k3] = a[k1][k2][k3]
    return new
# %%

# %% ============== Calculation (adaptation of crispy 0.7.3) ============= %% #
class Calculation(object):
    """Initialize a Crispy/Quanty calculation object.

    Please refer to brixs/examples/multiplet folder
    """

    def __init__(self, element    = None,
                       charge     = None,
                       symmetry   = None,
                       experiment = None,
                       edge       = None,
                       #
                       polarization = 'Isotropic',
                       nPsis        = None,
                       nPsisAuto    = 1,
                       k            = (0, 0, 1),  # quantization axis for calculating expected values (it does not influence the calculation!!)
                       #  nConfigurations = 1,
                       denseBorder = 2000,
                       #
                       magneticField = 0.002,
                       magneticFieldOrientation = [0, 0, 1],
                       temperature   = 10,
                       #
                       tth   = 130,
                       R     = [],
                       #
                       E        = None,  # only for RIXS, incident photon energy. If None, the resonance value will be selected
                       xMin     = None,
                       xMax     = None,
                       xNPoints = None,
                       #
                       gamma1 = None, #0.1, # core-hole lifetime in eV (For RIXS: incident energy core-hole lifetime in eV)
                       gamma2 = None, # For RIXS: energy loss core-hole lifetime in eV
                       ):
        # primary attributes (cannot be modified)
        self._set_primary_attributes(element, charge, symmetry, experiment, edge)

        branch = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]
        self._configurations  = branch['configurations']
        self._block           = self.configurations[0][1][:2]
        self._nElectrons      = int(self.configurations[0][1][2:])

        self._resonance = branch['axes'][0][4]

        if self.experiment == 'RIXS':
            self._xLabel = branch['axes'][1][0]
        else:
            self._xLabel = branch['axes'][0][0]
        self._yLabel = 'Intensity'
        

        # Calculation parameters
        self._nPsis           = branch['number of states']
        self._nPsisMax        = self.nPsis
        self.nPsisAuto        = nPsisAuto
        if self.nPsisAuto == 0:
            self.nPsis = nPsis
        self.polarization        = polarization       
        self._nConfigurationsMax = 1
        self.nConfigurations     = 1  # this is always 1 in Crispy if LMCT or MLCT is False
        self._hamiltonianTerms = branch['hamiltonian terms']
        
        # axis to calculate expected values
        self.k = k

        # geometry
        assert isinstance(R, Iterable), 'R must be a list of type [[R1, angle1], [R2, angle2]]'
        if len(R) > 0:
            for _ in R:
                assert isinstance(_, Iterable), 'R must be a list of type [[R1, angle1], [R2, angle2]]'
                assert isinstance(_[0], str), 'Rotations must be string'
                assert _[0] in ['x', 'y', 'z'], 'Rotations must be `x`, `y`, or `z`'
        assert tth >= 0 and tth <= 180, 'tth must be a number (in degrees) between 0 and 180'
        self._R   = R
        self._tth = tth

        self._LVin = (0, 0, 1)
        self._LHin = (0, 1, 0)
        self._kin  = (1, 0, 0)
        self._update_geometry()

        # spectrum attributes
        if self.experiment == 'RIXS':
            self.E = E
            self.xMin     = xMin
            self.xMax     = xMax
            self.xNPoints = xNPoints
            self.gamma1   = gamma1
            self.gamma2   = gamma2

            self._is_energy_map = False  # this is for internal use only (it is used at q.update_lua_script() and q.run_rixs_energy_map()
        else:
            self._E = None
            self.xMin     = xMin
            self.xMax     = xMax
            self.xNPoints = xNPoints
            self.gamma1   = gamma1
            self._gamma2  = None

        # set basename (for filepath_lua and filepath_spec)
        # if self.experiment in ['RIXS', 'XES']:
        #     _shortEdge = self.edge[-5:-1]
        # else:
        #     _shortEdge = self.edge[-3:-1]
        # self._baseName = '{}{}_{}_{}_{}'.format(
        #     self.element, self.charge, self.symmetry, _shortEdge,
        #     self.experiment)

        # get template filename
        self._templateName = branch['template name']
        self._templatePath = settings.TEMPLATES_FOLDERPATH/self.templateName
        with open(self.templatePath) as p:
            self.template = p.read()

        # other attributes
        self.verbosity   = '0x0000'  # in the future, verbosity could be a parameters
        self.denseBorder = str(denseBorder)  # Default 2000
        self._lua_script = None  # stores lua script

        # hamiltonianData and hamiltonianState
        self._hamiltonianData = odict()
        self._hamiltonianState = odict()
        self._fixedTermsParameters = odict()

        branch = settings.PARAMETERS['elements'][self.element]['charges'][self.charge] # noqa
        for label, configuration in self.configurations:
            label = '{} Hamiltonian'.format(label)
            terms = branch['configurations'][configuration]['terms']

            for term in self.hamiltonianTerms:
                if term in ('Atomic', 'Magnetic Field', 'Exchange Field'):
                    node = terms[term]
                else:
                    node = terms[term]['symmetries'][self.symmetry]

                parameters = node['parameters']['variable']
                for parameter in parameters:
                    if 'Atomic' in term or 'Hybridization' in term:
                        if parameter[0] in ('F', 'G'):
                            scaleFactor = 0.8
                            data = [parameters[parameter], scaleFactor]
                        elif parameter[0] == 'ζ':
                            scaleFactor = 1.0
                            data = [parameters[parameter], scaleFactor]
                        else:
                            data = parameters[parameter]
                    else:
                        data = parameters[parameter]

                    self.hamiltonianData[term][label][parameter] = data

                parameters = terms[term]['parameters']['fixed']
                for parameter in parameters:
                    value = parameters[parameter]
                    self._fixedTermsParameters[term][parameter] = value

                if term in ('Atomic', 'Crystal Field', 'Magnetic Field'):
                    self.hamiltonianState[term] = 2
                else:
                    self.hamiltonianState[term] = 0
        self._hamiltonianState = _hamiltonianState(self.hamiltonianState, parent=self)
        self._hamiltonianData  = _hamiltonianData(self.hamiltonianData)

        # experiment attributes
        self.temperature               = temperature
        self._magneticField            = magneticField 
        self._magneticFieldOrientation = magneticFieldOrientation 
        self.magneticField             = self.magneticField 
        self.magneticFieldOrientation  = self.magneticFieldOrientation 
        
        return
        
    def _set_primary_attributes(self, element, charge, symmetry, experiment, edge):
        """Validate and set primary attributes."""

        # validation
        if element is None:
            element = 'Ni'
        assert element in settings.ELEMENTS, f'Not a valid element.\nValid elements are: {settings.ELEMENTS}'

        valid_charges = tuple(settings.PARAMETERS['elements'][element]['charges'].keys())
        if charge is None:
            charge = valid_charges[0]
        assert charge in valid_charges, f'Not a valid oxidation state for {element}\nValid oxidation states are: {valid_charges}'

        valid_symmetries = tuple(settings.PARAMETERS['elements'][element]['charges'][charge]['symmetries'].keys())
        if symmetry is None:
            symmetry = valid_symmetries[0]
        assert symmetry in valid_symmetries, f'Not a valid symmetry for {element} {charge}\nValid symmetries are: {valid_symmetries}'

        valid_experiments = tuple(settings.PARAMETERS['elements'][element]['charges'][charge]['symmetries'][symmetry]['experiments'].keys())
        if experiment is None:
            experiment = valid_experiments[0]
        assert experiment in valid_experiments, f'Not a valid experiment for {element} {charge} {symmetry}\nValid experiments are: {valid_experiments}'

        valid_edges = tuple(settings.PARAMETERS['elements'][element]['charges'][charge]['symmetries'][symmetry]['experiments'][experiment]['edges'].keys())
        if edge is None:
            edge = valid_edges[0]
        assert edge in valid_edges, f'Not a valid edge for {element} {charge} {symmetry} {experiment}\nValid edges are: {valid_edges}'

        # initialize primary attributes
        self._element    = element
        self._charge     = charge
        self._symmetry   = symmetry
        self._experiment = experiment
        self._edge       = edge
        
        return

    # PRIMARY ATTRIBUTES
    @property
    def element(self):
        return self._element
    @element.setter
    def element(self, value):
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge). Please, start a new Calculation() object')
    @element.deleter
    def element(self):
        raise AttributeError('Cannot delete object.')

    @property
    def charge(self):
        return self._charge
    @charge.setter
    def charge(self, value):
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge). Please, start a new Calculation() object')
    @charge.deleter
    def charge(self):
        raise AttributeError('Cannot delete object.')

    @property
    def symmetry(self):
        return self._symmetry
    @symmetry.setter
    def symmetry(self, value):
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge). Please, start a new Calculation() object')
    @symmetry.deleter
    def symmetry(self):
        raise AttributeError('Cannot delete object.')

    @property
    def experiment(self):
        return self._experiment
    @experiment.setter
    def experiment(self, value):
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge). Please, start a new Calculation() object')
    @experiment.deleter
    def experiment(self):
        raise AttributeError('Cannot delete object.')

    @property
    def edge(self):
        return self._edge
    @edge.setter
    def edge(self, value):
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge). Please, start a new Calculation() object')
    @edge.deleter
    def edge(self):
        raise AttributeError('Cannot delete object.')
    

    # ATTRIBUTES SET AUTOMATICALY WITH PRIMARY ATTRIBUTES
    @property
    def xLabel(self):
        return self._xLabel
    @xLabel.setter
    def xLabel(self, value):
        raise AttributeError('`xLabel` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed. Please, start a new Calculation() object')
    @xLabel.deleter
    def xLabel(self):
        raise AttributeError('Cannot delete object.')

    @property
    def yLabel(self):
        return self._yLabel
    @yLabel.setter
    def yLabel(self, value):
        raise AttributeError('`yLabel` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed. Please, start a new Calculation() object')
    @yLabel.deleter
    def yLabel(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nElectrons(self):
        return self._nElectrons
    @nElectrons.setter
    def nElectrons(self, value):
        raise AttributeError('`nElectrons` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed. Please, start a new Calculation() object')
    @nElectrons.deleter
    def nElectrons(self):
        raise AttributeError('Cannot delete object.')

    @property
    def configurations(self):
        return self._configurations
    @configurations.setter
    def configurations(self, value):
        raise AttributeError('`configurations` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed. Please, start a new Calculation() object')
    @configurations.deleter
    def configurations(self):
        raise AttributeError('Cannot delete object.')

    @property
    def block(self):
        return self._block
    @block.setter
    def block(self, value):
        raise AttributeError('`block` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed. Please, start a new Calculation() object')
    @block.deleter
    def block(self):
        raise AttributeError('Cannot delete object.')

    @property
    def hamiltonianTerms(self):
        return self._hamiltonianTerms
    @hamiltonianTerms.setter
    def hamiltonianTerms(self, value):
        raise AttributeError('`hamiltonianTerms` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed. Please, start a new Calculation() object')
    @hamiltonianTerms.deleter
    def hamiltonianTerms(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def hamiltonianState(self):
        return self._hamiltonianState
    @hamiltonianState.setter
    def hamiltonianState(self, value):
        raise AttributeError('Cannot replace `hamiltonianState`. Please, edit one q.hamiltonianState[term] at a time')
    @hamiltonianState.deleter
    def hamiltonianState(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def hamiltonianData(self):
        return self._hamiltonianData
    @hamiltonianData.setter
    def hamiltonianData(self, value):
        raise AttributeError('Cannot replace `hamiltonianData`. Please, edit one q.hamiltonianData[term] at a time')
    @hamiltonianData.deleter
    def hamiltonianData(self):
        raise AttributeError('Cannot delete object.')


    # SPECTRA ATTRIBUTES
    @property
    def E(self):
        return self._E
    @E.setter
    def E(self, value):
        if self.experiment == 'RIXS':
            if value is None:
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][4]
            assert value > 0, 'E cannot be negative'
            self._E = value
        else:
            raise AttributeError(f'`E` (incidente photon energy) is only used for calculation of RIXS spectra and is not used for experiment of type: {self.experiment}') 
    @E.deleter
    def E(self):
        raise AttributeError('Cannot delete object.')


    @property
    def xMin(self):
        return self._xMin
    @xMin.setter
    def xMin(self, value):
        if value is None:
            if self.experiment == 'RIXS':
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][1][1]
            else:
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][1]
        if self.experiment == 'RIXS':
            pass
        else:
            assert value > 0, 'xMin cannot be negative'
        self._xMin = value
    @xMin.deleter
    def xMin(self):
        raise AttributeError('Cannot delete object.')

    @property
    def xMax(self):
        return self._xMax
    @xMax.setter
    def xMax(self, value):
        if value is None:
            if self.experiment == 'RIXS':
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][1][2]
            else:
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][2]
        if self.experiment == 'RIXS':
            pass
        else:
            assert value > 0, 'xMax cannot be negative'
        self._xMax = value
    @xMax.deleter
    def xMax(self):
        raise AttributeError('Cannot delete object.')

    @property
    def xNPoints(self):
        return self._xNPoints
    @xNPoints.setter
    def xNPoints(self, value):
        if value is None:
            if self.experiment == 'RIXS':
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][1][3]
            else:
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][3]
        assert value >= 2, 'xNPoints cannot be less than 2.\nThe CreateResonantSpectra() function from Quanty prevents xNPoints to be less than 2.'
        self._xNPoints = int(value)
    @xNPoints.deleter
    def xNPoints(self):
        raise AttributeError('Cannot delete object.')

    @property
    def resonance(self):
        return self._resonance
    @resonance.setter
    def resonance(self, value):
        raise AttributeError('`resonance` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed. Please, start a new Calculation() object')
        # if value is None:
        #     value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][4]
        # assert value > 0, 'xEdge cannot be negative'
        # self._xEdge = value
    @resonance.deleter
    def resonance(self):
        raise AttributeError('Cannot delete object.')


    @property
    def k(self):
        return self._k
    @k.setter
    def k(self, value):
        assert len(value) == 3, 'k must be a vector, like [0, 1, 0].'
        self._k = br.normalize_vector(value)
    @k.deleter
    def k(self):
        raise AttributeError('Cannot delete object.')


    # GEOMETRY ATTRIBUTES
    @property
    def tth(self):
        return self._tth
    @tth.setter
    def tth(self, value):
        assert value >= 0 and value <= 180, 'tth must be a number (in degrees) between 0 and 180'
        self._tth = value
        self._update_geometry()
    @tth.deleter
    def tth(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def R(self):
        return self._R
    @R.setter
    def R(self, value):
        assert isinstance(value, Iterable), 'R must be a list of type [[R1, angle1], [R2, angle2]]'
        if len(value) > 0:
            for _ in value:
                assert isinstance(_, Iterable), 'R must be a list of type [[R1, angle1], [R2, angle2]]'
                assert isinstance(_[0], str), 'Rotations must be string'
                assert _[0] in ['x', 'y', 'z'], 'Rotations must be `x`, `y`, or `z`'
        self._R = value
        self._update_geometry()
    @R.deleter
    def R(self):
        raise AttributeError('Cannot delete object.')

    @property
    def cf(self):
        return self._cf
    @cf.setter
    def cf(self, value):
        raise AttributeError('`cf` cannot be modified. It is defined via automatically and depends on q.R attribute')
    @cf.deleter
    def cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def kin(self):
        return self._kin
    @kin.setter
    def kin(self, value):
        assert len(value) == 3, 'kin must be a vector, like [0, 1, 0].'
        self._kin = br.normalize_vector(value)
        self._update_geometry()
    @kin.deleter
    def kin(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LVin(self):
        return self._LVin
    @LVin.setter
    def LVin(self, value):
        assert len(value) == 3, 'LVin must be a vector, like [0, 1, 0].'
        self._LVin = br.normalize_vector(value)
        self._update_geometry()
    @LVin.deleter
    def LVin(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LHin(self):
        return self._LHin
    @LHin.setter
    def LHin(self, value):
        assert len(value) == 3, 'LHin must be a vector, like [0, 1, 0].'
        self._LHin = br.normalize_vector(value)
        self._update_geometry()
    @LHin.deleter
    def LHin(self):
        raise AttributeError('Cannot delete object.')

    @property
    def kout(self):
        return self._kout
    @kout.setter
    def kout(self, value):
        raise AttributeError('`kout` cannot be modified. It is defined via automatically and depends on q.kin and q.tth')
        # assert len(value) == 3, 'kout must be a vector, like [0, 1, 0].'
        # self._kout = br.normalize_vector(value)
        # self._update_geometry()
    @kout.deleter
    def kout(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LVout(self):
        return self._LVout
    @LVout.setter
    def LVout(self, value):
        raise AttributeError('`LVout` cannot be modified. It is defined via automatically and depends on q.LVin and q.tth')
        # assert len(value) == 3, 'LVout must be a vector, like [0, 1, 0].'
        # self._LVout = br.normalize_vector(value)
        # self._update_geometry()
    @LVout.deleter
    def LVout(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LHout(self):
        return self._LHout
    @LHout.setter
    def LHout(self, value):
        raise AttributeError('`LHout` cannot be modified. It is defined via automatically and depends on q.LHin and q.tth')
        # assert len(value) == 3, 'LHout must be a vector, like [0, 1, 0].'
        # self._LHout = br.normalize_vector(value)
        # self._update_geometry()
    @LHout.deleter
    def LHout(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def kin_cf(self):
        return self._kin_cf
    @kin_cf.setter
    def kin_cf(self, value):
        raise AttributeError('`kin` cannot be modified. It is defined via automatically and depends on q.LVin and q.R')
        # assert len(value) == 3, 'k1 must be a vector, like [0, 1, 0].'
        # self._kin_cf = br.normalize_vector(value)
        # self._update_geometry()
    @kin_cf.deleter
    def kin_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LVin_cf(self):
        return self._LVin_cf
    @LVin_cf.setter
    def LVin_cf(self, value):
        raise AttributeError('`LVin_cf` cannot be modified. It is defined via automatically and depends on q.LVin and q.R')
        # assert len(value) == 3, 'LVin_cf must be a vector, like [0, 1, 0].'
        # self._LVin_cf = br.normalize_vector(value)
        # self._update_geometry()
    @LVin_cf.deleter
    def LVin_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LHin_cf(self):
        return self._LHin_cf
    @LHin_cf.setter
    def LHin_cf(self, value):
        raise AttributeError('`LHin_cf` cannot be modified. It is defined via automatically and depends on q.LHin and q.R')
        # assert len(value) == 3, 'LHin_cf must be a vector, like [0, 1, 0].'
        # self._LHin_cf = br.normalize_vector(value)
        # self._update_geometry()
    @LHin_cf.deleter
    def LHin_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def kout_cf(self):
        return self._kout_cf
    @kout_cf.setter
    def kout_cf(self, value):
        raise AttributeError('`kout_cf` cannot be modified. It is defined via automatically and depends on q.kout and q.R')
        # assert len(value) == 3, 'kout_cf must be a vector, like [0, 1, 0].'
        # self._kout_cf = br.normalize_vector(value)
        # self._update_geometry()
    @kout_cf.deleter
    def kout_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LVout_cf(self):
        return self._LVout_cf
    @LVout_cf.setter
    def LVout_cf(self, value):
        raise AttributeError('`LVout_cf` cannot be modified. It is defined via automatically and depends on q.LVout and q.R')
        # assert len(value) == 3, 'LVout_cf must be a vector, like [0, 1, 0].'
        # self._LVout_cf = br.normalize_vector(value)
        # self._update_geometry()
    @LVout_cf.deleter
    def LVout_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LHout_cf(self):
        return self._LHout_cf
    @LHout_cf.setter
    def LHout_cf(self, value):
        raise AttributeError('`LHout_cf` cannot be modified. It is defined via automatically and depends on q.LHout and q.R')
        # assert len(value) == 3, 'LHout_cf must be a vector, like [0, 1, 0].'
        # self._LHout_cf = br.normalize_vector(value)
        # self._update_geometry()
    @LHout_cf.deleter
    def LHout_cf(self):
        raise AttributeError('Cannot delete object.')

    # EXPERIMENT ATTRIBUTES
    @property
    def temperature(self):
        return self._temperature
    @temperature.setter
    def temperature(self, value):
        assert value >= 0, 'Temperature cannot be negative.'
        if value == 0:
            print('Temperature = 0\nnPsi set to 1.')
            self.nPsis = 1
            self.nPsisAuto = 0
        self._temperature = value
    @temperature.deleter
    def temperature(self):
        raise AttributeError('Cannot delete object.')

    @property
    def magneticField(self):
        return self._magneticField
    @magneticField.setter
    def magneticField(self, value):
        # print(value)
        # small = np.finfo(np.float32).eps  # ~1.19e-7 
        if value is None:
            raise ValueError('Magnetic field value must be number bigger than 0.002 T. To turn off magnetic field hamiltonian, please use hamiltonianState["Magnetic Field"] = False')     
        elif value < 0:
            raise ValueError('Magnetic field value must be positive.')
        elif value < 0.002:
            raise ValueError('Magnetic field cannot be smaller than 0.002 T. Please, turn off the magnetic field hamiltonian using hamiltonianState["Magnetic Field"] = False')     
        self._magneticField = value
        self._update_magnetic_field_hamiltonian_data()      
    @magneticField.deleter
    def magneticField(self):
        raise AttributeError('Cannot delete object.')

    @property
    def magneticFieldOrientation(self):
        return self._magneticFieldOrientation
    @magneticFieldOrientation.setter
    def magneticFieldOrientation(self, value):
        assert len(value) == 3, 'magneticFieldOrientation must be a vector, like [0, 1, 0].'
        self._magneticFieldOrientation = br.normalize_vector(value)
        self._update_magnetic_field_hamiltonian_data()      
    @magneticFieldOrientation.deleter
    def magneticFieldOrientation(self):
        raise AttributeError('Cannot delete object.')

    # CALCULATION ATTRIBUTES
    @property
    def polarization(self):
        return self._polarization
    @polarization.setter
    def polarization(self, value):
        error = "Invalid value for `polarization`\nValid values are:\nXAS: 'isotropic', 'linear', and 'circular'\nXES: 'isotropic'\nXPS: 'isotropic'\nRIXS: 'isotropic', 'linear' (not for 1s), and 'circular' (not for 1s)"
        assert type(value) == str, error
        if self.experiment == 'RIXS' or self.experiment == 'XAS':
            assert value.lower() in ['isotropic', 'circular', 'linear'], error
            if self.experiment == 'RIXS' and self.edge.startswith('K'):
                assert value.lower() in ['isotropic', ], error
        elif self.experiment == 'XPS' or self.experiment == 'XES':
            assert value.lower() in ['isotropic', ], error
        self._polarization = value
    @polarization.deleter
    def polarization(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nPsis(self):
        return self._nPsis
    @nPsis.setter
    def nPsis(self, value):
        if value is None:
            self.nPsisAuto = 1
            value = self._nPsisMax
        else:
            assert value > 0, 'The number of states must be larger than zero.'
            assert value <= self._nPsisMax, f'The selected number of states exceeds the maximum.\nMaximum {self._nPsisMax}'
        self._nPsis     = value
        # self._nPsisAuto = 0
    @nPsis.deleter
    def nPsis(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nPsisMax(self):
        return self._nPsisMax
    @nPsisMax.setter
    def nPsisMax(self, value):
        raise AttributeError('Cannot manually edit nPsisMax object.')
    @nPsisMax.deleter
    def nPsisMax(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nPsisAuto(self):
        return self._nPsisAuto
    @nPsisAuto.setter
    def nPsisAuto(self, value):
        if value == 0:
            self._nPsisAuto = 0
        elif value == 1:
            self._nPsisAuto = 1
        else:
            raise ValueError('nPsiAuto can only be 0 or 1')
    @nPsisAuto.deleter
    def nPsisAuto(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nConfigurationsMax(self):
        return self._nConfigurationsMax
    @nConfigurationsMax.setter
    def nConfigurationsMax(self, value):
        raise AttributeError('Cannot manually edit nConfigurationsMax object.')
    @nConfigurationsMax.deleter
    def nConfigurationsMax(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nConfigurations(self):
        return self._nConfigurations
    @nConfigurations.setter
    def nConfigurations(self, value):
        assert value <= self.nConfigurationsMax, f'The maximum number of configurations is {self.nConfigurationsMax}'
        assert value >= 1, f'The minimum number of configurations is 1'
        assert isinstance(value, int), f'nConfigurations must be a int, not type: {type(value)}'
        self._nConfigurations = value
    @nConfigurations.deleter
    def nConfigurations(self):
        raise AttributeError('Cannot delete object.')
    

    # OTHER
    @property
    def verbosity(self):
        return self._verbosity
    @verbosity.setter
    def verbosity(self, value):
        # -- For maximal verbosity
        # Verbosity(0xFFFF)
        # -- For minimal verbosity
        # Verbosity(0x0)
        # -- For standard verbosity
        # Verbosity(0xF)
        self._verbosity = value
    @verbosity.deleter
    def verbosity(self):
        raise AttributeError('Cannot delete object.')

    @property
    def denseBorder(self):
        return self._denseBorder
    @denseBorder.setter
    def denseBorder(self, value):
        self._denseBorder = value
    @denseBorder.deleter
    def denseBorder(self):
        raise AttributeError('Cannot delete object.')

    @property
    def lua_script(self):
        return self._lua_script
    @lua_script.setter
    def lua_script(self, value):
        self._lua_script = value
    @lua_script.deleter
    def lua_script(self):
        raise AttributeError('Cannot delete object.')

    @property
    def templatePath(self):
        return self._templatePath
    @templatePath.setter
    def templatePath(self, value):
        raise AttributeError('`templatePath` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed. Please, start a new Calculation() object')
    @templatePath.deleter
    def templatePath(self):
        raise AttributeError('Cannot delete object.')

    @property
    def templateName(self):
        return self._templateName
    @templateName.setter
    def templateName(self, value):
        raise AttributeError('`templateName` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed. Please, start a new Calculation() object')
    @templateName.deleter
    def templateName(self):
        raise AttributeError('Cannot delete object.')

    @property
    def template(self):
        return self._template
    @template.setter
    def template(self, value):
        self._template = value
    @template.deleter
    def template(self):
        raise AttributeError('Cannot delete object.')

    # BROADENING
    @property
    def gamma1(self):
        return self._gamma1
    @gamma1.setter
    def gamma1(self, value):
        if value is None:
            value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][5][0]
        assert value > 0, 'gamma1 cannot be negative'
        self._gamma1 = value
    @gamma1.deleter
    def gamma1(self):
        raise AttributeError('Cannot delete object.')

    @property
    def gamma2(self):
        return self._gamma2
    @gamma2.setter
    def gamma2(self, value):
        if self.experiment == 'RIXS':
            if value is None:
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][1][5][0]       
            assert value > 0, 'gamma2 cannot be negative'
            self._gamma2 = value
        else:
            raise AttributeError(f'`gamma2` is only used for calculation of RIXS spectra and is not used for experiment of type: {self.experiment}') 
    @gamma2.deleter
    def gamma2(self):
        raise AttributeError('Cannot delete object.')


    # SUPPORT METHODS
    def _update_geometry(self):
        """Converts lab coordinates to CF coordinates."""
        # geometry (lab coordinates)
        # lab coordinates: z points up and the beam points in x
        # LH, LV, kin, 
        # LHout, LVout, kout
        # LVin = (0, 0, 1)
        # LHin = (0, 1, 0)
        # kin  = (1, 0, 0)
        LVin = self.LVin
        LHin = self.LHin
        kin  = self.kin

        tth  = self.tth
        R    = self.R
        assert br.is_perpendicular(LVin, LHin) == True, 'LV and LH must be perpendicular'
        assert br.is_perpendicular(LVin, kin) == True, 'LV and Kin must be perpendicular'
        assert br.is_perpendicular(LHin, kin) == True, 'LV and Kin must be perpendicular'
        kout  = br.rotate_vector('z', tth, kin)  # out-going beam (lab coordinates)
        LHout = br.rotate_vector('z', tth, LHin)
        LVout = br.rotate_vector('z', tth, LVin)

        # defining CF (prime) coordinate system
        x_prime, y_prime, z_prime = (1, 0, 0), (0, 1, 0), (0, 0, 1)
        if len(R) > 0:
            for _rotation in R:
                _R = _rotation[0]
                _angle = _rotation[1]
                x_prime, y_prime, z_prime = br.rotate_system(_R, _angle, x_prime, y_prime, z_prime)
        # Lab coordinate system to Octahedra coordinate system
        # incoming beam
        kin_cf = br.change2prime(x_prime, y_prime, z_prime, kin)
        # incoming polarization
        LHin_cf = br.change2prime(x_prime, y_prime, z_prime, LHin)
        LVin_cf = br.change2prime(x_prime, y_prime, z_prime, LVin)
        # out-going beam
        kout_cf = br.change2prime(x_prime, y_prime, z_prime, kout)
        # out-going polarization
        LHout_cf = br.change2prime(x_prime, y_prime, z_prime, LHout)
        LVout_cf = br.change2prime(x_prime, y_prime, z_prime, LVout)
    
        self._kin   = [float(_) for _ in kin]
        self._LHin  = [float(_) for _ in LHin]
        self._LVin  = [float(_) for _ in LVin]
        self._kout  = [float(_) for _ in kout]
        self._LHout = [float(_) for _ in LHout]
        self._LVout = [float(_) for _ in LVout]

        self._kin_cf   = [float(_) for _ in kin_cf]
        self._LVin_cf  = [float(_) for _ in LVin_cf]
        self._LHin_cf  = [float(_) for _ in LHin_cf]
        self._kout_cf  = [float(_) for _ in kout_cf]
        self._LVout_cf = [float(_) for _ in LVout_cf]
        self._LHout_cf = [float(_) for _ in LHout_cf]
        self._cf    = [x_prime, y_prime, z_prime]
        return
    
    def _update_magnetic_field_hamiltonian_data(self):
        """Updates the hamiltanian data given Mag. field and Mag. field orientation"""
        value = self.magneticField
        k1    = np.array(self.magneticFieldOrientation)

        TESLA_TO_EV = 5.788e-05  # muB Bohr magneton in units of electron Volt 
        value = value * TESLA_TO_EV
        
        # updating hamiltonian data
        configurations = self.hamiltonianData['Magnetic Field']
        for configuration in configurations:
            parameters = self.hamiltonianData['Magnetic Field'][configuration]
            for i, parameter in enumerate(parameters):
                value2 = float(value * np.abs(k1[i]))
                # value2 = float(value * k1[i] * TESLA_TO_EV)

                # fix negative zero problem (same as in crispy)
                if abs(value2) == 0.0:
                    value2 = 0.0
                
                self.hamiltonianData['Magnetic Field'][configuration][parameter] = value2
        return 
    
    def _termSuffix(self, term):
        """sets up the i, m, f suffix (H_i, H_m, H_f, ...)"""
        term = term.lower()
        replacements = [(' ', '_'), ('-', '_'), ('(', ''), (')', '')]
        for replacement in replacements:
            term = term.replace(*replacement)
        return term

    # GET ATTRS AND METHODS
    def get_attrs(self):
        """returns attrs""" 
        # return [key for key in self.__dir__.keys() if key.startswith('_') == False] + [key[1:] for key in self.__dict__.keys() if key.startswith('_') == True] 
        return [key for key in self.__dir__() if key.startswith('_') == False and callable(getattr(self, key))==False]

    def _get_methods(self):
        """return a list of methods available (including hidden ones)"""

        methodList = []
        for method_name in dir(self):
            if method_name in self.get_attrs():
                pass
            else:
                try:
                    if callable(getattr(self, method_name)):
                        methodList.append(str(method_name))
                except Exception:
                    methodList.append(str(method_name))
        return methodList
    
    def get_methods(self):
        """returns a list of methods available"""
        return [_ for _ in self._get_methods() 
                if _.startswith('_') == False]
    
    # PLOT
    def plot_geometry(self, show_magnetic_field=True, show_legend=True):
        """Plot a 3D sketch of the experimental geometry

        Args:
            show_magnetic_field (bool, optional): If True, shows arrow indicating
                direction of the magnetic field. Default is True.
            show_legend (bool, optional): If True, shows legend. Default is True.

        Returns:
            None
        """
        fig = br.figure()
        ax1 = fig.add_subplot(121, projection='3d', proj_type = 'ortho')
        ax2 = fig.add_subplot(122, projection='3d', proj_type = 'ortho')


        # lab
        k1_l = [-4, 0, 0] + list(np.array(self.kin)*4)
        _ = ax1.quiver(*k1_l, color='black', arrow_length_ratio=0.1, lw=4)

        LH1_l = [-3, 0, 0] + list(self.LHin)
        _ = ax1.quiver(*LH1_l, color='black', arrow_length_ratio=0.1, lw=1)

        LV1_l = [-3, 0, 0] + list(self.LVin)
        _ = ax1.quiver(*LV1_l, color='red', arrow_length_ratio=0.1, lw=1)

        k2_l = [0, 0, 0] + list(np.array(self.kout)*4)
        a = ax1.quiver(*k2_l, color='green', arrow_length_ratio=0.1, lw=4)

        LH2_l = [np.array(self.kout)[0]*2, np.array(self.kout)[1]*2, 0] + list(self.LHout)
        _ = ax1.quiver(*LH2_l, color='black', arrow_length_ratio=0.1, lw=1)

        LV2_l = [np.array(self.kout)[0]*2, np.array(self.kout)[1]*2, 0] + list(self.LVout)
        _ = ax1.quiver(*LV2_l, color='red', arrow_length_ratio=0.1, lw=1)

        br.draw_octahedron(ax1, x=self.cf[0], y=self.cf[1], z=np.array(self.cf[2])*2)


        # CF
        k1_l = [0, 0, 0] + list(-np.array(self.kin_cf)*4)
        _ = ax2.quiver(*k1_l, color='black', arrow_length_ratio=0, lw=4, label='k$_{in}$')

        LH1_l = [0, 0, 0] + list(self.LHin_cf)
        LH1_l = [-np.array(self.kin_cf)[0]*2, -np.array(self.kin_cf)[1]*2, -np.array(self.kin_cf)[2]*2] + list(self.LHin_cf)
        _ = ax2.quiver(*LH1_l, color='black', arrow_length_ratio=0.1, lw=1)

        LV1_l = [-np.array(self.kin_cf)[0]*2, -np.array(self.kin_cf)[1]*2, -np.array(self.kin_cf)[2]*2] + list(self.LVin_cf)
        _ = ax2.quiver(*LV1_l, color='red', arrow_length_ratio=0.1, lw=1)

        k2_l = [0, 0, 0] + list(np.array(self.kout_cf)*4)
        a = ax2.quiver(*k2_l, color='green', arrow_length_ratio=0.1, lw=4, label='k$_{out}$')

        LH2_l = [np.array(self.kout_cf)[0]*2, np.array(self.kout_cf)[1]*2, np.array(self.kout_cf)[2]*2] + list(self.LHout_cf)
        _ = ax2.quiver(*LH2_l, color='black', arrow_length_ratio=0.1, lw=1, label='LH')

        LV2_l = [np.array(self.kout_cf)[0]*2, np.array(self.kout_cf)[1]*2, np.array(self.kout_cf)[2]*2] + list(self.LVout_cf)
        _ = ax2.quiver(*LV2_l, color='red', arrow_length_ratio=0.1, lw=1, label='LV')

        br.draw_octahedron(ax2, z=(0, 0, 2))

        # final
        ax1.set_title('lab coordinate system')
        ax2.set_title('CF coordinate system')
        for ax in (ax1, ax2):
            ax.set_xlim(-4, 4)
            ax.set_ylim(-4, 4)
            ax.set_zlim(-4, 4)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_zlabel('z')

        if show_magnetic_field:
            _x, _y, _z = self.magneticFieldOrientation
            magneticFieldOrientation = [-_x/2, -_y/2, -_z/2] + [2*_x/2, 2*_y/2, 2*_z/2]
            _ = ax2.quiver(*magneticFieldOrientation, color='magenta', arrow_length_ratio=0.3, lw=3, label='Mag. Field')
        
            R = self.R
            x_prime, y_prime, z_prime = (1, 0, 0), (0, 1, 0), (0, 0, 1)
            if len(R) > 0:
                for _rotation in R[::-1]:
                    _R = _rotation[0]
                    _angle = -_rotation[1]
                    x_prime, y_prime, z_prime = br.rotate_system(_R, _angle, x_prime, y_prime, z_prime)
        
            _x, _y, _z = br.change2prime(x_prime, y_prime, z_prime, [_x, _y, _z])
            magneticFieldOrientation = [-_x/2, -_y/2, -_z/2] + [2*_x/2, 2*_y/2, 2*_z/2]
            _ = ax1.quiver(*magneticFieldOrientation, color='magenta', arrow_length_ratio=0.3, lw=3)

        if show_legend:
            br.leg(ax=ax2)
        return

    # LUA SCRIPT
    def update_lua_script(self):
        """Creates lua script based on calculation parameters. Saves it to q.lua_script."""

        # initialization =====================================
        replacements = odict()

        # other ==============================================
        replacements['$DenseBorder'] = self.denseBorder
        replacements['$Verbosity']   = self.verbosity

        # calculation parameters ============================
        replacements['$NConfigurations'] = self.nConfigurations
        subshell = self.configurations[0][1][:2]
        subshell_occupation = self.configurations[0][1][2:]
        replacements['$NElectrons_{}'.format(subshell)] = subshell_occupation
        replacements['$NPsisAuto'] = self.nPsisAuto
        replacements['$NPsis'] = self.nPsis

        # temperature =======================================
        replacements['$T'] = self.temperature

        # Spectrum and gamma parameters =====================
        if self.experiment in ['RIXS', ]:
            if self._is_energy_map == False:
                replacements['$Emin1'] = self.E
                replacements['$Emax1'] = self.E + 10
                replacements['$NE1'] = 1
            replacements['$Eedge1'] = self.resonance
            replacements['$Gamma1'] = self.gamma1

            replacements['$Emin2']  = self.xMin
            replacements['$Emax2']  = self.xMax
            replacements['$NE2']    = self.xNPoints - 1
            replacements['$Gamma2'] = self.gamma2  # For RIXS: Broadening of the energy loss

        else:
            replacements['$Emin1']  = self.xMin
            replacements['$Emax1']  = self.xMax
            replacements['$NE1']    = self.xNPoints - 1
            replacements['$Eedge1'] = self.resonance
            replacements['$Gamma1'] = self.gamma1
            if self.experiment == 'XES':
                replacements['$Emin1'] = self.xMin + 20
                replacements['$Emax1'] = self.xMax + 20

            
        # remove artificial broadening =====================
        if self.experiment == 'XPS' or self.experiment == 'XAS':
            pattern = r"Gmin1 = $Gmin1 - Gamma"
            subst   = r"-- Gmin1 = $Gmin1 - Gamma"
            replacements[pattern] = subst

            pattern = r"Gmax1 = $Gmax1 - Gamma"
            subst   = r"-- Gmax1 = $Gmax1 - Gamma"
            replacements[pattern] = subst

            pattern = r"Egamma1 = ($Egamma1 - Eedge1) + DeltaE"
            subst   = r"-- Egamma1 = ($Egamma1 - Eedge1) + DeltaE"
            replacements[pattern] = subst

            pattern = r"G.Broaden(0, {{Emin, Gmin1}, {Egamma1, Gmin1}, {Egamma1, Gmax1}, {Emax, Gmax1}})"
            subst   = r"-- G.Broaden(0, {{Emin, Gmin1}, {Egamma1, Gmin1}, {Egamma1, Gmax1}, {Emax, Gmax1}})"
            replacements[pattern] = subst

        # geometry =========================================
        s = '{{{0:.8g}, {1:.8g}, {2:.8g}}}'

        k1 = np.array(self.kin_cf)
        k1 = k1 / np.linalg.norm(k1)
        replacements['$k1'] = s.format(k1[0], k1[1], k1[2])

        eps11 = np.array(self.LVin_cf)
        eps11 = eps11 / np.linalg.norm(eps11)
        replacements['$eps11'] = s.format(eps11[0], eps11[1], eps11[2])

        eps12 = np.array(self.LHin_cf)
        eps12 = eps12 / np.linalg.norm(eps12)
        replacements['$eps12'] = s.format(eps12[0], eps12[1], eps12[2])

        if self.experiment == 'RIXS':
            k2 = np.array(self.kout_cf)
            k2 = k2 / np.linalg.norm(k2)
            replacements['$k2'] = s.format(k2[0], k2[1], k2[2])

            eps21 = np.array(self.LVout_cf)
            eps21 = eps21 / np.linalg.norm(eps21)
            replacements['$eps21'] = s.format(eps21[0], eps21[1], eps21[2])

            eps22 = np.array(self.LHout_cf)
            eps22 = eps22 / np.linalg.norm(eps22)
            replacements['$eps22'] = s.format(eps22[0], eps22[1], eps22[2])

        k = np.array(self.k)
        k = k / np.linalg.norm(k)
        replacements['$k'] = s.format(k[0], k[1], k[2])

        # to calculate =====================================
        if self.polarization.lower() == 'isotropic':
            polarization = 'Isotropic'
        elif self.polarization.lower() == 'circular':
            polarization = 'Circular Dichroism'
        elif self.polarization.lower() == 'linear':
            polarization = 'Linear Dichroism'
        replacements['$spectra'] = f"'{polarization}'"
    
        # hamiltonian =======================================
        hamiltonianState = {key:(1 if value else 0) for key, value in self.hamiltonianState.items()}
        for term in self.hamiltonianData:
            configurations = self.hamiltonianData[term]
            for configuration, parameters in configurations.items():
                if 'Initial' in configuration:
                    suffix = 'i'
                elif 'Intermediate' in configuration:
                    suffix = 'm'
                elif 'Final' in configuration:
                    suffix = 'f'
                for parameter, data in parameters.items():
                    # Convert to parameters name from Greek letters.
                    parameter = parameter.replace('ζ', 'zeta')
                    parameter = parameter.replace('Δ', 'Delta')
                    parameter = parameter.replace('σ', 'sigma')
                    parameter = parameter.replace('τ', 'tau')
                    parameter = parameter.replace('μ', 'mu')
                    parameter = parameter.replace('ν', 'nu')

                    scaleFactor = None
                    try:
                        value, scaleFactor = data
                    except TypeError:
                        value = data

                    # if self.magneticField == 0:
                    #     small = np.finfo(np.float32).eps  # ~1.19e-7
                    #     if parameter == 'Bx':
                    #         value = k1[0] * small
                    #     elif parameter == 'By':
                    #         value = k1[1] * small
                    #     elif parameter == 'Bz':
                    #         value = k1[2] * small

                    key = '${}_{}_value'.format(parameter, suffix)
                    replacements[key] = '{}'.format(value)

                    if scaleFactor is not None:
                        key = '${}_{}_scale'.format(parameter, suffix)
                        replacements[key] = '{}'.format(scaleFactor)

            checkState = hamiltonianState[term]
            if checkState > 0:
                checkState = 1

            termSuffix = self._termSuffix(term)
            replacements['$H_{}'.format(termSuffix)] = checkState

            try:
                parameters = self._fixedTermsParameters[term]
            except KeyError:
                pass
            else:
                for parameter in parameters:
                    value = parameters[parameter]
                    replacements['${}'.format(parameter)] = value

        # experiment =======================================
        replacements['$Experiment'] = self.experiment

        # return spectra ====================================
        # this is a small fix so I don't have to change all templates
        # I changed for RIXS manually
        if self.experiment == 'XPS':
            _core = self.templateName.split('.')[0].split('_')[-1]
            pattern = r"SaveSpectrum(Giso / #T_" + _core + r", 'iso')"
            subst   = r"io.write('Here starts ISO spectrum:')" + '\n' +\
                      r"    print(-1 / math.pi * (Giso / #T_" + _core + "))" + '\n' +\
                      r"    io.write('Here ends ISO spectrum')" + '\n'
            replacements[pattern] = subst
        elif self.experiment == 'XES':
            pattern = r"Giso.Print({{'file', '$BaseName_iso.spec'}})"
            subst   = r"io.write('Here starts ISO spectrum:')" + '\n' +\
                      r"print(Giso)" + '\n' +\
                      r"io.write('Here ends ISO spectrum')" + '\n'
            replacements[pattern] = subst
        elif self.experiment == 'RIXS':
            pass
        elif self.experiment == 'XAS':
            pattern = r"SaveSpectrum(Giso, 'iso')"
            subst   = r"io.write('Here starts ISO spectrum:')" + '\n' +\
                      r"    print(-1 / math.pi * Giso)" + '\n' +\
                      r"    io.write('Here ends ISO spectrum')" + '\n'
            replacements[pattern] = subst

            pattern = r"SaveSpectrum(Gr, 'r')"
            subst   = r"io.write('Here starts CR spectrum:')" + '\n' +\
                      r"    print(-1 / math.pi * Gr)" + '\n' +\
                      r"    io.write('Here ends CR spectrum')" + '\n'
            replacements[pattern] = subst

            pattern = r"SaveSpectrum(Gl, 'l')"
            subst   = r"io.write('Here starts CL spectrum:')" + '\n' +\
                      r"    print(-1 / math.pi * Gl)" + '\n' +\
                      r"    io.write('Here ends CL spectrum')" + '\n'            
            replacements[pattern] = subst

            pattern = r"SaveSpectrum(Gr - Gl, 'cd')"
            subst   = r""           
            replacements[pattern] = subst

            pattern = r"SaveSpectrum(Gv, 'v')"
            subst   = r"io.write('Here starts LV spectrum:')" + '\n' +\
                      r"    print(-1 / math.pi * Gv)" + '\n' +\
                      r"    io.write('Here ends LV spectrum')" + '\n' 
            replacements[pattern] = subst

            pattern = r"SaveSpectrum(Gh, 'h')"
            subst   = r"io.write('Here starts LH spectrum:')" + '\n' +\
                      r"    print(-1 / math.pi * Gh)" + '\n' +\
                      r"    io.write('Here ends LH spectrum')" + '\n' 
            replacements[pattern] = subst

            pattern = r"SaveSpectrum(Gv - Gh, 'ld')"
            subst   = r""           
            replacements[pattern] = subst

        # replacement =====================================
        self.lua_script = copy.deepcopy(self.template)
        for replacement in replacements:
            self.lua_script = self.lua_script.replace(
                replacement, str(replacements[replacement]))
        return 
    
    def save_lua_script(self, filepath):
        """Save lua script to a file"""
        # self.update_lua_script()
        assert self.lua_script is not None, 'It seems like q.lua_script has not been created yet. Use q.update_lua_script()'
        with open(filepath, 'w') as f:
            f.write(self.lua_script)
        return

    def save_template(self, filepath):
        """Save lua template to a file"""
        # self.update_lua_script()
        with open(filepath, 'w') as f:
            f.write(self.template)
        return

    def replace(self, old, new):
        """replace strings on the template

        Args:
            old (str): text to be replaced. Must be contained in a single line.
            new (str): new text. Can be multiline.
        
        Returns:
            None
        """
        self.template = self.template.replace(old, new)
        return 
    
    # SAVE LOAD ATTRS
    def get_parameters(self):
        """Returns a dictionary with all calculation paramters."""
        p = dict(copy.deepcopy(self.hamiltonianData))
        for k1 in p:
            p[k1] = dict(p[k1])
        for k1 in p:
            for k2 in p[k1]:
                p[k1][k2] = dict(p[k1][k2])

        return dict(initial = dict(element    = self.element,
                                   charge     = self.charge,
                                   symmetry   = self.symmetry,
                                   experiment = self.experiment,
                                   edge       = self.edge),
                    geometry = dict(tth      = float(self.tth),
                                    R        = list(self.R), 
                                    kin      = list(self.kin),
                                    LHin     = list(self.LHin),
                                    LVin     = list(self.LVin),
                                    kout     = list(self.kout),
                                    LHout    = list(self.LHout),
                                    LVout    = list(self.LVout),
                                    kin_cf   = list(self.kin_cf),
                                    LVin_cf  = list(self.LVin_cf),
                                    LHin_cf  = list(self.LHin_cf),
                                    kout_cf  = list(self.kout_cf),
                                    LVout_cf = list(self.LVout_cf),
                                    LHout_cf = list(self.LHout_cf),
                                    cf       = [list(_) for _ in self.cf]),
                    experiment = dict(temperature      = float(self.temperature),
                                      magneticField    = float(self.magneticField),
                                      magneticFieldOrientation = list([float(_) for _ in self.magneticFieldOrientation])),
                    hamiltonian = dict(hamiltonianTerms = self.hamiltonianTerms,
                                       hamiltonianState = dict(self.hamiltonianState),
                                       hamiltonianData  = p),
                    calculation = dict(templateName   = self.templateName,
                                       templatePath   = str(self.templatePath),
                                       configurations = self.configurations,
                                       nConfigurationsMax = self.nConfigurationsMax,
                                       block          = self.block,
                                       nElectrons     = self.nElectrons,
                                       polarization   = self.polarization,
                                       nPsis          = self.nPsis,
                                       nPsisMax       = self.nPsisMax,
                                       nPsisAuto      = self.nPsisAuto,
                                       gamma1         = self.gamma1,
                                       xLabel         = self.xLabel,
                                       xMin           = self.xMin,
                                       xMax           = self.xMax,
                                       xNPoints       = self.xNPoints,
                                       resonance      = self.resonance,
                                       gamma2         = self.gamma2,
                                       yLabel         = self.yLabel,
                                       verbosity      = self.verbosity,
                                       denseBorder    = self.denseBorder)
            )

    def save_parameters(self, filepath):
        """Save calculation parameters to a text file
        
        This function is a bit finicky because it use json package and the 
        parameters have to be perfectly formated for it to work. """
        filepath = Path(filepath)
        pretty_print = True
        par = self.get_parameters()

        with open(str(filepath), 'w') as file:
            if pretty_print:
                file.write(json.dumps(par, indent=4, sort_keys=False))
            else:
                file.write(json.dumps(par))
        return

    # RUN
    def run(self, update=True):
        """Run Quanty.

        Args:
            update (bool, optional): if True, updates the lua script before 
                running calculation.

        Returns:
            XPS returns 1 spectrum and calculation output

            XES returns 1 spectrum and calculation output

            XAS returns 1 spectrum (iso) or 2 spectra and calculation output

            RIXS returns 1 spectra (iso) or one dictionary with 4 spectra
             and calculation output 

             
        Examples:

            XPS, isotropic
                >>> iso, out = q.run()
            
            XES, isotropic
                >>> iso, out = q.run()

            RIXS, isotropic
                >>> iso, out = q.run()
            RIXS, linear
                >>> {vv, vh, hv, hh, v, h}, out = q.run()
            RIXS, circular
                >>> {rv, rh, lv, lh, r, l}, out = q.run()
        
            XAS, isotropic
                >>> iso, out = q.run()
            XAS, linear
                >>> {LV, LH}, out = q.run()
            XAS, circular
                >>> {CR, CL}, out = q.run()
        """
        # update and run lua script
        if update:
            self.update_lua_script()
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            fp.write(self.lua_script.encode('utf-8'))
            fp.close()
            _out = quanty(fp.name)
        os.unlink(fp.name)
        
        # get spectra from output
        _x = np.linspace(self.xMin, self.xMax, self.xNPoints)  # _x = np.linspace(self.xMin, self.xMax, self.xNPoints + 1)
        if self.experiment == 'XPS':
            if self.polarization.lower() == 'isotropic':
                _out2 = _out.split('Here starts ISO spectrum:')[1].split('Here ends ISO spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                iso = br.Spectrum(_x, _y)
                out = _out.split('Here starts ISO spectrum:')[0]

                parameters = self.get_parameters()
                for _ss in (iso, ):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])

                return iso, out
        elif self.experiment == 'XES':
            if self.polarization.lower() == 'isotropic':
                _out2 = _out.split('Here starts ISO spectrum:')[1].split('Here ends ISO spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]][::-1]
                iso = br.Spectrum(_x, _y/np.abs(np.max(_y)))
                out = _out.split('Here starts ISO spectrum:')[0]

                parameters = self.get_parameters()
                for _ss in (iso, ):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])

                return iso, out
        elif self.experiment == 'RIXS':
            # When crispy plot the data, it considers xMax to be xMax-step
            # I still don't understand why this is like this
            # step = (self.E + 10 - self.E)/(2)
            # _x3 = np.linspace(self.E, self.E + 10-step, 2)
            # print(_x3)
            # _x2 = np.linspace(self.yMin, self.yMax, self.yNPoints)
            if self.polarization.lower() == 'isotropic':

                # =============== iso =============== #
                _out2 = _out.split('Here starts Giso spectrum:')[1].split('Here ends Giso spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                s = br.Spectrum(_x, _y)

                # calculation metadata
                parameters = self.get_parameters()
                for _name in parameters:
                    s.__setattr__(_name, parameters[_name])
                s.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(s.hamiltonian['hamiltonianData'])

                # output is everything before first spectrum
                out = _out.split('Here starts Giso spectrum:')[0]

                return s, out
            elif self.polarization.lower() == 'linear':

                # =============== vv =============== #
                _out2 = _out.split('Here starts Gvv spectrum:')[1].split('Here ends Gvv spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                s_vv = br.Spectrum(_x, _y)

                # =============== vh =============== #
                _out2 = _out.split('Here starts Gvh spectrum:')[1].split('Here ends Gvh spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                s_vh = br.Spectrum(_x, _y)

                # =============== hv =============== #
                _out2 = _out.split('Here starts Ghv spectrum:')[1].split('Here ends Ghv spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                s_hv = br.Spectrum(_x, _y)

                # =============== hh =============== #
                _out2 = _out.split('Here starts Ghh spectrum:')[1].split('Here ends Ghh spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                s_hh = br.Spectrum(_x, _y)

                # =============== v and h =============== #
                s_v = br.Spectra([s_vv, s_vh]).calculate_average()
                s_h = br.Spectra([s_hv, s_hh]).calculate_average()

                # calculation metadata
                parameters = self.get_parameters()
                for _s in (s_vv, s_vh, s_hv, s_hh, s_v, s_h):
                    for _name in parameters:
                        _s.__setattr__(_name, parameters[_name])
                    _s.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_s.hamiltonian['hamiltonianData'])
                    
                # output is everything before first spectrum
                out = _out.split('Here starts Gvv spectrum:')[0]

                return {'vv':s_vv, 'vh':s_vh, 'hv':s_hv, 'hh':s_hh, 'v':s_v, 'h':s_h}, out
            elif self.polarization.lower() == 'circular':

                # =============== rv =============== #
                _out2 = _out.split('Here starts Grv spectrum:')[1].split('Here ends Grv spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                s_rv = br.Spectrum(_x, _y)

                # =============== rh =============== #
                _out2 = _out.split('Here starts Grh spectrum:')[1].split('Here ends Grh spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                s_rh = br.Spectrum(_x, _y)

                # =============== lv =============== #
                _out2 = _out.split('Here starts Glv spectrum:')[1].split('Here ends Glv spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                s_lv = br.Spectrum(_x, _y)

                # =============== lh =============== #
                _out2 = _out.split('Here starts Glh spectrum:')[1].split('Here ends Glh spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                s_lh = br.Spectrum(_x, _y)

                # =============== v and h =============== #
                s_r = br.Spectra([s_rv, s_rh]).calculate_average()
                s_l = br.Spectra([s_lv, s_lh]).calculate_average()

                # calculation metadata
                parameters = self.get_parameters()
                for _s in (s_rv, s_rh, s_lv, s_lh, s_l, s_r):
                    for _name in parameters:
                        _s.__setattr__(_name, parameters[_name])
                    _s.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_s.hamiltonian['hamiltonianData'])
                
                # output is everything before first spectrum
                out = _out.split('Here starts Grv spectrum:')[0]

                return {'rv':s_rv, 'rh':s_rh, 'lv':s_lv, 'lh':s_lh, 'r':s_r, 'l':s_l}, out
        elif self.experiment == 'XAS':
            if self.polarization.lower() == 'isotropic':

                # =============== iso =============== #
                _out2 = _out.split('Here starts ISO spectrum:')[1].split('Here ends ISO spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0])+self.resonance for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                iso = br.Spectrum(_x, _y)
                
                # calculation metadata
                parameters = self.get_parameters()
                for _ss in (iso, ):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])
                
                # output is everything before first spectrum
                out = _out.split('Here starts ISO spectrum:')[0]

                return iso, out
            elif self.polarization.lower() == 'linear':

                # =============== v =============== #
                _out2 = _out.split('Here starts LV spectrum:')[1].split('Here ends LV spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0])+self.resonance for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                cr = br.Spectrum(_x, _y)

                # =============== h =============== #
                _out2 = _out.split('Here starts LH spectrum:')[1].split('Here ends LH spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0])+self.resonance for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                cl = br.Spectrum(_x, _y)

                # calculation metadata
                parameters = self.get_parameters()
                for _ss in (cr, cl):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])
                
                # output is everything before first spectrum
                out = _out.split('Here starts LV spectrum:')[0]

                return {'v':cr, 'h':cl}, out
            elif self.polarization.lower() == 'circular':
                # =============== r =============== #
                _out2 = _out.split('Here starts CR spectrum:')[1].split('Here ends CR spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0])+self.resonance for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                cr = br.Spectrum(_x, _y)

                # =============== l =============== #
                _out2 = _out.split('Here starts CL spectrum:')[1].split('Here ends CL spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0])+self.resonance for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                cl = br.Spectrum(_x, _y)

                # calculation metadata
                parameters = self.get_parameters()
                for _ss in (cr, cl):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])
                
                # output is everything before first spectrum
                out = _out.split('Here starts CR spectrum:')[0]

                return {'r':cr, 'l':cl}, out
        return 

    def run_rixs_energy_map(self, Emin, Emax, npoints, update=True, wsl=False):
        """same as q.run(), but returns rixs spectra for many incident energies

        Args:
            Emin, Emax (number): minimum and maximum photon incident energy.
            npoints (int): number of incident energies.
            update (bool, optional): if True, updates the lua script before 
                running calculation. Must be True. The current implementation 
                requires the lua script to be updated before runing quanty.

        Returns:
            RIXS returns br.Spectra and br.Image objects. If polarization is 
                linear or circular, returns dictionaries with spectra and images
                for each polarization. The integration of the rixs spectra is 
                also returned as a xas (pfy) reference. Calculation output is also returned.
             
        Examples:
            RIXS, isotropic
                >>> iso, im, xas, out = q.run()
            RIXS, linear
                >>> ss, ims, xas, out = q.run()
                >>> ss.keys()   # -> hh, vh, hv, hh, v, h
                >>> ims.keys()  # -> hh, vh, hv, hh, v, h
            RIXS, circular
                >>> ss, ims, xas, out = q.run()
                >>> ss.keys()   # -> rv, rh, lv, lh, r, l
                >>> ims.keys()  # -> rv, rh, lv, lh, r, l
        """
        assert self.experiment in ['RIXS', ], 'This only works for experiment=RIXS'
        assert update, 'the current implementation requires `update=True`. If you need `update=False`, this function must be revised'
        assert npoints > 1, f'invalid npoints ({npoints}). npoints must be higher than 1'

        # update lua script
        self._is_energy_map = True
        try:
            self.update_lua_script()
        except Exception as e:
            raise e
        finally:
            self._is_energy_map = False

        # adjust lua script for energy map
        replacements = odict()
        replacements['$Emin1'] = Emin
        replacements['$Emax1'] = Emax
        replacements['$NE1']   = npoints - 1

        for replacement in replacements:
            self.lua_script = self.lua_script.replace(
                replacement, str(replacements[replacement]))

        # save temporary file
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            fp.write(self.lua_script.encode('utf-8'))
            fp.close()
            _out = quanty(fp.name, wsl=wsl)
        os.unlink(fp.name)
        
        # get spectra from output
        _x = np.linspace(self.xMin, self.xMax, self.xNPoints)  # _x = np.linspace(self.xMin, self.xMax, self.xNPoints + 1)

        if self.polarization.lower() == 'isotropic':

            # =============== iso =============== #
            _out2 = _out.split('Here starts Giso spectrum:')[1].split('Here ends Giso spectrum')[0].split('\n')
            _out3 = [_ for _ in _out2 if _ != '']
            # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
            ss = br.Spectra()
            if wsl:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.strip().replace('  ', ' ').split(' ')[i]) for _ in _out3[5:]]
                    ss.append(br.Spectrum(_x, _y))
            else:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.split(' ')[i]) for _ in _out3[5:]]
                    ss.append(br.Spectrum(_x, _y))

            # calculation metadata
            parameters = self.get_parameters()
            for _ss in (ss, ):
                _ss.E = np.linspace(Emin, Emax, npoints)
                for _name in parameters:
                    _ss.__setattr__(_name, parameters[_name])
                _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])
            
            # create images
            im = ss.stack_spectra_as_columns()
            im.x_centers = ss.E

            # create xas
            xas = br.Spectrum(x=ss.E, y=ss.calculate_y_sum())

            # output is everything before first spectrum
            out = _out.split('Here starts Giso spectrum:')[0]

            return ss, im, xas, out
        elif self.polarization.lower() == 'linear':

            # =============== vv =============== #
            _out2 = _out.split('Here starts Gvv spectrum:')[1].split('Here ends Gvv spectrum')[0].split('\n')
            _out3 = [_ for _ in _out2 if _ != '']
            ss_vv = br.Spectra()
            if wsl:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.strip().replace('  ', ' ').split(' ')[i]) for _ in _out3[5:]]
                    ss_vv.append(br.Spectrum(_x, _y))
            else:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.split(' ')[i]) for _ in _out3[5:]]
                    ss_vv.append(br.Spectrum(_x, _y))

            # =============== vh =============== #
            _out2 = _out.split('Here starts Gvh spectrum:')[1].split('Here ends Gvh spectrum')[0].split('\n')
            ss_vh = br.Spectra()
            if wsl:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.strip().replace('  ', ' ').split(' ')[i]) for _ in _out3[5:]]
                    ss_vh.append(br.Spectrum(_x, _y))
            else:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.split(' ')[i]) for _ in _out3[5:]]
                    ss_vh.append(br.Spectrum(_x, _y))

            # =============== hv =============== #
            _out2 = _out.split('Here starts Ghv spectrum:')[1].split('Here ends Ghv spectrum')[0].split('\n')
            _out3 = [_ for _ in _out2 if _ != '']
            ss_hv = br.Spectra()
            if wsl:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.strip().replace('  ', ' ').split(' ')[i]) for _ in _out3[5:]]
                    ss_hv.append(br.Spectrum(_x, _y))
            else:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.split(' ')[i]) for _ in _out3[5:]]
                    ss_hv.append(br.Spectrum(_x, _y))

            # =============== hh =============== #
            _out2 = _out.split('Here starts Ghh spectrum:')[1].split('Here ends Ghh spectrum')[0].split('\n')
            _out3 = [_ for _ in _out2 if _ != '']
            ss_hh = br.Spectra()
            if wsl:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.strip().replace('  ', ' ').split(' ')[i]) for _ in _out3[5:]]
                    ss_hh.append(br.Spectrum(_x, _y))
            else:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.split(' ')[i]) for _ in _out3[5:]]
                    ss_hh.append(br.Spectrum(_x, _y))

            # =============== v and h =============== #
            ss_v = br.Spectra()
            ss_h = br.Spectra()
            for i in range(len(ss_hh)):
                _s = br.Spectra([ss_vv[i], ss_vh[i]]).calculate_average()
                ss_v.append(_s)

                _s = br.Spectra([ss_hv[i], ss_hh[i]]).calculate_average()
                ss_h.append(_s)

            # calculation metadata
            parameters = self.get_parameters()
            for _ss in (ss_vv, ss_vh, ss_hv, ss_hh, ss_v, ss_h):
                _ss.E = np.linspace(Emin, Emax, npoints)
                for _name in parameters:
                    _ss.__setattr__(_name, parameters[_name])
                _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])
            
            # create images
            ss  = {'vv':ss_vv, 'vh':ss_vh, 'hv':ss_hv, 'hh':ss_hh, 'v':ss_v, 'h':ss_h}
            ims = {}
            for _pol in ss:
                ims[_pol] = ss[_pol].stack_spectra_as_columns()
                ims[_pol].x_centers = ss[_pol].E

            # create xas
            xas = {}
            for _pol in ss:
                xas[_pol] = br.Spectrum(x=ss[_pol].E, y=ss[_pol].calculate_y_sum())

            # output is everything before first spectrum
            out = _out.split('Here starts Gvv spectrum:')[0]
                
            return ss, ims, xas, out
        elif self.polarization.lower() == 'circular':
            # =============== rv =============== #
            _out2 = _out.split('Here starts Grv spectrum:')[1].split('Here ends Grv spectrum')[0].split('\n')
            _out3 = [_ for _ in _out2 if _ != '']
            ss_rv = br.Spectra()
            if wsl:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.strip().replace('  ', ' ').split(' ')[i]) for _ in _out3[5:]]
                    ss_rv.append(br.Spectrum(_x, _y))
            else:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.split(' ')[i]) for _ in _out3[5:]]
                    ss_rv.append(br.Spectrum(_x, _y))

            # =============== rh =============== #
            _out2 = _out.split('Here starts Grh spectrum:')[1].split('Here ends Grh spectrum')[0].split('\n')
            _out3 = [_ for _ in _out2 if _ != '']
            # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
            ss_rh = br.Spectra()
            if wsl:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.strip().replace('  ', ' ').split(' ')[i]) for _ in _out3[5:]]
                    ss_rh.append(br.Spectrum(_x, _y))
            else:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.split(' ')[i]) for _ in _out3[5:]]
                    ss_rh.append(br.Spectrum(_x, _y))

            # =============== lv =============== #
            _out2 = _out.split('Here starts Glv spectrum:')[1].split('Here ends Glv spectrum')[0].split('\n')
            _out3 = [_ for _ in _out2 if _ != '']
            # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
            _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
            ss_lv = br.Spectra()
            if wsl:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.strip().replace('  ', ' ').split(' ')[i]) for _ in _out3[5:]]
                    ss_lv.append(br.Spectrum(_x, _y))
            else:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.split(' ')[i]) for _ in _out3[5:]]
                    ss_lv.append(br.Spectrum(_x, _y))

            # =============== lh =============== #
            _out2 = _out.split('Here starts Glh spectrum:')[1].split('Here ends Glh spectrum')[0].split('\n')
            _out3 = [_ for _ in _out2 if _ != '']
            # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
            ss_lh = br.Spectra()
            if wsl:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.strip().replace('  ', ' ').split(' ')[i]) for _ in _out3[5:]]
                    ss_lh.append(br.Spectrum(_x, _y))
            else:
                for i in np.arange(2, (npoints+1)*2, 2):
                    _y = [float(_.split(' ')[i]) for _ in _out3[5:]]
                    ss_lh.append(br.Spectrum(_x, _y))
            
            # =============== r and l =============== #
            ss_r = br.Spectra()
            ss_l = br.Spectra()
            for i in range(len(ss_lh)):
                _s = br.Spectra([ss_rv[i], ss_rh[i]]).calculate_average()
                ss_r.append(_s)

                _s = br.Spectra([ss_lv[i], ss_lh[i]]).calculate_average()
                ss_l.append(_s)

            # calculation metadata
            parameters = self.get_parameters()
            for _ss in (ss_rv, ss_rh, ss_lv, ss_lh, ss_r, ss_l):
                _ss.E = np.linspace(Emin, Emax, npoints)
                for _name in parameters:
                    _ss.__setattr__(_name, parameters[_name])
                _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])

            # create images
            ss  = {'rv':ss_rv, 'rh':ss_rh, 'lv':ss_lv, 'lh':ss_lh, 'r':ss_r, 'l':ss_l}
            ims = {}
            for _pol in ss:
                ims[_pol] = ss[_pol].stack_spectra_as_columns()
                ims[_pol].x_centers = ss[_pol].E

            # create xas
            xas = {}
            for _pol in ss:
                xas[_pol] = br.Spectrum(x=ss[_pol].E, y=ss[_pol].calculate_y_sum())

            # output is everything before first spectrum
            out = _out.split('Here starts Grv spectrum:')[0]

            return ss, ims, xas, out
        return 

    def get_initial_wavefunctions(self):
        """Returns a string with details of the initial havefunctions

        Warning:
            This method was only tested for Cu2+ in a Oh and D4h environment.

        Returns:
            string
        """
        # initialization =====================================
        replacements = odict()
        replacements['-- Normalize dZ to unity.'] = text_for_printing_initial_wavefunctions

        # replacement =====================================
        self.update_lua_script()
        for replacement in replacements:
            self.lua_script = self.lua_script.replace(
                replacement, str(replacements[replacement]))

        # create and run lua script
        with tempfile.NamedTemporaryFile(delete=False) as fp:
            fp.write(self.lua_script.encode('utf-8'))
            fp.close()
            _out = quanty(fp.name)
        os.unlink(fp.name)

        return _out.split('Here starts the initial wavefunctions ===')[1].split('Here ends the initial wavefunctions ===')[0]

# %% =========================== Experimental ============================ %% #
# The functions down below are used to print/save the initial wavefunctions 
# specifically for Cu2+ (d9). I didn't test it for other ions.
text_for_printing_initial_wavefunctions = r"""--------------------------------------------------------------------------------
-- base wavefunctions for printing/saving initial wavefunctions
--------------------------------------------------------------------------------
-- z2    = Y^0
-- x2-y2 = 1/sqrt(2) [Y-2 + Y2]
-- xy    = 1/i*sqrt(2) [Y-2 - Y2]
-- yz    = i/sqrt(2) [Y-1 + Y1]
-- xz    = 1/sqrt(2) [Y-1 - Y1]

function replace_char(pos, str, r)
    return str:sub(1, pos-1) .. r .. str:sub(pos+1)
end

base_psis = {}
base_dets = {}
base_sphe = {'Y^-2(-)', 'Y^-2(+)', 'Y^-1(-)', 'Y^-1(+)', 'Y^0(-)', 'Y^0(+)', 'Y^1(-)', 'Y^1(+)', 'Y^2(-)', 'Y^2(+)'}
if H_3d_ligands_hybridization_lmct == 1 or H_3d_ligands_hybridization_mlct == 1 then
    for i = 1, 20, 1 do
        _det  = replace_char(6+i, '11111111111111111111111111', '0') 
        _psi = NewWavefunction(NFermions, NBosons, {{_det, 1}})
        table.insert(base_psis, _psi)

        _det2 = replace_char(6, _det, '1 ') 
        if i < 11 then
            _det2 = replace_char(18, _det2, ' 1') 
        else
            _det2 = replace_char(17, _det2, '1 ') 
        end
        table.insert(base_dets, _det2)
    end
else
    for i = 1, 10, 1 do
        _det  = replace_char(6+i, '1111111111111111', '0') 
        _psi = NewWavefunction(NFermions, NBosons, {{_det, 1}})
        table.insert(base_psis, _psi)

        _det2 = replace_char(6, _det, '1 ') 
        table.insert(base_dets, _det2)
    end
end

-- z2    = Y^0
psi_z2_m    = base_psis[5]
psi_z2_p    = base_psis[6]

-- x2-y2 = 1/sqrt(2) [Y-2 + Y2]
psi_x2_y2_m = math.sqrt(1/2)*base_psis[1] + math.sqrt(1/2)*base_psis[9]
psi_x2_y2_p = math.sqrt(1/2)*base_psis[2] + math.sqrt(1/2)*base_psis[10]

-- xy    = 1/i*sqrt(2) [Y-2 - Y2]
psi_xy_m = (1/I)*math.sqrt(1/2)*base_psis[1] - (1/I)*math.sqrt(1/2)*base_psis[9]
psi_xy_p = (1/I)*math.sqrt(1/2)*base_psis[2] - (1/I)*math.sqrt(1/2)*base_psis[10]

-- yz    = i/sqrt(2) [Y-1 + Y1]
psi_yz_m = I*math.sqrt(1/2)*base_psis[3] + I*math.sqrt(1/2)*base_psis[7]
psi_yz_p = I*math.sqrt(1/2)*base_psis[4] + I*math.sqrt(1/2)*base_psis[8]

-- xz    = 1/sqrt(2) [Y-1 - Y1]
psi_xz_m = math.sqrt(1/2)*base_psis[3] - math.sqrt(1/2)*base_psis[7]
psi_xz_p = math.sqrt(1/2)*base_psis[4] - math.sqrt(1/2)*base_psis[8]


--------------------------------------------------------------------------------
-- functions for saving/printing initial wavefunctions
--------------------------------------------------------------------------------
-- z2    = Y^0
-- x2-y2 = 1/sqrt(2) [Y-2 + Y2]
-- xy    = 1/i*sqrt(2) [Y-2 - Y2]
-- yz    = i/sqrt(2) [Y-1 + Y1]
-- xz    = 1/sqrt(2) [Y-1 - Y1]


function print_psi(Psis_i, base_psis, base_dets, base_sphe, H_i, verbose)

    -- wavefunction name and energy
    for i, psi in ipairs(Psis_i) do
        number_of_eigenfunctions = i
    end
    
    -- initial header
    text = ''
    text = text .. '#########################\n'
    text = text .. '# INITIAL WAVEFUNCTIONS #\n'
    text = text .. '#########################\n'
    text = text .. '# Number of wavefunctions: ' .. number_of_eigenfunctions .. '\n'
    text = text .. '# number Energy wavefunction\n'
    for i, psi in ipairs(Psis_i) do
        text = text .. '    ' .. i .. '      <<E' .. i .. '>>  <<psi' .. i .. '>>' .. '\n'
    end
    text = text .. '# \n'

    for i, psi in ipairs(Psis_i) do
        -- calculate energy
        if H_i ~= nil then
            _E = Complex.Re(psi * H_i * psi) -- Taking the real part only because sometimes the precision error creates a small im component
        else
            _E =  0
        end

        -- wavefunction name and energy
        text = text .. '###################\n'
        text = text .. '# WAVEFUNCTION: ' .. i ..' # \n'
        text = text .. '###################\n'
        text = text .. '# E' .. i .. ' = ' .. _E.. '\n'
        -- text = text .. '# \n'

        -- M holes
        text = text .. '# =============== M holes ===============\n'
        if H_3d_ligands_hybridization_lmct == 1 then
            text = text .. '# 2p       3d         L     Y_2^ml(ms) factor\n'
        else
            text = text .. '# 2p       3d    Y_2^ml(ms) factor\n'
        end
        for j = 1, 10, 1 do
            factor = psi*base_psis[j]
            if math.abs(Complex.Im(factor)) > epsilon then
                _factor = '(' .. Complex.Re(factor) .. ', ' .. Complex.Im(factor) .. '*I)'
                text = text .. base_dets[j]  .. '  ' .. base_sphe[j] .. '    ' .. _factor .. '\n'
            else
                text = text .. base_dets[j]  .. '  ' .. base_sphe[j] .. '    ' .. Complex.Re(factor) .. '\n'
            end
        end

        -- find wavefunction with highest factor
        -- highest = psi*base_psis[1]
        factor = psi*base_psis[1]
        highest = Complex.Re(factor)^2 + Complex.Im(factor)^2
        for j = 2, 10, 1 do
            _factor = psi*base_psis[j]
            _temp = Complex.Re(_factor)^2 + Complex.Im(_factor)^2
            if _temp > highest then
                highest = _temp
            end
        end

        -- find all wavefunctions where contribution is larger than 10% of highest contribution
        _final = 'psi' .. i .. ' = '
        for j = 1, 10, 1 do 
            _factor = psi*base_psis[j]
            _temp = Complex.Re(_factor)^2 + Complex.Im(_factor)^2
            if math.abs(_temp) > math.abs(highest)*0.1 then
                factor = psi*base_psis[j]
                if math.abs(Complex.Im(_factor)) > epsilon then
                    _factor = '(' .. tonumber(string.format("%.3f", math.abs(Complex.Re(factor)))) .. ', ' .. tonumber(string.format("%.3f", math.abs(Complex.Im(factor))))  .. '*I)'
                    _final = _final .. ' + ' .. _factor .. '*' .. base_sphe[j]
                else
                    _factor = tonumber(string.format("%.3f", math.abs(Complex.Re(factor))))
                    if Complex.Re(factor) < 0 then
                        _final = _final .. ' - ' .. _factor .. '*' .. base_sphe[j]
                    else
                        _final = _final .. ' + ' .. _factor .. '*' .. base_sphe[j]
                    end
                end

            end
        end
        -- _final = string.sub(_final, 1, -3)

        -- print final psi_i
        text = text .. '# \n'
        text = text .. '# Estimated wavefunction in terms of spherical harmonics:\n'
        text = text .. _final .. '\n'
        text = text .. '# \n'

        -- L holes
        if H_3d_ligands_hybridization_lmct == 1 or H_3d_ligands_hybridization_mlct == 1 then
            text = text .. '# =============== L holes ===============\n'
            text = text .. '# 2p       3d         L    factor\n'
            for j = 11, 20, 1 do
                factor = psi*base_psis[j]
                if math.abs(Complex.Im(factor)) > epsilon then
                    _factor = '(' .. Complex.Re(factor) .. ', ' .. Complex.Im(factor) .. '*I)'
                    text = text .. base_dets[j]  .. '  ' .. _factor .. '\n'
                else
                    text = text .. base_dets[j]  .. '  '  .. Complex.Re(factor) .. '\n'
                end
                -- text = text .. base_dets[j] .. ' ' .. psi*base_psis[j] .. '\n'
            end
        end

        -- end of wavefunction table
        text = text .. '# \n'

        -- update header
        text = string.gsub(text, '<<psi' .. i .. '>>', _final)
        text = string.gsub(text, '<<E' .. i .. '>>', _E)

    end


    print('Here starts the initial wavefunctions ===')
    print(text)
    print('Here ends the initial wavefunctions ===')

end

print_psi(Psis_i, base_psis, base_dets, base_sphe, H_i, 1)

os.exit(0)


-- Normalize dZ to unity."""
