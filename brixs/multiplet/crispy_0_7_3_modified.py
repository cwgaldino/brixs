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

Usage:

    >>> import brixs as br
    >>> import brixs.multiplet as multiplet
    >>>
    >>> multiplet.settings.QUANTY_FILEPATH = r'C:\Users\gald\Documents\quanty\quanty_win\QuantyWin64.exe'
    >>>
    >>> todo


Notes:

    The assumed experimental geometry is given below:


    Top view (lab coordinates)
                           .
                         .   .
    LH y (010)         .       .
        |            .           .
    ────O────────  .   octahedron  . ────> x (100)
     LV z (001)      .           .
                       .       .
                         .   .
                           .

    The outgoing beam (kout) is then defined based on the tth angle. See below

          kout
            .
             .             
              .
               .     tth
                .        
    kin ─────── ||-------------  

    Finally, we can define the incoming and outgoing beam direction as

    kin  = 100
    LHin = 010
    LVin = 001

    kout  = kin*Rz(tth)
    LHout = LHin*Rz(tth)
    LHout = LHin*Rz(tth)

    Note that, the direction of the outgoing beam is given by kout, which is 
    defined by a rotation of kin by tth angle (same for the polarization vectors).

    For a octahedral local environment, the default Hamiltonian H (Crystal Field, Exchange direction, etc...) 
    is such that the planar ligands point to the x and y coordinates and the apical 
    ligands are laid along the z axis. Here, I am calling the coordinate system 
    of the Hamiltonian the cf coordinate system. If you want H to be rotated in a 
    different orientation in relation to the incoming beam, one can you the parameter R.
     
    R is a list of rotations one should do with the local environment before 
    running the calculations. For instance R = [['x', 90], ] will rotate the 
    a octahedra 90 degrees around the 'x' axis (100). This will make the apical 
    ligands to point along the 'y' axis (010). 

    One can double check the experimental geometry by using the method 
    q.plot_geometry()


TO DO:
    ( ) Check if broaden works with non-monotonic data.
    ( ) Not sure if broaden works fine. Compare with crispy.
    ( ) gaussian broaden
    ( ) lorentizian broaden

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
def quanty(filepath):
    """Run Quanty.

    Args:
        filepath (string or pathlib.Path): path to file.

    Returns:
        Calculation output (stdout).
    """
    quanty_exe = Path(settings.QUANTY_FILEPATH)
    if is_windows:
        if quanty_exe.is_absolute():
            quanty = subprocess.Popen([f"{quanty_exe}", f"{filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            quanty = subprocess.Popen([f"./{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif is_linux:
        if quanty_exe.is_absolute():
            quanty = subprocess.Popen([f"{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        else:
            quanty = subprocess.Popen([f"./{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elif is_mac:
        quanty = subprocess.Popen([f"./{quanty_exe} {filepath}"], shell=True, close_fds=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # print('here')
    # print(quanty)

    
    output = quanty.stdout.read().decode("utf-8")
    # print(output)

    error  = quanty.stderr.read().decode("utf-8")
    # print(error)

    # print('there')
    if error != '':
        raise RuntimeError(f"Error while reading file: {filepath}. \n {error}")

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
    q._xGamma         = calculation.pop('xGamma')
    q._xLabel         = calculation.pop('xLabel')
    q._xMin           = calculation.pop('xMin')
    q._xMax           = calculation.pop('xMax')
    q._xNPoints       = calculation.pop('xNPoints')
    q._xEdge          = calculation.pop('xEdge')
    q._yGamma         = calculation.pop('yGamma')
    q._yLabel         = calculation.pop('yLabel')
    q._yMin           = calculation.pop('yMin')
    q._yMax           = calculation.pop('yMax')
    q._yNPoints       = calculation.pop('yNPoints')
    q._yEdge          = calculation.pop('yEdge')
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

    Args:
        element (string, optional): transition metals, lanthanides and actinides.
            default is 'Ni'.
        charge (string, optional): suitable oxidation state of the element as
            '1+', '2+' and so on. default depends on the element.
        symmetry (string, optional): local symmetry. Possible values are 'Oh',
            'D4h', 'Td', 'C3v', and 'D3h'. default depends on the element.
        experiment (string, optional): experiment to simulate spectrum. Possible
            values are 'XAS', 'XES', 'XPS', and 'RIXS'. default is 'XAS'.
        edge (string, optional): investigated edge.
            * For 'XAS' and 'XPS', possible values are 'K (1s)', 'L1 (2s)', 'L2,3 (2p)', 'M1 (3s)', and 'M2,3 (3p)'.
            * For 'XES', possible values are 'Ka (1s2p)' and 'Kb (1s3p)'.
            * For 'RIXS', possible values are 'K-L2,3 (1s2p)', 'K-M2,3 (1s3p)', and 'L2,3-M4,5 (2p3d)'
            default is depends on the experiment type and element.

        polarization (list, optional): type of spectrum to calculate.

            * For 'XAS', possible values are `isotropic`, `linear`, `circular`.

            * For 'XES' and 'XPS', possible value is `isotropic`.

            * For 'RIXS', possible values are `isotropic`, `linear`, `circular`.

        nPsis (number, optional): number of states to calculate. If None, 
            nPsiAuto will be set to 1 and the maximum value will be used. 
            Default is None. The max value is defined by q.nPsisMax
        nPsisAuto (int, optional): If 1, a suitable value will be picked for 
            nPsis. Default is 1.
        k (vector, optional): quantization axis for calculating expected values. 
            It does not influence the calculation, just the table printed at the end.
            Default is (001).
        denseBorder (int, optional): Number of determinants in the bases where 
            quanty switch from dense methods to sparse methods. The higher this 
            number, the more precise are the calculations, but also more expensive. 
            Default is 2000, which should be enough. If in doubt, one can try 
            this to see how the final results depends on it. Final result 
            shouldn't improve much for values above 2000.

        magneticField (number, optional): Magnetic field value in Tesla. If zero
            or None, a very small magnetic field (0.002 T) will be used to
            have nice expected values for observables. To force a zero magnetic
            field you have to turn it off via the q.hamiltonianState. default
            is 0.002 T. This can be overwritten by changing the hamiltonian data
            directly.
        magneticFieldOrientation (vector, optional): Direction of the magnetic field.
            Default is (001). This can be overwritten by changing the hamiltonian data
            directly.
        temperature (number, optional): temperature in Kelvin. default is 10 K.
            If temperature is zero, nPsi is set to 1 and nPsiAuto to 0.

        tth (number): number from 0 to 180 indicating the 2theta scattering angle.
            Default is 130.
        R (list): list of rotations one should do with the local environment before 
            running the calculations. For instance R = [['x', 90], ] will rotate the 
            a octahedra 90 degrees around the 'x' axis (100). This will make the apical 
            ligands to point along the 'y' axis (010). Example: R=[['x', 90], ['y', 45]]. 
            One can double check the experimental geometry by using the method 
            q.plot_geometry(). Default is []. 
        
        xMin, xMax, xNPoints (number, optional): minimum, maximum, and number of points.

            * For 'XAS', incident photon (in eV).
            * For 'XES', emission photon energy (in eV)
            * For 'XPS', biding energy (in eV)
            * For 'RIXS', incident photon energy (in eV).

            If None, a suitable value depending on the element will
            be chosen. default is None.

        yMin, yMax, yNPoints (number, optional): minimum, maximum, and number of points.

            * For 'XAS', 'XES' and 'XPS', n.a.
            * For 'RIXS', energy transfer (in eV).
            If None, a suitable value depending on the element will
            be chosen. default is None.
        
        xGamma (number, optional): 
            * For 'XAS', 'XES', and 'XPS', core-hole lifetime (defines lorentzian broadening - FWHM)
            * For 'RIXS', correlation lenght (lorentzian broadening of the `incident energy`)
        yGamma (tuple or list, optional): core-hole lifetime. 
            * For 'XAS', 'XES', and 'XPS', n.a.
            * For 'RIXS', core-hole lifetime (defines lorentzian broadening of the `energy loss`)

    Attributes:
        All initial args are also attributes, plus we have extra attrs to be 
        edited after creating a Calculation object. Note that some of them are
        read-only and some are defined automatically.

        nConfigurations (int, optional): number of configurations. If None,
            a suitable value will be calculated. default is None. If LMCT or
            MLCT are False in hamiltonianState, then the nConfigurations is 
            always 1 (there is only one final configuration). If LMCT or MLCT 
            are True, then there could be multiple final states where the ligand
            lends the Metal a electron (or more) and vice-versa. The maximum number
            of final states depends on how many holes the Metal or the ligand has.
            To check the maximum number of final states allowed, check q.nConfigurationsMax.
            Example: for XAS Cu2+ (d9) the initial state is p6d9 and the final 
            state is p5d10. If LMCT is `on`, then we have two possible initial 
            states, p6 d9 L10 and p6 d10 L9. nConfigurations is basically, the 
            maximum number of electrons to be shared between M and L.
        nConfigurationsMax (int, read only): Maximum number of final configuration allowed.
        nPsisMax (int, read only): maximum possible number of wavefunctions

        hamiltonianState (dict, optional): a dictionary that turns on or off the
            different contributions to the hamiltonian. default is such that
            'Atomic', 'Crystal Field', and 'Magnetic Field' terms will be True.
        hamiltonianData (dict, optional): dictionary with values for the strength
            of each different contributions to the hamiltonian. The default value
            ise different for each element and charge stat. Values are in eV.

        templatePath, templateName: where template lus script is stored
        template (string): loaded template
        lua_script (string): lua script to be run
        
        verbosity: Verbosity for quanty. By default it is set to minimal verbosity.
        
        yLabel, xLabel: Axis labels.
        nElectrons: number of electrons
        configurations: Initial, Intermediate, and final configurations. Not always accurate I think.
        block: d block, f block, ...

        xEdge: Energy value of the resonance.
        yEdge: Not used.

        kin, LVin, LHin: incoming vectors in lab coordinates
        kout, LVout, LHout: outgoing vectors in lab coordinates
        kin_cf, LVin_cf, LHin_cf: incoming vectors in cf coordinates
        kout_cf, LVout_cf, LHout_cf: outgoing vectors in cf coordinates
        cf: cf cordinate system in terms of lab coordinate system.
        
    Methods:
        q.get_attrs(): returns a list of available attrs
        q.get_methods(): returns a list of methods available

        q.plot_geometry(): Plot a 3D sketch of the experimental geometry

        q.get_parameters(): Returns a dictionary with all calculation parameters.
        q.save_parameters(): Save calculation parameters to a file.

        q.run(): run quanty calculation and return spectra.
        q.get_initial_wavefunctions(): Only tested for Cu2+ Oh and d4h. Returns the 
            initial groundstate wavefunctions

    Extra functions from the module:
        out = multiplet.quanty(): Run quanty calculation from a filepath.
        q   = multiplet.load_calculation(): Load calculation from parameters file.
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
                       magneticField = 0,
                       magneticFieldOrientation = [0, 0, 1],
                       temperature   = 10,
                       #
                       tth   = 130,
                       R     = [],
                       #
                       xMin     = None,
                       xMax     = None,
                       xNPoints = None,
                       xEdge    = None, # resonance value
                       yMin     = None,
                       yMax     = None,
                       yNPoints = None,
                       yEdge    = None, # ????
                       #
                       xGamma = None, #0.1, # core-hole lifetime in eV (For RIXS: incident energy core-hole lifetime in eV)
                       yGamma = None, # For RIXS: energy loss core-hole lifetime in eV
                       ):
        # primary attributes
        self._set_primary_attributes(element, charge, symmetry, experiment, edge)

        # set secondary parameters
        branch = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]
        # print(branch)
        self._configurations  = branch['configurations']
        self._block           = self.configurations[0][1][:2]
        self._nElectrons      = int(self.configurations[0][1][2:])

        # Calculation parameters
        self._nPsis           = branch['number of states']
        self._nPsisMax        = self.nPsis
        self.nPsisAuto        = nPsisAuto
        if self.nPsisAuto == 0:
            self.nPsis = nPsis
        self.polarization        = polarization       
        self._nConfigurationsMax = 1
        self.nConfigurations     = 1  # this is always 1 in Crispy if LMCT or MLCT is False
        self.hamiltonianTerms = branch['hamiltonian terms']
        
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
        self._update_geometry()

        # spectrum attributes
        self._xLabel  = branch['axes'][0][0]
        self.xMin     = xMin
        self.xMax     = xMax
        self.xNPoints = xNPoints
        self.xEdge    = xEdge
        self.xGamma   = xGamma
        if self.experiment == 'RIXS':
            self._yLabel  = branch['axes'][1][0]
            self.yMin     = yMin
            self.yMax     = yMax
            self.yNPoints = yNPoints
            self.yEdge    = yEdge
            self.yGamma   = yGamma
        else:
            self._yLabel  = None
            self.yMin     = None
            self.yMax     = None
            self.yNPoints = None
            self.yEdge    = None
            self.yGamma   = None

        # set basename (for filepath_lua and filepath_spec)
        # if self.experiment in ['RIXS', 'XES']:
        #     _shortEdge = self.edge[-5:-1]
        # else:
        #     _shortEdge = self.edge[-3:-1]
        # self._baseName = '{}{}_{}_{}_{}'.format(
        #     self.element, self.charge, self.symmetry, _shortEdge,
        #     self.experiment)

        # get template filename
        self.templateName = branch['template name']
        self.templatePath = settings.TEMPLATES_FOLDERPATH/self.templateName
        with open(self.templatePath) as p:
            self.template = p.read()

        # other attributes
        self.verbosity   = '0x0000'  # in the future, verbosity could be a parameters
        self.denseBorder = str(denseBorder)  # Default 2000
        self._lua_script = None  # stores lua script

        # hamiltonianData and hamiltonianState
        self.hamiltonianData = odict()
        self.hamiltonianState = odict()
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
        self.hamiltonianState = _hamiltonianState(self.hamiltonianState, parent=self)
        self.hamiltonianData  = _hamiltonianData(self.hamiltonianData)

        # experiment attributes
        self.temperature   = temperature
        self._magneticField = magneticField 
        self._magneticFieldOrientation = magneticFieldOrientation 
        self._update_magnetic_field_hamiltonian_data()

        
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
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge).\nPlease, start a new Calculation() object')
    @element.deleter
    def element(self):
        raise AttributeError('Cannot delete object.')

    @property
    def charge(self):
        return self._charge
    @charge.setter
    def charge(self, value):
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge).\nPlease, start a new Calculation() object')
    @charge.deleter
    def charge(self):
        raise AttributeError('Cannot delete object.')

    @property
    def symmetry(self):
        return self._symmetry
    @symmetry.setter
    def symmetry(self, value):
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge).\nPlease, start a new Calculation() object')
    @symmetry.deleter
    def symmetry(self):
        raise AttributeError('Cannot delete object.')

    @property
    def experiment(self):
        return self._experiment
    @experiment.setter
    def experiment(self, value):
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge).\nPlease, start a new Calculation() object')
    @experiment.deleter
    def experiment(self):
        raise AttributeError('Cannot delete object.')

    @property
    def edge(self):
        return self._edge
    @edge.setter
    def edge(self, value):
        raise AttributeError('Primary attributes cannot be changed (element, charge, symmetry, experiment, edge).\nPlease, start a new Calculation() object')
    @edge.deleter
    def edge(self):
        raise AttributeError('Cannot delete object.')
    

    # ATTRIBUTES SET AUTOMATICALY WITH PRIMARY ATTRIBUTES
    @property
    def xLabel(self):
        return self._xLabel
    @xLabel.setter
    def xLabel(self, value):
        raise AttributeError('`xLabel` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed.\nPlease, start a new Calculation() object')
    @xLabel.deleter
    def xLabel(self):
        raise AttributeError('Cannot delete object.')

    @property
    def yLabel(self):
        return self._yLabel
    @yLabel.setter
    def yLabel(self, value):
        raise AttributeError('`yLabel` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed.\nPlease, start a new Calculation() object')
    @yLabel.deleter
    def yLabel(self):
        raise AttributeError('Cannot delete object.')

    @property
    def nElectrons(self):
        return self._nElectrons
    @nElectrons.setter
    def nElectrons(self, value):
        raise AttributeError('`nElectrons` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed.\nPlease, start a new Calculation() object')
    @nElectrons.deleter
    def nElectrons(self):
        raise AttributeError('Cannot delete object.')

    @property
    def configurations(self):
        return self._configurations
    @configurations.setter
    def configurations(self, value):
        raise AttributeError('`configurations` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed.\nPlease, start a new Calculation() object')
    @configurations.deleter
    def configurations(self):
        raise AttributeError('Cannot delete object.')

    @property
    def block(self):
        return self._block
    @block.setter
    def block(self, value):
        raise AttributeError('`block` is defined based primary attributes (element, charge, symmetry, experiment, edge) and cannot be changed.\nPlease, start a new Calculation() object')
    @block.deleter
    def block(self):
        raise AttributeError('Cannot delete object.')



    # SPECTRA ATTRIBUTES
    @property
    def xMin(self):
        return self._xMin
    @xMin.setter
    def xMin(self, value):
        if value is None:
            value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][1]
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
            value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][2]
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
            value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][3]
        assert value >= 2, 'xNPoints cannot be less than 2.\nThe CreateResonantSpectra() function from Quanty prevents xNPoints to be less than 2.'
        self._xNPoints = int(value)
    @xNPoints.deleter
    def xNPoints(self):
        raise AttributeError('Cannot delete object.')

    @property
    def xEdge(self):
        return self._xEdge
    @xEdge.setter
    def xEdge(self, value):
        if value is None:
            value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][4]
        assert value > 0, 'xEdge cannot be negative'
        self._xEdge = value
    @xEdge.deleter
    def xEdge(self):
        raise AttributeError('Cannot delete object.')


    @property
    def yMin(self):
        return self._yMin
    @yMin.setter
    def yMin(self, value):
        if value is None:
            if self.experiment != 'RIXS':
                value = None
            else:
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][1][1]
        self._yMin = value
    @yMin.deleter
    def yMin(self):
        raise AttributeError('Cannot delete object.')

    @property
    def yMax(self):
        return self._yMax
    @yMax.setter
    def yMax(self, value):
        if value is None:
            if self.experiment != 'RIXS':
                value = None
            else:
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][1][2]
        self._yMax = value
    @yMax.deleter
    def yMax(self):
        raise AttributeError('Cannot delete object.')

    @property
    def yNPoints(self):
        return self._yNPoints
        # try:
        #     return self._yNPoints + 1
        # except TypeError:
        #     return self._yNPoints
    @yNPoints.setter
    def yNPoints(self, value):
        if value is None:
            if self.experiment != 'RIXS':
                self._yNPoints = None
                return
            else:
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][1][3]
        assert value >= 2, 'yNPoints cannot be less than 2.\nThe CreateResonantSpectra() function from Quanty prevents yNPoints to be less than 2.'
    
        self._yNPoints = int(value)
    @yNPoints.deleter
    def yNPoints(self):
        raise AttributeError('Cannot delete object.')

    @property
    def yEdge(self):
        return self._yEdge
    @yEdge.setter
    def yEdge(self, value):
        if value is None:
            if self.experiment != 'RIXS':
                self._yEdge = None
                return
            else:
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][1][4]
        assert value > 0, 'yEdge cannot be negative'
        self._yEdge = value
    @yEdge.deleter
    def yEdge(self):
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
        raise AttributeError('`cf` is defined via _update_geometry() and depends on q.R attribute')
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
        assert len(value) == 3, 'kout must be a vector, like [0, 1, 0].'
        self._kout = br.normalize_vector(value)
        self._update_geometry()
    @kout.deleter
    def kout(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LVout(self):
        return self._LVout
    @LVout.setter
    def LVout(self, value):
        assert len(value) == 3, 'LVout must be a vector, like [0, 1, 0].'
        self._LVout = br.normalize_vector(value)
        self._update_geometry()
    @LVout.deleter
    def LVout(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LHout(self):
        return self._LHout
    @LHout.setter
    def LHout(self, value):
        assert len(value) == 3, 'LHout must be a vector, like [0, 1, 0].'
        self._LHout = br.normalize_vector(value)
        self._update_geometry()
    @LHout.deleter
    def LHout(self):
        raise AttributeError('Cannot delete object.')
    
    @property
    def kin_cf(self):
        return self._kin_cf
    @kin_cf.setter
    def kin_cf(self, value):
        assert len(value) == 3, 'k1 must be a vector, like [0, 1, 0].'
        self._kin_cf = br.normalize_vector(value)
        self._update_geometry()
    @kin_cf.deleter
    def kin_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LVin_cf(self):
        return self._LVin_cf
    @LVin_cf.setter
    def LVin_cf(self, value):
        assert len(value) == 3, 'LVin_cf must be a vector, like [0, 1, 0].'
        self._LVin_cf = br.normalize_vector(value)
        self._update_geometry()
    @LVin_cf.deleter
    def LVin_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LHin_cf(self):
        return self._LHin_cf
    @LHin_cf.setter
    def LHin_cf(self, value):
        assert len(value) == 3, 'LHin_cf must be a vector, like [0, 1, 0].'
        self._LHin_cf = br.normalize_vector(value)
        self._update_geometry()
    @LHin_cf.deleter
    def LHin_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def kout_cf(self):
        return self._kout_cf
    @kout_cf.setter
    def kout_cf(self, value):
        assert len(value) == 3, 'kout_cf must be a vector, like [0, 1, 0].'
        self._kout_cf = br.normalize_vector(value)
        self._update_geometry()
    @kout_cf.deleter
    def kout_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LVout_cf(self):
        return self._LVout_cf
    @LVout_cf.setter
    def LVout_cf(self, value):
        assert len(value) == 3, 'LVout_cf must be a vector, like [0, 1, 0].'
        self._LVout_cf = br.normalize_vector(value)
        self._update_geometry()
    @LVout_cf.deleter
    def LVout_cf(self):
        raise AttributeError('Cannot delete object.')

    @property
    def LHout_cf(self):
        return self._LHout_cf
    @LHout_cf.setter
    def LHout_cf(self, value):
        assert len(value) == 3, 'LHout_cf must be a vector, like [0, 1, 0].'
        self._LHout_cf = br.normalize_vector(value)
        self._update_geometry()
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
        if value is None or value == 0:
            value = 0.002
        elif value < 0:
            raise ValueError('Magnetic field value must be positive.')
        elif value < 0.002:
            raise ValueError('Magnetic field cannot be smaller than 0.002 T.\nTurn off the magnetic field hamiltonian using hamiltonianState["Magnetic Field"] = False')     
       
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
        self._nPsisAuto = 0
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
        self._templatePath = value
    @templatePath.deleter
    def templatePath(self):
        raise AttributeError('Cannot delete object.')

    @property
    def templateName(self):
        return self._templateName
    @templateName.setter
    def templateName(self, value):
        self._templateName = value
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
    def xGamma(self):
        return self._xGamma
    @xGamma.setter
    def xGamma(self, value):
        if value is None:
            value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][0][5][0]
        # assert value > 0, 'xGamma cannot be negative'
        self._xGamma = value
    @xGamma.deleter
    def xGamma(self):
        raise AttributeError('Cannot delete object.')

    @property
    def yGamma(self):
        return self._yGamma
    @yGamma.setter
    def yGamma(self, value):
        if value is None:
            if self.experiment == 'RIXS':
                value = settings.PARAMETERS['elements'][self.element]['charges'][self.charge]['symmetries'][self.symmetry]['experiments'][self.experiment]['edges'][self.edge]['axes'][1][5][0]       
            else:
                self._yGamma = None
                return
        # assert value > 0, 'yGamma cannot be negative'
        self._yGamma = value
    @yGamma.deleter
    def yGamma(self):
        raise AttributeError('Cannot delete object.')


    # SUPPORT METHODS
    def _update_geometry(self):
        """Converts lab coordinates to CF coordinates."""
        # geometry (lab coordinates)
        # lab coordinates: z points up and the beam points in x
        # LH, LV, kin, 
        # LHout, LVout, kout
        LVin = (0, 0, 1)
        LHin = (0, 1, 0)
        kin  = (1, 0, 0)
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
    
        self._kin   = kin
        self._LHin  = LHin
        self._LVin  = LVin
        self._kout  = kout
        self._LHout = LHout
        self._LVout = LVout

        self._kin_cf   = kin_cf
        self._LVin_cf  = LVin_cf
        self._LHin_cf  = LHin_cf
        self._kout_cf  = kout_cf
        self._LVout_cf = LVout_cf
        self._LHout_cf = LHout_cf
        self._cf    = [x_prime, y_prime, z_prime]
        return
    
    def _update_magnetic_field_hamiltonian_data(self):
        """Updates the hamiltanian data given Mag. field and Mag. field orientation"""
        value = self.magneticField
        k1    = np.array(self.magneticFieldOrientation)

        TESLA_TO_EV = 5.788e-05
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
        return [key for key in self.__dict__.keys() if key.startswith('_') == False] + [key[1:] for key in self.__dict__.keys() if key.startswith('_') == True] 

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

        if show_legend:
            br.leg(ax=ax2)


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
        replacements['$Emin1'] = self.xMin
        replacements['$Emax1'] = self.xMax
        replacements['$NE1'] = self.xNPoints - 1
        replacements['$Eedge1'] = self.xEdge
        replacements['$Gamma1'] = self.xGamma
        if self.experiment == 'XES':
            replacements['$Emin1'] = self.xMin + 20
            replacements['$Emax1'] = self.xMax + 20
        if self.experiment in ['RIXS', ]:
            # The Lorentzian broadening along the incident axis cannot be
            # changed in the interface, and must therefore be set to the
            # final value before the start of the calculation.
            # replacements['$Gamma1'] = self.xLorentzian
            replacements['$Emin2']  = self.yMin
            replacements['$Emax2']  = self.yMax
            replacements['$NE2']    = self.yNPoints - 1
            replacements['$Eedge2'] = self.yEdge
            replacements['$Gamma2'] = self.yGamma  # For RIXS: Broadening of the energy loss

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

                    if self.magneticField == 0:
                        small = np.finfo(np.float32).eps  # ~1.19e-7
                        if parameter == 'Bx':
                            value = k1[0] * small
                        elif parameter == 'By':
                            value = k1[1] * small
                        elif parameter == 'Bz':
                            value = k1[2] * small

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
                    geometry = dict(tth      = self.tth,
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
                    experiment = dict(temperature      = self.temperature,
                                      magneticField    = self.magneticField,
                                      magneticFieldOrientation = list(self.magneticFieldOrientation)),
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
                                       xGamma         = self.xGamma,
                                       xLabel         = self.xLabel,
                                       xMin           = self.xMin,
                                       xMax           = self.xMax,
                                       xNPoints       = self.xNPoints,
                                       xEdge          = self.xEdge,
                                       yGamma         = self.yGamma,
                                       yLabel         = self.yLabel,
                                       yMin           = self.yMin,
                                       yMax           = self.yMax,
                                       yNPoints       = self.yNPoints,
                                       yEdge          = self.yEdge,
                                       verbosity      = self.verbosity,
                                       denseBorder    = self.denseBorder)
            )

    def save_parameters(self, filepath):
        """Save calculation parameters to a text file"""
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
            update (bool, optional): if True, updates the lua script.

        Returns:
            XPS returns 1 spectrum and calculation output

            XES returns 1 spectrum and calculation output

            XAS returns 1 spectrum (iso) or 2 spectra and calculation output

            RIXS returns 1 spectra (iso) or 4 spectra and calculation output
            For RIXS, each spectrum inside a Spectra objects has different 
            incident energy (according to xMin, xMax, and xNpoints). 

            XPS, isotropic
                iso, out = q.run()
            
            XES, isotropic
                iso, out = q.run()

            RIXS, isotropic
                iso, out = q.run()
            RIXS, linear
                vv, vh, hv, hh, out = q.run()
            RIXS, circular
                rv, rh, lv, lh, out = q.run()
        
            XAS, isotropic
                iso, out = q.run()
            XAS, linear
                LV, LH, out = q.run()
            XAS, circular
                CR, CL, out = q.run()
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
        _x = np.linspace(self.xMin, self.xMax, self.xNPoints)
        # _x = np.linspace(self.xMin, self.xMax, self.xNPoints + 1)
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
            step = (self.xMax - self.xMin)/(self.xNPoints)
            _x = np.linspace(self.xMin, self.xMax-step, self.xNPoints)
            _x2 = np.linspace(self.yMin, self.yMax, self.yNPoints)
            if self.polarization.lower() == 'isotropic':
                _out2 = _out.split('Here starts Giso spectrum:')[1].split('Here ends Giso spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                ss = br.Spectra()
                for j in range(self.xNPoints):
                    # _x2 = [float(_.split(' ')[0]) for _ in _out3[5:]]
                    _y = [float(_.split(' ')[2+j*2]) for _ in _out3[5:]]
                    _s = br.Spectrum(_x2, _y)
                    _s.E = _x[j]
                    ss.append(_s)
                ss.E = _x
                out = _out.split('Here starts Giso spectrum:')[0]

                parameters = self.get_parameters()
                for _ss in (ss, ):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])

                return ss, out
            elif self.polarization.lower() == 'linear':
                _out2 = _out.split('Here starts Gvv spectrum:')[1].split('Here ends Gvv spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                ss_vv = br.Spectra()
                for j in range(self.xNPoints):
                    # _x2 = [float(_.split(' ')[0]) for _ in _out3[5:]]
                    _y = [float(_.split(' ')[2+j*2]) for _ in _out3[5:]]
                    _s = br.Spectrum(_x2, _y)
                    _s.E = _x[j]
                    ss_vv.append(_s)
                ss_vv.E = _x
                out = _out.split('Here starts Gvv spectrum:')[0]

                _out2 = _out.split('Here starts Gvh spectrum:')[1].split('Here ends Gvh spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                ss_vh = br.Spectra()
                for j in range(self.xNPoints):
                    # _x2 = [float(_.split(' ')[0]) for _ in _out3[5:]]
                    _y = [float(_.split(' ')[2+j*2]) for _ in _out3[5:]]
                    _s = br.Spectrum(_x2, _y)
                    _s.E = _x[j]
                    ss_vh.append(_s)
                ss_vh.E = _x

                _out2 = _out.split('Here starts Ghv spectrum:')[1].split('Here ends Ghv spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                ss_hv = br.Spectra()
                for j in range(self.xNPoints):
                    # _x2 = [float(_.split(' ')[0]) for _ in _out3[5:]]
                    _y = [float(_.split(' ')[2+j*2]) for _ in _out3[5:]]
                    _s = br.Spectrum(_x2, _y)
                    _s.E = _x[j]
                    ss_hv.append(_s)
                ss_hv.E = _x

                _out2 = _out.split('Here starts Ghh spectrum:')[1].split('Here ends Ghh spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                ss_hh = br.Spectra()
                for j in range(self.xNPoints):
                    # _x2 = [float(_.split(' ')[0]) for _ in _out3[5:]]
                    _y = [float(_.split(' ')[2+j*2]) for _ in _out3[5:]]
                    _s = br.Spectrum(_x2, _y)
                    _s.E = _x[j]
                    ss_hh.append(_s)
                ss_hh.E = _x

                parameters = self.get_parameters()
                for _ss in (ss_vv, ss_vh, ss_hv, ss_hh):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                        for _s in _ss:
                            _s.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])
                    
                return ss_vv, ss_vh, ss_hv, ss_hh, out
            elif self.polarization.lower() == 'circular':
                _out2 = _out.split('Here starts Grv spectrum:')[1].split('Here ends Grv spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                ss_rv = br.Spectra()
                for j in range(self.xNPoints):
                    # _x2 = [float(_.split(' ')[0]) for _ in _out3[5:]]
                    _y = [float(_.split(' ')[2+j*2]) for _ in _out3[5:]]
                    _s = br.Spectrum(_x2, _y)
                    _s.E = _x[j]
                    ss_rv.append(_s)
                ss_rv.E = _x
                out = _out.split('Here starts Grv spectrum:')[0]

                _out2 = _out.split('Here starts Grh spectrum:')[1].split('Here ends Grh spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                ss_rh = br.Spectra()
                for j in range(self.xNPoints):
                    # _x2 = [float(_.split(' ')[0]) for _ in _out3[5:]]
                    _y = [float(_.split(' ')[2+j*2]) for _ in _out3[5:]]
                    _s = br.Spectrum(_x2, _y)
                    _s.E = _x[j]
                    ss_rh.append(_s)
                ss_rh.E = _x

                _out2 = _out.split('Here starts Glv spectrum:')[1].split('Here ends Glv spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                ss_lv = br.Spectra()
                for j in range(self.xNPoints):
                    # _x2 = [float(_.split(' ')[0]) for _ in _out3[5:]]
                    _y = [float(_.split(' ')[2+j*2]) for _ in _out3[5:]]
                    _s = br.Spectrum(_x2, _y)
                    _s.E = _x[j]
                    ss_lv.append(_s)
                ss_lv.E = _x

                _out2 = _out.split('Here starts Glh spectrum:')[1].split('Here ends Glh spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                ss_lh = br.Spectra()
                for j in range(self.xNPoints):
                    # _x2 = [float(_.split(' ')[0]) for _ in _out3[5:]]
                    _y = [float(_.split(' ')[2+j*2]) for _ in _out3[5:]]
                    _s = br.Spectrum(_x2, _y)
                    _s.E = _x[j]
                    ss_lh.append(_s)
                ss_lh.E = _x

                parameters = self.get_parameters()
                for _ss in (ss_rv, ss_rh, ss_lv, ss_lh):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                        for _s in _ss:
                            _s.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])

                return ss_rv, ss_rh, ss_lv, ss_lh, out
        elif self.experiment == 'XAS':
            if self.polarization.lower() == 'isotropic':
                _out2 = _out.split('Here starts ISO spectrum:')[1].split('Here ends ISO spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5::]]
                iso = br.Spectrum(_x, _y)
                out = _out.split('Here starts ISO spectrum:')[0]

                parameters = self.get_parameters()
                for _ss in (iso, ):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])

                return iso, out
            elif self.polarization.lower() == 'linear':
                _out2 = _out.split('Here starts LV spectrum:')[1].split('Here ends LV spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                cr = br.Spectrum(_x, _y)

                _out2 = _out.split('Here starts LH spectrum:')[1].split('Here ends LH spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                cl = br.Spectrum(_x, _y)

                out = _out.split('Here starts LV spectrum:')[0]

                parameters = self.get_parameters()
                for _ss in (cr, cl):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])
                
                return cr, cl, out
            elif self.polarization.lower() == 'circular':
                _out2 = _out.split('Here starts CR spectrum:')[1].split('Here ends CR spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                cr = br.Spectrum(_x, _y)

                _out2 = _out.split('Here starts CL spectrum:')[1].split('Here ends CL spectrum')[0].split('\n')
                _out3 = [_ for _ in _out2 if _ != '']
                # _x = [float(_.split(' ')[0]) for _ in _out3[5:]]
                _y = [float(_.split(' ')[2]) for _ in _out3[5:]]
                cl = br.Spectrum(_x, _y)

                out = _out.split('Here starts CR spectrum:')[0]

                parameters = self.get_parameters()
                for _ss in (cr, cl):
                    for _name in parameters:
                        _ss.__setattr__(_name, parameters[_name])
                    _ss.hamiltonian['hamiltonianData'] = remove_greek_letters_from_hamiltonianData(_ss.hamiltonian['hamiltonianData'])

                return cr, cl, out
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
