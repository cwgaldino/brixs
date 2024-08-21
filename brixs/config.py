#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""BRIXS settings file


###################
# DEVELOPERS NOTE #
###################

Internal variables:

################################
# Reserved and forbidden words #
################################
1) These words will raise an error if the user try to define attrs with these words

2) This avoids errors when converting from one object to another (e.g. im -> s)

3) This also avoid that import methods or variables are overwritten

reserved words:  methods, vars, and pseudovars for a object, i.e., the reserved 
words of a Spectrum object are all the words representing methods, vars, and pseudovars
in a Spectrum object (pseudovars are vars defined via @property decorator)

forbidden words: all methods from a object + pseudovar from all other objects, i.e.,
the forbidden words of a Spectrum object are all methods defined for a spectrum 
object plus all methods and pseudovars of all other objects (as long the pseudovars
from other objects are not equal to pseudovars in Spectrum).

#########
# extra #
#########
`extra` allows for posterior addition of new methods and attrs after brixs has 
been imported. For example, the `model` module requires
 the `model` attr to be initialized inside Spectrum.__init__(). However, the
module `model` only exist if we import it, see below

import brixs as br
import brixs.model

to make the attr brixs.model.model() to be added inside Spectrum.__init__() we
add the model class to extra. See below.

br.settings._extra['Spectrum']['model'] = Model

Note that `self` is passed to Model.

#############
# Modifiers #
#############
settings._modifiers allow for adding new objects that respond accordingly when
set_shift, set_calib, ... are called. For example, when one calls s.set_shift(10)
them the spectrum s is shifted by 10, and also everything else inside 
settings._modifiers['shift']. If 'model' is in settings._modifiers['shift'], then
s.set_shift() also calls s.model.set_shift().

For now, settings._modifiers is only implemented for Spectrum and Spectra. But 
it should be trivial to expand it to Image and PhotonEvents. I just did not do 
it yet, because I see no need.

"""

# %% ---------------------------- help text ------------------------------- %% #
text = """**figure**
FIGURE_POSITION (tuple): default position on screen in px (x, y) for new figures.
    If None, the matplotlib default will be used. Default is None.

FIGURE_SIZE (tuple): default size in px (w, h) for new figures. If None, the 
    matplotlib default will be used. Default is None. This is overwritten by
    `figsize` argument when calling a new figure.

FIGURE_DPI (float): default DPI value for new figures. If None, the 
    matplotlib default will be used. Default is None.

FIGURE_FORCE_NEW_WINDOW (bool): calling `plot()` functions from within brixs,
    forces a new figure to open every time. Default is False. Do not affect 
    figures created by standard matplotlib functions like `plt.figure()`.

FIGURE_GRID (tuple): if not False, FIGURE_GRID will spread new figures
    from left to right and top to bottom on the screen. The number of rows
    and columns if defined by FIGURE_GRID = (rows, cols). Default is False.

FIGURE_GRID_OFFSET (tuple): if FIGURE_GRID is not False, FIGURE_GRID_OFFSET 
    defines the spacing in px between new figures in the grid. Default is (40, 0).
"""

# %% -------------------------- obsolete text ----------------------------- %% #
# FIGURE_FORCE_ON_TOP (bool): new figures are open on top of other windows.
#     Default is False.
# **check**
# MAX_ERROR_STEP_X (float): max allowed error between floats, when checking if
#     spectra have the same x axis or when checking if spectrum have uniform
#     x axis. MAX_ERROR_STEP_X is given in percentage of the average x step.
#     Default is 0.1%.
# %%

# %% =============================== settings ============================= %% #
class _settings():
    """BRIXS settings"""
    def __init__(self):
        ##########
        # figure #
        ##########
        self.FIGURE_POSITION          = None
        self.FIGURE_SIZE              = None
        self.FIGURE_DPI               = None

        self.FIGURE_FORCE_NEW_WINDOW  = False
        # self.FIGURE_FORCE_ON_TOP      = False

        self.FIGURE_GRID              = False  # (rows, cols) -> (1, 1) same as False 
        self.FIGURE_GRID_OFFSET       = (40, 0)

        ############
        # internal #
        ############
        self._forbidden_words = {'Spectrum':[], 'Spectra':[], 'Image':[], 'PhotonEvents':[], 'Dummy':[]}
        self._reserved_words  = {'Spectrum':    {'methods':[], 'vars':[], 'pseudovars': []}, 
                                 'Spectra':     {'methods':[], 'vars':[], 'pseudovars': []}, 
                                 'Image':       {'methods':[], 'vars':[], 'pseudovars': []}, 
                                 'PhotonEvents':{'methods':[], 'vars':[], 'pseudovars': []},
                                 'Dummy':       {'methods':[], 'vars':[], 'pseudovars': []}}
        self._extra    = {'Spectrum':{}, 'Spectra':{}, 'Image':{}, 'PhotonEvents':{}, 'Dummy':{}}
        self._modfiers = {'shift':[], 'offset':[], 'factor':[], 'calib':[]}

        #############
        # help text #
        #############
        self._help = text

    def __str__(self):
        text  =  '========== figure settings ============\n'
        text += f'FIGURE_POSITION:         {self.FIGURE_POSITION}\n' +\
                f'FIGURE_SIZE:             {self.FIGURE_SIZE}\n' +\
                f'FIGURE_DPI:              {self.FIGURE_DPI}\n' +\
                f'FIGURE_FORCE_NEW_WINDOW: {self.FIGURE_FORCE_NEW_WINDOW}\n' +\
                f'FIGURE_GRID:             {self.FIGURE_GRID}\n' +\
                f'FIGURE_GRID_OFFSET:      {self.FIGURE_GRID_OFFSET}'
        # text += '\n======== spectra settings ===========\n'  
        # text += f'MAX_ERROR_STEP_X:        {self.MAX_ERROR_STEP_X}\n'
        return text
                # f'FIGURE_FORCE_ON_TOP:     {self.FIGURE_FORCE_ON_TOP}\n' +\

    def help(self):
        print(self._help)
    
    def pretty_print(self):
        print(self.__str__())
# %%
        
# %% ============================ initialization ========================== %% #
global settings
settings = _settings()