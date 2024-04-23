#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""BRIXS settings file"""

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

**check**
MAX_ERROR_STEP_X (float): max allowed error between floats, when checking if
    spectra have the same x axis or when checking if spectrum have uniform
    x axis. MAX_ERROR_STEP_X is given in percentage of the average x step.
    Default is 0.1%.
"""

# %% -------------------------- obsolete text ----------------------------- %% #
# FIGURE_FORCE_ON_TOP (bool): new figures are open on top of other windows.
#     Default is False.
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

        #########
        # check #
        #########
        self.MAX_ERROR_STEP_X = 0.1

        ############
        # internal #
        ############
        self._figure_count   = 1
        self._forbidden_words = {'Spectrum':[], 'Spectra':[], 'Image':[], 'PhotonEvents':[]}
        self._reserved_words  = {'Spectrum':    {'methods':[], 'vars':[], 'pseudovars': []}, 
                                 'Spectra':     {'methods':[], 'vars':[], 'pseudovars': []}, 
                                 'Image':       {'methods':[], 'vars':[], 'pseudovars': []}, 
                                 'PhotonEvents':{'methods':[], 'vars':[], 'pseudovars': []}}

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
                f'FIGURE_GRID_OFFSET:      {self.FIGURE_GRID_OFFSET}\n'
        text += '\n======== spectra settings ===========\n'  
        text += f'MAX_ERROR_STEP_X:        {self.MAX_ERROR_STEP_X}\n'
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