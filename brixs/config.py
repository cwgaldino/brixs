#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""BRIXS settings file

################################################################################
######################### Reserved and forbidden words #########################
################################################################################
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

################################################################################
############################ brixs expansion modules ###########################
################################################################################
brixs expansion modules allows for posterior addition of new methods and attrs 
after brixs has been imported. 

For a module to be considered a `expansion module` it must satisfie two conditions:

1. brixs must not depend on it
2. the module must create objects that are attributes of brixs objects, where 
    objects are required to be initialized with __init__(), or objects are 
    required to react when s.copy(), s.set_shift() (and other modifiers) are 
    called.*
     
* The reason for an object to be initialized with __init__() is that the object
can have access to the data inside the brixs object. For example,

We have an expansion module that creates a object type Model() for each 
instance of Spectrum, i.e.,

s1.model = Model(s1)
s2.model = Model(s2)

note that s1 and s2 are an argument for initializing Model(). This way, 
model from s1 is different from model from s2. s1.model has access to data
from s1 and s2.model has access to data from s2.

In this example, we are manually initializing s.model. If we want to 
programmatically initialize s.model, we have to do it during s.__init__()
so `self` can be passed to Model(). 

To do that, we must include Model() in `br.settings._init['Spectrum']`.

Note that, when we do

import brixs as br

From now on, the definition of Spectrum is well defined. However, 
Spectrum.model does not exist, because Spectrum.model isn't defined at
brixs.py.

Then we import the `model module`

import brixs.model

Once this is imported, model will be included in the Spectrum.__init__() of 
Spectrum via br.settings._init['Spectrum'] 

print(br.settings._init['Spectrum']) ---> {'model': Model}

From now on, every new Spectrum object will come equipped with s.model -> Model()

##################
# settings._init #
##################

br.settings._init['brixs_object']  --->  {'attr_name': Object, }

for every entry inside settings._init, a attr will be created with name 
'attr_name' and it will be linked to an Object. The parent brixs object (`self`)
 is passed to the object. Example,

br.settings._init['Spectrum'] = {'model': Model, }

When a Spectrum s1 is initialized, s1 will have an attr 'model'

s1 = br.Spectrum()
print(s1.model) ---> Model()

Where Model was created by calling 

br.model.Model(parent=s1)

##################
# settings._copy #
##################

br.settings._copy['brixs_object']  --->  ['attr_name', ]

for every entry in settings._copy, the attr_name.copy() will be called together with
brixs_object.copy(). Example.

br.settings._copy['Spectrum'] = ['model', ]

When we call copy() for Spectrum s1 like, 

s2 = s1.copy()

we have that s1.model.copy() is also called. Internally, this is happening:

s2.model = s1.model.copy()

Therefore, s2 is a complete copy of s1, including s1.model.


####################################
# _shift, _offset, _factor, _calib #
####################################

br.settings._shift['brixs_object']  --->  ['attr_name', ]
br.settings._offset['brixs_object']  --->  ['attr_name', ]
br.settings._factor['brixs_object']  --->  ['attr_name', ]
br.settings._calib['brixs_object']  --->  ['attr_name', ]

for every entry in settings._shift, the attr_name.set_shift() will be called together with
brixs_object.set_shift(). Example.

let's say 

br.settings._shift['Spectrum'] = ['model', ]

then,

s1 = br.Spectrum()

s2 = s1.set_shift(10)  ---> this will also call s1.model.set_shift(10)

therefore, s2 will be a shifted copy of s1 where model is also shifted.

WARNING: For now, [_shift, _offset, _factor, _calib] is only implemented for 
Spectrum and Spectra. But it should be trivial to expand it to Image and 
PhotonEvents. I just did not do it yet, because I see no need.

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

FIGURE_ONCLICK (tuple): change the default `onclick` function for br.figure()
    figures. See br.figmanip._onclick() for an example.
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
        self.FIGURE_ONCLICK           = None 

        ############
        # internal #
        ############
        self._forbidden_words = {'Spectrum':[], 'Spectra':[], 'Image':[], 'PhotonEvents':[], 'Dummy':[]}
        self._reserved_words  = {'Spectrum':    {'methods':[], 'vars':[], 'pseudovars': []}, 
                                 'Spectra':     {'methods':[], 'vars':[], 'pseudovars': []}, 
                                 'Image':       {'methods':[], 'vars':[], 'pseudovars': []}, 
                                 'PhotonEvents':{'methods':[], 'vars':[], 'pseudovars': []},
                                 'Dummy':       {'methods':[], 'vars':[], 'pseudovars': []}}
        
        ###########################
        # brixs expansion modules #
        ###########################
        self._init   = {'Spectrum':{}, 'Spectra':{}, 'Image':{}, 'PhotonEvents':{}, 'Dummy':{}}
        self._copy   = {'Spectrum':[], 'Spectra':[], 'Image':[], 'PhotonEvents':[], 'Dummy':[]}
        self._shift  = {'Spectrum':[], 'Spectra':[], 'Image':[], 'PhotonEvents':[], 'Dummy':[]}
        self._offset = {'Spectrum':[], 'Spectra':[], 'Image':[], 'PhotonEvents':[], 'Dummy':[]}
        self._factor = {'Spectrum':[], 'Spectra':[], 'Image':[], 'PhotonEvents':[], 'Dummy':[]}
        self._calib  = {'Spectrum':[], 'Spectra':[], 'Image':[], 'PhotonEvents':[], 'Dummy':[]}

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