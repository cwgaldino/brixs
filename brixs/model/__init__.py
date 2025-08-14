
# %% ------------------------- special import ----------------------------- %% #
# import lmfit
# import scipy

# %% -------------------------- brixs import ------------------------------ %% #
import brixs as br

# %% -------------------------- model import ------------------------------ %% #
from .model import Model
# %%

# %% ======================= add model to brixs =========================== %% #
br.settings._init['Spectrum']['model'] = Model
br.settings._init['Spectra']['model']  = Model

for object in ('Spectrum', 'Spectra'):
    br.settings._copy[object].append('model')
    br.settings._shift[object].append('model')
    br.settings._offset[object].append('model')
    br.settings._factor[object].append('model')
    br.settings._calib[object].append('model')

# %% ========================= model functions ============================ %% #
from .model_functions import *

# %% ====================== update forbidden words ======================== %% #
for obj in ('Spectrum', 'Spectra'):
    br.settings._reserved_words[obj]['pseudovars'].append('model')
for obj in ('Image', 'PhotonEvents'):
    br.settings._forbidden_words[obj].append('model')
