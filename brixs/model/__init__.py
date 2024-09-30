
# %% ------------------------- special import ----------------------------- %% #
# import lmfit
# import scipy

# %% -------------------------- brixs import ------------------------------ %% #
import brixs as br

# %% -------------------------- model import ------------------------------ %% #
from .model import Model
# %%

# %% ======================= add model to brixs =========================== %% #
br.settings._extra['Spectrum']['model'] = Model
br.settings._extra['Spectra']['model']  = Model

br.settings._modifiers['shift'].append('model')
br.settings._modifiers['offset'].append('model')
br.settings._modifiers['factor'].append('model')
br.settings._modifiers['calib'].append('model')

# %% ========================= model functions ============================ %% #
from .model_functions import *

# %% ====================== update forbidden words ======================== %% #
for obj in ('Spectrum', 'Spectra'):
    br.settings._reserved_words[obj]['pseudovars'].append('model')
for obj in ('Image', 'PhotonEvents'):
    br.settings._forbidden_words[obj].append('model')
