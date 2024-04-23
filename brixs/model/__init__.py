
# %% ------------------------- special import ----------------------------- %% #
# import lmfit
# import scipy

# %% -------------------------- brixs import ------------------------------ %% #
import brixs as br

# %% -------------------------- model import ------------------------------ %% #
from .model import Model
# %%

# %% ====================== add model to objects ========================== %% #
br.Spectrum.model = Model()
br.Spectra.model  = Model()

# %% ========================= model functions ============================ %% #
from .model_functions import *

# %% ====================== update forbidden words ======================== %% #
for obj in ('Spectrum', 'Spectra'):
    br.settings._reserved_words[obj]['pseudovars'].append('model')
for obj in ('Spectrum', 'Spectra', 'Image', 'PhotonEvents'):
    br.settings._forbidden_words[obj].append('model')