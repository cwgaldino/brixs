# settings (must be imported first)
from .config import settings

# backpack
from .backpack import *

# brixs
from .brixs import Image
from .brixs import PhotonEvents
from .brixs import Spectrum
from .brixs import Spectra
from .brixs import Dummy
from .brixs import figure
from .brixs import subplots
from .brixs import get_functions

# %% ===================================================================== %% #
# %% ============== Reserved and forbidden attribute names =============== %% #
# %% ===================================================================== %% #
# DEVELOPER NOTE: Any changes to this section must also be applied to its copy in config.py

# Reserved words are the names of all public methods, variables, and
# pseudovariables (@property attributes) defined by each brixs object.
# Each object has its own set of reserved words. These are used internally
# to build the corresponding list of forbidden words.

# Forbidden words are defined for each brixs object and include:
#   - all method names defined across all brixs objects;
#   - pseudovariables (@property attributes) defined by other objects, unless 
#     also defined by the current object;
#   - variables defined by other objects, unless also defined by the current
#     object.

# Brixs prevents users from creating attributes with forbidden names. This
# avoids accidentally overwriting built-in methods or attributes and, more
# importantly, prevents naming conflicts when one brixs object generates an
# object of another type (e.g. Image.calculate_histogram() returns a Spectrum).
# In such cases, user-defined attributes transferred to the new object must
# not conflict with methods or attributes of the destination object.

import types
import inspect

# %% =========== Collect reserved words for each brixs object ============ %% #
for cls in (Spectrum, Spectra, Image, PhotonEvents, Dummy):
    name  = cls.__name__
    _s    = cls()

    for attr in dir(_s):
        # ignore special attributes (doublescore methods and variables)
        if attr.startswith('__'):
            continue
        value = inspect.getattr_static(_s, attr)

        if isinstance(value, property):
            settings._reserved_words[name]['pseudovars'].append(attr)
        elif isinstance(value, (types.FunctionType, staticmethod, classmethod)):
            settings._reserved_words[name]['methods'].append(attr)
    settings._reserved_words[name]['vars'] = list(_s.__dict__)
del _s, name  

# %% ====== Build the list of forbidden words for each brixs object ====== %% #
for name in settings._forbidden_words:
    for name2 in settings._forbidden_words:

        # Methods from every brixs object are always forbidden
        settings._forbidden_words[name] += settings._reserved_words[name2]['methods']

        if name != name2:
            # Pseudovariables from other objects are forbidden unless they are
            # also defined by the current object
            settings._forbidden_words[name] += [
                pseudovar
                for pseudovar in settings._reserved_words[name2]['pseudovars']
                if pseudovar not in settings._reserved_words[name]['pseudovars']
            ]

            # Variables from other objects are forbidden unless they are also
            # defined by the current object
            settings._forbidden_words[name] += [
                var
                for var in settings._reserved_words[name2]['vars']
                if var not in settings._reserved_words[name]['vars']
            ]

    # Keep each forbidden-word list unique
    settings._forbidden_words[name] = [str(_) for _ in remove_duplicates(settings._forbidden_words[name])]

# %% ============================== addons =============================== %% #
from .addons import labels
from .addons import finder
# import .addons


