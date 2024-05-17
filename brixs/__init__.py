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

# reserved and forbidden words
# these words will raise an error if the user try to define attrs with these words
# this avoids errors when converting from one object to another (e.g. im -> s)
# this also avoid that import methods or variables are overwritten
# reserved words:  methods, vars, and pseudovars for each object
# forbidden words: all methods from all objects + pseudovar from all other objects
import types
for object in (Spectrum, Spectra, Image, PhotonEvents):
    name  = object.__name__
    _s    = object()
    attrs = _s.__dir__()
    for attr in attrs:
        if attr.startswith('_'):
            pass
        elif type(getattr(_s, attr)) == types.MethodType:
            settings._reserved_words[name]['methods'] += [attr]
        else:
            settings._reserved_words[name]['pseudovars'] += [attr]
    settings._reserved_words[name]['vars'] = list(_s.__dict__)
del _s, attrs, name    
for name in settings._forbidden_words:
    for name2 in settings._forbidden_words:
        settings._forbidden_words[name] += settings._reserved_words[name2]['methods']
        if name != name2:
            settings._forbidden_words[name] += [pseudovar for pseudovar in settings._reserved_words[name2]['pseudovars'] if pseudovar not in settings._reserved_words[name]['pseudovars']]
    settings._forbidden_words[name] = remove_duplicates(settings._forbidden_words[name])

# addons
from .addons import labels
from .addons import finder
# import .addons
