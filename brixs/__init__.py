# settings (must be imported first)
from .config import settings

# backpack
from backpack import *

# reserved words
# these words will raise an error if the user try to define attrs with these words
# this avoids errors when converting from one object to another (e.g. im -> s)
# from .brixs import Image as _Image
# from .brixs import PhotonEvents as _PhotonEvents
# from .brixs import Spectrum as _Spectrum
# from .brixs import Spectra as _Spectra
# settings._reserved_words  = remove_duplicates(list(_Spectrum().__dict__.keys()) + list(_Spectra().__dict__.keys()) + list(_Image().__dict__.keys()) + list(_PhotonEvents().__dict__.keys()))
# settings._reserved_words += [_[1:] for _ in settings._reserved_words if _.startswith('_')]
# del _Image, _PhotonEvents, _Spectrum, _Spectra
# print(settings._reserved_words)

# brixs
from .brixs import Image
from .brixs import PhotonEvents
from .brixs import Spectrum
from .brixs import Spectra
from .brixs import figure
from .brixs import subplots

# reserved and forbidden words
# these words will raise an error if the user try to define attrs with these words
# this avoids errors when converting from one object to another (e.g. im -> s)
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

# model (not sure if we should import this here)
from .model import __init__
from .model import peaks  # temporary (just for backward compatibility)
from .model import voigt_fwhm

# support
from .finder import *
#from .finder import _get_function_args_and_default_values

# beamlines
from beamlines import *

# extra
from .extra.crystal import *
from .extra import labels