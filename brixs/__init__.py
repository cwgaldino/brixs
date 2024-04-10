# settings (must be imported first)
from .config import settings

# backpack
from backpack import *

# brixs
from .brixs import Image
from .brixs import PhotonEvents
from .brixs import Spectrum
from .brixs import Spectra
from .brixs import figure
from .brixs import subplots

# reserved words
# these words will raise an error if the user try to define attrs with these words
# this avoids errors when converting from one object to another (e.g. im -> s)
settings._reserved_words  = remove_duplicates(list(Spectrum().__dict__.keys()) + list(Spectra().__dict__.keys()) + list(Image().__dict__.keys()) + list(PhotonEvents().__dict__.keys()))
settings._reserved_words += [_[1:] for _ in settings._reserved_words if _.startswith('_')]

# model (not sure if we should import this here)
from .model import __init__
from .model import peaks  # temporary (just for backward compatibility)

# support
from .finder import *

# beamlines
from beamlines import *

# extra
from .extra.crystal import *
from .extra import labels