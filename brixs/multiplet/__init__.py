from pathlib import Path
import collections
import gzip
import json

class odict(collections.OrderedDict):
    def __missing__(self, key):
        value = self[key] = type(self)()
        return value

from .crispy_0_7_3_modified import settings
settings.TEMPLATES_FOLDERPATH = Path(__file__).parent/'templates'
settings.PARAMETERS_FILEPATH  = Path(__file__).parent/'parameters.json.gz'


# make parameters a class attribute. This speeds up the creation
# of a new calculation object; significantly.
with gzip.open(settings.PARAMETERS_FILEPATH, 'r') as f:
    settings._PARAMETERS = json.loads(f.read().decode('utf-8'), object_pairs_hook=odict)
settings._ELEMENTS = list(settings.PARAMETERS['elements'].keys())

# 
from .crispy_0_7_3_modified import Calculation, load_calculation, quanty, remove_greek_letters_from_hamiltonianData 
