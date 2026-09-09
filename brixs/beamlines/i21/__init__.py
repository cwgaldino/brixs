from pathlib import Path
from .support import _read_csv
from .core import *
from .advanced import *
from .advanced import _process
from .cosmic_rays import *
from .cosmic_rays import _detect_outliers_via_row_median

values = _read_csv(filepath=Path(__file__).parent/'attrs/attrsrixs.csv', delimiter=',', comments='#', strip=True)
settings.METADATA['rixs'] = values
values = _read_csv(filepath=Path(__file__).parent/'attrs/attrsxas.csv', delimiter=',', comments='#', strip=True)
settings.METADATA['xas'] = values
# values = _read_csv(filepath=Path(__file__).parent/'attrs/attrslinescans.csv', delimiter=',', comments='#', strip=True)
# settings.METADATA['linescans'] = values

