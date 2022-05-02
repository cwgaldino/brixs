# ==============================================================================
# %% EXAMPLE: Spectrum ====================================== 26/04/2022 =======
# ==============================================================================

# %% imports ===================================================================
import brixs as br
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

%matplotlib qt5
%load_ext autoreload
%autoreload 2

# %% initial definitions =======================================================
filepath = Path(r'../fixtures/ADRESS/Cu_0005_d1.h5')

# %% import photonn events =====================================================
s = br.read_ADRESS(filepath)
print(s.x)
print(s.y)

# %% plot  =====================================================================
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()


# %% other methods =============================================================

interp
extract
crop
floor
flip
