# ==============================================================================
# %% EXAMPLE: Spectrum ====================================== 26/08/2022 =======
# ==============================================================================

# %% imports ===================================================================
import brixs as br
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

%matplotlib qt5
%load_ext autoreload
%autoreload 2

# %% Create object ============++++++===========================================
v = np.array([[0, 0], [1, 1], [2, 2], 3, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]])


v = np.array([[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]])
v
s = br.Spectrum(data=[[0, 1, 2, 3, 4, 5, 6, 7, 8, 9], [0, 1, 2, 3, 4, 4, 3, 2, 1, 0]])
# im = br.Image('<path-to-Image>')  # must be a txt with the pixel matrix


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

# %% modifiers =================================================================
s.offset
s.calib
s.factor
s.shift
s.shift_roll
s.shift_interp

# offset
s.set_offset(20)
s.plot()

# modifiers are absolute quantities
s.set_offset(0)
s.plot()

# modifiers can be set via the attribute
s.offset = 30
s.plot()

# for relative modifications use assignment operators
s.offset -= 30
s.plot()

# %% hard shift modifier =======================================================
s = br.read_ADRESS(filepath)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()

# hard ('hard', 'x') shift
# note how the whole spectra moves to the right
# y shape is fully preserved
# x axis is not preserved
s.set_shift(500, mode='hard')
s.plot()

# %% soft shift modifier =======================================================
s = br.read_ADRESS(filepath)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()

# soft ('soft', 'y') shift
# y shape is interpolated in a different position
# x axis is fully preserved
# in this plot, data note that the left edge of the data is meaningless
s.set_shift(500, mode='soft')
s.plot()

# %% roll shift modifier =======================================================
s = br.read_ADRESS(filepath)

# plot
fig = br.backpack.figure()
br.backpack.set_window_position(2048, 232)
s.plot()

# roll ('roll', 'r') shift
# y shape is fully preserved
# x axis is fully preserved
# the data points that fall outside the x range are put on the other side of the
# spectrum. Therefore, the left edge of this plot is meaningless.
s.set_shift(500, mode='roll')
s.plot()
