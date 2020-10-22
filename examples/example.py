#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simulate image from RIXS detector.

TODO:
Remove outliers from offsets.
Evaluate necessity for calculating curvature again
"""

# standard libraries
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import copy

# backpack
import sys
sys.path.append('/home/galdino/github/py-backpack')
import backpack.filemanip as fmanip
import backpack.figmanip as figmanip
import backpack.arraymanip as manip
from backpack.arraymanip import index
import importlib
importlib.reload(fmanip)
importlib.reload(manip)
importlib.reload(figmanip)

# rixs_utils
sys.path.append('../')
import rixs_utils as rixs
rixs = importlib.reload(rixs)

%matplotlib qt5
# figmanip.set_default_window_position((351, 1137))
figmanip.set_default_window_position((1212, 30))
# figmanip.getWindowPosition()

# %% initialize folders ========================================================
def create_folders(root):
    root = Path(root)
    root.mkdir(exist_ok=True)

    # remove
    try:
        fmanip.rmdir(root/'photon_events')
        fmanip.rmdir(root/'photon_events_curvature')
        fmanip.rmdir(root/'spectra')
        fmanip.rmdir(root/'spectra_aligned')
        fmanip.rmdir(root/'spectra_final')
    except FileNotFoundError:
        pass

    Path(root/'photon_events').mkdir(exist_ok=True)
    Path(root/'photon_events_curvature').mkdir(exist_ok=True)
    Path(root/'spectra').mkdir(exist_ok=True)
    Path(root/'spectra_aligned').mkdir(exist_ok=True)
    Path(root/'spectra_final').mkdir(exist_ok=True)

    return root

# %% Spectrum =================================================================
rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)

c   = 1           # eV
w   = 20 * 10**-3 # eV
#             rel area, rel center, rel width],
excitations = [[0.2,       1.5,        2.5],
               [0.3,         2,          4],
               [0.3,         4,          4],
              ]

I = rixs.simulate_spectrum(c, w, excitations)
energy = np.linspace(-5, 10, 10000)  # eV

spectrum_simulated = rixs.spectrum(data=np.vstack([energy, I(energy)]).T)

# %% plot =============================
# rixs = importlib.reload(rixs)
# fmanip = importlib.reload(fmanip)
#
# spectrum_simulated.plot(marker=None)










# ==============================================================================
# %% ================================== TEST 0 =================================
# ==============================================================================
temp = fmanip.parsed_filelist(string='data')
number = int(list(temp.keys())[-1] +1)

# root = Path('data_27')
root = create_folders(root=f'data_{number}')

rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)

output_folderpath = Path(root/'photon_events')
prefix   = 'data'
detector = dict(dispersion = 8.45 * (10**-3 / 10**-6),  # eV/m
                x_max      = 52.22e-3, # m
                y_max      = 25.73e-3, # m
                angle      = 0         # deg
               )
background = 0
noise      = 0
x_spread   = 0#5e-6
y_spread   = 0.85e-6

exposure = 600e4
n_images = 10

j = 0
for i in range(0, n_images):
    y_zero         = -50 #+ np.random.uniform(-1, 1) # eV
    simulated_data = rixs.simulate_photon_events(I, background=background,
                                                    noise=noise,
                                                    exposure=exposure,
                                                    y_zero=y_zero,
                                                    **detector)
    photon_events  = rixs.photon_events(data=simulated_data, x_max=detector['x_max'], y_max=detector['y_max'])
    photon_events.apply_spread(x_spread=x_spread, y_spread=y_spread)
    photon_events.save(filepath=output_folderpath/f'{prefix}_{j}.dat')
    j += 1

# %% plot =============================
# input_folderpath = Path(root/'photon_events')
# filelist = fmanip.parsed_filelist(dirpath=input_folderpath, ref=-1)
#
# file = 0
# photon_events  = rixs.photon_events(filepath=filelist[file])
# photon_events.plot(pointsize=3)






# %% adjust binning for curvature ==============================================
# rixs = importlib.reload(rixs)
# fmanip = importlib.reload(fmanip)
#
# input_folderpath = Path(root/'photon_events')
# filelist = fmanip.parsed_filelist(dirpath=input_folderpath, ref=-1)
#
# file=0
# # first iteration
# binx_1 = 4
# biny_1 = 40
# photon_events = rixs.photon_events(filepath=filelist[file], x_max=detector['x_max'], y_max=detector['y_max'])
# photon_events.calculate_curvature(ref=0, type='max', poly_order=1, binx=binx_1, biny=biny_1, y_min=2e-3, y_max=10e-3)
#
# ax = photon_events.plot(pointsize=3, show_binx=True, show_biny=True, show_curvature=True, show_lim=True)
# plt.title("first iteraction")
#
# ax = photon_events.plot_offsets()
# photon_events.plot_curvature(ax=ax)
# plt.title("first iteraction")
#
# photon_events.plot_columns()
# plt.title("first iteraction")
#
# # %%
# photon_events.apply_curvature()
# photon_events.plot_columns()
# plt.title("first iteraction: after fix")
#
# # %% second iteration
# binx_2 = 6
# biny_2 = 500
# photon_events.calculate_curvature(ref=0, type='max', poly_order=1, binx=binx_2, biny=biny_2, y_min=4e-3, y_max=9e-3)
#
# ax = photon_events.plot(pointsize=3, show_binx=True, show_curvature=True, show_lim=True)
# plt.title("second iteraction")
#
# ax = photon_events.plot_offsets()
# photon_events.plot_curvature(ax=ax)
# plt.title("second iteraction")
#
# photon_events.plot_columns()
# plt.title("second iteraction")
#
# # %%
# photon_events.apply_curvature()
# photon_events.plot_columns()
# plt.title("second iteraction: after fix")
#
# # %% third iteration
# binx_3 = 10
# biny_3 = 20000
# elastic_center = photon_events.y_centers[np.argmax(photon_events.columns[0])]
# photon_events.calculate_curvature(ref=0, poly_order=1, binx=binx_3, biny=biny_3, y_min=elastic_center-0.08e-3, y_max=elastic_center+0.08e-3)
#
# ax = photon_events.plot(pointsize=3, show_binx=True, show_curvature=True, show_lim=True)
# plt.title("third iteraction")
#
# ax = photon_events.plot_offsets()
# photon_events.plot_curvature(ax=ax)
# plt.title("third iteraction")
#
# photon_events.plot_columns()
# plt.title("third iteraction")
#
# # %%
# photon_events.apply_curvature()
# photon_events.plot_columns()
# plt.title("third iteraction: after fix")
#
# # %% final verification
# photon_events.calculate_curvature(ref=0, poly_order=1, binx=6, biny=8000)
#
# ax = photon_events.plot(pointsize=3, show_curvature=True)
# plt.title("final")
#
# ax = photon_events.plot_offsets()
# photon_events.plot_curvature(ax=ax)
# plt.title("final")
#
# photon_events.plot_columns()
# plt.title("final")




# %% load photon_events and fix curvature ======================================
rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)

input_folderpath = Path(root/'photon_events')
output_folderpath = Path(root/'photon_events_curvature')
prefix = 'data'
filelist = fmanip.parsed_filelist(dirpath=input_folderpath, ref=-1)

for file in filelist:
    photon_events = rixs.photon_events(filepath=filelist[file], x_max=detector['x_max'], y_max=detector['y_max'])
    # photon_events.calculate_curvature(ref=0, type='max', poly_order=1, binx=binx_1, biny=biny_1, y_min=2e-3, y_max=10e-3)
    # photon_events.apply_curvature()
    # photon_events.calculate_curvature(ref=0, type='max', poly_order=1, binx=binx_2, biny=biny_2, y_min=4e-3, y_max=9e-3)
    # photon_events.apply_curvature()
    # elastic_center = photon_events.y_centers[np.argmax(photon_events.columns[0])]
    # photon_events.calculate_curvature(ref=0, poly_order=1, binx=binx_3, biny=biny_3, y_min=elastic_center-0.08e-3, y_max=elastic_center+0.08e-3)
    # photon_events.apply_curvature()
    photon_events.save(filepath=output_folderpath/f'{prefix}_{file}.dat')

# %% plot =============================
# input_folderpath = Path(root/'photon_events_curvature')
# filelist = fmanip.parsed_filelist(dirpath=input_folderpath, ref=-1)
#
# temp = list(filelist.keys())
# temp.sort()
# file_iter = iter(temp)
# # %%
# file = next(file_iter)
# # file = 8
# photon_events  = rixs.photon_events(filepath=filelist[file], x_max=detector['x_max'], y_max=detector['y_max'])
# photon_events.plot()
# plt.ylim(6, 9)
#
# # %%
# plt.close('all')








# %% calculate spectrum =============================================================
rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)

input_folderpath = Path(root/'photon_events_curvature')
output_folderpath = Path(root/'spectra')
prefix = 'spectrum'
filelist = fmanip.parsed_filelist(dirpath=input_folderpath, ref=-1)

biny = 90000  # Smaller biny decreases resolution. Large biny increases noise
biny = 90000  # Smaller biny decreases resolution. Large biny increases noise

for file in filelist:
    photon_events  = rixs.photon_events(filepath=filelist[file], x_max=detector['x_max'], y_max=detector['y_max'])
    photon_events.calculate_spectrum(biny=biny)
    photon_events.spectrum.save(filepath=output_folderpath/f'{prefix}_{file}.dat')

# %% plot =============================
rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)

input_folderpath = Path(root/'spectra')
spectra = rixs.spectra(input_folderpath)
spectra.plot(vertical_increment=0, marker='x')
figmanip.setWindowPosition()
figmanip.zoom(6.6e-3, 7.4e-3)


# %% calculate and apply shifts =============================================================
rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)
manip = importlib.reload(manip)

input_folderpath = Path(root/'spectra')
output_folderpath = Path(root/'spectra_aligned')
filelist = fmanip.parsed_filelist(dirpath=input_folderpath, ref=-1)

spectra = rixs.spectra()
for file in filelist:
    spectra.append(filepath=filelist[file])

# spectra.calculate_shifts()
# spectra.apply_shifts()
# spectra.crop()
#
# spectra.calculate_shifts()
# spectra.apply_shifts()
# spectra.interpolate(num=len(spectra.spectrum[0].data[:, 0]))
#
# spectra.calculate_shifts(type='fit', start=5.8e-3, stop=6.2e-3)
# spectra.apply_shifts()
# spectra.interpolate(num=len(spectra.spectrum[0].data[:, 0]))


spectra.save(output_folderpath)

# %% plot =============================
rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)

input_folderpath = Path(root/'spectra_aligned')
spectra = rixs.spectra(input_folderpath)
spectra.plot(factor=1e3)
figmanip.setWindowPosition()


# %% sum spectra ===============================================================
rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)

input_folderpath = Path(root/'spectra_aligned')
output_folderpath = Path(root/'spectra_final')

spectra = rixs.spectra(input_folderpath)
spectra.calculate_sum()

spectra.sum.save(output_folderpath/'final.dat')

# %% plot =============================
# rixs = importlib.reload(rixs)
# fmanip = importlib.reload(fmanip)
#
# input_folderpath = Path(root/'spectra_final')
# spectrum = rixs.spectrum(filepath=input_folderpath/'final.dat')
# ax = spectrum.plot(normalized=False)


# %% calibrate spectra ===============================================================
rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)

input_folderpath = Path(root/'spectra_final')
output_folderpath = Path(root/'spectra_final')

spectrum = rixs.spectrum(filepath=input_folderpath/'final.dat')
spectrum_calibrated = spectrum.calibrate(dispersion=detector['dispersion'])

spectrum_calibrated.save(output_folderpath/'final_calibrated.dat')

# %% plot =============================
# rixs = importlib.reload(rixs)
# fmanip = importlib.reload(fmanip)
#
# input_folderpath = Path(root/'spectra_final')
# spectrum = rixs.spectrum(filepath=input_folderpath/'final_calibrated.dat')
# shift = c - spectrum.data[np.argmax(spectrum.data[:, 1]), 0]
# ax = spectrum.plot(normalized=True, shift=shift)
#
# spectrum_simulated.plot(ax, marker=None)
# figmanip.zoom(0, 6)






# ==============================================================================
# %% ============================= Result 0, 1 =================================
# ==============================================================================
rixs = importlib.reload(rixs)
fmanip = importlib.reload(fmanip)

input_folderpath = Path(root/'spectra_final')
spectrum = rixs.spectrum(filepath=input_folderpath/'final_calibrated.dat')
shift = c - spectrum.data[np.argmax(spectrum.data[:, 1]), 0]
spectrum_shifted = spectrum.shift(shift)

ax = spectrum_shifted.plot(normalized=True)#, shift=shift)
fit, popt = spectrum_shifted.peak_fit(0.9, 1.1)
print(popt[2])

spectrum_simulated.plot(ax, marker=None)
plt.plot(fit[:, 0], fit[:, 1])
figmanip.setWindowPosition()
# figmanip.zoom(0.96, 1.05)

# save width to file
f = open(str(root/'spectra_final/result.txt'), 'a')
f.write(str(popt[2])+'\n')
f.close()
