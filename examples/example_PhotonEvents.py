# ==============================================================================
# %% EXAMPLE: photon events ================================= 22/02/2022 =======
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
folderpath = Path(r'../fixtures/ADRESS')

# %% import photonn events =====================================================
pe  = br.read_ADRESS_pe(folderpath/'Cu_0005_d1.h5')

# %% PhotonEvents basic attributes =============================================
pe.shape  # size of the detector
pe.data
pe.x
pe.y
pe.I
len(pe)  # number of photon events

# %% detector data (exclusive of ADRESS data) ==================================
pe.nd

# %% setting new parameters ====================================================
pe.temperature = 10
print(pe.temperature)

# %% plot ======================================================================
plt.figure()
pe.plot()

plt.figure()
plt.scatter(pe.x, pe.y, s=1)

# %% binning ===================================================================
pe.binning(15, 10)

plt.figure()
pe.reduced.plot(colorbar=True)

plt.figure()
pe.reduced.columns.plot()

plt.figure()
pe.reduced.rows.plot()

plt.figure()
pe.plot()
plt.vlines(pe.reduced.x+pe.bins_size[1]/2, 0, pe.shape[0], color='red')
plt.hlines(pe.reduced.y+pe.bins_size[0]/2, 0, pe.shape[1], color='red')

# %% calculate and fit shifts ==================================================
pe.binning(3000, 10)
pe.calculate_shifts()
p, f, s_fit = pe.calculated_shifts.polyfit(deg=2)

plt.figure()
pe.calculated_shifts.plot(marker='o', lw=0, color='black')
s_fit.plot(color='red')

plt.figure()
pe.plot()
s_fit.plot(offset=730, factor=-1, color='red')

# %% set shifts ================================================================
pe.set_shifts(p=p)

plt.figure()
pe.plot()

# %% spectrum ==================================================================
s = pe.calculate_spectrum(6000)

plt.figure()
s.plot()

# %% fix curvature =============================================================
pe  = br.read_ADRESS_pe(folderpath/'Cu_0005_d1.h5')
pe.binning(3000, 10)
pe.fix_curvature()

plt.figure()
pe.plot()

pe.p  # this can be used to apply the same correction in other files
pe.f
pe.shifts

# %% save and load data ========================================================
pe.save(r'test.txt')
pe.temperature = 10

pe2 = br.PhotonEvents()
pe2.load(r'test.txt')
pe2.temperature

# %%

























# calculate calibration factor =================
disp, sss = br.calculate_calib_ADRESS(adress_example, 'Cu', 5, 15)
print(np.mean(disp))

for ccd in range(3):
    plt.figure()
    plt.title(f'ccd {ccd}: calib={round(disp[ccd], 7)}')
    for s in sss[ccd]:
        s.plot(label=s.scan)
        s.fit.plot(color='black')
    plt.xlabel('Energy loss (eV)')
    plt.ylabel('Intensity (arb. units)')
    plt.legend()

plt.figure()
for ss in sss:
    x = ss.calib_calculated['values']
    y = ss.calib_calculated['centers']
    y_fit = ss.calib_calculated['func'](x)
    plt.scatter(x, y, label=f'ccd {ss.ccd}')
    plt.plot(x, y_fit, color='black')
plt.xlabel('Photon energy (eV)')
plt.ylabel('Position of the elastic line (subpixel)')
plt.legend()

plt.figure()
for ss in sss:
    x = ss.calib_calculated['values']
    y = ss.calib_calculated['fwhms']*1000
    plt.scatter(x, y, label=f'ccd {ss.ccd}')
plt.xlabel('Photon energy (eV)')
plt.ylabel('FWHM of the elastic line (meV)')
plt.legend()

# %% data reduction ==========================
def data_reduction(self, calib):
    # align
    self.calculate_shifts()
    self.set_shift()
    s = self.calculate_sum()

    # calibrate
    s.calib = calib

    # shift to zero
    s.find_peaks()
    s.fit_peaks()
    s.shift = -s.peaks[-1]['c']

    return s

# attaching it to spectra
br.Spectra.data_reduction = data_reduction

# loading data
ss = br.read_ADRESS(adress_example2, 'Cu', 40)

# get final spectra
s = ss.data_reduction(calib=0.012796)

# plot
plt.figure()
s.plot()
plt.xlabel('Energy loss (eV)')
plt.ylabel('Intensity (arb. units)')

# %% map ========================================
scanlist = [36, 37, 38, 39, 40]

sss = []
ss  = br.Spectra()
for scan in scanlist:
    ss_temp = br.read_ADRESS(adress_example2, 'Cu', scan)
    sss.append(ss_temp)

    s = ss_temp.data_reduction(calib=0.012796)
    ss.append(s)

# plot
plt.figure()
ss.plot()
plt.xlabel('Energy loss (eV)')
plt.ylabel('Intensity (arb. units)')

# plot each ccd data for each scan for verification
for ss_temp in sss:
    plt.figure()
    plt.title(f'scan: {ss.scan}')
    for s_temp in ss_temp:
        s_temp.plot(label=f'ccd {s.ccd}')
    plt.legend()
    plt.xlabel('Energy loss (eV)')
    plt.ylabel('Intensity (arb. units)')
plt.close('all')

# plot map
ss.interp()
ss.plot_map()
