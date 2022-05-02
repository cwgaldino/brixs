# ==============================================================================
# %% BRIXS quickstart ====================================== 22/02/2022 ========
# ==============================================================================

# %% ADRESS beamline ===========================================================
import brixs as br
from pathlib import Path
import matplotlib.pyplot as plt

# initial definitions ==========================
adress_example1 = Path('<path-to-brixs>/fixtures/ADRESS/calib')
adress_example2 = Path('<path-to-brixs>/fixtures/ADRESS/map')


# reading data ==================================
s  = br.read_ADRESS(adress_example/'Cu_0005_d1.h5')
ss = br.read_ADRESS(adress_example, 'Cu', 5)

s  = br.read_ADRESS(filepath=adress_example/'Cu_0005_d1.h5')
ss = br.read_ADRESS(folderpath=adress_example, prefix='Cu', n=5)

pe  = br.read_pe_ADRESS(adress_example/'Cu_0005_d1.h5')
pes = br.read_pe_ADRESS(adress_example, 'Cu', 5)

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
