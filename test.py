import os
os.chdir('/home/galdino/github/brixs/')
import brixs
import matplotlib.pyplot as plt
import numpy as np
from brixs import arraymanip as am
from brixs.arraymanip import index
from brixs import filemanip as fm
from brixs import figmanip as figm
plt.ion()
# simulating a generic spectrum
I = brixs.dummy_spectrum(0, 0.2, excitations=[[0.5, 2, 2], [0.5, 4, 2]])
# simulating the photon_event list(where we're using energy in eV's and length in meters)
photon_events = brixs.dummy_photon_events(I, background=0.02,
                                                noise=0.05,
                                                exposure=50e4,
                                                dispersion= 8.45 * (10**-3 / 10**-6),
                                                x_max=52.22e-3,
                                                y_max=25.73e-3,
                                                y_zero_energy=-20,
                                                angle=2,
                                                psf_fwhm=(3e-6, 1e-6))

p = brixs.photon_events(data=photon_events)

p.apply_correction(lambda x, y: (x*10**3, y*10**3))
p.set_binning((10, 1000))
p.calculate_offsets(ranges=[[0.5666516987271235,  5.161124931649747]])
p.fit_offsets()
p.plot(show_offsets=True, show_offsets_fit=True)
p.plot_columns(columns=(0, 2, 5), vertical_increment=10, show_ranges=True)

figm.set_onclick(format='png', resolution=300, folder='/home/galdino/github/brixs/docs/source/_figs')

import os
os.chdir('/home/galdino/github/brixs/')
import brixs
import matplotlib.pyplot as plt
import numpy as np
from brixs import figmanip as figm
figm.set_onclick(format='png', resolution=300, folder='/home/galdino/github/brixs/docs/source/_figs')

plt.ion()



import brixs
import matplotlib.pyplot as plt
import numpy as np
# simulating a generic spectrum
x = np.linspace(-5, 10, 400)
data = []
for i in range(12):
    I = brixs.dummy_spectrum(0+i*0.1, 0.2, excitations=[[0.5, 2, 2], [0.5, 4, 2]])
    noise = np.random.default_rng().normal(-0.05, 0.05, size=400)
    data.append(brixs.spectrum(data=np.vstack((x, I(x)+noise)).transpose()))
# data[-1].plot()
s = brixs.spectra(data=data)
s.calculate_shifts(ref=0, mode='cross-correlation', ranges=None, check_x=True, verbose=True)
s.plot(show_ranges=True)


ax = data[0].plot()
data[0].data[:, 1]
step = data[0].x[1] - data[0].x[0]
data[0].y
data[0].apply_shift(10, mode='roll')
data[0].plot(ax)
data[0].apply_shift(-10, mode='roll')
# data[0].plot(ax)
data[0].apply_shift(10*step, mode='hard')
data[0].plot(ax)
data[0].apply_shift(10*step, mode='soft')
data[0].plot(ax)


s.calculate_shifts(ref=0, mode='fit', ranges=[[-2, 2]], check_x=True, verbose=True)

d = np.diff(s.spectrum[0].data[:, 0])
(max(d) - min(d))*100/np.mean(d))


# simulating the photon_event list
data = brixs.dummy_photon_events(I, noise=0.02, background=0.01, y_zero_energy=-20)
# initializing photon_events object
p = brixs.photon_events(data=data)
p.

# set binning
p.set_binning((10, 50))
p.plot(show_bins=True)
print(p.hist)
p.plot_columns(columns='all', vertical_increment=100)
# fitting offsets
p.set_binning((10, 1000))
p.calculate_offsets(ranges=[[0,  0.005]])
p.fit_offsets()
p.plot(show_offsets=True, show_offsets_fit=True)
# remove offsets
p.offsets_correction()
p.plot()
