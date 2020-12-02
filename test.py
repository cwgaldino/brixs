import os
os.chdir('/home/galdino/github/brixs/')

# brixs
import brixs
from brixs import figmanip as figm
figm.set_onclick(format='png', resolution=300, folder='/home/galdino/github/brixs/docs/source/_figs')

# standard libraries
import matplotlib.pyplot as plt
plt.ion()

# simulating a generic spectrum
I = brixs.dummy_spectrum(0, 0.2, excitations=[[0.5, 2, 2], [0.5, 4, 2]])
# simulating the photon_event list(where we're using energy in eV's and length in meters)
data = brixs.dummy_photon_events(I, background=0.02,
                                    noise=0.05,
                                    exposure=50e4,
                                    dispersion= 8.45 * (10**-3 / 10**-6),
                                    x_max=52.22e-3,
                                    y_max=25.73e-3,
                                    y_zero_energy=-20,
                                    angle=2,
                                    psf_fwhm=(3e-6, 1e-6))

# initializing photon_events object
p = brixs.photon_events(data=data)
# change leght unit from m to mm
p.apply_correction(lambda x, y: (x*10**3, y*10**3))

# ploting data
ax = p.plot()

# set binning and calculating offsets
p.set_binning((10, 1000))
p.calculate_offsets(ranges=[[0.5666516987271235,  5.161124931649747]])
p.fit_offsets()
ax = p.plot(show_bins=(True, False), show_offsets=True, show_offsets_fit=True)

# plot all columns
ax = p.plot_columns(vertical_increment=10, show_ranges=True)
figm.zoom(0, 7)
ax.set_xlabel('y position (mm)')
ax.set_ylabel('intensity')

# apply offset correction
p.offsets_correction()

# plot all columns after correction
ax = p.plot_columns(vertical_increment=10, show_ranges=True)
figm.zoom(0, 7)
ax.set_xlabel('y position (mm)')
ax.set_ylabel('intensity')

# calculate spectrum
p.calculate_spectrum(y_bins=2000)

# plot final spectrum
ax = p.spectrum.plot()
figm.zoom(2.5, 4.5)
ax.set_xlabel('y position (mm)')
ax.set_ylabel('intensity')





class car():

    def __init__(self):
        self.color = 'green'

    def change_color(self):
        self.color = 'red'

car1 = car()
car2 = car()
car_list = [car1, car2]
class cars():

    def __init__(self, car_list):
        self.data = car_list

    def change_color(self):
        for car in self.data:
            car.change_color()

many_cars = cars(car_list)

many_cars.data[0].color
many_cars.data[1].color
many_cars.change_color()
many_cars.data[0].color
many_cars.data[1].color


import os
os.chdir('/home/galdino/github/brixs/')

# brixs
import brixs
from brixs import figmanip as figm
figm.set_onclick(format='png', resolution=300, folder='/home/galdino/github/brixs/docs/source/_figs')

# standard libraries
import numpy as np
import matplotlib.pyplot as plt
plt.ion()


# creating a dummy list of spectra with a misalignment between them
x = np.linspace(-5, 10, 400) # energy (eV)
data = []
for i in range(12):
    I = brixs.dummy_spectrum(0+i*0.1, 0.2, excitations=[[0.5, 2, 2], [0.5, 4, 2]])
    spectrum = brixs.spectrum(data=np.column_stack((x, I(x))))
    data.append(spectrum)

# data is list of 12 spectrum objects
# print(data)
# [<brixs.brixs.spectrum at 0x7f2fbc78e390>,
#  <brixs.brixs.spectrum at 0x7f2fbc78ee80>,
#  <brixs.brixs.spectrum at 0x7f2fbc78e2e8>,
#  <brixs.brixs.spectrum at 0x7f2fbc78e048>,
#  <brixs.brixs.spectrum at 0x7f2fbc78e6d8>,
#  <brixs.brixs.spectrum at 0x7f2fbc78ef28>,
#  <brixs.brixs.spectrum at 0x7f2fbc78e438>,
#  <brixs.brixs.spectrum at 0x7f2fbc78e160>,
#  <brixs.brixs.spectrum at 0x7f2fbc78e240>,
#  <brixs.brixs.spectrum at 0x7f2fbc78e0b8>,
#  <brixs.brixs.spectrum at 0x7f2fbc78e1d0>,
#  <brixs.brixs.spectrum at 0x7f2fbc78e518>]

# plot
ax = data[0].plot()
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Intensity')

# add random noise to the data
for i in range(len(data)):
    noise = np.random.normal(-0.05, 0.05, size=len(x))
    f = lambda x, y: (x, y+noise)
    data[i].apply_correction(f)

# plot
ax = data[0].plot()
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Intensity')

# initialize spectra object
s = brixs.spectra(data=data)
ax = s.plot(vertical_increment=0.5, show_ranges=True)
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Intensity')

# calculate shifts
s.calculate_shifts(ref=0, mode='cross-correlation', ranges=[[-1, 6]])
s.shifts_correction()

# plot after correction
ax = s.plot(vertical_increment=0.5, show_ranges=True)

# calculate sum
s.calculate_sum()
ax = s.sum.plot()
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Intensity')
