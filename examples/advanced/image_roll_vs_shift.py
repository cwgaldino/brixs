"""Quick example to illustrate the difference between Image.set_vertical_roll 
and Image.set_vertical_shift.

- Roll does not care about the x and y centers. It just rolls the data in the image.

- Shift will roll the data in relation to the x and y centers. 
"""

# %%
import matplotlib.pyplot as plt
plt.ion()
import brixs as br
import numpy as np


# % generic example data
x = np.linspace(1, 6, 100)

ss = br.Spectra()
ss.values = []
for i in range(14, 30):
    ss.append(br.Spectrum(x=x, y=np.sin(i*0.1*x)))
    ss.values.append(i)

br.figure()
_ = ss.plot()

# %% Create intensity map with vertical roll = 10
# Note that rolled out pixels are re-introduced on the other side of the image
# It is up to you to decide if rolled out pixels make sense in your analysis.
# In same cases, the image needs to be cropped to remove rolled out pixels.
br.figure()
ss.create_intensity_map(ss.values).set_vertical_roll(10).plot()

## %% Create intensity map with vertical roll = 90
br.figure()
ss.create_intensity_map(ss.values).set_vertical_roll(90).plot()

# %% Create intensity map with vertical roll = 110
# note that the roll is larger than the number of rows (which is the number of spectra)
# This roll will be equivalent to roll = 10
br.figure()
ss.create_intensity_map(ss.values).set_vertical_roll(110).plot()

# %% Create intensity map with vertical shift = 10
# Note that the y axis (y_centers) are shifted by 10 units, but the data is not 
# rolled.
br.figure()
ss.create_intensity_map(ss.values).set_vertical_shift(10).plot()

# %% Create intensity map with vertical shift = np.linspace(1, 2, len(ss))
# Note how the first column has no roll, only the axis is shifted by 1
# The second column is rolled by 1 pixel and the third column by 3 pixels
# All columns are rolled in relation to the first column.
# Again, it is up to you to decide if rolled out pixels make sense in your analysis.
shift = np.linspace(1, 2, len(ss))
br.figure()
ss.create_intensity_map(ss.values).set_vertical_shift(shift).plot()

# %% Create intensity map with vertical shift = np.linspace(10, 12, len(ss))
# Note that this will give the same result as shift = np.linspace(1, 2, len(ss))
# with the exception that the y axis (y_centers) are shifted by 10 units.
shift = np.linspace(10, 12, len(ss))
br.figure()
ss.create_intensity_map(ss.values).set_vertical_shift(shift).plot()

# %% Create intensity map with vertical shift = np.linspace(1, 8, len(ss))
# Note that the y axis goes from 1 to 6, therefore, a shift of 8 will roll
# some columns by more than the axis range which would mean that some 
# columns maybe shouldn't even be in the image. As the pixels are
# reintroduced on the other side of the image, this may create non-sensical 
# results. 
# In the past, this would result in a ERROR, but for now, the code will allow
# this to happen, but it is up to you to decide if rolled out pixels make sense
# in your analysis.
shift = np.linspace(1, 8, len(ss))
br.figure()
ss.create_intensity_map(ss.values).set_vertical_shift(shift).plot()

# %% An alternative to the rolled out pixels is to extend the y axis (y_centers)

# extend the lower end with zeros
ss[0].check_step()  # trigger the calculation of the step size for spectrum 0
step = ss[0].step

xleft = np.arange(-2, 1, step)
left = br.Spectrum(x=xleft, y=np.zeros(len(xleft)))

xright = np.arange(6 + step, 15, step)
right = br.Spectrum(x=xright, y=np.zeros(len(xright)))

# concatenate the left and right spectra to the original spectra
ss2 = br.Spectra()  # create new spectra to not modify the original spectra
for i in range(len(ss)):
    ss2.append(br.Spectra([left, ss[i], right]).concatenate())

# plot for verification
br.figure()
_ = ss2.plot()  

# plot image for verification
br.figure()
ss2.check_monotonicity()
ss2.create_intensity_map(ss.values).plot()

# now, we apply a large shift
shift = np.linspace(1, 8, len(ss))
x = np.arange(-2, 15, step)
br.figure()
ss2.interp(x=x).create_intensity_map(ss.values).set_vertical_shift(shift).plot()
# Note how we may have to quicly interpolate the x axis 
# transforming spectra to image requires stacking spectra as columns. 
# Therefore, they must have the same y axis to be able to "fit" the image
# Also, shiftting requires that the y axis is uniform (unvariable step)
# because otherwise rolling pixels would have no-sense


