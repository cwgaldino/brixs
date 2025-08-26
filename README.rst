##########################################################
BRIXS
##########################################################

BRIXS is an object-oriented (OO) python package for processing/analysis of XAS and RIXS spectra.

.. Click here `https://cwgaldino.github.io/brixs/ <https://cwgaldino.github.io/brixs/>`_ for brixs documentation.

##########################################################
Installation
##########################################################

Currently, the recommended method for installing BRIXS is via pip command:

.. code-block:: bash

   pip install brixs

But brixs can also be used by downloading the ZIP file directly from the GitHub repository.
For that, download the repository as a ZIP file. Extract all and remove the contents from
the folder. Rename the folder from `brixs-main` to `brixs`. In the end, you should
have a folder named `brixs` with some files inside (README.rst, LICENCE, ...) and 
another folder named `brixs` inside. Then add brixs to the "path" using the code below at the top of you script,

   .. code-block:: python

      import sys
      sys.path.append('<path-to-brixs>')
      import brixs as br

   note that `<path-to-brixs>` must point to the folder that has the folder `brixs` and 
   the file README.rst. 


##########################################################
Overview
##########################################################

BRIXS is based on four major objects:

.. code-block:: python

   im = br.Image()
   pe = br.PhotonEvents()
   s  = br.Spectrum()
   ss = br.Spectra()


Once a BRIXS object has been created, one can use the methods,

.. code-block:: python

   s  = br.Spectrum()
   s.get_methods()
   s.get_attrs()

to print all methods (functions) and attrs related to that object (in this case `s`).

The description of methods can be accessed via the python `help()` function. See 
example below,

.. code-block:: python

   s  = br.Spectrum()
   help(s.set_shift)


##########################################################
Introduction
##########################################################

Object-oriented is a programming model that organizes software design 
around objects rather than functions. An OO approach makes sense for data 
processing/analysis of XAS and RIXS because it allows more intuitive syntax
and make the project easier to maintain and upgrade as it grows.

At the lab (beamline), the output of an experiment is often a data file which is
saved on a computer. We wish to open this file, perform some operations (processing),
and display it in a meaningful way. Below we can see some typical examples of 
everyday data processing and how OO can facilitate data processing. 

Suppose that our dataset is composed of two data arrays (xy-pair), where y=f(x).
Each dataset also have a temperature (T) and pressure (P) associated to it. A 
typical script for reading, operation, and displaying one dataset from a file would
look like this: 

.. code-block:: python

   # import packages
   from mypackage import read_data
   import matplotlib.pyplot as plt
   import numpy as np
   import scipy

   # import data
   x, y, T, P = read_data(...)  

   # initial processing
   x = x + 10          # shift the x-axis
   y = y * 2.1         # apply a multiplicative factor

   # interpolating data
   x_new = np.arange(0, 100, 1000)  # define a more suitable x-axis
   x     = np.interp(x_new, x, y)   # interpolate data to the new x-axis

   # fit data with a gaussian peak
   gaussian = lambda x, mu, sig: 1.0/(np.sqrt(2.0*np.pi)*sig)*np.exp(-np.power((x-mu)/sig, 2.0)/2)
   guess    = [max(y), x[argmax(y)], x[argmax(y)]*0.1]         
   popt, _ = scipy.optimize.curve_fit(gaussian, x, y, p0=guess)

   # display data
   plt.figure()                           # open new figure
   plt.plot(x, y, label=f"T={T}, P={P}")  # plot data
   plt.plot(x, gaussian(x, popt*))        # plot fitting 
   plt.legend()                           # legend
   plt.show()                             # show figure

In an OO approach, the same processing would look like this,

.. code-block:: python

   # import packages
   import brixs as br

   # import data
   s = br.Spectrum(...)   

   # initial processing
   s = s.set_shift(10)           # shift the x-axis
   s = s.set_factor(2.1)         # apply a multiplicative factor
   s = s.interp(0, 10, 1000)     # interpolate data
   fit, popt, sigma, model = s.fit_peak()  # fit data with a gaussian peak

   # display data
   br.figure()                        # open new figure
   s.plot(label=f"T={s.T}, P={s.P}")  # plot data
   fit.plot()                     # plot fitting 
   br.leg()                           # legend
   plt.show()                         # show figure

or we can use a one-liner:

   .. code-block:: python

   # import packages
   import brixs as br

   # data processing
   s = br.Spectrum(...).set_shift(10).set_factor(2.1).interp(0, 10, 1000)
   fit, popt, sigma, model = s.fit_peak()

We can argue that keeping track and labeling data is more intuitive in an OO approach as the 
number of variables is drastically reduced. For instance, if one tries to load 2 different 
datasets we have 2 variables with OO vs 8 variables using a functional approach. See below: 

.. code-block:: python

   # dataset number 1 and 2
   s1 = br.Spectrum(...)
   s2 = br.Spectrum(...)

   # x and y arrays as well values of temperature and pressure are store inside the object
   print('x1', s1.x)
   print('T1', s1.T)

where in a functional approach, the number of variables can easily start to became overwhelming

.. code-block:: python
   
   # dataset number 1 and 2
   x1, y1, T1, P1 = read_data(...)  
   x2, y2, T2, P2 = read_data(...)  

At the same time, OO does not limit the most experienced users, because you 
can always simulate a `functional approach` by extracting the x and y data from the 
object like,

.. code-block:: python

   s1 = br.Spectrum(...)
   x = s.x
   y = s.y

also new metadata can be added on the fly,

.. code-block:: python

   s.angle = 12.53

Just like metadata, repetitive tasks can be added to the object,

.. code-block:: python

   # define new method
   def common_processing(s):
      s = s.shift = 10
      s = s.factor = 2.1
      s = s.interp(0, 10, 1000)
      return s.fit_peak()

   # add new method to all Spectrum objects
   br.Spectrum.processing = common_processing

   # from now on, one can use new method on every Spectrum object
   s1.processing()
   s2.processing()

##########################################################
Core
##########################################################

BRIXS is based on four major objects:

.. code-block:: python

   im = br.Image()
   pe = br.PhotonEvents()
   s  = br.Spectrum()
   ss = br.Spectra()

The **Image** object is used for handling 2D arrays, like detector images. For detectors
capable of single photon count, one can use a centroid algorithm to get a sub-pixel
resolution. The output a a centroid algorithm is a photon events list, which is 
handled by the **PhotonEvents** object. Either way, the detector data is eventually 
turned into a spectrum which is handled by the **Spectrum** object. This is the most 
rich object of the BRIXS package so far. Batch operation, data alignment, or any
data manipulation that requires comparison between many spectra can be done via 
the **Spectra** object. Having only four classes makes the code easy to maintain. 
A better description of each object will be given later in this readme. 

##########################################################
Support modules
##########################################################

BRIXS also comes with additional functionally from supporting modules.

================================
backpack
================================

*Backpack* is a module with quality-of-life (QOL) functions. This module is completely
independent from brixs. As for the time of writing, these are the submodules: 

.. code-block:: python

   brixs.figmanip          # matplotlib QOL functions
   brixs.filemanip         # file reading and saving QOL functions
   brixs.arraymanip        # array manipulation QOL functions
   brixs.numanip           # float/int manipulation QOL functions
   brixs.interact          # user interaction QOL functions

See the documentation for a description of the functions available. 
All functions within *backpack* are readily available when brixs is imported. For
instance, the function *brixs.arraymanip.check_monotonicity* which checks the 
monotonicity of an array can be called directly from brixs:

.. code-block:: python

      # import brixs
      import brixs as br

      # array
      a = [1, 2, 3, 4, 5, 6]

      # check monotonicity 
      br.check_monotonicity(a)

================================
finder
================================

Module for quickly saving and recovering processed spectra so you can avoid 
running functions multiple types with same input parameters. Finder is imported
with brixs

.. code-block:: python

   # import brixs
   import brixs as br

   # set finder folderpath
   br.finder.folderpath = '<folderpath>'

   # apply the decorator to you function
   @br.finder.track
   def processing_function(a, b, c):
      s = <does something with a, b and c and returns s>
      return s

   # processing may take a while if it is the first time you run
   s = processing_function(a=1, b=2, c=3)  

   # if you run processing with same parameters, it runs instantly because 
   # finder recovers already processed spectrum
   s = processing_function(a=1, b=2, c=3)  

Full description of finder functionally can be found inside brixs.addons.finder.py file.

================================
labels
================================

Module with common x and y labels for xas and rixs plots.

.. code-block:: python

   # this
   br.labels.xas()

   # is the same as this
   plt.xlabel('Photon Energy (eV)')
   plt.ylabel('Intensity')

Full description of labels can be found inside brixs.addons.labels.py file.



================================
model
================================

Module for data fitting. For enabling fitting functionally do

.. code-block:: python

   # enable fitting functionality
   import brixs.model

   # model functions are then available
   br.model.gaussian()

This module is fully implemented, but improvements are often implemented. 
More information can be found inside brixs.model.model.py file.

================================
beamlines
================================

BRIXS objects and modules are independent of `data collection methods` and 
`data file types`. All "file reading" functions specific for each lab or beamline
which reads a file and converts it to one of the 4 major BRIXS objects can be
found in the `beamlines` folder. 

For example, data collected at I21 beamline of Diamond Light Source, can use 
imported as a **Image** object using the code below,

.. code-block:: python

   # method 1
   import brixs.beamlines.I21 as I21
   im = I21.read(<filepath>)

   # method 2
   from brixs.beamlines.I21 import read
   im = read(<filepath>)

Please, refer to the folder "beamlines" to see if code has been implemented for 
the beamline of interest and refer to the beamline's .py files for more information.

================================
crystal
================================

Module with functions for calculating momentum transfer in single crystals. 
It is assumed that the photon hits the crystal surface at a angle th and is 
scattered in a 2th angle. See drawing inside brixs.crystal.crystal.py file for 
more information. This module can be used like this

.. code-block:: python

      # import 
      import brixs.crystal

      # functions available
      br.ev2angstrom()
      br.calculate_q_transfer()
      br.momentum2rlu()

This module has some limitations. Please refer to the file brixs.crystal.crystal.py
for more information.




##########################################################
Requirements
##########################################################

================================
Base (required)
================================

- numpy
- matplotlib

Some modules require additional imports:

================================
brixs.model
================================

- scipy
- lmfit

================================
brixs.crystal
================================

- pbcpy

================================
brixs.beamlines
================================

Some modules here might require 

- h5py

