#######
BRIXS
#######

BRIXS is an object-oriented (OO) python package for processing/analysis of XAS and RIXS spectra.

Click here `https://cwgaldino.github.io/brixs/ <https://cwgaldino.github.io/brixs/>`_ for brixs documentation.

.. contents:: Table of Contents


Introduction
##########################################################

Object-oriented is a programming model that organizes software design 
around objects, rather than functions. An OO approach makes sense for data 
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
   s.shift  = 10          # shift the x-axis
   s.factor = 2.1         # apply a multiplicative factor
   s.interp(0, 10, 1000)  # interpolate data
   s.fit_peak()           # fit data with a gaussian peak

   # display data
   br.figure()                        # open new figure
   s.plot(label=f"T={s.T}, P={s.P}")  # plot data
   s.model.plot()                     # plot fitting 
   br.leg()                           # legend
   plt.show()                         # show figure

We can argue that keeping track and labeling data is more intuitive in an OO approach as the 
number of variables is drastically reduced. For instance, if one tries to load 2 different 
datasets we have 2 variables with OO vs 8 variables using a functional approach: 

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
      s.shift = 10
      s.factor = 2.1
      s.interp(0, 10, 1000)
      s.fit_peak()

   # add new method to all Spectrum objects
   br.Spectrum.processing = common_processing

   # from now on, one can use new method on every Spectrum object
   s1.processing()
   s2.processing()


Basics
#############

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

BRIXS also comes with additional functionally from supporting modules. See
below: 

backpack
==========

*Backpack* is a module with quality-of-life (QOL) functions. This module is completely
independent from brixs. As for the time of writing, these are the submodules: 

.. code-block:: python

   brixs.figmanip          # matplotlib QOL functions
   brixs.filemanip         # file reading and saving QOL functions
   brixs.arraymanip        # array manipulation QOL functions
   brixs.numanip           # float/int manipulation QOL functions
   brixs.interact          # user interaction QOL functions

See the documentation (PUT LINK HERE) for a description of the functions available. 
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

finder
==========



labels
==========


model 
==========


crystal
==========

Module with function for calculating momentum transfer in single crystals. 
It is assumed that the photon hits the crystal surface at a angle th and is 
scattered in a 2th angle as the drawing below::
     \      /.
      \    /   .
       \  /     .
   th ( \/       .
   ┌──────────┐  . 2th
   ├ crystal  ┤  .
   └──────────┘ .
           \   .
            \.        
            
This module can be used like this

.. code-block:: python

      # import 
      import brixs.crystal

      # functions available
      br.ev2angstrom()
      br.calculate_q_transfer()
      br.momentum2rlu()

The description of each function can be accessed via the python help() function or 
by reading the documentation (PUT LINK HERE).





BRIXS
depends on them, but they do not depend on BRIXS. Therefore, they can be used independently.
Their functionality ranges from, array operations to spreadsheet manipulation.
More specific and focused modules are found in the support folder. 

All BRIXS classes and modules are independent of `data collection methods` and 
`data file types`. All the "file reading" functions, which reads a file and converts it
to one of the 4 major BRIXS objects, are stored in the `file_reading` folder. For
example, for data collected at ADRESS beamline of PSI, one can use the following
line of code to get the spectrum from scan number 56,

.. code-block:: python

   from brixs.file_reading import ADRESS
   ADRESS.read(folderpath, 'Cu_', 56)

The modules inside the `file_reading` folder, may be used as an example for 
writing new modules for other file types.


Installation
#############

There are two recommended methods:

1. Using pip 

.. code-block::    
   
   pip install git+https://github.com/cwgaldino/brixs

or

2. Cloning (or downloading) the GitHub repository then adding brixs to the "path":

.. code-block:: python

    import sys
    sys.path.append('<path-to-brixs>')
    import brixs as br



Requirements
############

Base (required):

- numpy



- matplotlib


- scipy
.. - lmfit >= 1.2.2

Reciprocal space calculations:

- pbcpy

Other packages may be required for file reading functions.


Usage
#############

Refer to the examples folder on GitHub for code examples.


Documentation
#############

Click here `https://cwgaldino.github.io/brixs/ <https://cwgaldino.github.io/brixs/>`_ for brixs documentation.
