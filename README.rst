#######
BRIXS
#######

BRIXS is a python package for processing and analysis of XAS and RIXS spectra.

Click here `https://cwgaldino.github.io/brixs/ <https://cwgaldino.github.io/brixs/>`_ for brixs documentation.

.. contents:: Table of Contents


Introduction (might be outdate, I'm rewriting it)
#################################################

BRIXS makes use of object oriented (OO) programming to easy the processing and 
analysis of data.
This way, instead of `functions` one would have `classes` and `methods`. 
This is makes sense because we often see ourselves doing the same operations over 
and over again for multiple datasets.

For instance, an example of spectrum processing would look something like this 

.. code-block:: python

   # import data
   x, y = read_data(...)  

   # initial processing
   x = x + 10          # shift the x-axis
   y = y * 2.1         # apply a multiplicative factor

   # interpolating data
   x_new = np.arange(0, 100, 1000)  # define a more suitable x-axis
   x     = np.interp(x_new, x, y)   # interpolate data to the new x-axis

   # fit data with a gaussian peak (assuming you have a gaussian function defined elsewhere)
   guess = [max(y),            # amplitude guess
            x[argmax(y)],      # peak center guess
            x[argmax(y)]*0.1]  # with guess
   popt, _ = scipy.optimize.curve_fit(gaussian, x, y, p0=guess)

   # plotting data
   plt.figure()                     # open new figure
   plt.plot(x, y)                   # plot data
   plt.plot(x, gaussian(x, popt*))  # plot fitting 

In a OO approach, the same processing would look like this,

.. code-block:: python

   s = br.Spectrum(...)   # import data

   s.shift  = 10          # shift the x-axis
   s.factor = 2.1         # apply a multiplicative factor
   s.interp(0, 10, 1000)  # interpolate data
   s.fit_peak()           # fit data with a gaussian peak

   plt.figure()    # open new figure
   s.plot()        # plot data
   s.peaks.plot()  # plot fitting 

In this simple case, we can see that code readability is better in the OO case. 
At the same time, this does not limit the most experienced users, because you 
can always apply a `functional approach` by extracting the x and y data from the 
object like,

.. code-block:: python

   x = s.x
   y = s.y


Keeping track and labeling data is also more intuitive in a OO approach. One can
load different datasets and easily label them, 

.. code-block:: python

   s1 = br.Spectrum(...)
   s2 = br.Spectrum(...)
   data_collected_yesterday = br.Spectrum(...)
   data_collected_today = br.Spectrum(...)


and metadata can be stored directly inside objects,

.. code-block:: python

   s = br.Spectrum(...)      # import data
   print(s.H)                # Metadata: Magnetic Field     
   print(s.T)                # Metadata: Temperature           

moreover, new metadata can be added on the fly,

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

The Image object is used for handling 2D arrays, like detector images. For detectors
capable of single photon count, one can use a centroid algorithm to get a sub-pixel
resolution. The output a a centroid algorithm is a photon events list, which is 
handled by the PhotonEvents object. Either way, the detector data is eventually 
turned into a spectrum which is handled by the Spectrum object. This is the most 
rich object of the BRIXS package so far. Batch operation, data alignment, or any
data manipulation that requires comparison between many spectra can be done via 
the Spectra object.

Having only four classes makes the code easy to maintain. 
.. Despite of that, BRIXS also have a secondary class

.. .. code-block:: python

..    peaks = br.Peaks()

.. which are used for peak finding and fitting. These are nothing more than `LMFIT
.. Parameters object <https://lmfit.github.io/lmfit-py/parameters.html#>`_ with some extra features.

.. BRIXS also have additional smaller modules with everyday functions, which we 
.. call `backpack functions`.

.. .. code-block:: python

..    br.figmanip
..    br.filemanip
..    br.arraymanip
..    br.interact
..    br.model_functions
..    br.xlsx

.. These modules are self-contained and independent of the rest of the project. BRIXS
.. depends on them, but they do not depend on BRIXS. Therefore, they can be used independently.
.. Their functionality ranges from, array operations to spreadsheet manipulation.
.. More specific and focused modules are found in the support folder. 

.. All BRIXS classes and modules are independent of `data collection methods` and 
.. `data file types`. All the "file reading" functions, which reads a file and converts it
.. to one of the 4 major BRIXS objects, are stored in the `file_reading` folder. For
.. example, for data collected at ADRESS beamline of PSI, one can use the following
.. line of code to get the spectrum from scan number 56,

.. .. code-block:: python

..    from brixs.file_reading import ADRESS
..    ADRESS.read(folderpath, 'Cu_', 56)

.. The modules inside the `file_reading` folder, may be used as an example for 
.. writing new modules for other file types.


Installation
#############

.. There are two recommended methods:

.. 1. Using pip (This method isn't recommended for now

.. .. code-block::    
   
..    pip install git+https://github.com/cwgaldino/brixs


.. or

.. 2. 
Cloning (or downloading) the GitHub repository then adding brixs to the "path":

.. code-block:: python

    import sys
    sys.path.append('<path-to-brixs>')
    import brixs as br

Soon, this package will be available via pip.


Requirements
############

Base (required):

- numpy

- scipy

- matplotlib

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
