BRIXS' core
=================================

BRIXS has 4 core objects:

.. toctree::
   :maxdepth: 1

   Spectrum
   Spectra
   Image
   PhotonEvents

**Image**: used for handling 2D arrays, like detector images. 

**PhotonEvents**: used for handling x, y arrays with uncorrelated x, y, i.e., 
y is not a function of x. Usually, this is the case for detectors capable of 
single photon count (where one can use a centroid algorithm to get a sub-pixel
resolution). The output a a centroid algorithm is a photon events list, which is 
handled by this object. 

**Spectrum**: used for handling x, y arrays where y=f(x). This is the most 
rich object of the BRIXS package so far. 

**Spectra**: used for batch operation, data alignment, or any
data manipulation that requires comparison between many spectra

Having only four classes makes the code easy to maintain and improves reusability.


