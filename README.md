# BRIXS

Python package for analysis of XAS/RIXS spectra.

Documentation: https://cwgaldino.github.io/brixs/

WARNING: The code has been in active development this past few days and some documentation might not be up to date.

## Installation

There are two recommended methods:

1. Using pip

```python
    pip install git+https://github.com/cwgaldino/brixs
```

or

2. Cloning (or downloading) the GitHub repository then adding brixs to the "path":

```python
    import sys
    sys.path.append('<path-to-brixs>')
    import brixs as br
```

## Requirements

Base (required):
- numpy
- scipy
- matplotlib

Reading files:
- detect_delimiter
- h5py (only for reading .hdf files)

Reciprocal space calculations:
- pbcpy

## Usage (not an extensive description)

Refer to the examples folder on GitHub for code examples.

BRIXS is based on four objects:

```python
    import brixs as br

    im = br.Image()
    pe = br.PhotonEvents()
    s  = brixs.Spectrum()
    ss = brixs.Spectra()
```

Also two support objects:

```python
    import brixs as br

    p  = br.Peak()
    ps = br.Peaks()
```



### Image attributes and methods

```python
    # basic
    im.data              # np.array (float)
    im.vmin              # float (read only)
    im.vmax              # float (read only)
    im.shape             # tuple (read only)
    im.x_centers         # np.array (float)
    im.y_centers         # np.array (float)
    im.x_edges           # np.array (float)
    im.y_edges           # np.array (float)

    # binning
    im.nbins             # np.array [runs Image.binning()]
    im.bins_size         # np.array [runs Image.binning()]
    im.reduced           # br.Image (read only)

    # shifts
    im.shifts_v          # np.array
    im.shifts_h          # np.array
    im.p                 # np.array      (read only)
    im.f                 # function f(x) (read only)
    im.calculated_shifts  # brixs.Spectrum (read only)

    # spectrum (Computed Attributes)
    im.histogram         # brixs.Spectrum
    im.spectrum          # brixs.Spectrum
    im.spectrum_v        # brixs.Spectrum
    im.spectrum_h        # brixs.Spectrum
    im.columns           # brixs.Spectra
    im.rows              # brixs.Spectra


    # methods
    im.save()
    im.load()

    im.floor()
    im.crop()

    im.pcolormesh()
    im.imshow()
    im.plot()

    im.binning()
    im.calculate_histogram()
    im.calculate_spectrum()
    im.calculate_shifts()
    im.set_shifts()
    im.fix_curvature()
```


### PhotonEvents attributes and methods

```python
    # basic
    pe.data          # np.array
    pe.x             # float (read only)
    pe.y             # float (read only)
    pe.I             # tuple (read only)
    pe.shape         # np.array

    # binning
    pe.nbins             # np.array [runs Image.binning()]
    pe.bins_size         # np.array [runs Image.binning()]
    pe.reduced           # br.Image (read only)

    # shifts
    pe.shifts_v          # np.array
    pe.shifts_h          # np.array
    pe.p                 # np.array      (read only)
    pe.f                 # function f(x) (read only)
    pe.calculated_shift  # brixs.Spectrum (read only)


    # methods
    pe.save()
    pe.load()

    pe.floor()
    pe.crop()

    pe.pcolormesh()
    pe.imshow()
    pe.plot()

    pe.binning()
    pe.calculate_histogram()
    pe.calculate_spectrum()
    pe.calculate_shift()
    pe.set_shift()
    pe.fix_curvature()
```

### Spectrum attributes and methods

```python
    # basic
    s.data              # np.array (2 column array)
    s.x                 # np.array (array)
    s.y                 # np.array (array)
    s.area              # float

    # modifiers
    s.offset        # float [runs Spectrum.set_offset()]
    s.factor        # float [runs Spectrum.set_factor()]
    s.calib         # float [runs Spectrum.set_calib()]
    s.shift         # float [runs Spectrum.set_shift()]
    s.shift_roll    # int   [runs Spectrum.set_shift()]
    s.shift_interp  # float [runs Spectrum.set_shift()]

    # check
    s.step          # float
    s.monotonicity  # string

    # fit
    s.fit     # brixs.Spectrum
    s.residue # brixs.Spectrum
    s.guess   # brixs.Spectrum
    s.R2      # float
    s.pcov    # np.array

    # peaks
    s.peaks   # brixs.peaks


    # methods
    s.save()
    s.load()

    s.plot()

    s.check_step_x()
    s.check_monotonicity()
    s.fix_monotonicity()

    s.set_calib()
    s.set_shift()
    s.set_factor()
    s.set_offset()

    s.interp()
    s.extract()
    s.crop()
    s.floor()
    s.flip()
    s.normalize()
    s.zero()
    s.calculate_area()

    s.find_peaks()
    s.fit_peak()
    s.fit_peaks()
    s.polyfit()
    s.apply_correction()
```


### Spectra attributes and methods

```python
    # basic
    ss.data         # list of brixs.Spectrum
    ss.area         # list (computed attribute)

    # modifiers
    ss.offset        # list (computed attribute)
    ss.factor        # list (computed attribute)
    ss.calib         # list (computed attribute)
    ss.shift         # list (computed attribute)
    ss.shift_roll    # list (computed attribute)
    ss.shift_interp  # list (computed attribute)

    # calculated modifiers
    ss.calculated_offset  # brixs.Spectrum
    ss.calculated_factor  # brixs.Spectrum
    ss.calculated_calib   # brixs.Spectrum
    ss.calculated_shift   # brixs.Spectrum

    # check
    ss.step          # float
    ss.length        # float
    ss.x             # list
    ss.monotonicity  # string

    # fit
    ss.fit     # brixs.Spectrum
    ss.residue # brixs.Spectrum
    ss.guess   # brixs.Spectrum
    ss.R2      # float

    # peaks
    ss.peaks   # list  [runs Spectra.get_peaks()]
    ss.error   # list  [runs Spectra.get_errors()]

    # map
    ss.map     # brixs.Image (computed attribute)


    # methods
    ss.save()
    ss.load()

    ss.plot()

    ss.check_length()
    ss.check_step_x()
    ss.check_same_x()
    ss.check_monotonicity()
    ss.fix_monotonicity()

    ss.set_calib()
    ss.set_shift()
    ss.set_factor()
    ss.set_offset()

    ss.interp()
    ss.extract()
    ss.crop()
    ss.floor()
    ss.flip()
    ss.concatenate()
    ss.align()
    ss.normalize()

    ss.calculate_shift()
    ss.calculate_factor()
    ss.calculate_offset()
    ss.calculate_calib()

    ss.calculate_sum()
    ss.calculate_map()

    ss.get_peaks()
    ss.plot_peaks()
    ss.find_peaks()
    ss.fit_peak()
    ss.fit_peaks()
    ss.polyfit()
```





## Information for developers

I will list here some practices that I have been trying to follow while developing brixs. This should facilitate the maintenance and expansion of the package.

- the metaclass `_Meta()` helps other classes creating read-only and non-removable attributes. The variables `_read_only = []` and `_non_removable = []` should be defined at the beginning of a class before __init__(). Attributes will be created based on the names listed inside these lists.
- __init__() methods should initialize *ALL* attributes, except for computed attributes.
- __init__() should call _sort_args().
- _sort_args() should sort the arguments and return a data type or a filepath to __init__(), and nothing else
- every time a new attribute is added to an object, one has to worry about if and how it is going to be saved/loaded from a file, if it should be reseted when other attributes change, and if it should be passed to child objects.
- every time an attribute is changed inside an object, one has to worry about which "check attributes" must be reseted.
