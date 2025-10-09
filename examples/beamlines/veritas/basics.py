#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""The brixs.beamlines.veritas module offers a python solution for data processing
and data analysis for data collected at the VERITAS beamline of MAX-IV.

Example data for running this example can be downloaded at this Onedrive link
https://1drv.ms/f/c/6666810d266a3cc9/EjNzJssVTWhFqAOx0-nDCxoBSNt2cvvF8zqOmNdx_zzslQ?e=BOAgx4
Note that some less important files have been removed from this example data to
make the folder lighter and easier to download.

Besides brixs requirements (numpy and matplotlib), brixs.beamlines.veritas 
module also requires h5py.

Onsite, the path of the top (main) directory is 

'/data/visitors/veritas/<proposal number>/<date>/'

for example:

'/data/visitors/veritas/20230650/2023112108/'

Inside this folder there are two other folders:

    process  --> python scripts and data analysis files will be saved here
    raw      --> data will be saved here

data is saved on hdf5 files and the user must must define a dataset name. For 
each dataset, there will be two hdf5 files

    <dataset>.h5      --> for xas, linescans, and meshscans
    <dataset>_DLD.h5  --> for rixs


One can use this code directly on MAX-IV JupyterHub through the
following link https://jupyterhub.maxiv.lu.se
Note that to access it you need to be onsite or connected to a MAX-IV VPN.

As for today, I do not think one can install brixs via pip directly on your 
MAX-IV JupyterHub workspace. Therefore, you have to upload a zipped version of brixs
to your folder.

For that, download brixs directly from github

https://github.com/cwgaldino/brixs

click on Code/Download ZIP. The downloaded folder will be named brixs-main.zip.
Click and drag this zipped folder to jupyter slurm directory
To unzip it, open a new tap on jupyter and click on Terminal.
In the terminal, type  

unzip brixs-main.zip

The unzipped folder named `brixs-main` will appear. That's all. 

In you script, make sure to add this folder to your path using the two lines 
below
>>> import sys
>>> sys.path.insert(0, 'brixs-main/')

Note that brixs requires only numpy and matplotlib to work. 
The brixs.beamlines.veritas modules requires h5py.

See examples below.
"""

# %% ========================== standard imports ========================= %% #
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
from pprint import PrettyPrinter
pretty_print = PrettyPrinter(depth=4).pprint

# %% =========================== brixs imports =========================== %% #
#### The following two lines are only necessary if you are importing brixs 
#### from a local file
# import sys
# sys.path.insert(0, 'brixs-main/')
import brixs as br
import brixs.beamlines.veritas as veritas

# %% ============================== settings ============================= %% #
# matplotlib (optional)
get_ipython().run_line_magic('matplotlib', 'qt5')
plt.ion()

# brixs (optional)
# br.settings.FIGURE_POSITION = (283, 567)
# br.get_window_position()
# %%

# %% ============================= folderpaths =========================== %% #
# the path to the top/main folder is going to be used constantly, so we are better
# off defining a variable for it
TOP = Path(r'C:\Users\galdin_c\github\brixsexampledata\beamlines\veritas')
RAW = TOP/'raw'

DLD = RAW/'dataset1_DLD.h5'
XAS = RAW/'dataset1.h5'
# %%


# %  ===================================================================== %% #
# %  ======================== supporting functions ======================= %% #
# %% ===================================================================== %% #

# help can be used for any function
help(veritas.scanlist)

# get list of scans: veritas.scanlist(filepath)
veritas.scanlist(DLD)
veritas.scanlist(XAS)

# get available data from scan: veritas.tree(filepath, scan)
scan = 140
print(veritas.tree(scan, XAS))

# get metadata for quick inspection: veritas.get_metadata(filepath, scan, entry)
scan = 140
veritas.get_metadata(scan, XAS, 'measurement/pre_scan_snapshot/a_slit1_hz')

# metadata being imported when scan is loaded
pretty_print(veritas.xas_attrs)
pretty_print(veritas.rixs_attrs)
# veritas.xas_attrs[<processing>][<name>] = <address in the hdf5 file (must match exactly)>
# name can be anything (something easy to remember)
# processing is what time of processing that attrs needs, when in doubt, just use 'raw'
     
# adding new metadata to be loaded
veritas.xas_attrs.keys()
veritas.xas_attrs['round2']['name_of_new_attr'] = 'measurement/pre_scan_snapshot/beamline_energy'

# deleting metadata
del veritas.xas_attrs['round2']['name_of_new_attr']
# %%

# %  ===================================================================== %% #
# %  =============================== basic =============================== %% #
# %% ===================================================================== %% #

# %% xas
scan = 197
TEY, MCP, TFY, RMU = veritas.read(scan, XAS)

# metadata is saved inside each scan
TEY.get_attrs()
print(TEY.sample_x)
print(TEY.scanned_motors)

# plot
br.figure()
TEY.plot()
# %%

# %% plot multiple xas divided by I0, floored, and normalized
scans = (197, 198)

fig, axes = br.subplots(1, 3, layout='constrained', sharex=True, sharey=True)
for scan in scans:
    TEY, MCP, TFY, RMU = veritas.read(scan, XAS)
    TFY = TFY/RMU
    TEY = TEY/RMU

    MCP.floor(limits=(565, 570)).normalize(1, limits=(591, 595)).plot(ax=axes[0], label=f'{scan}')
    TFY.floor(limits=(565, 570)).normalize(1, limits=(591, 595)).plot(ax=axes[1])
    TEY.floor(limits=(565, 570)).normalize(1, limits=(591, 595)).plot(ax=axes[2])

br.leg(ax=axes[0])
plt.suptitle('xas')
br.labels.xas('MCP', ax=axes[0])
br.labels.xas('TFY', ax=axes[1])
br.labels.xas('TEY', ax=axes[2])
# %%

# %% linescan
scan = 145
TEY, MCP, TFY, RMU = veritas.read(scan, XAS)
     
# plot
br.figure()
TEY.plot()

plt.grid()
plt.ylabel('TEY')
plt.xlabel(TEY.scanned_motors[0])
# %%

# %% mesh scan
# Mesh scan have two profiles, stairs and snake
# This code only works for stairs
#          ┌──── b       /\        a
#          |            /  \      /
#     ┌────┘           /    \    /
#     |               /      \  /
# ────┘              /        \/
scan = 211
TEY, MCP, TFY, RMU = veritas.read(scan, XAS)

# this scan was done with this command
print(TEY.command)

# use the following lines to flip the image if necessary
# TEY = TEY.flipx()
# TEY = TEY.flipy()

# plot
fig, axes = br.subplots(1, 3, layout='constrained')

TEY.pcolormesh(ax=axes[0], colorbar=True)
TFY.pcolormesh(ax=axes[1], colorbar=True)
MCP.pcolormesh(ax=axes[2], colorbar=True)

axes[0].set_title('TEY')
axes[1].set_title('TFY')
axes[2].set_title('MCP')
# %%

# %% rixs
scan = 302
pe = veritas.read(scan, DLD)

# photon events x and y photon hit positions
pe.x
pe.y

# plot
br.figure()
pe.plot()
# %%
