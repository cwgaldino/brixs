#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Basic plotting examples"""

# standard imports
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

# default brixs import
import brixs as br

# creating and plotting a spectrum
x = np.linspace(-np.pi, np.pi, 100)
s = br.Spectrum(x=x, y=np.sin(x))

# s.plot() is the same as plt.plot(s.x, s.y)
br.figure()
s.plot()
plt.plot(s.x, s.y)

# br.figure() is the same as plt.plot()
# but br.figure() has mouse callbacks implemented
# mouse left click copies the data cursor x position to the clipboard
# mouse right click copies the data cursor y position to the clipboard
# mouse middle click copies the FIGURE cursor (x, y) position to the clipboard
# figure cursor position from 0 to 1 where (0, 0) is the bottom left corner

# TODO


