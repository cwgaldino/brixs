#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""BRIXS configuration file

FIGURE_FORCE_NEW_WINDOW

FIGURE_FORCE_ON_TOP

FIGURE_POSITION
Default figure position

FIGURE_SIZE
Default figure size

FIGURE_DPI
Default figure resolution




"""

# FIGURE_FIX_RESOLUTION
# if True, figure resolution won't change when moving figure to a monitor with 
# different screen resolution. Figure size won't change.

class _settings():
    def __init__(self):
        self.FIGURE_FORCE_NEW_WINDOW  = False
        self.FIGURE_FORCE_ON_TOP      = False

        self.FIGURE_POSITION = None
        self.FIGURE_SIZE     = None
        self.FIGURE_DPI      = None

        # not implemmented yet
        # self.FIGURE_UPDATE_POSITION   = False  # if True, open next figure with same position as the current figure
        # self.FIGURE_UPDATE_SIZE       = False  # if True, open next figure with same size as the current figure
        
        # self.FIGURE_FIX_RESOLUTION    = False 

        # self.DEFAULT_TAIL_FOR_PEAKS = None
        # self.DEFAULT_SPECTRUM_CALIB = 1
        self.MAX_ERROR_STEP_X = 0.1

    def save(self):
        f = open("demofile2.txt", "a")
        f.write("Now the file has more content!")
        f.close()

    def __str__(self):
        text = f'FIGURE_FORCE_NEW_WINDOW: {self.FIGURE_FORCE_NEW_WINDOW}\n' +\
               f'FIGURE_POSITION: {self.FIGURE_POSITION}\n' +\
               f'MAX_ERROR_STEP_X: {self.MAX_ERROR_STEP_X}\n' #+\
               # f'DEFAULT_TAIL_FOR_PEAKS: {self.DEFAULT_TAIL_FOR_PEAKS}\n' #+\
               # f'DEFAULT_SPECTRUM_CALIB: {self.DEFAULT_SPECTRUM_CALIB}\n' #+\
        return text

    def __repr__(self):
        return self.__str__()

global settings
settings = _settings()
