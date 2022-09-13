#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""BRIXS configuration file"""

class _settings():
    def __init__(self):
        self.ALWAYS_PLOT_NEW_WINDOW = False
        self.FIGURE_POSITION        = None
        # self.DEFAULT_TAIL_FOR_PEAKS = None
        # self.DEFAULT_SPECTRUM_CALIB = 1
        self.MAX_ERROR_STEP_X = 0.1


    def __str__(self):
        text = f'ALWAYS_PLOT_NEW_WINDOW: {self.ALWAYS_PLOT_NEW_WINDOW}\n' +\
               f'FIGURE_POSITION: {self.FIGURE_POSITION}\n' +\
               f'MAX_ERROR_STEP_X: {self.MAX_ERROR_STEP_X}\n' #+\
               # f'DEFAULT_TAIL_FOR_PEAKS: {self.DEFAULT_TAIL_FOR_PEAKS}\n' #+\
               # f'DEFAULT_SPECTRUM_CALIB: {self.DEFAULT_SPECTRUM_CALIB}\n' #+\
        return text

    def __repr__(self):
        return self.__str__()

settings = _settings()
