#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""ASCAN quick viwer for IPE beamline

TODO:
[ ] display cursor position for plot in all tabs
[ ] online plot (Thread), it is crashing all the time
"""

# %% ============================= USER PROPOSAL ========================== %% #
proposal = '20221532'
# %%

# %% ==========================  STANDARDS IMPORTS ======================== %% #
from collections.abc import Iterable
import matplotlib.pyplot as plt
from pathlib import Path
import numpy as np
import sys
import os
# %%

# %% =============================== BASE ================================= %% #
if os.name == 'nt':  # Windows
    BASE = Path('Z:/')
elif os.name == 'posix':  # Linux
    BASE = Path('/ibira')
# %%
    
# %% =========================== FOLDERPATHS ============================== %% #
IBIRA = BASE/'lnls/beamlines/ipe'

IPE   = IBIRA/'apps/IPE'

USER  = IBIRA/'proposals'/str(proposal)
PROC  = USER/'proc'
DATA  = PROC/'DATA'

BRIXS  = PROC/'brixs'
# %%

# %% =========================== BRIXS IMPORT ============================= %% #
sys.path.insert(0, '/usr/local/scripts/apps/brixs')
import brixs as br

# from brixs.beamlines import IPE as _ipe
import brixs.beamlines.ipe as ipe
import brixs.addons.fitting
# %%

# %% ========================= PYQTGRAPH IMPORT =========================== %% #
import pyqtgraph as pg
from pyqtgraph.console import ConsoleWidget
from pyqtgraph.dockarea.Dock import Dock
from pyqtgraph.dockarea.DockArea import DockArea

from pyqtgraph.Qt import QtWidgets
from pyqtgraph.Qt import QtGui
# %%

# %% ===================== THREAD IMPORT (AUTOUPDATE) ===================== %% #
from PyQt5.QtCore import QObject, QThread, pyqtSignal
import time

# Create a worker class (tentative worker class for online plotting - Thread)
class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)

    def run(self):
        """Long-running task."""
        print('auto update on')
        while r1.isChecked():

            # time
            time.sleep(2)

            # update list if number of files change
            fl = br.filelist(vars['folderpath'])
            if len(fl) != int(l1.count()):
                update_filelist()
            
            # if last item is selected, update plot
            if len(fl)-1 in [item[0] for item in vars['selected']]:
                update_plot()

            # testing
            # update_filelist()
            # print('gg')
            # for i in range(5):
            #     sleep(1)
                # self.progress.emit(i + 1)

        print('auto update off')
        self.finished.emit()
# %%

# ============================================================================ #
# %% ============================ WINDOW ================================== %% #
# ============================================================================ #

# %% Initial definitions
name                = "Quick viewer"   # app name doesn't matter
window_title        = 'BRIXS: quick viewer'
initial_window_size = (1000, 500)
# default_folderpath  = Path('./').absolute()
default_folderpath  = str(IBIRA)+'/proposals/'
width = 2
vars = {'folderpath': default_folderpath,
        'selected':[],
        'thread': None,
        'worker': None,
        'floor': ['', ''],
        'norm':  ['', '']
        }
# %%

# %% Initialization
# app = QtGui.QApplication([])
# win = QtGui.QMainWindow()
app = QtWidgets.QApplication([])
win = QtWidgets.QMainWindow()
area = DockArea()
win.setCentralWidget(area)
# %%

# %% initial setup
win.setCentralWidget(area)
win.resize(*initial_window_size)
win.setWindowTitle(window_title)
# %%

# %% support functions
def get_folderpath():
    """get text from editline, call update_folderpath"""
    print('get_folderpath')
    update_folderpath(e1.text())
    return

def select_folderpath():
    """opens file dialog window, call update_folderpath"""
    print('select_folderpath')
    dialog     = QtWidgets.QFileDialog()
    folderpath = dialog.getExistingDirectory(None, "Select Folder")
    update_folderpath(folderpath)
    return

def go_back_directory():
    """goes up one directory"""
    print('go_back_directory')
    update_folderpath(vars['folderpath'].parent)
    return

def get_selection():
    """puts selected files in vars, call update plot"""
    print('get_selection')

    # get selection
    selected  = [(index.row(), item.text()) for index, item in zip(l1.selectedIndexes(), l1.selectedItems())]
    # print(l1.selectedIndexes()[0].__dir__())

    # folderpath
    folderpath = vars['folderpath']

    # check if selection is dir
    if len(selected) == 0:
        pass
    if len(selected) == 1:
        if (folderpath/selected[0][1]).is_dir():
            update_folderpath(folderpath/selected[0][1])
            selected = []

    # save selection
    vars['selected'] = selected

    # update plot
    update_plot()

def get_normalization():    
    """puts normalization ranges in vars, update plot"""
    print('get_normalization')
    
    floor = [e2.text(), e3.text()]
    norm  = [e4.text(), e5.text()]
    
    if floor[0] == '': 
        floor[0] = None
    else:
        floor[0] = float(floor[0])
    if floor[1] == '': 
        floor[1] = None
    else:
        floor[1] = float(floor[1])
    if norm[0] == '':  
        norm[0] = None
    else:
        norm[0] = float(norm[0])
    if norm[1] == '':
        norm[1] = None
    else:
        norm[1] = float(norm[1])

    # TODO: make sure input is valid
    #
    #
    
    vars['floor'] = floor
    vars['norm']  = norm
    
# %% update
def update_folderpath(folderpath):
    """puts folderpath in vars, update text, call update filelist"""
    print('update_folderpath')
    vars['folderpath'] = Path(folderpath)
    e1.setText(str(folderpath))  
    update_filelist()
    return

def update_filelist():
    """updates filelist given the folderpath in vars"""
    print('update_filelist')
    
    # selected items
    selected = vars['selected']

    # update list
    fl = br.filelist(vars['folderpath'])
    l1.clear()
    l1.insertItems(0, [p.name for p in fl])

    # select items back
    for item in selected:
        l1.item(item[0]).setSelected(True)
    
    # scroll to bottom if last item is selected
    if len(fl)-1 in [item[0] for item in vars['selected']]:
        l1.scrollToBottom()
    # l1.scrollToBottom()
    return 

def update_plot():
    print('update_plot')
    
    # get normalization
    # get_normalization()

    # color cycler
    cycler = ('w', 'b', 'g', 'r', 'c', 'm', 'y')

    # clear
    for ax in axes:
        ax.clear()
    w3.clear()
    w4.clear()
    w5.clear()
    w6.clear()
    legend.clear()
    w7.clear()

    # plot
    for i, item in enumerate(vars['selected']):
        # get data
        filename = item[1]
        ss = ipe.read(fpath=vars['folderpath']/filename)

        detectors = ['_tey', '_fy', '_pd', '_i0']
        xps_detectors = ['XPS' + suffix for suffix in detectors]
        rixs_detectors = ['RIXS' + suffix for suffix in detectors]
        if all(element in ss.detectors for element in xps_detectors):
            detectors=xps_detectors

        elif all(element in ss.detectors for element in rixs_detectors):
            detectors=rixs_detectors   

        out = list()
        for detector in detectors:
            out.append(ss.get_by_attr('detector', detector, closest=False, verbose=False))
        TEY, TFY, I0, PD = out

        #################
        # plot complete #
        #################
        # first row
        ref = axes[0].plot(TEY.x, TEY.y, pen={'color':cycler[i%len(cycler)], 'width':width}, name=TEY.scan)
        axes[1].plot(TFY.x, TFY.y,       pen={'color':cycler[i%len(cycler)], 'width':width})
        axes[2].plot(I0.x, I0.y,         pen={'color':cycler[i%len(cycler)], 'width':width})
        # pen={'color':'w', 'width':1.5}
        
        # second row
        axes[3].plot((TEY/I0).x, (TEY/I0).y, pen={'color':cycler[i%len(cycler)], 'width':width})
        axes[4].plot((TFY/I0).x, (TFY/I0).y, pen={'color':cycler[i%len(cycler)], 'width':width})
        axes[5].plot((I0/I0).x, (I0/I0).y,   pen={'color':cycler[i%len(cycler)], 'width':width})

        # third row
        try:
            s = (TEY/I0).floor(limits=vars['floor']).normalize(value=1, limits=vars['norm'])
        except Exception as e:
            print(e)
            s = TEY/I0
        axes[6].plot(s.x, s.y, pen={'color':cycler[i%len(cycler)], 'width':width})
        
        try:
            s = (TFY/I0).floor(limits=vars['floor']).normalize(value=1, limits=vars['norm'])
        except Exception as e:
            print(e)
            s = TEY/I0
        axes[7].plot(s.x, s.y, pen={'color':cycler[i%len(cycler)], 'width':width})
        # axes[8].plot((TEY/I0).x, (TEY/I0).y, pen={'color':cycler[i%len(cycler)], 'width':width})
        # axes[8].plot((TFY/I0).x, (TFY/I0).y, pen={'color':cycler[i%len(cycler)], 'width':width})

        # plot single
        w3.plot(TEY.x, TEY.y, pen={'color':cycler[i%len(cycler)], 'width':width})
        w4.plot(TFY.x, TFY.y, pen={'color':cycler[i%len(cycler)], 'width':width})
        w5.plot(I0.x, I0.y, pen={'color':cycler[i%len(cycler)], 'width':width})
        w6.plot(PD.x, PD.y, pen={'color':cycler[i%len(cycler)], 'width':width})


        # legend
        legend.addItem(ref, ref.name())

        # metadata
        # a = pg.TextItem(text=TEY.scan, color=(200, 200, 200), html=None, anchor=(0, 0), border=None, fill=None, angle=0, rotateAxis=None)
        text = ''
        for attr in TEY.get_attrs():
            if attr in ('EPOCH', 'SETPOINT_RIXS_X', 'SETPOINT_RIXS_Y', 'PHASE', 'DVF', 'ENERGY_SP', 'SP_ENERGY', 'SP_PHASE', 'RING_CURRENT', 'TIMESTAMP'):
                pass
            else:
                text += attr + ': ' + str(getattr(TEY, attr)) + '\n'
        a = pg.TextItem(text=text, color=cycler[i%len(cycler)], anchor=(0, 0))
        a.setPos(8*(i%3), -10*(i//3))
        w7.addItem(a)

    # autorange reset
    for ax in axes:
        ax.autoRange()

    pass
# %%

# %% AUTO UPDATE
def auto_update():
    """Tentative function for online plotting (Thread)"""
    if r1.isChecked():
        # Step 2: Create a QThread object
        vars['thread'] = QThread()

        # Step 3: Create a worker object
        vars['worker'] = Worker()

        # Step 4: Move worker to the thread
        vars['worker'].moveToThread(vars['thread'] )

        # Step 5: Connect signals and slots
        vars['thread'].started.connect(vars['worker'].run)
        vars['worker'].finished.connect(vars['thread'].quit)
        vars['worker'].finished.connect(vars['worker'].deleteLater)
        vars['thread'].finished.connect(vars['thread'].deleteLater)
        # vars['worker'].progress.connect(self.reportProgress)

        # Step 6: Start the thread
        vars['thread'].start()

    else:
        pass
# %%

# %% cursor position
def mouse_moved(event):
    """get and display cursor position"""
    pos = event[0]  ## using signal proxy turns original arguments into a tuple
    
    ax = False
    if axes[0].sceneBoundingRect().contains(pos):
        ax = axes[0]
    elif axes[1].sceneBoundingRect().contains(pos):
        ax = axes[1]
    elif axes[2].sceneBoundingRect().contains(pos):
        ax = axes[2]
    elif axes[3].sceneBoundingRect().contains(pos):
        ax = axes[3]
    elif axes[4].sceneBoundingRect().contains(pos):
        ax = axes[4]
    elif axes[5].sceneBoundingRect().contains(pos):
        ax = axes[5]
    elif axes[6].sceneBoundingRect().contains(pos):
        ax = axes[6]
    elif axes[7].sceneBoundingRect().contains(pos):
        ax = axes[7]
    
    if ax:
        if hasattr(ax, 'vb'):
            mousePoint = ax.vb.mapSceneToView(pos)
        
            x = mousePoint.x()
            y = mousePoint.y()
            
            # round x
            if abs(x) < 0:
                if abs(x) < 0.01:
                    if abs(x) < 0.0001:
                        n = False
                    else:
                        n=6
                else:
                    n = 4
            elif abs(x) < 10:
                n = 6
            elif abs(x) < 1000:
                n = 4
            elif abs(x) < 10000:
                n = 2
            else:
                n = 1
                
            if n:
                x = round(x, n)
                
            # round x
            if abs(y) < 0:
                if abs(y) < 0.01:
                    if abs(y) < 0.0001:
                        n = False
                    else:
                        n=6
                else:
                    n = 4
            elif abs(y) < 10:
                n = 6
            elif abs(y) < 1000:
                n = 4
            elif abs(y) < 10000:
                n = 2
            else:
                n = 1
                
            if n:
                y = round(y, n)
                
            t4.setText('x: ' + str(x))
            t5.setText('y: ' + str(y))
# %%

# ============================================================================ #
# %% ============================= DOCKS ================================== %% #
# ============================================================================ #


# %% ============================== dock 1 ================================ %% #
################
# initializing #
################
d1 = Dock('Directory selection', size=(200, 1))    
w1 = pg.LayoutWidget()
d1.addWidget(w1)
area.addDock(d1, 'left')

###########
# widgets #
###########
ICON = os.path.dirname(brixs.__file__)+'/beamlines/ipe/quickviewer/'
e1 = QtWidgets.QLineEdit()
e1.textChanged.connect(get_folderpath)
w1.addWidget(e1, row=0, col=0, colspan=3)

b1 = QtWidgets.QPushButton()
b1.setIcon(QtGui.QIcon(ICON+'back.png'))
b1.clicked.connect(go_back_directory)
w1.addWidget(b1, row=1, col=0)

b2 = QtWidgets.QPushButton()
b2.setIcon(QtGui.QIcon(ICON+'folder.png'))
b2.clicked.connect(select_folderpath)
w1.addWidget(b2, row=1, col=1)

b3 = QtWidgets.QPushButton()
b3.setIcon(QtGui.QIcon(ICON+'update.png'))
b3.clicked.connect(update_filelist)
w1.addWidget(b3, row=1, col=2)

l1 = QtWidgets.QListWidget()
l1.setSelectionMode(QtWidgets.QListWidget.ExtendedSelection)
l1.itemSelectionChanged.connect(get_selection)
w1.addWidget(l1, row=2, col=0, colspan=3)

r1 = QtWidgets.QRadioButton()
r1.setIcon(QtGui.QIcon(ICON+'update.png'))
r1.clicked.connect(auto_update)
w1.addWidget(r1, row=3, col=0, colspan=3)

# # normalization
# e2 = QtWidgets.QLineEdit()
# e2.textChanged.connect(update_plot)
# w1.addWidget(e2, row=3, col=1)

# e3 = QtWidgets.QLineEdit()
# e3.textChanged.connect(update_plot)
# w1.addWidget(e3, row=3, col=2)

# e4 = QtWidgets.QLineEdit()
# e4.textChanged.connect(update_plot)
# w1.addWidget(e4, row=3, col=3)

# e5 = QtWidgets.QLineEdit()
# e5.textChanged.connect(update_plot)
# w1.addWidget(e5, row=3, col=4)
# %%

# %% ============================== dock 2 ================================ %% #
# initializing
d2 = Dock('Complete', size=(600, 1))
area.addDock(d2, 'right')

# widgets 
axes = []
w2 = pg.GraphicsLayoutWidget(show=True)
axes.append(w2.addPlot())
axes.append(w2.addPlot())
axes.append(w2.addPlot())

w2.nextRow()
axes.append(w2.addPlot())
axes.append(w2.addPlot())
axes.append(w2.addPlot())

w2.nextRow()
axes.append(w2.addPlot())
axes.append(w2.addPlot())

d2.addWidget(w2)

# mouse mode
for ax in axes:
    ax.vb.setMouseMode(ax.vb.RectMode)

# connect axes
for i, ax in enumerate(axes):
    if i < len(axes)-1:
        axes[i].setXLink(axes[i+1])

# labels
axes[0].setLabel('left', 'TEY')
axes[1].setLabel('left', 'TFY')
axes[2].setLabel('left', 'I0')

axes[3].setLabel('left', 'TEY/I0')
axes[4].setLabel('left', 'TFY/I0')
axes[5].setLabel('left', 'I0/I0')

axes[6].setLabel('left', 'TEY/I0 norm.')
axes[7].setLabel('left', 'TFY/I0 norm.')
# axes[8].setLabel('left', 'set normalization')

# connect cursor moved
for ax in axes:
    proxy = pg.SignalProxy(ax.scene().sigMouseMoved, rateLimit=60, slot=mouse_moved)
# %%

# %% ============================== dock 3 ================================ %% #
# initializing
d3 = Dock('TEY', size=(600, 1))
area.addDock(d3, 'below', d2)

# widgets 
w3 = pg.PlotWidget()
v3 = w3.getViewBox()
v3.setMouseMode(v3.RectMode)
d3.addWidget(w3)
# %%

# %% ============================== dock 4 ================================ %% #
# initializing
d4 = Dock('TFY', size=(600, 1))
area.addDock(d4, 'below', d3)

# widgets 
w4 = pg.PlotWidget()
v4 = w4.getViewBox()
v4.setMouseMode(v4.RectMode)
d4.addWidget(w4)
# %%

# %% ============================== dock 5 ================================ %% #
# initializing
d5 = Dock('I0', size=(600, 1))
area.addDock(d5, 'below', d4)

# widgets 
w5 = pg.PlotWidget()
v5 = w5.getViewBox()
v5.setMouseMode(v5.RectMode)
d5.addWidget(w5)
# %%

# %% ============================== dock 6 ================================ %% #
# initializing
d6 = Dock('PD', size=(600, 1))
area.addDock(d6, 'below', d5)

# widgets 
w6 = pg.PlotWidget()
v6 = w6.getViewBox()
v6.setMouseMode(v6.RectMode)
d6.addWidget(w6)
# %%

# %% ============================== dock 7 ================================ %% #
# initializing
d7 = Dock('Metadata', size=(600, 1))
area.addDock(d7, 'below', d6)

# widgets 
w7 = pg.PlotWidget()
w7.setXRange(0, 20, padding=0)
w7.setYRange(0, -20, padding=0)
# v7 = w7.getViewBox()
# v7.setMouseMode(v7.RectMode)
d7.addWidget(w7)
# %%

# %% ============================== dock 8 ================================ %% #
# initializing
d8 = Dock('info', size=(100, 1))
area.addDock(d8, 'right', d7)

# widgets 
w8 = pg.PlotWidget()
w8.getPlotItem().hideAxis('bottom')
w8.getPlotItem().hideAxis('left')
legend = w8.addLegend()
d8.addWidget(w8)
# %%

# %% ============================== dock 9 ================================ %% #
################
# initializing #
################
d9 = Dock('Normalization', size=(100, 0.2))    
w9 = pg.LayoutWidget()
d9.addWidget(w9)
area.addDock(d9, 'bottom', d8)

###########
# widgets #
###########
e2 = QtWidgets.QLineEdit()
e2.textChanged.connect(get_normalization)
w9.addWidget(e2, row=0, col=1)

e3 = QtWidgets.QLineEdit()
e3.textChanged.connect(get_normalization)
w9.addWidget(e3, row=1, col=1)

e4 = QtWidgets.QLineEdit()
e4.textChanged.connect(get_normalization)
w9.addWidget(e4, row=2, col=1)

e5 = QtWidgets.QLineEdit()
e5.textChanged.connect(get_normalization)
w9.addWidget(e5, row=3, col=1)

t0 = QtWidgets.QLabel()
t0.setText('pre start')
w9.addWidget(t0, row=0, col=0)

t1 = QtWidgets.QLabel()
t1.setText('pre stop')
w9.addWidget(t1, row=1, col=0)

t2 = QtWidgets.QLabel()
t2.setText('post start')
w9.addWidget(t2, row=2, col=0)

t3 = QtWidgets.QLabel()
t3.setText('post stop')
w9.addWidget(t3, row=3, col=0)

# %%

# %% ============================== dock 10 ================================ %% #
################
# initializing #
################
d10 = Dock('Cursor', size=(100, 0.1))    
w10 = pg.LayoutWidget()
d10.addWidget(w10)
area.addDock(d10, 'bottom', d8)

###########
# widgets #
###########
t4 = QtWidgets.QLabel()
t4.setText('x: -')
w10.addWidget(t4, row=0, col=0)

t5 = QtWidgets.QLabel()
t5.setText('y: -')
w10.addWidget(t5, row=1, col=0)

# %%

# %% ================================ run ================================= %% #
if __name__ == '__main__':
    update_folderpath(default_folderpath)
    win.show()
    app.exec()
    sys.exit()
