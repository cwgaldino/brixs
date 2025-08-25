#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fit."""

# standard libraries
import copy
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import collections

# backpack
import sys
sys.path.append('/home/galdino/github/py-backpack')
# sys.path.append(r'C:\Users\carlo\github\py-backpack')
import backpack.filemanip as fm
import backpack.figmanip as figm
import backpack.arraymanip as am
from backpack.libremanip2 import soffice
import backpack.datafit
from backpack.libremanip2 import get_cell_value_from_sheets

import importlib
import backpack.libremanip2
importlib.reload(backpack.libremanip2)
from backpack.libremanip2 import soffice
importlib.reload(backpack.datafit)
importlib.reload(am)
importlib.reload(fm)
importlib.reload(figm)

%matplotlib qt5
# figm.set_default_window_position((1037, 1137))
figm.set_default_window_position((1280, 31))
# figm.getWindowPosition()


# %% ======================== ALIGN DATA =======================================

# %% ========================= filelist =======================================
exp_folder = Path('dataTidy/8_allData')
filelist_misaligned = fm.parsed_filelist(exp_folder, '*.dat', ref=0)

# Temperature list
T_list = {}
for i in filelist_misaligned:
    T_list[i] = int(filelist_misaligned[i].name.split('_')[-1].split('K')[0])

# %% exp data =================================================================
# data = fm.load_data(filelist[6], delimiter=' ')

nested_dict = lambda: collections.defaultdict(nested_dict)
data_misaligned = nested_dict()
for i in filelist_misaligned:
    data_misaligned[i] = fm.load_data(filelist_misaligned[i], delimiter=' ')
    # data_misaligned[i]['filename'] = filelist[T]
    # data_misaligned[T]['name']     = filelist[T].name
    # data_misaligned[T]['data']     = fm.load_data(filelist_misaligned[T], delimiter=' ')

# %% plot =====================================================================
fig = figm.figure()
figm.setWindowPosition()

i = 65
plt.plot(data_misaligned[i]['rshift'], data_misaligned[i]['intensity'], label=i, marker='o', ms=3)
plt.legend()
plt.xlabel('Raman shift (cm$^{-1}$)')
plt.ylabel('Intensity')
plt.tight_layout()

# %% plot second der ==========================================================
fig = figm.figure()
figm.setWindowPosition()
importlib.reload(manip)

i = 1
factor = 1e-1

x = data_misaligned[i]['rshift']
y = data_misaligned[i]['intensity']
plt.plot(x, y, label=i)

x_diff, y_diff = am.derivative(x, y, order=1)
plt.plot(x_diff, y_diff*factor, label=r'2$^{nd}$ der $\times$ '+str(factor))

x_diff, y_diff = am.derivative(x, y, order=1, window_size=2)
plt.plot(x_diff, y_diff*factor, label=r'smoothed 2$^{nd}$ der $\times$ '+str(factor))

# figm.zoom(7682.4, 7793.4)
plt.legend()
plt.xlabel('Raman shift (cm$^{-1}$)')
plt.ylabel('Intensity')

# %% support functions ===================================================
from backpack.model_functions import square

def find_square(x, y, test_width=2, plot=False):

    # i = 40
    # x = data[i]['rshift']
    # y = data[i]['intensity']
    # x, y = am.extract(x, y, [[455, 525]])
    # test_width=2

    # initial parameters
    offset = np.mean(y)
    step = np.mean(np.diff(xt))

    # square function
    xs = np.arange(x[0]-len(x)*step, xt[-1]+len(x)*step, step)
    amp = 0.5
    c = 470
    ys = square(xs, amp, c, test_width) + offset

    # cross correlation
    y_cc = np.correlate(ys, y, mode='full')
    ref_shift = (x[-1]-x[0]) + (x[0]-xs[0])
    shift = -(np.argmax(y_cc)*step - ref_shift)

    # derivative
    x_cc = np.arange(0, len(y_cc), 1)
    dx_cc, dy_cc = am.derivative(x_cc, y_cc, window_size=2)

    # trimming edges
    y_cc_2 = y_cc[np.argmax(y_cc)-50:np.argmax(y_cc)+50]
    x_cc_2 = np.arange(np.argmax(y_cc)-50, np.argmax(y_cc)+50, 1)
    dx_cc_2, dy_cc_2 = am.derivative(x_cc_2, y_cc_2, window_size=5)

    # calculating shift
    shift_0 = -(dx_cc_2[np.argmax(dy_cc_2)]*step - ref_shift)
    shift_1 = -(dx_cc_2[np.argmin(dy_cc_2)]*step - ref_shift)

    if plot:
        # plot derivative
        fig = figm.figure()
        figm.setWindowPosition()
        plt.plot(y_cc, marker='o', label='cross correlation')
        plt.plot(x_cc_2, y_cc_2)

        offset_2 = np.mean(y_cc_2)
        factor = 10
        plt.plot(dx_cc,   dy_cc*factor + offset_2-5, marker='o', label='derivative')
        plt.plot(dx_cc_2, dy_cc_2*factor + offset_2-5, marker='o')

        plt.vlines(dx_cc_2[np.argmax(dy_cc_2)], offset_2-10, offset_2+10, )
        plt.vlines(dx_cc_2[np.argmin(dy_cc_2)], offset_2-10, offset_2+10, )
        plt.legend()

        # plot data
        fig = figm.figure()
        figm.setWindowPosition()
        plt.plot(x, y, marker='o', ms=3, label='data')
        plt.vlines(c+shift_0, min(y), max(y))
        plt.vlines(c+shift_1, min(y), max(y))
        plt.legend()

    # return c+shift_1, c+shift_0
    # return c+shift_1-int(test_width/2), c+shift_0+int(test_width/2)
    return c+shift_1-(shift_0 - shift_1)/4/2, c+shift_0+(shift_0 - shift_1)/4/2



# %% test ==========================================================
i = 2
x = data_misaligned[i]['rshift']
y = data_misaligned[i]['intensity']
xt, yt = am.extract(x, y, [[455, 525]])
a, b = find_square(xt, yt, test_width=2, plot=True)

# plot data
fig = figm.figure()
figm.setWindowPosition()
plt.plot(xt, yt, marker='o', ms=3, label='data')
plt.vlines(a, min(yt), max(yt), color='red')
plt.vlines(b, min(yt), max(yt), color='green')
plt.legend()

# %% align data =====================================================
for i in data_misaligned:
    x = data_misaligned[i]['rshift']
    y = data_misaligned[i]['intensity']

    # find square
    xt, yt = am.extract(x, y, [[455, 525]])
    a, b = find_square(xt, yt, test_width=2, plot=False)
    c = (a+b)/2
    shift = 485-c

    datatemp = {'rshift': np.append(x[:am.index(x, a)], x[am.index(x, b):])+shift,
                'intensity':np.append(y[:am.index(x, a)], y[am.index(x, b):])}
    fm.save_data(datatemp, f'dataTidy/9_allData_aligned/{i:02d}_Co3O2BO3_{T_list[i]}K.dat')
    # plt.plot(datatemp['rshift'], datatemp['intensity'], marker='o', ms=3)

















# %% ========================= FIT =======================================

# %% ========================= filelist =======================================
exp_folder = Path('dataTidy/9_allData_aligned')
filelist = fm.parsed_filelist(exp_folder, '*.dat', ref=0)

# Temperature list
T_list = {}
for i in filelist:
    T_list[i] = int(filelist[i].name.split('_')[-1].split('K')[0])

# %% exp data =================================================================
nested_dict = lambda: collections.defaultdict(nested_dict)
data = nested_dict()
for i in filelist:
    data[i] = fm.load_data(filelist[i])

# %% plot =====================================================================
fig = figm.figure()
figm.setWindowPosition()

i = 65
plt.plot(data[i]['rshift'], data[i]['intensity'], marker='o', ms=3, label=i)
plt.legend()
plt.xlabel('Raman shift (cm$^{-1}$)')
plt.ylabel('Intensity')
plt.tight_layout()


# %% connect to calc ==========================================================
import sys
sys.path.append('/home/galdino/github/py-backpack')
import importlib
import backpack.libremanip2
importlib.reload(backpack.libremanip2)
from backpack.libremanip2 import soffice
import backpack.datafit
importlib.reload(backpack.datafit)


# filepath = 'fit_parameters.ods'
# try:
#     libreoffice.terminate(ask=False)
# except: pass
libreoffice = soffice(norestore=True)
# libreoffice.send('print(local)', True)
# calcObject = libreoffice.openCalc(filepath=filepath)
calcObject = libreoffice.openCalc()
# calcObject.filepath
# %% define background functions ===============================================
bkg1 = lambda x, d0, d1: d0*x + d1
bkg2 = lambda x, d0: d0 + 0*x

# %% check fit (NEW) ======
to_exclude = [86, 84, 81, 70, 69,  67, 63, 61,58,  57, 55, 53, 52, 50, 49, 48, 47, 45, 44, 43, 42, 41, 40, 39, 37, 36, 35, 34, 32, 31, 30, 29, 27, 26, 25, 24, 23, 22] + \
                [1, 2, 3, 7, 8, 9, 10, 11, 13, 14, 15, 16 , 17,19,20,21, 4, 6]+\
                [84, 79, 78, 77, 28, 38, 33, 59, 62, 66, 68, 73]
to_exclude = list(dict.fromkeys(to_exclude))  # remove duplicates
to_exclude.sort()  # sort
to_exclude_str = [str(x) for x in to_exclude]
len(to_exclude_str)
names = calcObject.get_sheets_name()[1:]
len(names)

names_filtered = [name for name in names if name.split('_')[0] not in to_exclude_str]
len(names_filtered)

sheetObject_list = calcObject.get_sheets(names_filtered)
len(sheetObject_list)
T_s = [int(name.split('_')[1]) for name in names_filtered]
i_s = [int(name.split('_')[0]) for name in names_filtered]


fig = figm.figure()
plt.errorbar(T_s, get_cell_value_from_sheets(sheetObject_list, row=24, col='I', format='number'),
              get_cell_value_from_sheets(sheetObject_list, row=24, col='J', format='number'), marker='o')
# fig = figm.figure()
# plt.plot(i_s, get_cell_value_from_sheets(sheetObject_list, row=24, col='I', format='number'), marker='o')

# %% RESULTS ===================================================================
data = dict(i=i_s,
            T=T_s,
            amp=get_cell_value_from_sheets(sheetObject_list, row=2, col='I', format='number'),
            error_amp=get_cell_value_from_sheets(sheetObject_list, row=2, col='J', format='number'),
            c=get_cell_value_from_sheets(sheetObject_list, row=3, col='I', format='number'),
            error_c=get_cell_value_from_sheets(sheetObject_list, row=3, col='J', format='number'),
            w=get_cell_value_from_sheets(sheetObject_list, row=4, col='I', format='number'),
            error_w=get_cell_value_from_sheets(sheetObject_list, row=4, col='J', format='number'),
            )
fm.save_data(data, 'fit_results/fwhmVoigt#0.txt', data_format='%.3e')

data = dict(i=i_s,
            T=T_s,
            amp=get_cell_value_from_sheets(sheetObject_list, row=18, col='I', format='number'),
            error_amp=get_cell_value_from_sheets(sheetObject_list, row=18, col='J', format='number'),
            c=get_cell_value_from_sheets(sheetObject_list, row=19, col='I', format='number'),
            error_c=get_cell_value_from_sheets(sheetObject_list, row=19, col='J', format='number'),
            w=get_cell_value_from_sheets(sheetObject_list, row=20, col='I', format='number'),
            error_w=get_cell_value_from_sheets(sheetObject_list, row=20, col='J', format='number'),
            )
fm.save_data(data, 'fit_results/fwhmVoigt#3.txt', data_format='%.3e')

data = dict(i=i_s,
            T=T_s,
            amp=get_cell_value_from_sheets(sheetObject_list, row=22, col='I', format='number'),
            error_amp=get_cell_value_from_sheets(sheetObject_list, row=22, col='J', format='number'),
            c=get_cell_value_from_sheets(sheetObject_list, row=23, col='I', format='number'),
            error_c=get_cell_value_from_sheets(sheetObject_list, row=23, col='J', format='number'),
            w=get_cell_value_from_sheets(sheetObject_list, row=24, col='I', format='number'),
            error_w=get_cell_value_from_sheets(sheetObject_list, row=24, col='J', format='number'),
            )
fm.save_data(data, 'fit_results/fwhmVoigt#4.txt', data_format='%.3e')

data = dict(i=i_s,
            T=T_s,
            amp=get_cell_value_from_sheets(sheetObject_list, row=26, col='I', format='number'),
            error_amp=get_cell_value_from_sheets(sheetObject_list, row=26, col='J', format='number'),
            c=get_cell_value_from_sheets(sheetObject_list, row=27, col='I', format='number'),
            error_c=get_cell_value_from_sheets(sheetObject_list, row=27, col='J', format='number'),
            w=get_cell_value_from_sheets(sheetObject_list, row=28, col='I', format='number'),
            error_w=get_cell_value_from_sheets(sheetObject_list, row=28, col='J', format='number'),
            )
fm.save_data(data, 'fit_results/fwhmVoigt#5.txt', data_format='%.3e')

# %%
i = 0
sheet = sheetObject_list[i_s.index(i)]
x = data[i]['rshift']
y = data[i]['intensity']
sheet.update_model()
sheet.model(x, *sheet.p_fitted)
# %% Fit: first try ============================================================
# sheet = calcObject.get_sheets(1)
# calcObject.get_sheets_name()
# calcObject.save()
# sheet.get_name()
# sheet.name
# sheet.set_name('Sheet1')
# sheet.get_last_row()
# sheet.get_last_col()
# libreoffice.send("print(sheet.getRowDescriptions()[-1].split(' ')[-1])", True)
# sheet.get_cell_value(row=19, col='H', format='formula')
# sheet.get_cell_value(row=19, col='H', format='string')
# sheet.get_cells_value(row_start=2, col_start='H',row_stop=10, col_stop='J', format='number')
# libreoffice.send('print(2+2)', True)
# write(libreoffice.python, 'print(3+2)')
i = 82
sheet = sheetObject_list[i_s.index(i)]
x = data[i]['rshift']
y = data[i]['intensity']

x2fit, y2fit = am.extract(x, y, [[300, 700]])#476]], [486, 500]])
# ties = [[426, 434.57, 12], [523.8, 538.8, 12], [555.2, 565.4, 12], [609.4, 622, 12], [661.4, 673.2, 12]]
ties = [[420, 440.59, 12], [449.82, 462.8, 12], [555.2, 565.4, 12], [609.4, 622, 12], [657.76, 673.2, 12]]
# ties = None
sheet.fit(x2fit, y2fit, ties=ties)

# % plot
sheet = sheetObject_list[i_s.index(i)]
x = data[i]['rshift']
y = data[i]['intensity']
ax = sheet.plot_fit(x, y, ax=None, show_exp=True, show_derivative=False, show_submodels=True, smoothing=10, ties=ties, submodels_bkg='bkg2', derivative_order=1, derivative_offset=None, derivative_factor=None, derivative_window_size=1)

try:
    plt.title((sheet.name, np.round(sheet.residue, 3)))
except:
    plt.title(sheet.name)
figm.zoom(290, 700)
plt.xlabel('Raman shift (cm$^{-1}$)')
plt.ylabel('Intensity')

figm.zoom(0, 1000)



# %% Delete all sheets but the first ===========================================
# r = calcObject.get_sheets_name()
# calcObject.remove_sheets(r[1:])
# calcObject.save()

# %% create one sheet for each temperature =====================================
# calcObject.insert_sheets([f'{key}_{T_list[key]}' for key in T_list])

# %% copy values ===============================================================
sheet, = calcObject.get_sheets(1)
values = sheet.get_cells_value()

# values
sheet_names = calcObject.get_sheets_name()
for s in sheet_names[1:]:
    print(s)
    sheet2fill, = calcObject.get_sheets_by_name(s)
    sheet2fill.set_cells_value(data=values)

# %% copy formatting ===========================================================
sheet, = calcObject.get_sheets(1)
values = sheet.get_cells_value()

# formating
sheet_names = calcObject.get_sheets_name()
for s in sheet_names[1:]:
    formating = sheet.get_cells_formatting()
    sheet2fill, = calcObject.get_sheets_by_name(s)
    sheet2fill.set_cells_formatting(cell_formatting_list=formating)

# col Width
sheet_names = calcObject.get_sheets_name()
for s in sheet_names[1:]:
    sheet2fill, = calcObject.get_sheets_by_name(s)
    for col in range(1, sheet.get_last_col()):
        sheet2fill.set_col_width(width=sheet.get_col_width(col)[0], col=col)

# %% copy merged ===============================================================
sheet, = calcObject.get_sheets(1)
toMerge = sheet.get_merged()

# formating
sheet_names = calcObject.get_sheets_name()
for s in sheet_names[1:3]:
    sheet2fill, = calcObject.get_sheets_by_name(s)
    for range in ranges:
        sheet2fill.merge(*toMerge)









# %% Fit many ==================================================================
sheet_list = calcObject.get_sheets_name()[1:]
# sheets = {int(key.split('_')[0]): next(calcObject.get_sheets_by_name(key)) for key in sheet_list}
sheets = [next(calcObject.get_sheets_by_name(key)) for key in sheet_list]

# for i in range(0, 3):
for i in range(67, 72):#len(sheets)):
    print(i)

    # exp data
    x = data[i]['rshift']
    y = data[i]['intensity']

    # fit parameters
    x2fit, y2fit = am.extract(x, y, [[300, 700]])#476]], [486, 500]])
    ties = None#[[400, 451, 2]]

    # fit
    sheet = sheets[i]
    sheet.fit(x2fit, y2fit, ties=ties, save=False)


# %% plot ==========================================================
# plt.close('all')
i = 75

sheet = sheets[i]
x = data[i]['rshift']
y = data[i]['intensity']

ax = sheet.plot_fit(x, y, ax=None, show_exp=True, show_derivative=False, show_submodels=True, smoothing=10, ties=None, submodels_bkg='bkg2', derivative_order=1, derivative_offset=None, derivative_factor=None, derivative_window_size=1)
figm.setWindowPosition()

try:
    plt.title(f'{i, np.round(sheet.residue, 3)}')
except AttributeError:
    plt.title(f'{i}')
figm.zoom(290, 700)
plt.xlabel('Raman shift (cm$^{-1}$)')
plt.ylabel('Intensity')

# %%
to_exclude = [86, 84, 81, 70, 69,  67, 63, 61,58,  57, 55, 53, 52, 50, 49, 48, 47, 45, 44, 43, 42, 41, 40, 39, 37, 36, 35, 34, 32, 31, 30, 29, 27, 26, 25, 24, 23, 22] + \
                [1, 2, 3, 9, 7, 8, 9, 10, 11, 13, 14, 15, 16 , 17,]+\
                [84, 79, 78, 77]


# %% fit individual ============================================================
i = 38

# exp data
x = data[i]['rshift']
y = data[i]['intensity']

# fit parameters
x2fit, y2fit = am.extract(x, y, [[300, 700]])#476]], [486, 500]])
ties = None#[[400, 451, 2]]

# fit
sheet = sheets[i]
sheet.fit(x2fit, y2fit, ties=ties, save=False)

# plot
ax = sheet.plot_fit(x, y, ax=None, show_exp=True, show_derivative=False, show_submodels=True, smoothing=10, ties=None, submodels_bkg='bkg2', derivative_order=1, derivative_offset=None, derivative_factor=None, derivative_window_size=1)
figm.setWindowPosition()

plt.title(f'{i, np.round(sheet.residue, 3)}')
figm.zoom(250, 700)
plt.xlabel('Raman shift (cm$^{-1}$)')
plt.ylabel('Intensity')

# %% adjust parameter ==========================================================
value = 6
for i in range(0, len(sheets)):
    sheets[i].set_cell_value(value, 16, 'F')


# %% Residue ===================================================================
residue = []
for i in range(len(sheets)):
    name =  sheets[i].object.getName()
    T = int(name.split('_')[-1])
    residue.append((i, T, sheets[i].residue))
residue = np.array(residue)

fig = figm.figure()
figm.setWindowPosition()
plt.plot(residue[:, 1], residue[:, 2], marker='o')
plt.title(np.mean(residue[:, 2]))
plt.xlabel('Temperature (K)')
plt.ylabel('Residue')

# %% result ====================================================================
# sheet_list = calcObject.get_sheets_name()[1:]
# sheets = [next(calcObject.get_sheets_by_name(key)) for key in sheet_list]
#
# for  sheet in sheets:
#     sheet.update_model()
#     sheet.get_parameters()

submodel = 'fwhmVoigt#0'
p = 'amp'
result = []
for i in range(len(sheets)):
    if i not in to_exclude: # 73
        name =  sheets[i].object.getName()
        try:
            if sheets[i].parameters[submodel][p]['fitted'][0] != '':
                T = int(name.split('_')[-1])
                result.append((i, T, sheets[i].parameters[submodel][p]['fitted'][0], sheets[i].parameters[submodel][p]['error'][0]))
        except KeyError: pass
result = np.array(result)

# plot
fig = figm.figure()
figm.setWindowPosition()
# plt.plot(result[:, 1], result[:, 2], marker='o')
plt.errorbar(result[:, 1], result[:, 2], yerr=result[:, 3]*2, linewidth=0.3, marker='o')
plt.title(submodel)
plt.xlabel('Temperature (K)')
plt.ylabel(p)


# %% parameter as function of fit ==============================================
for submodel in ['fwhmVoigt#0', 'fwhmVoigt#3', 'fwhmVoigt#4', 'fwhmVoigt#5']:
    result = []
    for i in range(len(sheets)):
        if i not in to_exclude: # 73
            name =  sheets[i].get_name()
            T = int(name.split('_')[-1])
            try:
                temp = []
                for p in ['amp', 'c', 'w']:
                    # if sheets[i].parameters[submodel][p]['fitted'][0] != '':
                    temp.append(sheets[i].parameters[submodel][p]['fitted'][0])
                    temp.append(sheets[i].parameters[submodel][p]['error'][0])
                result.append(np.hstack((i, T, temp)))
            except KeyError: pass
    result = np.array(result)
    fm.save_data(result, filepath=f'fit_results/{submodel}.txt', header='i, T, amp, error_amp, c, error_c, w, error_w', data_format=['%d', '%d', '%f', '%f', '%f',  '%f', '%f', '%f'])

# %% save fitted curves ========================================================
sheet_list = calcObject.get_sheets_name()[1:]
sheets = calcObject.get_sheets(sheet_list)
#
for  sheet in sheets:
    sheet.update_model()
    sheet.get_parameters()

libreoffice = soffice(norestore=True)
calcObject = libreoffice.openCalc()

for i in range(48, len(sheets)):
    # i=0
    sheets[i].update_model()
    # sheet.get_parameters()
    print(i)
    name =  sheets[i].get_name()
    T = int(name.split('_')[-1])

    # exp data
    x = data[i]['rshift']

    temp = {'rshift':[], 'intensity':[]}
    temp['rshift'] = np.linspace(min(x), max(x), len(x))
    temp['intensity'] = sheets[i].model(temp['rshift'], *sheets[i].p_fitted)
    fm.save_data(temp, filepath=f'fit_results/{i:02d}_Co3O2BO3_{T}K.dat', data_format='%f')
