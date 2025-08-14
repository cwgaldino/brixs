#! /usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fit.

TODO:
-method to plot parameter for various temperatures
-add fit button to libreoffice
-fit various data (like first and second der)


New atributes added to sheet object:
parameters
var_string
model_string
p_min
p_max
p_guess
p_fitted
p_error
linked_parameters
id_list
submodel
residue
p_cov

new methods:
get_parameters
update_model
update_submodels
fit

functions:
fake_sigma
"""

# standard libraries
import copy
import numpy as np
from pathlib import Path
import inspect
import importlib

# matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# fit
from scipy.optimize import curve_fit
from scipy.integrate import trapz

# model functions
from model_functions import *

import sys
sys.path.append('../../')
import libreoffice_wrapper as lw

plt.ion()

def index(x, value):
    """Returns the index of the element in array which is closest to value.

    Args:
        x (list or array): 1D array.
        value (float or int): value.

    Returns:
        index (int)
    """
    return np.argmin(np.abs(np.array(x)-value))

# %%

import importlib
importlib.reload(lw)

sys.path.append('/home/galdino/github/py-backpack')
import backpack.figmanip as figm
import backpack.filemanip as fm
import backpack.arraymanip as am
import time

pid = lw.start_soffice()
time.sleep(5)
soffice = lw.soffice()
calc = soffice.Calc('fit.ods')
calc.save()



# %%
def load_data(self, x, y, sigma=None):
    if len(x) == len(y):
        self.set_column(column='P', value=x, row_start='2')
        self.set_column(column='Q', value=y, row_start='2')
    else:
        raise ValueError('x and y must have the same length')
    if sigma is not None:
        if len(x) == len(sigma):
            self.set_column(column='R', value=sigma, row_start='2')
        else:
            raise ValueError('x and sigma must have the same length')
    else:
        sigma = [1]*len(x)
        self.set_column(column='R', value=sigma, row_start='2')

def get_x(self):
    return self.get_column('P', row_start='2')

def get_y(self):
    return self.get_column('Q', row_start='2')

def get_sigma(self):
    return self.get_column('R', row_start='2')

def get_data(self):
    length = self.get_column_length('P')
    return np.array(self.get_value('P', '2', 'R', length-1, format='number'))

def set_range2fit(self, value=None, global_sigma=None):
    """Build a sigma array which determines the uncertainty in ydata.

    Adaptaded from the `scipy.optimize.curve_fit() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html>`_ documentation:

        If we define residuals as ``r = ydata - model(xdata, *popt)``, then sigma
        for a 1-D data should contain values of standard deviations of errors in
        ydata. In this case, the optimized function is ``chisq = sum((r / sigma) ** 2)``.

    Args:
        x (list): x array.
        sigma (float, optional): sigma value to be used for all points in ``x``.
        sigma_specific (list, optional): list of triples specfing new sigma for specific ranges, e.g.,
            ``sigma_specific = [[x_init, x_final, sigma], [x_init2, x_final2, sigma2], ]``.

    Returns:
        array.
    """
    if global_sigma is not None:
        self.set_value('N1', value=global_sigma, format='number')

    length = self.get_column_length('L')
    if length > 6:
        self.clear('L', '7', 'N', length-1)

    if value is not None:
        try:
            for i, sigma in enumerate(value):
                if len(sigma) == 2:
                    sigma.append(1)
                if len(sigma) != 3:
                    raise ValueError(f'{sigma} is not a valid data range.')
            self.set_value(row='7', value=value, column='L')
        except TypeError:
            if len(value) == 2:
                value.append(1)
            if len(value) != 3:
                raise ValueError(f'{value} is not a valid data range.')
            self.set_value(row='7', value=[value], column='L')

def get_range2fit(self):
    length = self.get_column_length('L')
    if length == 6:
        return None
    else:
        return self.get_value('L', '7', 'N', length-1)

def set_global_sigma(self, value):
    self.set_value('N3', value=value)

def get_global_sigma(self):
    return self.get_value('N3', format='number')

def get_submodels(self):
    length = self.get_column_length('B')
    data = self.get_value('A', 1, 'I', length-1)

    submodels = {}
    for i, row in enumerate(data):
        name, arg, use, vary, bound_min, guess, bound_max, fitted, error = row[0:9]

        if use == 'y':
            if name not in submodels:
                submodels[name] = dict()
            if arg not in submodels[name]:
                submodels[name][arg] = dict()

            submodels[name][arg]['use'] = use
            submodels[name][arg]['vary'] = vary

            if bound_min != '':
                submodels[name][arg]['min'] = float(bound_min)
            else:
                submodels[name][arg]['min'] = -np.inf

            if bound_min != '':
                submodels[name][arg]['max'] = float(bound_max)
            else:
                submodels[name][arg]['max'] = np.inf

            try:
                submodels[name][arg]['guess'] = float(guess)
            except ValueError:
                submodels[name][arg]['guess'] = guess

            try:
                submodels[name][arg]['fitted'] = float(fitted)
                submodels[name][arg]['error'] = float(error)
            except ValueError:
                submodels[name][arg]['fitted'] = [0]*length
                submodels[name][arg]['error'] = [0]*length
            submodels[name][arg]['row'] = i

    return submodels

def get_model(self):
    submodels = self.get_submodels(self)
    var_string = ''
    model_string = ''
    p_guess = []
    p_min = []
    p_max = []
    p = 0
    for name in submodels:
        model_string += f"{name.split('#')[0]}(x, "
        for arg in submodels[name]:
            if submodels[name][arg]['vary'] == 'y':
                if type(submodels[name][arg]['guess']) is str:
                    name2 = submodels[name][arg]['guess'].split('_')[0]
                    arg2 = submodels[name][arg]['guess'].split('_')[1]
                    p2 = submodels[name2][arg2]['p']
                    submodels[name][arg]['p'] = p2
                    model_string += f'{arg}=p{p2}, '
                else:
                    var_string   += f'p{p}, '
                    model_string += f'{arg}=p{p}, '
                    p_guess.append(submodels[name][arg]['guess'])
                    p_min.append(submodels[name][arg]['min'])
                    p_max.append(submodels[name][arg]['max'])
                    submodels[name][arg]['p'] = p
                    p += 1
            else:
                value = submodels[name][arg]['guess']
                model_string += f'{arg}={value}, '
        model_string += ') + '
    model_string = f'lambda x, {var_string[:-2]}: {model_string[:-3]}'
    model = eval(model_string)
    return model, model_string, submodels, p_guess, p_min, p_max

def fit(self, save_final=True, save_partial_guess=False, save_partial_fit=False, print_time=True):
    if print_time: start_time = time.time()
    length = self.get_column_length('A')

    # data
    data = self.get_data(self)
    r_raw = self.get_range2fit(self)
    if r_raw is not None:
        r = [(x, y) for x, y, s in r_raw]
        x, y     = am.extract(data[:, 0], data[:, 1], ranges=r)
        _, sigma = am.extract(data[:, 0], data[:, 2], ranges=r)
    else:
        x = data[:, 0]
        y = data[:, 1]
        sigma = data[:, 2]

    # remove zero sigmas
    x2fit = np.array([x1 for x1, s in zip(x, sigma) if s != 0])
    y2fit = np.array([y1 for y1, s in zip(y, sigma) if s != 0])
    sigma = np.array([s for s in sigma if s != 0])

    # sigma
    global_sigma = self.get_global_sigma(self)
    sigma = global_sigma*sigma
    if global_sigma != 0:
        for s in r_raw:
            init = index(x, s[0])
            final = index(x, s[1])
            sigma[init:final+1] = sigma[init:final+1]*s[2]

    # get model, guess, bounds
    model, _, submodels, p_guess, p_min, p_max = self.get_model(self)
    y_guess = model(x2fit, *p_guess)

    # fit
    if print_time: start_time2 = time.time()
    p_fitted, p_cov = curve_fit(model, x2fit, y2fit, p_guess, sigma=sigma, bounds=[p_min, p_max])
    p_error = np.sqrt(np.diag(p_cov))  # One standard deviation errors on the parameters
    if print_time: elapsed_time2 = time.time() - start_time2

    # fitted values back to the spreadsheet
    p_fitted2 = [['', '']]*length
    # j = 0
    for name in submodels:
        for arg in submodels[name]:
            if submodels[name][arg]['vary'] == 'y':
                if type(submodels[name][arg]['guess']) is str:
                    name2 = submodels[name][arg]['guess'].split('_')[0]
                    arg2 = submodels[name][arg]['guess'].split('_')[1]
                    p2 = submodels[name2][arg2]['p']
                    p_fitted2[submodels[name][arg]['row']] = [p_fitted[p2], p_error[p2]]
                else:
                    p = submodels[name][arg]['p']
                    p_fitted2[submodels[name][arg]['row']] = [p_fitted[p], p_error[p]]
                    # j += 1
            else:
                p_fitted2[submodels[name][arg]['row']] = [submodels[name][arg]['guess'], 0]
    self.set_value('H', '2', 'I', length, value=p_fitted2)
    y_fitted = model(x2fit, *p_fitted)

    # get residue
    residue = trapz(abs(y2fit - model(x2fit, *p_fitted)), x2fit)
    self.set_value('N1', residue)

    # submodels fitted curves
    sub_fitted = []
    for name in submodels:
        string = 'lambda x: '+ name.split('#')[0] + '(x, '
        for arg in submodels[name]:
            p = p = submodels[name2][arg2]['p']
            value = p_fitted[p]
            string += arg + '=' + str(value) + ', '
        string += ')'
        s = eval(string)
        sub_fitted.append(s(x2fit))

    # submodels guessed curves
    sub_guess = []
    for name in submodels:
        string = 'lambda x: '+ name.split('#')[0] + '(x, '
        for arg in submodels[name]:
            p = p = submodels[name2][arg2]['p']
            value = p_guess[p]
            string += arg + '=' + str(value) + ', '
        string += ')'
        s = eval(string)
        sub_guess.append(s(x2fit))

    # clear old curves
    l = self.get_row_length('1')
    if l < 19: l = 19
    self.clear('T', '1', l, self.get_column_length('T'))

    # curves
    data_final = []
    data_final.append(x2fit)
    data_final.append(y2fit)
    data_final.append(sigma)
    data_final.append(y_guess)
    data_final.append(y_fitted)

    partial_guess = []
    partial_guess.append(x2fit)
    for s in sub_guess:
        partial_guess.append(s)

    partial_fitted = []
    partial_fitted.append(x2fit)
    for s in sub_fitted:
        partial_fitted.append(s)

    # save curves
    column = 19
    if save_final or save_partial_guess or save_partial_fit:
        self.set_value(column, '1', column+5-1, '1', value=[['x2fit', 'y2fit', 'sigma', 'y guess', 'y fit']])
        data2save = copy.copy(data_final)

    column = 24
    if save_partial_guess:
        self.set_value(column, '1', column+len(submodels)-1, '1', value=[[key+'_guess' for key in submodels.keys()]])
        column += len(partial_guess)-1
        for d in partial_guess[1:]:
            data2save.append(d)

    if save_partial_fit:
        self.set_value(column, '1', column+len(submodels)-1, '1', value=[[key+'_fit' for key in submodels.keys()]])
        for d in partial_fitted[1:]:
            data2save.append(d)

    if print_time: start_time3 = time.time()
    if save_final or save_partial_guess or save_partial_fit:
        self.set_value(19, '2', 19+len(data2save)-1, len(data2save[0]), value=lw.transpose(data2save))
    if print_time: elapsed_time3 = time.time() - start_time3

    if print_time:
        elapsed_time = time.time() - start_time
        print(f'fit elapsed time: {elapsed_time2}')
        print(f'save data elapsed time: {elapsed_time3}')
        print(f'other elapsed time: {elapsed_time-elapsed_time2-elapsed_time3}')
        print(f'total elapsed time: {elapsed_time}')
    return model, p_fitted, p_error, np.array(data_final), np.array(partial_guess), np.array(partial_fitted)


# %%
def refresh():
    importlib.reload(sys.modules[__name__])

sheet.load_data = load_data
sheet.get_x = get_x
sheet.get_y = get_y
sheet.get_sigma = get_sigma
sheet.get_data = get_data
sheet.set_range2fit = set_range2fit
sheet.get_range2fit = get_range2fit
sheet.set_global_sigma = set_global_sigma
sheet.get_global_sigma = get_global_sigma
sheet.fit = fit
sheet.get_submodels = get_submodels
sheet.get_model = get_model

bkg2 = lambda x, d0: d0 + 0*x

# sheet.set_global_sigma(sheet, value=10**-10)
# sheet.set_range2fit(sheet, value=[300, 700, 1])
# sheet.set_range2fit(sheet)
# sheet.get_submodels(sheet)
# sheet.get_model(sheet)
# model, p_fitted, p_error, data_final, partial_guess, partial_fitted = sheet.fit(sheet, save_final=False)
# model, p_fitted, p_error, data_final, partial_guess, partial_fitted = sheet.fit(sheet)
# model, p_fitted, p_error, data_final, partial_guess, partial_fitted = sheet.fit(sheet, save_partial_guess=True)
# model, p_fitted, p_error, data_final, partial_guess, partial_fitted = sheet.fit(sheet, save_partial_fit=True)
model, p_fitted, p_error, data_final, partial_guess, partial_fitted = sheet.fit(sheet, save_partial_guess=True, save_partial_fit=True)





# %%
start_time = time.time()
f = sheet.get_value('A1:AF6000')
elapsed_time = time.time() - start_time
print(elapsed_time)

start_time = time.time()
f = sheet.get_value('A1:AF1')
f = sheet.get_value('A1:AF1')
f = sheet.get_value('A1:AF1')
elapsed_time = time.time() - start_time
print(elapsed_time)



# %%
importlib.reload(lw)
pid = lw.start_soffice()
time.sleep(5)
soffice = lw.soffice()

calc = soffice.Calc('fit.ods')
sheet = calc.get_sheet_by_position(1)
sheet.load_data = load_data

data = fm.load_data('example/data/00_Co3O2BO3_15K.dat')
x, y = data['rshift'], data['intensity']
sheet.load_data(sheet, x, y)

calc = soffice.Calc('fit.ods')
sheet = calc.get_sheet_by_position(1)
sheet.set_sigma = set_sigma
sheet.set_sigma(sheet, global_sigma=10**-12, sigma_specific=[[14, 16, 2]])



# x = data['rshift']
# x = np.array([1,2,3])
type(x)
sheet.set_column(column='R', value=x[0:373], row_start='3')
sheet.set_column(column='R', value=x, row_start='3')
sheet.set_row(row='13', value=x[0:500])



x2 = lw.transpose(x[0:374])
sheet.set_value(value=x2, column_start=17, column_stop=17, row_start=2, row_stop=2 + len(x2)- 1, format='formula')




# %%
last_col = 'L'
header_row = 1
hashtag_col = 1

exp_def     = dict(linewidth=2, markersize=8, color='black')
guess_def   = dict(linewidth=2, linestyle='--', color='green')
fit_def     = dict(linewidth=2, color='red')
ties_def    = dict(markersize=8, color='orange')
der_def     = dict(linewidth=0, marker='o', color='black')
der_fit_def = dict(linewidth=2, color='red')
der_guess_def = dict(linewidth=2, linestyle='--', color='green')
sub_def     = dict(linewidth=1)


def update_model(self):
    refresh()

    var_string = ''
    model_string = ''
    self.p_min = []
    self.p_max = []
    self.p_guess = []
    self.p_fitted = []
    self.p_error = []

    # get parameters
    self.get_parameters()

    # get header and submodels
    header = self.get_row_values(header_row, col_stop=last_col)
    guess_col = header.index('guess')+1
    fitted_col = header.index('fitted')+1
    error_col = header.index('error')+1
    id_col = header.index('id')+1
    self.set_col_values(data=['' for i in range(self.get_last_row()-1)], row_start=header_row+1, col=id_col)

    p = 0
    x = 0
    self.linked_parameters = {}
    for submodel in self.parameters:

        # check if this submodel should be used
        if 'y' in [use for sublist in [self.parameters[submodel][arg]['use'] for arg in self.parameters[submodel]] for use in sublist]:

            # get tag
            try:
                submodel_tag = submodel.split('#')[-1]
            except IndexError:
                submodel_tag = None
            submodel_name = submodel.split('#')[0]

            # get arguments from function
            args_expected = list(inspect.signature(eval(submodel_name)).parameters)

            # initialize model
            model_string += f"{submodel_name}(x, "

            # build min, max, guess, model
            for arg in args_expected[1: ]:
                # print(submodel, arg)

                # check if submodel has active argument
                missing_arg = False
                try:
                    to_use = self.parameters[submodel][arg]['use'].index('y')
                except ValueError:
                    missing_arg = True
                if missing_arg: raise MissingArgument(submodel, arg)

                # check if parameter must vary =========================================================
                vary = list(self.parameters[submodel][arg]['vary'])[to_use]
                hashtag = list(self.parameters[submodel][arg]['#'])[to_use]
                # linked parameter ===================================
                if vary != 'y' and vary != 'n':
                    submodel2link = vary.split(',')[0]
                    arg2link = vary.split(',')[-1]
                    if submodel2link in self.parameters and arg2link in self.parameters[submodel]:
                        to_use_linked = self.parameters[submodel2link][arg2link]['use'].index('y')
                        vary2 = list(self.parameters[submodel2link][arg2link]['vary'])[to_use_linked]
                    else:
                        raise ValueError(f"Cannot find submodel '{submodel2link}' with arg '{arg2link}'.")
                    while vary2 != 'y' and vary2 != 'n':
                        # print(vary2)
                        submodel2link = vary2.split(',')[0]
                        arg2link = vary2.split(',')[-1]
                        # print(submodel2link, arg2link)
                        if submodel2link in self.parameters and arg2link in self.parameters[submodel]:
                            to_use_linked = self.parameters[submodel2link][arg2link]['use'].index('y')
                            vary2 = list(self.parameters[submodel2link][arg2link]['vary'])[to_use_linked]
                        else:
                            raise ValueError(f"Cannot find submodel '{submodel2link}' with arg '{arg2link}'.")

                    if vary2 == 'n': #
                        v = list(self.parameters[submodel2link][arg2link]['guess'])[to_use_linked]
                        model_string += f'{v}, '
                        self.set_cell_value(value='-', row=hashtag+header_row+1, col=id_col)
                        self.set_cell_value(value=v, row=hashtag+header_row+1, col=guess_col)
                        self.set_cell_value(value=v, row=hashtag+header_row+1, col=fitted_col)
                        self.set_cell_value(value=0, row=hashtag+header_row+1, col=error_col)
                    else:
                        if submodel2link+','+arg2link in self.linked_parameters:
                            x_temp = self.linked_parameters[submodel2link+','+arg2link]
                            self.set_cell_value(value='x' + str(x_temp), row=hashtag+header_row+1, col=id_col)
                            self.parameters[submodel][arg]['id'][to_use] = 'x' + str(x_temp)
                            # var_string += f'x{x_temp}, '
                            model_string += f'x{x_temp}, '
                        else:
                            self.linked_parameters[submodel2link+','+arg2link] = x
                            self.set_cell_value(value='x' + str(x), row=hashtag+header_row+1, col=id_col)
                            self.parameters[submodel][arg]['id'][to_use] = 'x' + str(x)
                            # var_string += f'x{x}, '
                            model_string += f'x{x}, '
                            x += 1


                # fixed parameter ================================
                elif vary == 'n':
                    v = list(self.parameters[submodel][arg]['guess'])[to_use]
                    model_string += f'{v}, '
                    self.set_cell_value(value='-', row=hashtag+header_row+1, col=id_col)
                    self.parameters[submodel][arg]['id'][to_use] = '-'
                    self.set_cell_value(value=v, row=hashtag+header_row+1, col=fitted_col)
                    self.parameters[submodel][arg]['fitted'][to_use] = v
                    self.set_cell_value(value=0, row=hashtag+header_row+1, col=error_col)
                    self.parameters[submodel][arg]['error'][to_use] = 0


                # variable parameter =============================
                else:
                    self.p_min.append(list(self.parameters[submodel][arg]['min'])[to_use])
                    self.p_max.append(list(self.parameters[submodel][arg]['max'])[to_use])
                    self.p_guess.append(list(self.parameters[submodel][arg]['guess'])[to_use])
                    self.p_fitted.append(list(self.parameters[submodel][arg]['fitted'])[to_use])
                    self.p_error.append(list(self.parameters[submodel][arg]['error'])[to_use])

                    try:
                        if submodel+','+arg in self.linked_parameters:
                            x_temp = self.linked_parameters[submodel+','+arg]
                            self.set_cell_value(value='x' + str(x_temp), row=hashtag+header_row+1, col=id_col)
                            self.parameters[submodel][arg]['id'][to_use] = 'x' + str(x_temp)
                            var_string += f'x{x_temp}, '
                            model_string += f'x{x_temp}, '
                        else:
                            self.set_cell_value(value='p' + str(p), row=hashtag+header_row+1, col=id_col)
                            self.parameters[submodel][arg]['id'][to_use] = 'p' + str(p)
                            var_string += f'p{p}, '
                            model_string += f'p{p}, '
                            p += 1
                    except UnboundLocalError:
                        var_string += f'p{p}, '
                        model_string += f'p{p}, '
                        p += 1

            model_string += ') + '

    # finish model
    self.id_list = [s.strip() for s in eval('["' + var_string[:-2].replace(',', '","') + '"]')]

    self.model_string = f'lambda x, {var_string[:-2]}: {model_string[:-3]}'
    self.model = eval(self.model_string)

    # check guess, min, max ============================
    if '' in self.p_guess:
        guess_missing = [self.id_list[i] for i, x in enumerate(self.p_guess) if x == '']
        raise ValueError(f'Parameters with id {guess_missing} do not have a guess value.')

    if '' in self.p_min:
        self.p_min = [-np.inf if x == '' else x for x in self.p_min]
    if '' in self.p_max:
        self.p_max = [np.inf if x == '' else x for x in self.p_max]
    if '' in self.p_fitted:
        self.p_fitted = [0 if x == '' else x for x in self.p_fitted]
    if '' in self.p_error:
        self.p_error = [0 if x == '' else x for x in self.p_error]

    # submodel
    self.update_submodels()


def update_submodels(self):

    self.submodel = {}

    for submodel in self.parameters:

        # check if submodel should be used
        if 'y' in [use for sublist in [self.parameters[submodel][arg]['use'] for arg in self.parameters[submodel]] for use in sublist]:
            self.submodel[submodel] = {'guess_string': '', 'fit_string':''}

            # get tag
            try:
                submodel_tag = submodel.split('#')[-1]
            except IndexError:
                submodel_tag = None
            submodel_name = submodel.split('#')[0]

            # get arguments from function
            import __main__
            try:
                args_expected = list(inspect.signature(eval(f'__main__.{submodel_name}')).parameters)
            except AttributeError:
                args_expected = list(inspect.signature(eval(submodel_name)).parameters)

            # initialize submodel
            self.submodel[submodel]['guess_string'] += f'{submodel_name}(x, '
            self.submodel[submodel]['fit_string'] += f'{submodel_name}(x, '


            for arg in args_expected[1: ]:

                # check if submodel has active argument
                missing_arg = False
                try:
                    to_use = self.parameters[submodel][arg]['use'].index('y')
                except ValueError:
                    missing_arg = True
                if missing_arg: raise MissingArgument(submodel, arg)

                # build min, max, guess, model
                if self.parameters[submodel][arg]['id'][to_use] != '-':
                    id = self.parameters[submodel][arg]['id'][to_use]
                    self.submodel[submodel]['guess_string'] += str(self.p_guess[self.id_list.index(id)]) + ', '
                    self.submodel[submodel]['fit_string']   += str(self.p_fitted[self.id_list.index(id)]) + ', '
                else:
                    self.submodel[submodel]['guess_string'] += str(self.parameters[submodel][arg]['guess'][to_use]) + ', '
                    self.submodel[submodel]['fit_string'] += str(self.parameters[submodel][arg]['fitted'][to_use]) + ', '

            self.submodel[submodel]['guess_string'] = self.submodel[submodel]['guess_string'][:-2] + ')'
            self.submodel[submodel]['fit_string'] = self.submodel[submodel]['fit_string'][:-2] + ')'

            self.submodel[submodel]['guess'] = eval(f'lambda x:' + self.submodel[submodel]['guess_string'])
            self.submodel[submodel]['fit'] = eval(f'lambda x:' + self.submodel[submodel]['fit_string'])


def fit(self, x, y, ties=None, global_sigma=1e-13, save=True):

    self.update_model()

    # fit
    if global_sigma is not None:
        sigma = fake_sigma(x, global_sigma=global_sigma, sigma_specific=ties)
    self.p_fitted, self.p_cov = curve_fit(self.model, x, y, self.p_guess, sigma=sigma, bounds=[self.p_min, self.p_max])
    self.p_error = np.sqrt(np.diag(self.p_cov))  # One standard deviation errors on the parameters

    # get residue
    self.residue = trapz(abs(y - self.model(x, *self.p_fitted)), x)

    # save to sheet and self.parameter =====================================================
    # get header and submodels
    header = self.get_row_values(header_row, col_stop=last_col)
    fitted_col = header.index('fitted')+1
    error_col = header.index('error')+1

    for submodel in self.parameters:
        # check if submodel should be used
        if 'y' in [use for sublist in [self.parameters[submodel][arg]['use'] for arg in self.parameters[submodel]] for use in sublist]:

            # get tag
            try:
                submodel_tag = submodel.split('#')[-1]
            except IndexError:
                submodel_tag = None
            submodel_name = submodel.split('#')[0]

            # get arguments from function
            args_expected = list(inspect.signature(eval(submodel_name)).parameters)

            # build min, max, guess, model
            for arg in args_expected[1: ]:

                # check if submodel has active argument
                missing_arg = False
                try:
                    to_use = self.parameters[submodel][arg]['use'].index('y')
                except ValueError:
                    missing_arg = True
                if missing_arg: raise MissingArgument(submodel, arg)

                hashtag = list(self.parameters[submodel][arg]['#'])[to_use]
                if self.parameters[submodel][arg]['id'][to_use] != '-':
                    id = self.parameters[submodel][arg]['id'][to_use]
                    v1 = self.p_fitted[self.id_list.index(id)]
                    v2 = self.p_error[self.id_list.index(id)]
                    self.set_cell_value(value=v1, row=hashtag+header_row+1, col=fitted_col)
                    self.set_cell_value(value=v2, row=hashtag+header_row+1, col=error_col)
                    self.parameters[submodel][arg]['fitted'][to_use] = v1
                    self.parameters[submodel][arg]['error'][to_use] = v2

    self.update_submodels()

    if save:
        self.calc.save()


def fake_sigma(x, global_sigma=10**-10, sigma_specific=None):
    """Build a fake sigma array which determines the uncertainty in ydata.

    Adaptaded from the `scipy.optimize.curve_fit() <https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html>`_ documentation:

        If we define residuals as ``r = ydata - model(xdata, *popt)``, then sigma
        for a 1-D data should contain values of standard deviations of errors in
        ydata. In this case, the optimized function is ``chisq = sum((r / sigma) ** 2)``.

    Args:
        x (list): x array.
        sigma (float, optional): sigma value to be used for all points in ``x``.
        sigma_specific (list, optional): list of triples specfing new sigma for specific ranges, e.g.,
            ``sigma_specific = [[x_init, x_final, sigma], [x_init2, x_final2, sigma2], ]``.

    Returns:
        array.
    """
    p_sigma = np.ones(len(x))*global_sigma

    if sigma_specific is not None:
        for sigma in sigma_specific:
            init = index(x, sigma[0])
            final = index(x, sigma[1])
            p_sigma[init:final] = global_sigma/sigma[2]

    return p_sigma





def plot_fit(self, x, y, ax=None, show_exp=True, show_derivative=False, show_submodels=False, smoothing=10, ties=None, submodels_bkg=None, derivative_order=1, derivative_offset=None, derivative_factor=None, derivative_window_size=1):

    if smoothing == 0:
        x_smooth = x
    elif smoothing < 1:
        smoothing = 1
        x_smooth = np.linspace(min(x), max(x), len(x)*smoothing)
    else:
        x_smooth = np.linspace(min(x), max(x), len(x)*smoothing)

    if ax is None:
        fig = figm.figure()
        ax = fig.add_subplot(111)

    # exp
    if show_exp:
        ax.plot(x, y, **exp_def, label='exp')

    # ties
    if ties is not None:
        for pair in ties:
            ax.plot(x[am.index(x, pair[0]):am.index(x, pair[1])], y[am.index(x, pair[0]):am.index(x, pair[1])], **ties_def)

    # derivative
    if show_derivative:
        x_der, y_der = am.derivative(x, y, order=derivative_order, window_size=derivative_window_size)

        if derivative_factor is None:
            derivative_factor = (max(y) - np.mean(y))/(max(y_der) - np.mean(y_der))
        if derivative_offset is None:
            derivative_offset = -(max(y_der*derivative_factor) - np.mean(y))

        ax.plot(x_der, y_der*derivative_factor+derivative_offset, **der_def)

    # fit
    self.update_model()
    y_fit = self.model(x_smooth, *self.p_fitted)
    ax.plot(x_smooth, y_fit, **fit_def)

    if show_derivative:
        x_fir_der, y_fit_der = am.derivative(x_smooth, y_fit, order=derivative_order)
        ax.plot(x_fir_der, y_fit_der*derivative_factor+derivative_offset, **der_fit_def)

    # submodels
    if show_submodels:
        sub_colors = cm.get_cmap('Set1').colors
        if submodels_bkg is not None:
            bkg = self.submodel[submodels_bkg]['fit'](x_smooth)
        else:
            bkg=0
        for submodel, color in zip(self.submodel, sub_colors):
            if submodel != submodels_bkg:
                ax.plot(x_smooth, self.submodel[submodel]['fit'](x_smooth) + bkg, **sub_def, color=color, label=submodel)

    plt.legend()
    return ax


def plot_guess(self, x, y, ax=None, show_exp=True, show_derivative=False, show_submodels=False, smoothing=10, ties=None, submodels_bkg=None, derivative_order=1, derivative_offset=None, derivative_factor=None, derivative_window_size=1):

    if smoothing == 0:
        x_smooth = x
    elif smoothing < 1:
        smoothing = 1
        x_smooth = np.linspace(min(x), max(x), len(x)*smoothing)
    else:
        x_smooth = np.linspace(min(x), max(x), len(x)*smoothing)

    if ax is None:
        fig = figm.figure()
        ax = fig.add_subplot(111)

    # exp
    if show_exp:
        ax.plot(x, y, **exp_def, label='exp')

    # ties
    if ties is not None:
        for pair in ties:
            ax.plot(x[am.index(x, pair[0]):am.index(x, pair[1])], y[am.index(x, pair[0]):am.index(x, pair[1])], **ties_def)

    # derivative
    if show_derivative:
        x_der, y_der = am.derivative(x, y, order=derivative_order, window_size=derivative_window_size)

        if derivative_factor is None:
            derivative_factor = (max(y) - np.mean(y))/(max(y_der) - np.mean(y_der))
        if derivative_offset is None:
            derivative_offset = -(max(y_der*derivative_factor) - np.mean(y))

        ax.plot(x_der, y_der*derivative_factor+derivative_offset, **der_def)

    # fit
    y_guess = self.model(x_smooth, *self.p_guess)
    ax.plot(x_smooth, y_guess, **guess_def)

    if show_derivative:
        x_fir_der, y_fit_der = am.derivative(x_smooth, y_guess, order=derivative_order)
        ax.plot(x_fir_der, y_fit_der*derivative_factor+derivative_offset, **der_guess_def)

    # submodels
    if show_submodels:
        sub_colors = cm.get_cmap('Set1').colors
        if submodels_bkg is not None:
            bkg = self.submodel[submodels_bkg]['guess'](x_smooth)
        else:
            bkg=0
        for submodel, color in zip(self.submodel, sub_colors):
            if submodel != submodels_bkg:
                ax.plot(x_smooth, self.submodel[submodel]['guess'](x_smooth) + bkg, **sub_def, color=color, label=submodel)
    plt.legend()
    return ax
