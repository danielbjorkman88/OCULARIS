# -*- coding: utf-8 -*-


# import os

# from pathlib import Path
import matplotlib.pyplot as plt
# import copy
# from scipy.optimize import minimize
import numpy as np
# from sklearn.metrics import mean_squared_error
from scipy.interpolate import interp1d
# from scipy.optimize import curve_fit
# from matplotlib import cm


# import numpy.polynomial.polynomial as poly

from utils.scan import Scan
# from configuration import config, data
# import dose_engine


def find_lateral_inflections(xes, yes):
    return find_inflection_L(xes, yes), find_inflection_R(xes, yes)


def find_lateral_sharpness(xes, yes):
    L80, L20 = find_80_20_L(xes, yes)
    R80, R20 = find_80_20_R(xes, yes)
    return L20, L80, R80, R20


def find_inflection_R(xes, yes):

    f = interp1d(xes, yes)
    f_rev = interp1d(yes, xes)

    subx = xes[xes > 0.99*f_rev(max(yes))]
    suby = f(subx)

    try:
        f_sub = interp1d(suby, subx)
    except:
        raise ValueError

    R50 = f_sub(max(suby)*0.5).max()
    return R50


def find_inflection_L(xes, yes):

    f = interp1d(xes, yes)
    f_rev = interp1d(yes, xes)

    subx = xes[xes < 0.99*f_rev(max(yes))]
    suby = f(subx)

    f_sub = interp1d(suby, subx)

    R50 = f_sub(max(suby)*0.5).max()
    return R50


def find_80_20_R(xes, yes):

    f = interp1d(xes, yes)
    f_rev = interp1d(yes, xes)

    subx = xes[xes > 0.99*f_rev(max(yes))]
    suby = f(subx)

    f_sub = interp1d(suby, subx)

    R80 = f_sub(max(suby)*0.8).max()
    R20 = f_sub(max(suby)*0.2).max()
    return R80, R20


def find_80_20_L(xes, yes):

    f = interp1d(xes, yes)
    f_rev = interp1d(yes, xes)

    subx = xes[xes < 0.99*f_rev(max(yes))]
    suby = f(subx)

    f_sub = interp1d(suby, subx)

    L80 = f_sub(max(suby)*0.8).max()
    L20 = f_sub(max(suby)*0.2).max()
    return L80, L20


def find_R90(xes, yes):

    f = interp1d(xes, yes)

    f_rev = interp1d(yes, xes)

    # P90 = f_rev(max(yes)*0.9).max()
    # P80 = f_rev(max(yes)*0.8).max()

    subx = xes[xes > 0.99*f_rev(max(yes))]
    suby = f(subx)

    f_sub = interp1d(suby, subx)

    # R100 = f_sub(max(suby)*1).max()
    R90 = f_sub(max(suby)*0.9).max()
    return R90



def find_R90_from_left(xes, yes):
    
    # xs, ys = xes[0:split], yes[0:split]
    
    tolerance = 0.01  # adjust this as needed
    for i in range(len(yes) - 1):
        
        xs = [xes[i], xes[i + 1] ]
        
        ys = [yes[i], yes[i + 1] ]
        
        f_rev = interp1d( ys, xs)
        
        try:
            return f_rev(0.9)
        except:
            pass
        
def find_R90_from_right(xes, yes):

    f = interp1d(xes, yes)
    f_rev = interp1d(yes, xes)

    subx = xes[xes > 0.9*f_rev(max(yes))]
    suby = f(subx)

    try:
        f_sub = interp1d(suby, subx)
    except:
        pass

    R90_R = f_sub(max(suby)*0.9).max()
    return R90_R
        
def find_R10_from_left(xes, yes):
    
    # xs, ys = xes[0:split], yes[0:split]
    
    tolerance = 0.01  # adjust this as needed
    for i in range(len(yes) - 1):
        
        xs = [xes[i], xes[i + 1] ]
        
        ys = [yes[i], yes[i + 1] ]
        
        f_rev = interp1d( ys, xs)
        
        try:
            return f_rev(0.1)
        except:
            pass
        
def find_R10_from_right(xes, yes):

    f = interp1d(xes, yes)
    f_rev = interp1d(yes, xes)

    subx = xes[xes > 0.1*f_rev(max(yes))]
    suby = f(subx)

    try:
        f_sub = interp1d(suby, subx)
    except:
        pass

    R10_R = f_sub(max(suby)*0.1).max()
    return R10_R
        

def find_R90_increasing(xes, yes):

    f = interp1d(xes, yes)

    f_rev = interp1d(yes, xes)

    # P90 = f_rev(max(yes)*0.9).max()
    # P80 = f_rev(max(yes)*0.8).max()

    subx = xes[xes < 1*f_rev(max(yes))]
    suby = f(subx)

    f_sub = interp1d(suby, subx)

    R100 = f_sub(max(suby)*1).max()
    R90 = f_sub(max(suby)*0.9).max()
    return R90




def config_scan(filename, path):
    scan = Scan(filename, path)
    scan.simple_config()
    scan.zes = np.flip(scan.zes)
    scan.zes -= find_R90(scan.zes, scan.val)
    return scan


# Bragg peak depth curve
def find_xpos_and_norm_factor(xes, yes):
    idx_split = -1000


    
    prev = -1
    for idx in np.flip(range(len(xes))):

        val = yes[idx]
        # print(val)
        
        if np.isnan(val):
            continue
        
        if val < 0.1:
            prev = val
            continue

        if val >= prev:
            prev = val
            continue
        
        # print(val)
        idx_split = idx
        # print(idx)
        break
    
    # print("idx_split" ,idx_split)
    
    # fig = plt.figure()
    
    # plt.plot(xes,yes)
    
    # plt.scatter(xes[idx_split], yes[idx_split])

    
    # plt.show()
    

    sub_xes = xes[idx_split:]
    sub_yes = yes[idx_split:]
    normfactor = np.max(sub_yes)
   # sub_yes = sub_yes/max(sub_yes)

    # R90 = find_R90(sub_xes, sub_yes)
    x_pos = xes[idx_split + 1]
    #sub_xes = sub_xes - R90
    return x_pos, normfactor


def find_DX(xes, yes, threshold):

    # f = interp1d(xes, yes)
    
    non_zero_count = np.count_nonzero(yes) + 1
    f_rev = interp1d(yes[:non_zero_count], xes[:non_zero_count], fill_value = "extrapolate")
    
    # f_rev = interp1d(yes, xes, fill_value = "extrapolate")
    # f_rev = interp1d(yes, xes)
    
    return f_rev(threshold)