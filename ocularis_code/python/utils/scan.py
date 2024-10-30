# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


# import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
# import matplotlib.patches as patches
from pathlib import Path
import math
# from src.BroadBeam import BroadBeam
# from src.CCD import CCD
from statistics import mean
# from src.DepthProfile import DepthProfile
# from src.BroadBeam import collimation_modifer
# from pathlib import Path
# import dose_engine

from scipy.optimize import curve_fit
from scipy.special import erf


def collimation_modifer3(r,  beta, r_half):
    return 0.5 - 0.5*erf(beta* (r - r_half)/(math.sqrt(2)))



class Scan():
    def __init__(self, filename, path, start = 58, title = ""):
        self.filename = filename
        self.path = Path(path)
        self.title = title
        self.start = start
        
        
    def simple_config(self):
        data = np.loadtxt(self.path / self.filename, skiprows = 7)

        self.val = data[0:,5]
        self.val_raw = data[0:,5]
        xes = data[0:,0] - min(data[0:,0])
        self.xes = xes - max(xes)/2
        # yes = data[0:,1] - min(data[0:,1])
        self.yes = data[0:,1] # - max(yes)/2
        self.zes = data[0:,2]
        
        self.zpos = round(-(self.zes[0] - self.start) + 2.51, 5)       
        
        
        self.val = self.val - min(self.val)
        normFactor = max(self.val)
        
        self.val = self.val/normFactor  
        
        
    def config_horizontal(self):           
        data = np.loadtxt(self.path / self.filename, skiprows = 7)

        self.val = data[0:,5]
        self.val_raw = data[0:,5]
        xes = data[0:,0] - min(data[0:,0])
        self.xes = xes - max(xes)/2
        # yes = data[0:,1] - min(data[0:,1])
        self.yes = data[0:,1] # - max(yes)/2
        self.zes = data[0:,2]
        
        self.zpos = round(-(self.zes[0] - self.start) + 2.51, 5)
        
        
        # if self.yes[0] < 0:
        #     self.yes = np.flip(self.yes)
        #     self.val = np.flip(self.val)
        #     self.val_raw = np.flip(self.val_raw)
        
        # if min(self.val[0:6]) < 0:
        self.val = self.val - self.val[0]
        
        
        # a = np.where( self.yes < -4 )[0]
        # b = np.where(self.yes > -10)[0]
        # #print(a, b)
        # start = max(min(a), min(b))

        # end = min(max(a), max(b))
        
        # self.normRangeBins = [start, end]
        

        
        
        normFactor = max(self.val)
        
        self.val = self.val/normFactor        
        
        # self.normRange = [self.yes[self.normRangeBins[0]], self.yes[self.normRangeBins[1]]]        
        #normFactor = max(self.val)
        

        
        
        f1 = interp1d(self.val[0:int(len(self.val)/2)], self.xes[0:int(len(self.val)/2)])  
        f2 = interp1d(self.val[int(len(self.val)/2):], self.xes[int(len(self.val)/2):])  
        
        self.intersects = [f1(0.5), f2(0.5)]
        
        self.shift = -(max(self.intersects) - (abs(self.intersects[0]) + abs(self.intersects[1]))/2)
        self.xes_centered = self.xes + self.shift
        
        
        

        # yes = self.yes
        
        
        
        
        # f1 = interp1d(self.val[0:int(len(self.val)/2)], yes[0:int(len(self.val)/2)])  
        # f2 = interp1d(self.val[int(len(self.val)/2):], yes[int(len(self.val)/2):])  
        
        # intersects = [f1(0.5), f2(0.5)]
        
        # shift = -(max(intersects) - (abs(intersects[0]) + abs(intersects[1]))/2)
        # yes_centered = yes +shift        
        
        
        

        f1c = interp1d(self.val[0:int(len(self.val)/2)], self.xes_centered[0:int(len(self.val)/2)])  
        f2c = interp1d(self.val[int(len(self.val)/2):], self.xes_centered[int(len(self.val)/2):])  
        
        
        self.intersects_centered = [f1c(0.5), f2c(0.5)]        
        self.mid = -(max(self.intersects_centered) - (abs(self.intersects_centered[0]) + abs(self.intersects_centered[1]))/2)

        
        
        
        
    def config_vertical(self):        
        data = np.loadtxt(self.path / self.filename, skiprows = 7)

        self.val = data[0:,5]
        self.val_raw = data[0:,5]
        self.xes = data[0:,0]
        yes = data[0:,1] - min(data[0:,1])
        self.yes = yes - max(yes)/2
        self.zes = data[0:,2]
        
        self.zpos = round(-(self.zes[0] - self.start) + 2.51, 5)
        
        
        if self.yes[0] < 0:
            self.yes = np.flip(self.yes)
            self.val = np.flip(self.val)
            self.val_raw = np.flip(self.val_raw)
        
        # if min(self.val[0:6]) < 0:
        self.val = self.val - self.val[0]
        
        
        a = np.where( self.yes < -4 )[0]
        b = np.where(self.yes > -10)[0]
        #print(a, b)
        start = max(min(a), min(b))

        end = min(max(a), max(b))
        
        self.normRangeBins = [start, end]
        

        
        
        normFactor = max(self.val) #mean(self.val[self.normRangeBins[0]:self.normRangeBins[1]])
        
        self.val = self.val/normFactor        
        
        self.normRange = [self.yes[self.normRangeBins[0]], self.yes[self.normRangeBins[1]]]        
        #normFactor = max(self.val)
        

        
        
        f1 = interp1d(self.val[0:int(len(self.val)/2)], self.yes[0:int(len(self.val)/2)])  
        f2 = interp1d(self.val[int(len(self.val)/2):], self.yes[int(len(self.val)/2):])  
        
        self.intersects = [f1(0.5), f2(0.5)]
        
        self.shift = -(max(self.intersects) - (abs(self.intersects[0]) + abs(self.intersects[1]))/2)
        self.yes_centered = self.yes + self.shift
        
        
        

        # yes = self.yes
        
        
        
        
        # f1 = interp1d(self.val[0:int(len(self.val)/2)], yes[0:int(len(self.val)/2)])  
        # f2 = interp1d(self.val[int(len(self.val)/2):], yes[int(len(self.val)/2):])  
        
        # intersects = [f1(0.5), f2(0.5)]
        
        # shift = -(max(intersects) - (abs(intersects[0]) + abs(intersects[1]))/2)
        # yes_centered = yes +shift        
        
        
        

        f1c = interp1d(self.val[0:int(len(self.val)/2)], self.yes_centered[0:int(len(self.val)/2)])  
        f2c = interp1d(self.val[int(len(self.val)/2):], self.yes_centered[int(len(self.val)/2):])  
        
        
        self.intersects_centered = [f1c(0.5), f2c(0.5)]        
        self.mid = -(max(self.intersects_centered) - (abs(self.intersects_centered[0]) + abs(self.intersects_centered[1]))/2)

        
        




    def curve_fit(self):
        # xes = self.yes_centered - min(self.intersects_centered)

        # yes = self.val
        # yes = yes[xes < 15]
        # xes = xes[xes < 15]
        
        # popt, pcov = curve_fit(collimation_modifer2, xes, yes)
        # self.beta = popt[0]
        
        xes = self.yes_centered
    
        yes = self.val_raw
        
        yes = yes - yes[0]
        
        
        yes1 = yes[xes < 0]
        self.xes1 = xes[xes < 0]
        
        yes2 = yes[xes > 0]
        self.xes2 = xes[xes > 0]    
        
        self.normFactor = max(yes1)
        self.normFactor2 = max(yes2)
        
        self.yes1 = yes1/self.normFactor
        self.yes2 = yes2/self.normFactor2
        
        
        
        
        
        f_sub = interp1d(self.yes1, self.xes1)
        self.R80 = f_sub(max(self.yes1)*0.8).max()
        self.R20 = f_sub(max(self.yes1)*0.2).max()
        
        f_sub_2 = interp1d(self.yes2, self.xes2)
        self.R80_2 = f_sub_2(max(self.yes2)*0.8).max()
        self.R20_2 = f_sub_2(max(self.yes2)*0.2).max()    
        
        
        self.popt1, self.pcov1 = curve_fit(collimation_modifer3, self.xes1, self.yes1, p0 = [-1, min(self.intersects_centered)])
        self.beta1 = self.popt1[0]
        
        # plt.plot(xes1, collimation_modifer3(xes1, popt[0], popt[1]), 'g--',
        #       label=f'fit: beta={round(popt[0],5)} r_half={round(popt[1],5)}, sd = {np.sqrt(np.diag(pcov))}' )
        
        
        self.popt2, self.pcov2 = curve_fit(collimation_modifer3, self.xes2, self.yes2, p0 = [1, max(self.intersects_centered)])     
        self.beta2 = self.popt2[0]
    
    
        self.beta = mean([abs(self.beta1), abs(self.beta2)])
        self.R80_R20 = mean([abs(self.R80 - self.R20),  abs(self.R80_2 - self.R20_2)])
    
    
    
    def show_fit(self):
    
        fig = plt.figure()    
    
        # xes = self.yes_centered - min(self.intersects_centered)

        # yes = self.val
        # yes = yes[xes < 15]
        # xes = xes[xes < 15]
    
        # normRange = self.normRange - min(self.intersects_centered)
    
    
        # plt.figure()
        
        # plt.plot(xes, yes, 'b-', label='data')  
        

        # plt.plot(xes, collimation_modifer2(xes, self.beta), 'g--',
        #           label='fit: beta=%5.3f' % self.beta)
        
        
        # plt.axvline(x= normRange[0] , color = "C1", linewidth = 1, linestyle = "--", label = f"Norm range from {round(normRange[0],3)} to {round(normRange[1],3)}")
        # plt.axvline(x= normRange[1] , color = "C1", linewidth = 1, linestyle = "--")
        
        # plt.ylim(-0.05, 1.05)
        
        # plt.title(f"{self.title} at z = {self.zpos} with beta = {self.beta} ")
        
        # plt.legend()
        
        # plt.show()
        
        plt.plot(self.xes1, collimation_modifer3(self.xes1, self.popt1[0], self.popt1[1]), 'g--',
              label=f'fit: beta={round(self.popt1[0],5)} r_half={round(self.popt1[1],5)}, sd = {np.sqrt(np.diag(self.pcov1))}' )
        
        
        # popt, pcov = curve_fit(collimation_modifer3, self.xes2, self.yes2, p0 = [1, max(scan.intersects_centered)])    
        plt.plot(self.xes2,self.normFactor2/self.normFactor *collimation_modifer3(self.xes2, self.popt2[0], self.popt2[1]), 'C1--',
              label=f'fit: beta={round(self.popt2[0],5)} r_half={round(self.popt2[1],5)}, sd = {np.sqrt(np.diag(self.pcov2))}' )    
        
        
        yes = self.val_raw
        yes = yes - yes[0]
        
        
        plt.plot(self.yes_centered, yes/self.normFactor, color = "C0")        
        
        # plt.plott(self.xes1,yes1, color = "C0")
        
        
        plt.axvline(x = self.R80, label = f'Left R80 - R20 = {self.R80 - self.R20}', color = "C3", linestyle = "--")
        plt.axvline(x = self.R20, color = "C3", linestyle = "--")    
        
        plt.axvline(x = self.R80_2, label = f'Right R80 - R20 = {abs(self.R80_2 - self.R20_2)}', color = "C4", linestyle = "--")
        plt.axvline(x = self.R20_2, color = "C4", linestyle = "--")        
        
    
        plt.title(f"{self.title} , z = {self.zpos}")
       
        plt.legend()
        plt.xlabel("y [mm]", fontsize = 14)
        
        plt.grid(linewidth = 0.3)
        
        xlength = 14
        fig.set_size_inches(xlength, xlength/1.61803398875)
        plt.show()
                

    def __repr__(self):
        return self.title + " . " + self.filename

