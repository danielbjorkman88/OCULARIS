#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

# import os, sys
# from matplotlib.colors import LogNorm
import numpy as np
# from matplotlib import cm
# import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy.special import erf
import pandas as pd
from scipy.interpolate import interp1d
from copy import deepcopy
from matplotlib.colors import LogNorm
import math
import scipy.io
from PIL import Image
from scipy import ndimage
import pytiff
import copy

from skimage import io


class CCD:
    # pass
    def __init__(self, filename_raw, filename_bkg, filename_flatfield):
        
        try:
            self.raw = scipy.io.loadmat(filename_raw)["rawCCD"]
        except:
            self.raw = io.imread(filename_raw)
            
            
        # print("raw loaded")
        try:
            self.bkg = scipy.io.loadmat(filename_bkg)["bkgCCD"]
        except:
            self.bkg = io.imread(filename_bkg)
        
                  
        # print("bkg loaded")
        # try:
        #     with pytiff.Tiff(filename_flatfield) as handle:
        #         self.im = handle[0:, 0:]
        # except:
        im = io.imread(filename_flatfield)
        
        self.im = im/np.max(im)

        # print("flatfield loaded")        
        #self.im = Image.open(filename_flatfield)
        #print("image")
        
        image = self.raw*self.im - self.bkg
        
        image = (self.raw - self.bkg) / self.im 
        
        self.inter_image = image
        result = ndimage.median_filter(image, size=5)
        self.inter_image2 = copy.deepcopy(result)
        result = result - result[int(result.shape[0])-1, int(result.shape[1]/2) ]

        self.image = result/result[int(result.shape[0]/2), int(result.shape[1]/2) ]
        
        
        #self.pixelsize = 0.0475*2 # mm/pixel
        self.pixelsize = 0.0939 # mm/pixel
        self.Npixels = self.raw.shape[0]
        self.minX = -self.Npixels/2 *self.pixelsize
        self.maxX = self.Npixels/2 *self.pixelsize
        
        
        
        
        self.trans_X, self.trans_Y = np.meshgrid(np.arange(self.minX + self.pixelsize,self.maxX + self.pixelsize, self.pixelsize),np.arange(self.minX + self.pixelsize,self.maxX + self.pixelsize, self.pixelsize))
        
        
        print("CCD image loaded")
        
        #self.trans_X2, self.trans_Y2 = np.meshgrid(np.arange(self.minX ,self.maxX, self.pixelsize),np.arange(self.minX,self.maxX, self.pixelsize))


# import os

# curr_path = os.getcwd()
# ccd_data_path = os.getcwd()[0:-11] + "\Data\profileinair"
# plot_path = os.getcwd()[0:-11] + "\plots"
# data_path = os.getcwd()[0:-11] + "\Data"


# os.chdir(data_path)



# raw = scipy.io.loadmat('raw.mat')["rawCCD"]
# print("raw done")
# bkg = scipy.io.loadmat('bkg.mat')["bkgCCD"]
# print("bkg done")


# with pytiff.Tiff('flatfield_gantry2_extrapolated_50.50.800_uint16.tif') as handle:
#     im = handle[0:, 0:]
    
    
    
# im = Image.open('flatfield_gantry2_extrapolated_50.50.800_uint16.tif')
# im2 = plt.imread('flatfield_gantry2_extrapolated_50.50.800_uint16.tif')
# print("image")

# image = raw*im - bkg
# result = ndimage.median_filter(image, size=5)
# result = result - result[int(result.shape[0])-1, int(result.shape[1]/2) ]

# image = result/result[int(result.shape[0]/2), int(result.shape[1]/2) ]


# pixelsize = 0.0475*2 # mm/pixel
# Npixels = raw.shape[0]
# minX = -Npixels/2 *pixelsize
# maxX = Npixels/2 *pixelsize




# trans_X, trans_Y = np.meshgrid(np.arange(minX + pixelsize,maxX + pixelsize, pixelsize),np.arange(minX + pixelsize,maxX + pixelsize, pixelsize))


#ccd = CCD('raw.mat', 'bkg.mat', 'flatfield_gantry2_extrapolated_50.50.800_uint16.tif')








