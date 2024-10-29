# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 15:08:35 2022

@author: bjoerk_c
"""


from matplotlib.pyplot import cm
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib import cm
import matplotlib.colors as mcolors
import copy

my_colors = ["C0", "C1", "C2", "C3", "C4", "C5", "C6", "C7",
             "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16"]


N = 10
rainbow_colors = cm.rainbow(np.linspace(0, 1,N))


def gamma_index_colormap(image):
    
    N_colors = 1000
    
    if np.min(image) < 0:
        image[image < 0] = 0
    
    n_colors_blue = int(N_colors/np.max(image))
    
    if np.max(image) < 1 or n_colors_blue >= N_colors or n_colors_blue < 0:
        color = plt.cm.Blues(np.linspace(0., 1, N_colors))
        return mcolors.LinearSegmentedColormap.from_list('my_colormap', color)
    
    colors1 = plt.cm.Blues(np.linspace(0., 1, n_colors_blue))
    
    n_colors_red = N_colors - n_colors_blue
    colors2 = plt.cm.autumn(np.linspace(0., 1, n_colors_red ) )
    
    colorIdx = 0
    for rowIdx in range(len(colors2)):
        colors2[rowIdx, 0:] = colors2[colorIdx, 0:]
        colorIdx +=1
        # if rowIdx % 150 == 0:
        #     colorIdx += 149
            
    color_list = np.vstack((colors1, colors2))
    new_cmap = mcolors.LinearSegmentedColormap.from_list('my_colormap', color_list)
    return new_cmap
