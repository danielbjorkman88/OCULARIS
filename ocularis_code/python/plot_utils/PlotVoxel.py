# -*- coding: utf-8 -*-
"""
Created on Thu Jan  6 13:53:12 2022

@author: bjoerk_c
"""

from plot_utils.Plotter import Plotter

from scipy.special import erf
from scipy.interpolate import interp1d
from pathlib import Path
from configuration import config, data
from matplotlib.pyplot import cm

import matplotlib.pyplot as plt

import numpy as np



from scipy.ndimage import zoom

from matplotlib import cm
import matplotlib.colors as mcolors




class PlotVoxel(Plotter):

    def __init__(self, algo, downsampling, Ncolors):
        super().__init__(algo)    
        self.downsampling = downsampling
        self.Ncolors = Ncolors

    def plot(self):

        
        # Ncolors = 100
        # downsampling = [0.25, 0.25, 1]
        Ncolors = self.Ncolors
        downsampling = self.downsampling
        
        dose_sub = zoom(self.algo.dose, (downsampling[0], downsampling[1], downsampling[2]))
        dose_sub = dose_sub/np.max(dose_sub)
        
        first_quadrant_sub = zoom(self.algo.first_quadrant_bev, (downsampling[0], downsampling[1], downsampling[2]))

        
        bool_array = first_quadrant_sub > 0.999999
        
        dose_sub[bool_array] = 0
        
        # dose_sub[0:int(dose_sub.shape[0]/2), 0:int(dose_sub.shape[1]/2), 0:] = 0
        
        
        x, y, z = np.indices((dose_sub.shape[0], dose_sub.shape[1], dose_sub.shape[2]))
        
        
        
        
        
        # draw cuboids in the top left and bottom right corners, and a link between
        # them
        # cube1 = (x < 3) & (y < 3) & (z < 3)
        # cube2 = (x >= 5) & (y >= 5) & (z >= 5)
        # link = abs(x - y) + abs(y - z) + abs(z - x) <= 2
        
        
        
        cmap = cm.get_cmap('jet', Ncolors)   # 21 = 3 colors per decade times 7 decades  
        
        # fig = plt.figure()
        # plt.pcolor( dose[0:, 0:, 0], cmap = cmap, vmin=0, vmax=1)
        # cbar = plt.colorbar()
        # plt.show()
        
        
        
        # colorIdx = 0
        # for rowIdx in range(len(colors2)):
        #     colors2[rowIdx, 0:] = colors2[colorIdx, 0:]
        #     if rowIdx % 150 == 0:
        #         colorIdx += 149
                
        # colors = np.vstack((colors1, colors2))
        
        colors1 = plt.cm.jet(np.linspace(0., 1, Ncolors))
        
        display_cube = dose_sub > 0.001
        
        
        
        colors = np.empty(dose_sub.shape, dtype=object)
        
        for i in range(Ncolors):
            colors[dose_sub > i/Ncolors] = mcolors.to_hex(colors1[i][0:3])
            
            
            
            
        
        # hot = dose_sub > 0.9
        # cool =  (dose_sub < 0.9)  & (dose_sub > 0.1)
        
        # colors = np.empty(dose_sub.shape, dtype=object)
        # colors[cool] = 'blue'
        # colors[hot] = 'red'
        
        # # combine the objects into a single boolean array
        # voxelarray = cube1 | cube2 | link
        # # set the colors of each object
        # colors = np.empty(voxelarray.shape, dtype=object)
        # colors[link] = 'red'
        # colors[link] = 'red'
        # colors[cube1] = 'blue'
        # colors[cube2] = 'green'
        
        # # and plot everything
        # ax = plt.figure().add_subplot(projection='3d')
        # ax.voxels(voxelarray, facecolors=colors, edgecolor='k')
        
        
        
        fig = plt.figure()
        
        ax = fig.add_subplot(projection='3d')
        ax.voxels(display_cube, facecolors=colors, alpha = 1)
        
        
        ax.set_xlabel('X', fontsize = 12)
        # ax1.set_xlim(-1, 1)
        ax.set_ylabel('Y', fontsize = 12)
        # ax1.set_ylim(-1, 1)
        ax.set_zlabel('Z', fontsize = 12)
        
        norm = mcolors.Normalize(vmin=0, vmax=1)
        m = cm.ScalarMappable(cmap=cmap, norm=norm)
        m.set_array([])
        cbar = plt.colorbar(m)
        cbar.ax.set_ylabel('Dose', rotation=270, size='xx-large')
        
        if self.algo.wedge.traced:
            plt.suptitle(f' Theta = {self.algo.config["Gantry rot theta"]} degrees , Phi = {self.algo.config["Gantry rot phi"]} degrees, Wedge cover = {self.algo.config["wedge_cover"]}', fontsize = 16)
        else:
            plt.suptitle(f' Theta = {self.algo.config["Gantry rot theta"]} degrees , Phi = {self.algo.config["Gantry rot phi"]} degrees', fontsize = 16)
        
        
        xlength = 12
        fig.set_size_inches(xlength, xlength/1.61803398875)
        #fig.set_size_inches(xlength, xlength)
        plt.show()
        
        

        try:
            plt.savefig(self.plot_path / "PlotVoxel.pdf",  bbox_inches = 'tight', pad_inches = 0.1)
        except:
            pass;
        
        
