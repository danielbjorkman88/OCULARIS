# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

from plot_utils.Plotter import Plotter

from scipy.special import erf
from scipy.interpolate import interp1d
from pathlib import Path

from matplotlib.pyplot import cm
import numpy as np
import matplotlib.pyplot as plt


import dose_engine
import plot_utils
import reference_frames


from sympy import symbols
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


# class Arrow3D(FancyArrowPatch):
#     def __init__(self, xs, ys, zs, *args, **kwargs):
#         FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
#         self._verts3d = xs, ys, zs

#     def draw(self, renderer):
#         xs3d, ys3d, zs3d = self._verts3d
#         xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M )
#         self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
#         FancyArrowPatch.draw(self, renderer)

#     def do_3d_projection(self, renderer):
#         xs3d, ys3d, zs3d = self._verts3d
#         xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
#         self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))


import matplotlib.pyplot as plt


class PlotReferenceFrames(Plotter):


    def plot(self):
        
            
        ref = reference_frames.Reference_Frame(self.algo.config, 'TreatmentRoom')
        bev = reference_frames.Reference_Frame(self.algo.config, 'Beams Eye View')
        

        q1,q2,q3 = symbols('x y z')
        
        rays = []
        
        rays.append(dose_engine.Ray(ref.central_axis.S, ref.central_axis.T))
        
        for x in np.linspace(-1, 1, num = 10):
            for y in np.linspace(-1, 1, num = 10):
                rays.append(dose_engine.Ray(ref.central_axis.S, ref.central_axis.T + x*ref.xvec_in_ref + y*ref.yvec_in_ref))
        
        
        fig = plt.figure()
        
        #------------------------------------------
        ax1 = fig.add_subplot(121, projection='3d')
        ax1.view_init(-70,90) 
        plt.title("Treatment Room", fontsize = 18, y=1.15)
        
        ref.add_frame_to_plot(ax1, ref)
        
        
        # arrow_prop_dict = dict(mutation_scale=20, arrowstyle='->', shrinkA=0, shrinkB=0)
        
        
        for ray in rays:
            plt.plot(ray.xes, ray.yes, ray.zes , color="C1")
            # a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color="C1")
            # ax1.add_artist(a)
            
        
        ax1.text(ref.origin_ref[0], ref.origin_ref[1], ref.origin_ref[2] -0.1, r'$0_{iso}$')
        
        
        ax1.set_xlabel('X')
        ax1.set_xlim(-1, 1)
        ax1.set_ylabel('Y')
        ax1.set_ylim(-1, 1)
        ax1.set_zlabel('Z')
        ax1.set_zlim(-1, 1)
        
        
        #------------------------------------------
        ax2 = fig.add_subplot(122, projection='3d')
        plt.title("BEV in Treatment Room frame", fontsize = 18, y=1.15)
        
        
        ax2.view_init(-70,90) #elev=0., azim=0
        
        bev.add_frame_to_plot(ax2, ref)
        
        
        for ray in rays:
            plt.plot(ray.xes, ray.yes, ray.zes , color="C1")
            # a = Arrow3D(ray.xes, ray.yes, ray.zes , **arrow_prop_dict, color="C1")
            # ax2.add_artist(a)
        
        
        
        ax2.text(bev.origin_ref[0], bev.origin_ref[1], bev.origin_ref[2] + 0.1, r'$0_{Source}$')
        
        ax2.set_xlabel('X')
        ax2.set_xlim(bev.origin_ref[0] + 1, bev.origin_ref[0] - 1)
        ax2.set_ylabel('Y')
        ax2.set_ylim(bev.origin_ref[1] + 1, bev.origin_ref[1] - 1)
        ax2.set_zlabel('Z')
        ax2.set_zlim(bev.origin_ref[2] + 1, bev.origin_ref[2] - 1)
        
         
        
        
        plt.suptitle(f"Reference frames. phi = {self.algo.config['Gantry rot phi']}, theta = {self.algo.config['Gantry rot theta']}", fontsize = 24)


        xlength = 12
        fig.set_size_inches(xlength, xlength/1.61803398875)
        #fig.set_size_inches(xlength, xlength)
        plt.show()
        try:
            plt.savefig(self.plot_path / "Reference systems_theta=_phi=_.pdf",  bbox_inches = 'tight', pad_inches = 0.1)
        except:
            pass;
            

