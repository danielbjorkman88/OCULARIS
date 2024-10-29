# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 16:14:26 2022

@author: bjoerk_c
"""

from sympy import symbols
from sympy.physics.vector import ReferenceFrame, express
import sympy as sym
from mpmath import radians
import numpy as np
import math
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from sympy.vector import CoordSys3D
import os, sys
import inspect
from configuration import constants

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

import dose_engine

import math as m

import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt
# from matplotlib.patches import FancyArrowPatch
# from mpl_toolkits.mplot3d import proj3d



# class Arrow3D(FancyArrowPatch):
#     def __init__(self, xs, ys, zs, *args, **kwargs):
#         FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
#         self._verts3d = xs, ys, zs

#     def draw(self, renderer):
#         xs3d, ys3d, zs3d = self._verts3d
#         xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M )
#         self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
#         FancyArrowPatch.draw(self, renderer)



def R_spherical_to_cartesian(azimuth, polar):
  return np.matrix([[ m.sin(azimuth)*m.cos(polar), m.cos(azimuth)*m.cos(polar), -m.sin(polar) ],
                   [ m.sin(azimuth)*m.sin(polar), m.cos(azimuth)*m.sin(polar) , m.cos(polar) ],
                   [ m.cos(azimuth)           , - m.sin(azimuth)            , 0 ]])

def find_axis_direction(theta, phi):
    
    if theta != 0 and phi != 0:
        v = np.asarray([math.sin(theta), math.sin(phi), math.cos(theta)])
        v = v/ np.linalg.norm(v)
        return v
    
    return np.asarray([math.sin(theta), math.sin(phi), math.cos(theta)*math.cos(phi)])



class Reference_Frame:
    def __init__(self, config, name_of_frame):

        if not type(config) == dict:
            raise ValueError
        
        self.name_of_frame = name_of_frame
        
        assert (name_of_frame in ['TreatmentRoom' , 'Beams Eye View' , 'Gaze'])
        
        
        self.config = config
        
        self.reference = ReferenceFrame('TreatmentRoom')
        self.coordsys_ref = CoordSys3D('TreatmentRoom')
        
        q1,q2,q3 = symbols('x y z')
        
        self.frame = ReferenceFrame(name_of_frame)
        self.coordsys = CoordSys3D(name_of_frame)
        
        

        S = [0,0, abs(constants["VPS"])]
        T = [0,0,-100]
        
        d_S_iso = math.dist(S, [0,0,0])
        d_T_iso = math.dist(T, [0,0,0])
        
        # This redefinition is used when moving the source for intra-fractional motion
        if "new_S" in self.config and "new_T" in self.config:
            
            self.new_S = self.config["new_S"]
            self.new_T = self.config["new_T"]

        else:
            axis_direction = find_axis_direction(self.config["Gantry rot theta"]*math.pi/180, self.config["Gantry rot phi"]*math.pi/180)
        
            self.new_S = axis_direction*d_S_iso
            self.new_T = -axis_direction*d_T_iso
            
        if not math.isclose(math.dist(S, T) , math.dist(self.new_S, self.new_T), rel_tol = 0.001):
            raise ValueError
        
        
        self.central_axis = dose_engine.Ray(self.new_S, self.new_T)
        
        

     
        
        if self.name_of_frame == 'Beams Eye View':
            
           
            self.coordsys = self.coordsys_ref.locate_new(self.name_of_frame, self.central_axis.S[0]*self.coordsys_ref.i + self.central_axis.S[1]*self.coordsys_ref.j + self.central_axis.S[2]*self.coordsys_ref.k) 
           
            self.frame.orient_body_fixed(self.reference, (radians(- self.config["Gantry rot phi"] ), radians(180 - self.config["Gantry rot theta"]),0), 'XYX') 
           
            if self.config["Gantry rot phi"] != 0 and  self.config["Gantry rot theta"] != 0:

                #TODO:
                # Doesnt rotate the frame as intended
                # This if-else block solves the definiton of base axis in this frame
                # however no frame rotation is taking place. Therefor

                
                # Defines base x and y vectors in Gantry reference system
                # x = np.random.randn(3)  # take a random vector
                x = np.asarray([1,-1,1e-16])
                x -= x.dot(self.central_axis.vec_norm) * self.central_axis.vec_norm  # make it orthogonal to central axis
                x /= np.linalg.norm(x)  # normalize vector
                y = np.cross(self.central_axis.vec_norm, x)      # cross product with central axis
                
                
                # Base axis in respective frame
                self.zvec_in_ref = self.central_axis.vec_norm
                self.xvec_in_ref = x
                self.yvec_in_ref = y
                
            else:
                # Base axis in respective frame
                self.xvec_in_ref = np.asarray(self.express_vector_in_ref([1,0,0]))
                self.yvec_in_ref = np.asarray(self.express_vector_in_ref([0,1,0]))
                self.zvec_in_ref = np.asarray(self.express_vector_in_ref([0,0,1]))
        
    

            #self.frame = self.reference.orientnew('Beams Eye View', 'Space', (-radians(self.config["Gantry rot phi"]), -radians(self.config["Gantry rot theta"]), 0), "123")
            # R = R_spherical_to_cartesian(-radians(self.config["Gantry rot theta"]), - radians(self.config["Gantry rot phi"]))
            # R = sym.Matrix(R)
            # self.frame = self.reference.orientnew('Beams Eye View', 'DCM', R)
            #self.coordsys = self.coordsys_ref.locate_new(name_of_frame, 0*self.coordsys_ref.i + 0*self.coordsys_ref.j + 0*self.coordsys_ref.k)
            
        
        else:
        
            
            # Base axis in respective frame
            self.xvec_in_ref = np.asarray(self.express_vector_in_ref([1,0,0]))
            self.yvec_in_ref = np.asarray(self.express_vector_in_ref([0,1,0]))
            self.zvec_in_ref = np.asarray(self.express_vector_in_ref([0,0,1]))
        
    
            
    
        self.origin_sym_ref = self.coordsys_ref.position_wrt(self.coordsys)           
        self.origin_ref = self.get_origin_in_ref()
        
        
        
        
    def __repr__(self):
        return f"{self.xvec_in_ref} , {self.yvec_in_ref} , {self.zvec_in_ref} "
    
    
    def express_vector_in_ref(self, vector):
        
        if self.frame.name == 'TreatmentRoom':
            return vector
        
        sympy_vector = vector[0]*self.reference.x + vector[1]*self.reference.y + vector[2]*self.reference.z
        return np.asarray([self.frame.x.dot(sympy_vector), self.frame.y.dot(sympy_vector), self.frame.z.dot(sympy_vector)], dtype='float64')
        
    
    def get_origin_in_ref(self):
        
        return np.asarray(self.coordsys.origin.express_coordinates(self.coordsys_ref), dtype='float64')
    
    
    def add_frame_to_plot(self, ax, ref_self, all_longer = 1):

        vec_x = self.xvec_in_ref * all_longer 
        vec_y = self.yvec_in_ref * all_longer
        vec_z = self.zvec_in_ref * all_longer
        
        
        # arrow_prop_dict = dict(mutation_scale=20, arrowstyle='->', shrinkA=0, shrinkB=0)
        
        vecs = [vec_x, vec_y , vec_z]
        colors = ['r', 'g', 'b']
        symbols = [r'$x$', r'$y$', r'$z$']
    
        
        for vec, color, symbol in zip(vecs, colors, symbols): 
            xes = [self.origin_ref[0], self.origin_ref[0] + vec[0]]
            yes = [self.origin_ref[1], self.origin_ref[1] + vec[1]]
            zes = [self.origin_ref[2], self.origin_ref[2] + vec[2]]
            
            # a = Arrow3D(xes, yes, zes , **arrow_prop_dict, color=color)
            # ax.add_artist(a)
            ax.plot(xes, yes, zes , color=color)
            
            text_pos = [self.origin_ref[0] + vec[0]*1.1, self.origin_ref[1] + vec[1]*1.1, self.origin_ref[2] + vec[2]*1.1 ]
            ax.text( text_pos[0], text_pos[1], text_pos[2], symbol, fontsize = 22)
    
    

        
        return ax