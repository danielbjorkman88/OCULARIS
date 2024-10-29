#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 15:35:06 2021

@author: bjoerk_c
"""



from scipy.special import erf
import math
from dose_engine.algo import Algo
from decorators import timer
import numpy as np
import jax.numpy as jnp
import jax

# def collimation_modifer(r, r_half, beta):
#     return 0.5 - 0.5*erf(beta* (r - r_half)/(math.sqrt(2)))


def collimation_modifer(dist, beta):
    return 0.5 + 0.5*jax.scipy.erf(beta* dist/(jnp.sqrt(2)))


class BroadBeam(Algo):
    
    
    

    
    @timer
    def calc(self):
        
        
        if self.config["Image"] == "xslice":
            self.logger.info("Calculating 2D distribution for Broad Beam...")
        else:
            self.logger.info("Calculating 3D distribution for Broad Beam...")
            

        
        
        # x_mid = self.config["Mesh dimensions"][0]/2
        # y_mid = self.config["Mesh dimensions"][1]/2
        
        
        indices = jnp.array(self.raytracer.traced_indices)
        x, y, z = indices[:, 0], indices[:, 1], indices[:, 2]
        z_rad = self.raytracer.depth[x,y,z]
        
        # if z_rad != 0:
            
        # x_pos = self.medium.resolution[0]*(x + 0.5)
        # y_pos = self.medium.resolution[1]*(y + 0.5)
        
        # r = math.sqrt( (x_pos - x_mid)**2 + (y_pos - y_mid)**2 ) 
        # r = self.radii[x,y,z]
        
        dist = self.dist_aperture[x,y,z]

        # r_half = self.f_half(z_rad)
        beta = self.f_beta(z_rad)                    
        
        
        non_nans = jnp.logical_not(jnp.isnan(beta))
        z_rad = z_rad[non_nans]
        beta = beta[non_nans]
        dist = dist[non_nans]
        x = x[non_nans]
        y = y[non_nans]
        z = z[non_nans]
        self.dose[x,y,z] = self.SOBP_f(z_rad) * collimation_modifer(dist, beta)

        
        self.dose_calculated = True
        self.logger.info(f"Calculation finished for {len(self.raytracer.traced_indices)} voxels")




            











