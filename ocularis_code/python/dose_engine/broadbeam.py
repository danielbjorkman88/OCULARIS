#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  4 15:35:06 2021

@author: bjoerk_c
"""


from scipy.special import erf
import numpy as np

from dose_engine.algo import Algo


def collimation_modifer(dist, beta):
    return 0.5 + 0.5*erf(beta * dist/(np.sqrt(2)))


class BroadBeam(Algo):

    def calc(self):
        
        if self.raytracer_configured == False:
            print("Configure raytracer first")
            return
        
        if self.config["Image"] == "xslice":
            self.logger.info("Calculating 2D distribution for Broad Beam...")
        else:
            self.logger.info("Calculating 3D distribution for Broad Beam...")

        indices = self.raytracer.traced_indices
        x, y, z = indices[:, 0], indices[:, 1], indices[:, 2]
        z_rad = self.raytracer.depth[x, y, z]

        dist = self.dist_aperture[x, y, z]

        beta = self.f_beta(z_rad)/self.penumbra_factor

        non_nans = np.logical_not(np.isnan(beta))
        z_rad = z_rad[non_nans]
        beta = beta[non_nans] #* (1/1.875) # *0.75
        dist = dist[non_nans]
        x = x[non_nans]
        y = y[non_nans]
        z = z[non_nans]
        self.dose[x, y, z] = self.SOBP_f(
            z_rad) * collimation_modifer(dist, beta)
        
        self.dose[self.dose < 0] = 0

        self.dose_calculated = True
        self.logger.info(
            f"Calculation finished for {len(self.raytracer.traced_indices)} voxels")
