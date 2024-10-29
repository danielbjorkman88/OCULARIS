# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 09:25:28 2022

@author: bjoerk_c
"""

import reference_frames
import geometry
import os
import sys
import inspect
import numpy as np
import math
import dose_engine

currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)


class Compensator(dose_engine.RayTracer):
    """
    https://diegoinacio.github.io/computer-vision-notebooks-page/pages/ray-intersection_sphere.html
    
    """

    def __init__(self, config, sphere_centre, radius, density=1):

        self.config = config

        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        geo = geometry.OptisGeometry(self.config)

        self.central_axis = ref.central_axis

        self.density = density

        self.traced = False

        self.rays = []
        self.intersects_xes = []
        self.intersects_yes = []
        self.intersects_zes = []

        geo = geometry.OptisGeometry(self.config)

        self.aperture_point = geo.aperture_point
        self.aperture_plane = geo.aperture_plane

        self.sphere_centre = sphere_centre
        self.radius = radius

    def trace_rays(self):

        O = self.central_axis.S
        Cs = self.sphere_centre
        OC_ = Cs - O

        for ray in self.rays:

            e_ = ray.vec_norm

            t = np.dot(OC_, ray.vec_norm)

            Pe = O + e_*t

            d = np.linalg.norm(Pe - Cs)

            sphere_intersect = None

            if(d > self.radius):
                pass
                #print("No intersection!")

            elif(d == self.radius):
                # Ps = Pe
                sphere_intersect = Pe
                # print(f'Intersection at {Ps}')

            else:
                i = (self.radius**2 - d**2)**0.5
                sphere_intersect = O + e_*(t - i)

            aperture_plane_intersect = ray.find_intersect_with_plane(
                self.aperture_plane, self.aperture_point)

            sphere_plane_intersect = ray.find_intersect_with_plane(
                self.aperture_plane, Cs)

            if sphere_intersect is not None:
                ray_exit_point = sphere_intersect
            else:
                ray_exit_point = sphere_plane_intersect

            ray.ray_exit_point = ray_exit_point
            ray.aperture_plane_intersect = aperture_plane_intersect

            d1, d2 = math.dist(aperture_plane_intersect,
                               O), math.dist(ray_exit_point, O)

            if d1 < d2:
                ray.add_depth(self.density*(d2 - d1))

        self.traced = True
