# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 09:25:28 2022

@author: bjoerk_c
"""
from dose_engine.raytracer import Point
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


class SkinPlane(dose_engine.RayTracer):
    """

    
    
    """

    def __init__(self, config, medium):

        self.config = config
        self.medium = medium

        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        geo = geometry.OptisGeometry(self.config)

        self.central_axis = ref.central_axis

        self.density = 1

        self.traced = False

        self.rays = []
        self.intersects_xes = []
        self.intersects_yes = []
        self.intersects_zes = []

        geo = geometry.OptisGeometry(self.config)

        self.skin_plane_point = geo.skin_plane_point
        self.skin_plane = geo.skin_plane

        self.mesh_plane = self.central_axis.vec_norm
        self.mesh_point = self.medium.mesh_apex


    def trace_rays(self):

        pair1 = [self.skin_plane, self.skin_plane_point]
        pair2 = [self.mesh_plane, self.mesh_point]

        for ray in self.rays:
            rayDirection = ray.vec_norm
            rayPoint = ray.S

            intersects = []

            pairs = [pair1, pair2]

            for pair in pairs:
                planeNormal = pair[0]
                planePoint = pair[1]

                # epsilon=1e-6

                ndotu = planeNormal.dot(rayDirection)

                if ndotu != 0:
                    w = rayPoint - planePoint
                    si = -planeNormal.dot(w) / ndotu
                    intersect = w + si * rayDirection + planePoint

                    distance_to_source = math.sqrt(
                        sum(list(map(lambda x, y: (x - y)**2, rayPoint, intersect))))
                    intersects.append(
                        Point(intersect, planeNormal, distance_to_source))
                    self.intersects_xes.append(intersect[0])
                    self.intersects_yes.append(intersect[1])
                    self.intersects_zes.append(intersect[2])

            intersect_aperture = intersects[0]
            intersect_wedge = intersects[1]

            d1, d2 = intersect_aperture.distance_to_source, intersect_wedge.distance_to_source

            if d1 < d2:
                ray.add_depth(self.density*(d2 - d1))

        self.traced = True
