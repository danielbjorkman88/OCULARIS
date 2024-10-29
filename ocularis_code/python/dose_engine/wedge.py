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


class Wedge(dose_engine.RayTracer):
    """
    wedge_insert angle = angle of which wedge is inserted towards central axis [degrees]
    wedge_r = wedge tip distance to central axis [mm]
    wedge_angle = angle between wedge plane and aperture plane [degrees]
    
    
    
    """

    def __init__(self, config, density = 1):

        self.config = config

        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        geo = geometry.OptisGeometry(self.config)

        self.central_axis = ref.central_axis

        self.wedge_insert_angle = self.config["wedge_insert_angle"]
        self.wedge_cover = self.config["wedge_cover"]
        self.wedge_angle = self.config["wedge_angle"]

        self.density = density 

        self.traced = False

        self.rays = []
        self.intersects_xes = []
        self.intersects_yes = []
        self.intersects_zes = []

        geo = geometry.OptisGeometry(self.config)

        self.wedge_plane = geo.wedge_plane
        self.wedge_apex_point = geo.wedge_apex_point
        self.aperture_point = geo.aperture_point
        self.aperture_plane = geo.aperture_plane

    def trace_rays(self):

        pair1 = [self.aperture_plane, self.aperture_point]
        pair2 = [self.wedge_plane, self.wedge_apex_point]

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

                    # print("Current")
                    # print("planepoint", planePoint)
                    # print("planeNormal", planeNormal)
                    # print("rayPoint", rayPoint)
                    # print("rayDirection", rayDirection)

                    distance_to_source = math.sqrt(
                        sum(list(map(lambda x, y: (x - y)**2, rayPoint, intersect))))
                    intersects.append(
                        Point(intersect, planeNormal, distance_to_source))
                    self.intersects_xes.append(intersect[0])
                    self.intersects_yes.append(intersect[1])
                    self.intersects_zes.append(intersect[2])

            # intersects1 = self.find_intersects(ray, self.aperture_plane, [self.aperture_point])
            # intersects2 = self.find_intersects(ray, self.wedge_plane, [self.wedge_apex_point])
            # intersects1.extend(intersects2)
            # print("Normal",intersects)
            # print("New", intersects1)

            intersect_aperture = intersects[0]
            intersect_wedge = intersects[1]

            d1, d2 = intersect_aperture.distance_to_source, intersect_wedge.distance_to_source

            if d1 < d2:
                ray.add_depth(self.density*(d2 - d1))

        self.traced = True
