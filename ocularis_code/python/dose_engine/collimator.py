# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 09:25:28 2022

@author: Daniel Björkman 2022
"""


import os
import sys
import inspect
import numpy as np
import math
import matplotlib.pyplot as plt
import warnings

from scipy.optimize import linprog
from utils.concave_hull import ConcaveHull
from shapely.geometry.polygon import Point as Point
from utils.vector_utils import angle
from scipy.spatial import Delaunay
from scipy.interpolate import interp1d

import reference_frames
import geometry
import dose_engine

currentdir = os.path.dirname(os.path.abspath(
    inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)

# colors = ["C0", "C1", "C2", "C3", "C4", "C5",
#           "C6", "C7", "C8", "C9", "C10", "C11"]


def p_is_on_axis(config, p_on_axis):

    ref = reference_frames.Reference_Frame(config, 'TreatmentRoom')

    suggested_vec = p_on_axis - ref.central_axis.S
    suggested_vec = suggested_vec / np.linalg.norm(suggested_vec)

    # print(suggested_vec, ref.central_axis.vec_norm)
    # print(abs(np.dot(suggested_vec, ref.central_axis.vec_norm)))

    epsilon = 1e-6

    return abs(np.dot(suggested_vec, ref.central_axis.vec_norm)) > 1 - epsilon


def p_belongs_to_ray(ray, p):
    suggested_vec = p - ray.S
    suggested_vec = suggested_vec / np.linalg.norm(suggested_vec)

    epsilon = 1e-8

    return abs(np.dot(suggested_vec, ray.vec_norm)) > 1 - epsilon


def perpendicular( a ) :
    b = np.empty_like(a)
    b[0] = -a[1]
    b[1] = a[0]
    return b


class Collimator(dose_engine.RayTracer):
    """
    Defines convex hull at aperture plane representing the collimator shape
    Expands aperture with specified radius
    Projects expanded aperture
    """

    def __init__(self, config):

        self.config = config
        self.target_rays = []

        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        geo = geometry.OptisGeometry(self.config)

        self.central_axis = ref.central_axis
        
        # The ruler length is central to making a concave hull.
        # Too long ruler and the shape will be convex
        self._ruler_length = 3.5


        self.target_rays_intersects_xes = []
        self.target_rays_intersects_yes = []
        self.target_rays_intersects_zes = []

        self.target_rays_intersects_in_plane_xes = []
        self.target_rays_intersects_in_plane_yes = []
        self.target_rays_intersects_in_plane_zes = []

        self.aperture_point = geo.aperture_point
        self.aperture_plane = geo.aperture_plane

        self.is_configured = False
        self.radius = None
        self.expansion = 0
        
        
        self.hull = []
        self.rays = []
        self.aperture_rays = []
        self.aperture_rays_expanded = []

        self.target_rays_aperture_intersect_points_in_plane = []
        self.target_rays_aperture_intersect_points_global = []

        self.points_expanded = []
        self.expanded_hull_points_global = []

    @property
    def ruler_length(self):
        return self._ruler_length

    @ruler_length.setter
    def ruler_length(self, length):
        self._ruler_length = length

    def project_aperture(self, p_on_axis):

        if not p_is_on_axis(self.config, p_on_axis):
            raise ValueError

        bev = reference_frames.Reference_Frame(self.config, 'Beams Eye View')

        points = []
        # ray_points = []
        for ray in self.aperture_rays_expanded:
            p = self.find_intersects(ray, self.aperture_plane, p_on_axis)[0].P
            # ray_points.append(p)
            
            vector = p - p_on_axis
            vector = vector / np.linalg.norm(vector)
            magnitude = math.dist(p, p_on_axis)

            x = - np.dot(bev.xvec_in_ref, magnitude*vector)
            y = np.dot(bev.yvec_in_ref, magnitude*vector)
            points.append((x, y))

        
        # for i in np.linspace(0,10,30):

        hull_in_plane = ConcaveHull(self.config)
        hull_in_plane.loadpoints(points)
        hull_in_plane.calculatehull(self.ruler_length) #
        hull_in_plane.configure(p_on_axis)
        
        # fig = plt.figure()
        
        # for p in points:
        #     plt.scatter(p[0], p[1], color = "C0")
            
        # for p in ray_points:
        #     plt.scatter(p[0], p[1], color = "C1")
        
        # for edge in self.expanded_hull.boundary_edges:
        #     xes = edge[0]
        #     yes = edge[1]
        #     zes = edge[2]
        #     plt.plot(xes, yes, zes, color = "C1", linestyle = "--")
            
        # for p in self.expanded_hull.points:
        #     plt.scatter(p[0], p[1], color = "k")
        
        # for p in self.aperture_hull.boundary_edges:
        #     xes = edge[0]
        #     yes = edge[1]
        #     zes = edge[2]
        #     plt.scatter(p[0], p[1], color = "grey")
        
        
        
        # plt.show()
        
        
        
        
        # print(i, len(points), len(hull_in_plane.boundary_edges))

        return hull_in_plane

    def add_ray(self, target_rays):

        if isinstance(target_rays, list):
            for ray in target_rays:
                self.target_rays.append(ray)
        else:
            self.target_rays.append(target_rays)

    def expand_aperture(self, radius):
        
        self.expansion = round(radius, 4)

        geo = geometry.OptisGeometry(self.config)
        bev = reference_frames.Reference_Frame(self.config, 'Beams Eye View')

        points_expanded = []
        points_interior = []

        xes = self.aperture_hull.boundary.boundary.xy[0]
        yes = self.aperture_hull.boundary.boundary.xy[1]
        
        
        # fig = plt.figure()
        
        # plt.scatter(xes, yes)

        for point_index in range(len(xes)):

            x = xes[point_index]
            y = yes[point_index]
            
            
            # Adds extra hull points in the horizontal plane
            if y == 0:
                fudge_factor = 1.02 # Fudgy factor necessary to make algorithm distinguish the added points
                p1 = np.asarray([x + radius*fudge_factor, 0])
                p2 = np.asarray([x - radius*fudge_factor, 0])
                points_expanded, points_interior, found = _add_points_to_expanded_hull(self, points_expanded, points_interior, p1, p2)
                continue


            p = np.asarray([x, y])

            next_index = point_index + 1

            if next_index >= len(xes):
                next_index -= len(xes)

            x_next = xes[next_index]
            y_next = yes[next_index]
            
            
            
            p_next = np.asarray([x_next , y_next])
            
            vec1 = p_next - p #np.asarray([x_next - x, y_next - y])
            vec1 = vec1/ np.linalg.norm(vec1)
            
            
            vec2 = np.asarray([0, 0])
            subtraction_of_index = 0
            while sum(vec2) == 0.0:
                prev_index = point_index - 1 - subtraction_of_index
                x_prev = xes[prev_index]
                y_prev = yes[prev_index]
                p_prev = np.asarray([x_prev, y_prev])
                vec2 = p_prev - p
                subtraction_of_index += 1
            
            vec2 = vec2/ np.linalg.norm(vec2)
            
            # angle_between = angle(vec1, vec2)

            angle_vec1_to_xaxis = angle(np.asarray([1, 0]), vec1)
            angle_vec2_to_xaxis = angle(np.asarray([1, 0]), vec2)
            
            # full_angle = abs(angle_vec1_to_xaxis) + abs(angle_vec2_to_xaxis)
            
            # if angle_vec2_to_xaxis > angle_vec1_to_xaxis:
            #     new_angle1 = angle_vec2_to_xaxis - full_angle/2
                
            
            
            # # if angle_vec1_to_xaxis > math.pi:
            # #     angle_vec1_to_xaxis = 2*math.pi - angle_vec1_to_xaxis

            new_angle1 = (angle_vec1_to_xaxis + angle_vec2_to_xaxis)/2

            new_angle2 = new_angle1 + math.pi

            p1 = p + radius * \
                np.asarray([math.cos(new_angle1), math.sin(new_angle1)])
            p2 = p + radius * \
                np.asarray([math.cos(new_angle2), math.sin(new_angle2)])
            
            p1_near = p + 0.15 * \
                np.asarray([math.cos(new_angle1), math.sin(new_angle1)])
            p2_near = p + 0.15 * \
                np.asarray([math.cos(new_angle2), math.sin(new_angle2)])
            
            
            
            
            perp_vector1 = perpendicular(vec1)
            perp_vector1 = perp_vector1/np.linalg.norm(perp_vector1)
            
            perp_vector2 = perpendicular(vec2)
            perp_vector2 = perp_vector2/np.linalg.norm(perp_vector2)            
            
            p3 = p + radius * perp_vector1
            p4 = p - radius * perp_vector1
            
            p5 = p + radius * perp_vector2
            p6 = p - radius * perp_vector2            
            
            points = [p1,p2] # [ p3, p4, p5, p6] #p1,p2,
            
            

            
            # found = None
            # for point in points:
                
                
                
            #     shapely_point = Point(point[0], point[1])
            
            
                # if not (np.isnan(point[0]) or np.isnan(point[1])):
                #     found = True
                #     if self.aperture_hull.boundary.contains(shapely_point):
                #         points_interior.append((point[0], point[1]))
                #     else:
                #         points_expanded.append((point[0], point[1]))            
                        
                # plt.scatter(point[0], point[1], color = "C1")
            
            # p1_tmp = p + vec1 + vec2
            # vec = p1_tmp - p
            # vec = vec/ np.linalg.norm(vec)
            # p1 = p + radius *vec
            # p2 = p - radius *vec
            
            # p1 = np.asarray(p1)
            # p2 = np.asarray(p2)
            
            # plt.quiver(p[0], p[1], vec1[0], vec1[1])
            # plt.quiver(p[0], p[1], vec2[0], vec2[1])
            
            # plt.scatter(p1[0], p1[1])
            # plt.scatter(p2[0], p2[1])
            
            indices_outside = _outside_expanded_hull(self, p1_near, p2_near)
            
            
            points_array = np.array(points)
            
            # result = np.delete(points_array, indices_outside, axis=0)
            
            points_outside = points_array[indices_outside]
            
            for p in points_outside:
                points_expanded.append(p)
                
                
            points_inside = points_array[~np.array(indices_outside)]
            
            for p in points_inside:
                points_interior.append(p)
            
            # if index_out > 1:
            #     a = 0
            
            #     fig = plt.figure()
                
            #     plt.scatter(xes, yes)
            #     plt.scatter(p1[0], p1[1], marker = "d")
            #     plt.scatter(p2[0], p2[1], marker = "d")
            #     plt.scatter(p[0],p[1], color = "k")
                
            #     plt.scatter(points_expanded[-1][0],points_expanded[-1][1], color = "C3")
                
            #     plt.scatter(p1_near[0], p1_near[1], color = "r")
            #     plt.scatter(p2_near[0], p2_near[1], color = "r")
                
            #     plt.show()
            
            # if index_out == 2:
            #     warnings.warn("Unexpected value.")
            
            # points_expanded.append(points[index_out])
            # points_interior.append(points[int(not index_out)])
            
            
            
            # points_expanded, points_interior, found = _add_points_to_expanded_hull(self, points_expanded, points_interior, p1, p2)
            
            
            # if p[0] < -0.92 and p[1] > 4.7:
                
            
            
            
            # shapely_point1 = Point(p1[0], p1[1])
            # shapely_point2 = Point(p2[0], p2[1])
            
            
            
            # if not (np.isnan(p1[0]) or np.isnan(p1[1])):
            #     found = True
            #     if self.aperture_hull.boundary.contains(shapely_point1):
            #         points_interior.append((p1[0], p1[1]))
            #     else:
            #         points_expanded.append((p1[0], p1[1]))
            #         # self.point_ray_pair_expanded_aperture.append(
            #         #     (p, self.target_rays[list(self.aperture_hull.boundary_vertices)[point_index]]))

            # if not (np.isnan(p2[0]) or np.isnan(p2[1])):
            #     found = True
            #     if self.aperture_hull.boundary.contains(shapely_point2):
            #         points_interior.append((p2[0], p2[1]))
            #     else:
            #         points_expanded.append((p2[0], p2[1]))
            #         # self.point_ray_pair_expanded_aperture.append(
            #         #     (p, self.target_rays[list(self.aperture_hull.boundary_vertices)[point_index]]))
                    
                    
            # if not found:
            #     # print(p1, p2)
            #     print(vec1, vec2)
            #     print(x,y)
            #     # print("vec",p1_tmp - p, vec)
            #     # # print(p_next - p)
            #     # print(p1_tmp , p , vec1 , vec2)
            #     # print(p + radius *vec)
            #     # print(p , radius ,vec)
                
                
            #     # print(p1, p2)
            #     print("")
            #     raise AssertionError
                
        
        
        # for p in points_expanded:
        #     plt.scatter(p[0], p[1], color = "C4")
        
        # plt.show()
        
        
        
        
        self.points_expanded = points_expanded

        self.expanded_hull = ConcaveHull(self.config)
        self.expanded_hull.loadpoints(self.points_expanded)
        self.expanded_hull.calculatehull(self.ruler_length)
        self.expanded_hull.configure(geo.aperture_point)

        for point, ray in zip(self.points_expanded, self.aperture_rays):
            p_global = geo.aperture_point - bev.xvec_in_ref * \
                point[0] + bev.yvec_in_ref*point[1]
            self.expanded_hull_points_global.append(p_global)
            self.aperture_rays_expanded.append(
                dose_engine.Ray(p_global, p_global + ray.vec_norm * 70.))

        # fig = plt.figure()

        # for ray in self.aperture_rays:
        #     plt.plot(ray.xes, ray.yes, color="C0")

        # for ray in self.aperture_rays_expanded:
        #     plt.plot(ray.xes, ray.yes, color="C1")

        # plt.show()

        # fig = plt.figure()

        # plt.plot(xes, yes, color="C0", marker="*")

        # # xes2 = self.expanded_hull.boundary.boundary.xy[0]
        # # yes2 = self.expanded_hull.boundary.boundary.xy[1]

        # # plt.plot(xes2, yes2, color="C1", marker=">")
        
        # plt.title("Eyeplan collimator shape")

        # plt.show()
        if not len(self.aperture_rays) == len(self.aperture_rays_expanded):
            raise AssertionError
    
        # assert len(self.aperture_rays) == len(self.aperture_rays_expanded)

    def sort_aperture_rays(self):
        """
        Method sorts aperture rays into sequence of aperture_hull
        
        """

        geo = geometry.OptisGeometry(self.config)
        bev = reference_frames.Reference_Frame(self.config, 'Beams Eye View')

        xes = self.aperture_hull.boundary.boundary.xy[0]
        yes = self.aperture_hull.boundary.boundary.xy[1]
        founds = np.zeros(len(xes))

        self.pair = {}

        for idx in range(len(xes)):
            found = False

            x = xes[idx]
            y = yes[idx]
            p = geo.aperture_point - bev.xvec_in_ref*x + bev.yvec_in_ref*y

            max_distance = 100000000
            ray_candidate = []
            for ray in self.aperture_rays:

                if (np.isclose(ray.aperture_intersect , p, rtol = 1e-3)).all():
                    if type(ray.idx) == list:
                        ray.idx = idx
                        self.pair[idx] = ray
                        # print(p_belongs_to_ray(ray,p))
                        found = True
                        founds[idx] = 1
                        break
                # else:
                #     pass
                    # #print("Hello")
                    # distance = math.dist(ray.aperture_intersect, p)

                    # if distance < max_distance:
                    #     print(p_belongs_to_ray(ray,p))
                    #     max_distance = distance
                    #     ray_candidate = ray

            # if type(ray_candidate) == list:
            #     Warning("No match for aperture ray")
            #     continue

            # # Protects from overwriting previous value
            # if type(ray_candidate.idx) == list:
            #     ray_candidate.idx = idx
            #     self.pair[idx] = ray_candidate
            #     found = True

            # if found == False:
            #     Warning("No match for aperture ray")
            #     self.pair[idx] = []

        sorted_aperture_rays = sorted(
            self.aperture_rays, key=lambda ray: ray.idx)

        self.aperture_rays = sorted_aperture_rays

        # fig = plt.figure()

        # for ray, color in zip(self.aperture_rays, colors):
        #     plt.plot(ray.xes, ray.yes, color = color)

        # for idx, color in zip(range(len(xes)), colors):
        #     plt.scatter(xes[idx], yes[idx], color = color)

        # plt.show()

    def define_aperture_hull_from_points(self, points):

        geo = geometry.OptisGeometry(self.config)
        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')
        bev = reference_frames.Reference_Frame(self.config, 'Beams Eye View')

        self.aperture_hull = ConcaveHull(self.config)
        self.aperture_hull.loadpoints(points)
        self.aperture_hull.calculatehull(self.ruler_length)
        self.aperture_hull.configure(geo.aperture_point)

        for x, y in points:
            T = geo.aperture_point - bev.xvec_in_ref*x + bev.yvec_in_ref*y
            ray = dose_engine.Ray(ref.central_axis.S, T)
            extended_T = T + ray.vec_norm*(math.dist(geo.aperture_point, geo.central_axis.T))
            ray = dose_engine.Ray(ref.central_axis.S, extended_T)
            ray.aperture_intersect = T
            self.aperture_rays.append(ray)

        self.sort_aperture_rays()

        self.is_configured = True

    def define_aperture_hull(self):

        self.intersects = []

        geo = geometry.OptisGeometry(self.config)
        bev = reference_frames.Reference_Frame(self.config, 'Beams Eye View')

        for ray in self.target_rays:

            tmp = self.find_intersects(
                ray, self.aperture_plane, self.aperture_point)
            ray.aperture_intersect = np.asarray(
                [tmp[0].P[0], tmp[0].P[1], tmp[0].P[2]])

            self.intersects.append(tmp)

        self.hull_points = np.ones([len(self.intersects), 3])
        for i in range(len(self.intersects)):
            val = self.hull_points[i, 0:]*list(self.intersects[i][0].P)

            self.hull_points[i, 0:] = val

        for i in range(len(self.intersects)):
            self.target_rays_intersects_xes.append(self.intersects[i][0].P[0])
            self.target_rays_intersects_yes.append(self.intersects[i][0].P[1])
            self.target_rays_intersects_zes.append(self.intersects[i][0].P[2])

        self.target_rays_aperture_intersect_points_global = list(zip(
            self.target_rays_intersects_xes, self.target_rays_intersects_yes, self.target_rays_intersects_zes))

        for point in self.target_rays_aperture_intersect_points_global:
            p = np.asarray(point)
            vec = p - geo.aperture_point

            p_in_plane = [- vec.dot(bev.xvec_in_ref), vec.dot(bev.yvec_in_ref)]
            self.target_rays_intersects_in_plane_xes.append(p_in_plane[0])
            self.target_rays_intersects_in_plane_yes.append(p_in_plane[1])

            self.target_rays_aperture_intersect_points_in_plane.append(
                p_in_plane)

        assert len(self.target_rays_aperture_intersect_points_in_plane) > 1
        
        
        # This block adds target rays in the horizontal and vertical planes for better comparisons for these planes
        additional_intersects, extra_rays = _find_extra_rays(self.config, self.ruler_length,  self.target_rays_aperture_intersect_points_in_plane)
        self.target_rays.extend(extra_rays)
        self.target_rays_aperture_intersect_points_in_plane.extend(additional_intersects)
        
        self.aperture_hull = ConcaveHull(self.config)
        self.aperture_hull.loadpoints(
            self.target_rays_aperture_intersect_points_in_plane)
        self.aperture_hull.calculatehull(self.ruler_length)
        self.aperture_hull.configure(geo.aperture_point)

        for idx in self.aperture_hull.boundary_vertices:
            self.aperture_rays.append(self.target_rays[idx])

        self.sort_aperture_rays()

        self.is_configured = True




def _add_points_to_expanded_hull(self, points_expanded, points_interior, p1, p2):
    
        found = False
    
        shapely_point1 = Point(p1[0], p1[1])
        shapely_point2 = Point(p2[0], p2[1])
        
        in_first = False
        in_second = False
        
        if not (np.isnan(p1[0]) or np.isnan(p1[1])):
            found = True
            if self.aperture_hull.boundary.contains(shapely_point1):
                points_interior.append((p1[0], p1[1]))
            else:
                points_expanded.append((p1[0], p1[1]))
                in_first = True
                # self.point_ray_pair_expanded_aperture.append(
                #     (p, self.target_rays[list(self.aperture_hull.boundary_vertices)[point_index]]))
            
            
            
        if not (np.isnan(p2[0]) or np.isnan(p2[1])):
            found = True
            if self.aperture_hull.boundary.contains(shapely_point2):
                points_interior.append((p2[0], p2[1]))
            else:
                points_expanded.append((p2[0], p2[1]))
                in_second = True
                # self.point_ray_pair_expanded_aperture.append(
                #     (p, self.target_rays[list(self.aperture_hull.boundary_vertices)[point_index]]))
            
            
            
        if in_first and in_second:
            print("In both")
            
        return points_expanded, points_interior, found

def _outside_expanded_hull(self, pA, pB):
    
        found = False
    
        shapely_point1 = Point(pA[0], pA[1])
        shapely_point2 = Point(pB[0], pB[1])
        
        first_outside = False
        second_outside = False
        
        if not (np.isnan(pA[0]) or np.isnan(pA[1])):
            found = True
            if self.aperture_hull.boundary.contains(shapely_point1):
                # points_interior.append((pA[0], pA[1]))
                pass
            else:
                # points_expanded.append((pA[0], pA[1]))
                first_outside = True
                # self.point_ray_pair_expanded_aperture.append(
                #     (p, self.target_rays[list(self.aperture_hull.boundary_vertices)[point_index]]))
            
            
            
        if not (np.isnan(pB[0]) or np.isnan(pB[1])):
            found = True
            if self.aperture_hull.boundary.contains(shapely_point2):
                pass
                # points_interior.append((pB[0], pB[1]))
            else:
                # points_expanded.append((pB[0], pB[1]))
                second_outside = True
                # self.point_ray_pair_expanded_aperture.append(
                #     (p, self.target_rays[list(self.aperture_hull.boundary_vertices)[point_index]]))
            
        # print(first_outside , second_outside)
        
        
        # indices_outside = np.zeros(2)
        # indices_outside[0] = int(first_outside)
        # indices_outside[1] = int(second_outside)
        # return indices_outside.astype(int)
    
    
        return [first_outside, second_outside]
        
        # if first_outside and not second_outside:
        #     return 0
        # elif second_outside and not first_outside:
        #     return 1
        # else:
        #     return 2
            
        # if in_first and in_second:
        #     print("In both")
            
        # return points_expanded, points_interior, found


def _find_extra_rays(config, ruler_length,  target_rays_aperture_intersect_points_in_plane):
    
    
    points = np.asarray(target_rays_aperture_intersect_points_in_plane)
    
    
    hull = ConcaveHull(config)
    hull.loadpoints(points)
    hull.calculatehull(ruler_length)
    
    xes, yes = hull.boundary.boundary.xy
    
    
    additional_intersects = []
    for i in range(len(xes) - 1):
        
        xs = [xes[i], xes[i+1]]
        ys = [yes[i], yes[i+1]]
        
        f = interp1d(xs, ys)
        f_rev = interp1d(ys, xs)
        try:
            # intersects_horizontal.append([f_rev(0), 0])
            additional_intersects.append([f_rev(0).cumsum()[0], 0])
        except:
            pass
        try:
            # intersects_vertical.append([0, f(0)])
            additional_intersects.append([0, f(0).cumsum()[0]])
        except:
            pass    
    
        


    geo = geometry.OptisGeometry(config)
    bev = reference_frames.Reference_Frame(config, 'Beams Eye View')
    
    
    extra_rays = []
    for p in additional_intersects:
        p_global = geo.aperture_point - bev.xvec_in_ref * \
            p[0] + bev.yvec_in_ref*p[1]
        
        ray = dose_engine.Ray(geo.central_axis.S, p_global)
        ray.aperture_intersect = ray.T
        extra_rays.append(ray)
    
    return additional_intersects, extra_rays
        