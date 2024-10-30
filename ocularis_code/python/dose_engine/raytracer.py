# -*- coding: utf-8 -*-
"""
Created on Fri Oct 15 11:39:39 2021

@author: bjoerk_c

This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


import numpy as np
import time
import math
# Ray Tracer


class Point():  # Point of intersect

    """
    Memory efficency improvement expected from making this class slotted
    
    """

    def __init__(self, P, plane, distance_to_source):
        self.P = P
        self.plane = plane  # plane of intersect
        self.distance_to_source = distance_to_source

    def __repr__(self):
        return "{} , {} , {}".format(self.P, self.plane, self.distance_to_source)



def _vectorized_find_intersects_with_planes(S, vecs, planeNormal, planePoints):

    # N = len(planePoints)
    
    #ndotu = planeNormal.dot(vecs)
    ndotu_s = np.inner(planeNormal, vecs)
    
    # mat = np.zeros([len(vecs), len(planePoints)])

    # if ndotu == 0:
    #     return  np.empty(shape=[0, 3])
        
    w_s = S - planePoints
    p = - (planeNormal[0]*w_s[0:,0] + planeNormal[1]*w_s[0:,1] + planeNormal[2]*w_s[0:,2])
    
    p_s = np.tile(p, (len(vecs),1))
    
    si_s = p_s /ndotu_s[:,None]
    
    w_s_planepoints = w_s + planePoints
    
    
    ray_intersects = vecs*si_s[0,0:][0:,None][:,np.newaxis]
    
    ray_intersects += w_s_planepoints[0]
    
    return ray_intersects


class RayTracer():
    """
    Ray Tracer 
    
    Locates intersections of ray with the planes seperating the voxels in depth cube. 
    These intersections are sorted by distance to source and then followed to assign relevant voxels with accumulated depth.
    
    Could be uppgraded into a Siddon raytracer for improved performance. 
    This would be done by finding the ray's intersections with all planes in a single step.
    
    Potential performance improvement could be achieved from integrating pyvista:
        https://docs.pyvista.org/examples/01-filter/poly-ray-trace.html

    
    
    """

    def __init__(self, config, medium):

        self.rays = []
        self.config = config

        self.depth = np.zeros(
            [self.config["Nvoxels"][0], self.config["Nvoxels"][1], self.config["Nvoxels"][2]])

        self.medium = medium

        self.checked_values = np.zeros(
            [self.config["Nvoxels"][0], self.config["Nvoxels"][1], self.config["Nvoxels"][2]])

        self.traced_indices = None

        self.resolution = [self.config["Mesh dimensions"][0]/self.config["Nvoxels"][0], config["Mesh dimensions"]
                           [1]/self.config["Nvoxels"][1], config["Mesh dimensions"][2]/self.config["Nvoxels"][2]]

        self.intersects_xes = []
        self.intersects_yes = []
        self.intersects_zes = []
        
        
        self.intersects_found = False

    def add_ray(self, rays):

        if isinstance(rays, list):
            for ray in rays:
                self.rays.append(ray)
        else:
            self.rays.append(rays)

    def find_intersects(self, ray, planeNormal, planePoints):

        intersects = []

        if isinstance(planePoints, list):
            planePoint = planePoints[0]
        else:
            planePoint = planePoints

        ndotu = planeNormal.dot(ray.vec_norm)

        if ndotu != 0:
            w = ray.S - planePoint
            si = -planeNormal.dot(w) / ndotu
            intersect = w + si * ray.vec_norm + planePoint

            # print("planepoint", planePoint)
            # print("planeNormal", planeNormal)
            # print("rayPoint", ray.S)
            # print("rayDirection", ray.vec_norm)
            # print("intersect", intersect)

            distance_to_source = np.sqrt(
                sum(list(map(lambda x, y: (x - y)**2, ray.S, intersect))))
            intersects.append(
                Point(intersect, planeNormal, distance_to_source))
            # self.intersects_xes.append(intersect[0])
            # self.intersects_yes.append(intersect[1])
            # self.intersects_zes.append(intersect[2])

        # print(intersects)

        return intersects








    # def find_intersects2(self, ray, planeNormal, planePoints):

    #     intersects = []

    #     planePoint = planePoints[0:2]

    #     print("New")

    #     # print("planeNormal", planeNormal)
    #     # print("planePoints",planePoints)

    #     ndotu = planeNormal.dot(ray.vec_norm)

    #     print("ndotu", ndotu)

    #     if ndotu != 0:
    #         print(ray.S, planePoint)
    #         w = ray.S - planePoint
    #         w = list(map(lambda x, y: x - y, ray.S, planePoints))
    #         print("w", w)
    #         si = -planeNormal.dot(w) / ndotu
    #         intersect = w + si * ray.vec_norm + planePoint

    #         # print("New")
    #         # print("planepoint", planePoint)
    #         # print("planeNormal", planeNormal)
    #         # print("rayPoint", ray.S)
    #         # print("rayDirection", ray.vec_norm)

    #         distance_to_source = np.sqrt(
    #             sum(list(map(lambda x, y: (x - y)**2, ray.S, intersect))))
    #         intersects.append(
    #             Point(intersect, planeNormal, distance_to_source))
    #         self.intersects_xes.append(intersect[0])
    #         self.intersects_yes.append(intersect[1])
    #         self.intersects_zes.append(intersect[2])

    #     return intersects

    def trace_rays(self):
        
        traced_indices = []
        
        xPlane = np.array([1, 0, 0])
        yPlane = np.array([0, 1, 0])
        zPlane = np.array([0, 0, 1])


        # planePoints_X = np.linspace(
        #     self.medium.mesh_origin[0], self.medium.mesh_apex[0], num=self.depth.shape[0] + 1)
        # planePoints_Y = np.linspace(
        #     self.medium.mesh_origin[1], self.medium.mesh_apex[1], num=self.depth.shape[1] + 1)
        # planePoints_Z = np.linspace(
        #     self.medium.mesh_origin[2], self.medium.mesh_apex[2], num=self.depth.shape[2] + 1)
        
        
        planePoints_x = np.zeros((self.depth.shape[0] + 1, 3))
        planePoints_x[0:, 0] = np.linspace(
            self.medium.mesh_origin[0], self.medium.mesh_apex[0], num=self.depth.shape[0] + 1)
        
        planePoints_y = np.zeros((self.depth.shape[1] + 1, 3))
        planePoints_y[0:, 1] = np.linspace(
            self.medium.mesh_origin[1], self.medium.mesh_apex[1], num=self.depth.shape[1] + 1)

        planePoints_z = np.zeros((self.depth.shape[2] + 1, 3))
        planePoints_z[0:, 2] = np.linspace(
            self.medium.mesh_origin[2], self.medium.mesh_apex[2], num=self.depth.shape[2] + 1)        
        
        
        skin_plane_z = self.config["skin_plane_point"][2]
        

        for ray in self.rays:
            
            # S = ray.S
            # T = ray.T

            intersects = []

            # rayDirection = ray.vec_norm
            # rayPoint = ray.S
            
            # start1 = time.time()
            z_intersects = _find_intersects_with_planes(
                ray, np.asarray([0, 0, 1]), planePoints_z)
            x_intersects = _find_intersects_with_planes(
                ray, np.asarray([1, 0, 0]), planePoints_x)
            y_intersects = _find_intersects_with_planes(
                ray, np.asarray([0, 1, 0]), planePoints_y)
            
            ray.z_intersects = z_intersects

            x_planes = np.tile(xPlane, (len(x_intersects), 1))
            y_planes = np.tile(yPlane, (len(y_intersects), 1))
            z_planes = np.tile(zPlane, (len(z_intersects), 1))

            con = np.concatenate((x_intersects, y_intersects, z_intersects))
            planes = np.concatenate((x_planes, y_planes, z_planes))

            indices = np.where((self.medium.minX <= con[0:, 0]) & (self.medium.maxX >= con[0:, 0]) & (self.medium.minY <= con[0:, 1]) & (
                self.medium.maxY >= con[0:, 1]) & (self.medium.minZ <= con[0:, 2]) & (self.medium.maxZ >= con[0:, 2]))
            
            con = con[indices]
            planes = planes[indices]
            
            
            # Temp
            self.my_intersects = con

            
            
            

            intersects, planes_sorted = _sort_by_distance_to_source(
                con, planes, ray.S)




            if len(intersects) == 0:
                continue
            
            
            self.intersects_found = True
                

            a0 = intersects[0]


            idx = np.floor(
                (a0[0] - self.medium.mesh_origin[0])/self.medium.resolution[0]).astype(int)
            idy = np.floor(
                (a0[1] - self.medium.mesh_origin[1])/self.medium.resolution[1]).astype(int)
            idz = np.floor(
                (a0[2] - self.medium.mesh_origin[2])/self.medium.resolution[2]).astype(int)
            


            # Adresses the case of first intersect is at the mesh faces associated with the mesh apex point
            if np.fmod(a0[0], self.resolution[0]) == 0 and ray.vec_norm[0] < 0:
                idx -= 1
            if np.fmod(a0[1], self.resolution[1]) == 0 and ray.vec_norm[1] < 0:
                idy -= 1
            if np.fmod(a0[2], self.resolution[2]) == 0 and ray.vec_norm[2] < 0:
                idz -= 1
            

            
            for intersect, plane in zip(intersects[1:], planes_sorted[1:]):

                # # Ends tracing for current ray in case the accumulated depth exceeds the max depth
                # if ray.acc_depth > 36:
                #     break

                a1 = intersect
                # plane = intersect.plane
        

                if 0 <= idx < self.depth.shape[0] and 0 <= idy < self.depth.shape[1] and 0 <= idz < self.depth.shape[2]:
                    
                    # print(self.config["Image"] == "3D")
                    
                    traced_indices.append([idx, idy, idz])
                    
                    # if self.config["Image"] == "3D": # 3D, xslice, yslice, zslice
                    #     print([idx, idy, idz])
                    #     traced_indices.append([idx, idy, idz])        
                    # elif self.config["Image"] == "xslice":
                    #     if idy == self.config["Slice"][1]:
                    #         traced_indices.append([idx, idy, idz])
                            
                    # D = np.sqrt(  (self.resolution[0]*(a1[0]-a0[0]))**2
                    #               + (self.resolution[1]*(a1[1]-a0[1]))**2
                    #               + (self.resolution[2]*(a1[2]-a0[2]))**2 )

                    D = np.sqrt(((a1[0]-a0[0]))**2 +
                                  ((a1[1]-a0[1]))**2 + ((a1[2]-a0[2]))**2)
                    
                    water_eq_density = self.medium.medium[idx, idy, idz]
                    
                    if ray.acc_depth > 2:
                        apa = 1
                    
                    if 5 < ray.T[0] < 6 and -2 < ray.T[1] < -1 and idz == 59 - 5:
                        apa = 1
                        
                    dist1 = math.dist(ray.surface_interaction_point, ray.S)
                    dist2 = math.dist(a1, ray.S)
                    
                    if ray.surface_behind_mesh:
                        # print("ray.surface_behind_mesh", ray.surface_behind_mesh)
                        # self.depth[idx, idy, idz] = 0
                        behind_distance = ray.surface_behind_skin_plane_by
                        if (ray.acc_depth + D/2) * water_eq_density >= behind_distance:
                            self.depth[idx, idy, idz] = (
                                ray.acc_depth + D/2) * water_eq_density - behind_distance
                            
                        ray.add_depth(D * water_eq_density)
                    
                    # Proper tracing starts after the skin plane
                    elif a0[2] < skin_plane_z - self.resolution[2]/2:
                    
                   
                        if not self.checked_values[idx, idy, idz]:
                            self.depth[idx, idy, idz] = (
                                ray.acc_depth + D/2) * water_eq_density
                            self.checked_values[idx, idy, idz] = 1
                        elif self.checked_values[idx, idy, idz] and ray.acc_depth + D/2 < self.depth[idx, idy, idz]:
                            self.depth[idx, idy, idz] = (
                                ray.acc_depth + D/2) * water_eq_density
                        

                        ray.add_depth(D * water_eq_density)
                        
                    elif a0[2] > skin_plane_z - self.resolution[2]/2:
                        self.depth[idx, idy, idz] = (
                                ray.acc_depth_before_skin_plane + D/2) * water_eq_density            
                            
                        ray.add_depth_before_skin_plane(D * water_eq_density)
                    else:
                        print("warning")
                        
                        
                    # dist3 = math.dist(a1, ray.surface_interaction_point)

                    # if dist2 < dist1:
                    #     # dist1 = math.dist(ray.surface_interaction_point, ray.S)
                    #     # dist2 = math.dist(a1, ray.S)
                    #     # dist3 = math.dist(a1, ray.surface_interaction_point)
                    #     # if dist2 < dist1:
                    #     #     continue
                    
                    #     self.depth[idx, idy, idz] = (
                    #         ray.acc_depth_before_skin_plane + D/2) * water_eq_density            
                        
                    #     ray.add_depth_before_skin_plane(D * water_eq_density)
                    
                        
                    #     # behind_distance = ray.surface_behind_mesh_by
                    #     # if (ray.acc_depth_before_skin_plane + D/2) * water_eq_density < behind_distance:
                    #     #     self.depth[idx, idy, idz] = (
                    #     #         ray.acc_depth_before_skin_plane + D/2) * water_eq_density - behind_distance
                    # else:
                    #     if ray.surface_behind_skin_plane:
                    #         # self.depth[idx, idy, idz] = 0
                    #         behind_distance = ray.surface_behind_skin_plane_by
                    #         if (ray.acc_depth + D/2) * water_eq_density >= behind_distance:
                    #             self.depth[idx, idy, idz] = (
                    #                 ray.acc_depth + D/2) * water_eq_density - behind_distance
                    #         #ray.add_depth(D * water_eq_density)
                            
                    #     else:
                            
                    #         if not self.checked_values[idx, idy, idz]:
                    #             self.depth[idx, idy, idz] = (
                    #                 ray.acc_depth + D/2) * water_eq_density
                    #             self.checked_values[idx, idy, idz] = 1
                                
                    #         elif self.checked_values[idx, idy, idz] and ray.acc_depth + D/2 < self.depth[idx, idy, idz]:
                    #             self.depth[idx, idy, idz] = (
                    #                 ray.acc_depth + D/2) * water_eq_density
                            
    
                        # ray.add_depth(D * water_eq_density)

                if (plane == xPlane).all() and idx <= self.depth.shape[0]:
                    if ray.vec[0] > 0:
                        idx += 1
                    elif ray.vec[0] < 0:
                        idx -= 1
                if (plane == yPlane).all() and idy <= self.depth.shape[1]:
                    if ray.vec[1] > 0:
                        idy += 1
                    elif ray.vec[1] < 0:
                        idy -= 1
                if (plane == zPlane).all() and idz <= self.depth.shape[2]:
                    if ray.vec[2] > 0:
                        idz += 1
                    elif ray.vec[2] < 0:
                        idz -= 1

                a0 = a1

        self.traced_indices = np.unique(np.array(traced_indices), axis = 0)


    def trace_rays_new(self):
        
        traced_indices = []
        
        xPlane = np.array([1, 0, 0])
        yPlane = np.array([0, 1, 0])
        zPlane = np.array([0, 0, 1])


        # planePoints_X = np.linspace(
        #     self.medium.mesh_origin[0], self.medium.mesh_apex[0], num=self.depth.shape[0] + 1)
        # planePoints_Y = np.linspace(
        #     self.medium.mesh_origin[1], self.medium.mesh_apex[1], num=self.depth.shape[1] + 1)
        # planePoints_Z = np.linspace(
        #     self.medium.mesh_origin[2], self.medium.mesh_apex[2], num=self.depth.shape[2] + 1)
        
        
        planePoints_x = np.zeros((self.depth.shape[0] + 1, 3))
        planePoints_x[0:, 0] = np.linspace(
            self.medium.mesh_origin[0], self.medium.mesh_apex[0], num=self.depth.shape[0] + 1)
        
        planePoints_y = np.zeros((self.depth.shape[1] + 1, 3))
        planePoints_y[0:, 1] = np.linspace(
            self.medium.mesh_origin[1], self.medium.mesh_apex[1], num=self.depth.shape[1] + 1)

        planePoints_z = np.zeros((self.depth.shape[2] + 1, 3))
        planePoints_z[0:, 2] = np.linspace(
            self.medium.mesh_origin[2], self.medium.mesh_apex[2], num=self.depth.shape[2] + 1)        
        
        
        S = self.rays[0].S
        vecs = np.asarray([ray.vec_norm for ray in self.rays])
        
        self.planePoints_x = planePoints_x
        self.planePoints_y = planePoints_y
        self.planePoints_z = planePoints_z
        
        intersect_array_x = _vectorized_find_intersects_with_planes(S, vecs, xPlane, planePoints_x)
        intersect_array_y = _vectorized_find_intersects_with_planes(S, vecs, yPlane, planePoints_y)
        intersect_array_z = _vectorized_find_intersects_with_planes(S, vecs, zPlane, planePoints_z)
        
        self.intersect_array_x = intersect_array_x
        self.intersect_array_y = intersect_array_y
        self.intersect_array_z = intersect_array_z
        
        
        for ray in self.rays:

            # S = ray.S
            # T = ray.T

            intersects = []

            # rayDirection = ray.vec_norm
            # rayPoint = ray.S
            
            # start1 = time.time()
            z_intersects = _find_intersects_with_planes(
                ray, np.asarray([0, 0, 1]), planePoints_z)
            x_intersects = _find_intersects_with_planes(
                ray, np.asarray([1, 0, 0]), planePoints_x)
            y_intersects = _find_intersects_with_planes(
                ray, np.asarray([0, 1, 0]), planePoints_y)

            x_planes = np.tile(xPlane, (len(x_intersects), 1))
            y_planes = np.tile(yPlane, (len(y_intersects), 1))
            z_planes = np.tile(zPlane, (len(z_intersects), 1))
            
            ray.x_intersects = x_intersects
            ray.y_intersects = y_intersects
            ray.z_intersects = z_intersects

            con = np.concatenate((x_intersects, y_intersects, z_intersects))
            planes = np.concatenate((x_planes, y_planes, z_planes))

            indices = np.where((self.medium.minX <= con[0:, 0]) & (self.medium.maxX >= con[0:, 0]) & (self.medium.minY <= con[0:, 1]) & (
                self.medium.maxY >= con[0:, 1]) & (self.medium.minZ <= con[0:, 2]) & (self.medium.maxZ >= con[0:, 2]))
            
            con = con[indices]
            planes = planes[indices]
            
            
            # Temp
            self.my_intersects = con

            
            
            

            intersects, planes_sorted = _sort_by_distance_to_source(
                con, planes, ray.S)




            if len(intersects) == 0:
                continue
            
            
            self.intersects_found = True
                

            a0 = intersects[0]


            idx = np.floor(
                (a0[0] - self.medium.mesh_origin[0])/self.medium.resolution[0]).astype(int)
            idy = np.floor(
                (a0[1] - self.medium.mesh_origin[1])/self.medium.resolution[1]).astype(int)
            idz = np.floor(
                (a0[2] - self.medium.mesh_origin[2])/self.medium.resolution[2]).astype(int)
            


            # Adresses the case of first intersect is at the mesh faces associated with the mesh apex point
            if np.fmod(a0[0], self.resolution[0]) == 0 and ray.vec_norm[0] < 0:
                idx -= 1
            if np.fmod(a0[1], self.resolution[1]) == 0 and ray.vec_norm[1] < 0:
                idy -= 1
            if np.fmod(a0[2], self.resolution[2]) == 0 and ray.vec_norm[2] < 0:
                idz -= 1
            

            
            for intersect, plane in zip(intersects[1:], planes_sorted[1:]):

                # # Ends tracing for current ray in case the accumulated depth exceeds the max depth
                # if ray.acc_depth > 36:
                #     break

                a1 = intersect
                # plane = intersect.plane
        

                if 0 <= idx < self.depth.shape[0] and 0 <= idy < self.depth.shape[1] and 0 <= idz < self.depth.shape[2]:
                    
                    # print(self.config["Image"] == "3D")
                    
                    traced_indices.append([idx, idy, idz])
                    
                    # if self.config["Image"] == "3D": # 3D, xslice, yslice, zslice
                    #     print([idx, idy, idz])
                    #     traced_indices.append([idx, idy, idz])        
                    # elif self.config["Image"] == "xslice":
                    #     if idy == self.config["Slice"][1]:
                    #         traced_indices.append([idx, idy, idz])
                            
                    # D = np.sqrt(  (self.resolution[0]*(a1[0]-a0[0]))**2
                    #               + (self.resolution[1]*(a1[1]-a0[1]))**2
                    #               + (self.resolution[2]*(a1[2]-a0[2]))**2 )

                    D = np.sqrt(((a1[0]-a0[0]))**2 +
                                  ((a1[1]-a0[1]))**2 + ((a1[2]-a0[2]))**2)
                    
                    water_eq_density = self.medium.medium[idx, idy, idz]

                    if not self.checked_values[idx, idy, idz]:
                        self.depth[idx, idy, idz] = (
                            ray.acc_depth + D/2) * water_eq_density
                        self.checked_values[idx, idy, idz] = 1
                    elif self.checked_values[idx, idy, idz] and ray.acc_depth + D/2 < self.depth[idx, idy, idz]:
                        self.depth[idx, idy, idz] = (
                            ray.acc_depth + D/2) * water_eq_density

                    ray.add_depth(D * water_eq_density)

                if (plane == xPlane).all() and idx <= self.depth.shape[0]:
                    if ray.vec[0] > 0:
                        idx += 1
                    elif ray.vec[0] < 0:
                        idx -= 1
                if (plane == yPlane).all() and idy <= self.depth.shape[1]:
                    if ray.vec[1] > 0:
                        idy += 1
                    elif ray.vec[1] < 0:
                        idy -= 1
                if (plane == zPlane).all() and idz <= self.depth.shape[2]:
                    if ray.vec[2] > 0:
                        idz += 1
                    elif ray.vec[2] < 0:
                        idz -= 1

                a0 = a1

        self.traced_indices = np.array(traced_indices)

#     def trace_rays_old(self):
        
#         traced_indices = []
        
#         xPlane = np.array([1, 0, 0])
#         yPlane = np.array([0, 1, 0])
#         zPlane = np.array([0, 0, 1])

#         for ray in self.rays:

#             S = ray.S
#             T = ray.T

#             intersects = []

#             rayDirection = ray.vec_norm
#             rayPoint = ray.S

#             planePoints_X = np.linspace(
#                 self.medium.mesh_origin[0], self.medium.mesh_apex[0], num=self.depth.shape[0] + 1)
#             planePoints_Y = np.linspace(
#                 self.medium.mesh_origin[1], self.medium.mesh_apex[1], num=self.depth.shape[1] + 1)
#             planePoints_Z = np.linspace(
#                 self.medium.mesh_origin[2], self.medium.mesh_apex[2], num=self.depth.shape[2] + 1)

#             for planeNormal in [xPlane, yPlane, zPlane]:

#                 if (planeNormal == xPlane).all():
#                     planePoints = planePoints_X
#                 elif (planeNormal == yPlane).all():
#                     planePoints = planePoints_Y
#                 elif (planeNormal == zPlane).all():
#                     planePoints = planePoints_Z

#                 for planePoint in planePoints:

#                     # epsilon=1e-6

#                     ndotu = planeNormal.dot(rayDirection)

#                     # print("ndotu", ndotu)

#                     if ndotu != 0:
#                         w = rayPoint - planePoint
#                         si = -planeNormal.dot(w) / ndotu
#                         intersect = w + si * rayDirection + planePoint

#                         if self.medium.mesh_origin[0] <= intersect[0] <= self.medium.mesh_apex[0] and self.medium.mesh_origin[1] <= intersect[1] <= self.medium.mesh_apex[1] and self.medium.mesh_origin[2] <= intersect[2] <= self.medium.mesh_apex[2]:
#                             distance_to_source = np.sqrt(
#                                 sum(list(map(lambda x, y: (x - y)**2, S, intersect))))
#                             intersects.append(
#                                 Point(intersect, planeNormal, distance_to_source))
#                             self.intersects_xes.append(intersect[0])
#                             self.intersects_yes.append(intersect[1])
#                             self.intersects_zes.append(intersect[2])

            

#             if len(intersects) == 0:
#                 continue
            
#             self.intersects_found = True
            
#             intersects = sorted(intersects, key=lambda x: x.distance_to_source)
            
#             self.my_intersects = intersects

#             a0 = intersects[0].P


#             idx = np.floor(
#                 (a0[0] - self.medium.mesh_origin[0])/self.medium.resolution[0]).astype(int)
#             idy = np.floor(
#                 (a0[1] - self.medium.mesh_origin[1])/self.medium.resolution[1]).astype(int)
#             idz = np.floor(
#                 (a0[2] - self.medium.mesh_origin[2])/self.medium.resolution[2]).astype(int)


#             # Adresses the case of first intersect is at the mesh faces associated with the mesh apex point
#             if np.fmod(a0[0], self.resolution[0]) == 0 and ray.vec_norm[0] < 0:
#                 idx -= 1
#             if np.fmod(a0[1], self.resolution[1]) == 0 and ray.vec_norm[1] < 0:
#                 idy -= 1
#             if np.fmod(a0[2], self.resolution[2]) == 0 and ray.vec_norm[2] < 0:
#                 idz -= 1
            
            
            
#             for intersect in intersects[1:]:

#                 # # Ends tracing for current ray in case the accumulated depth exceeds the max depth
#                 # if ray.acc_depth > 36:
#                 #     break

#                 a1 = intersect.P
#                 plane = intersect.plane

#                 if 0 <= idx < self.depth.shape[0] and 0 <= idy < self.depth.shape[1] and 0 <= idz < self.depth.shape[2]:
                    
                    
#                     # if self.config["Image"] == "3D": # 3D, xslice, yslice, zslice
#                     #     traced_indices.append([idx, idy, idz])        
#                     # elif self.config["Image"] == "xslice":
#                     #     # if idy == self.config["Slice"][1]:
#                     traced_indices.append([idx, idy, idz])
                            
#                     # D = np.sqrt(  (self.resolution[0]*(a1[0]-a0[0]))**2
#                     #               + (self.resolution[1]*(a1[1]-a0[1]))**2
#                     #               + (self.resolution[2]*(a1[2]-a0[2]))**2 )

#                     D = np.sqrt(((a1[0]-a0[0]))**2 +
#                                   ((a1[1]-a0[1]))**2 + ((a1[2]-a0[2]))**2)
                    

#                     water_eq_density = self.medium.medium[idx, idy, idz]

#                     if not self.checked_values[idx, idy, idz]:
#                         self.depth[idx, idy, idz] = (
#                             ray.acc_depth + D/2) * water_eq_density
#                         self.checked_values[idx, idy, idz] = 1
#                     elif self.checked_values[idx, idy, idz] and ray.acc_depth + D/2 < self.depth[idx, idy, idz]:
#                         self.depth[idx, idy, idz] = (
#                             ray.acc_depth + D/2) * water_eq_density

#                     ray.add_depth(D * water_eq_density)

#                 if (plane == xPlane).all() and idx <= self.depth.shape[0]:
#                     if ray.vec[0] > 0:
#                         idx += 1
#                     elif ray.vec[0] < 0:
#                         idx -= 1
#                 if (plane == yPlane).all() and idy <= self.depth.shape[1]:
#                     if ray.vec[1] > 0:
#                         idy += 1
#                     elif ray.vec[1] < 0:
#                         idy -= 1
#                 if (plane == zPlane).all() and idz <= self.depth.shape[2]:
#                     if ray.vec[2] > 0:
#                         idz += 1
#                     elif ray.vec[2] < 0:
#                         idz -= 1

#                 a0 = a1
#         self.traced_indices = np.array(traced_indices)


def _sort_by_distance_to_source(con, planes, S):
    func = lambda a, b: np.sqrt(np.sum((b - a)**2, axis=1))
    
    
    dist = func(con, S)

    intersects = [x for _, x in sorted(zip( dist, con), key=lambda x: x[0])]
    
    
    planes = [x for _, x in sorted(zip( dist, planes), key=lambda x: x[0])]

    return intersects, planes


def _find_intersects_with_planes(ray, planeNormal, planePoints):

    N = len(planePoints)
    
    ndotu = planeNormal.dot(ray.vec_norm)

    if ndotu == 0:
        return  np.empty(shape=[0, 3])
        
    w = ray.S - planePoints
    si = - (planeNormal[0]*w[0:,0] + planeNormal[1]*w[0:,1] + planeNormal[2]*w[0:,2]) /ndotu
    intersect_points = w + si[0:,None]*np.tile(ray.vec_norm, (N,1)) + planePoints

    return intersect_points



def _find_intersect(ray, planeNormal, planePoints):

    intersects = []

    if isinstance(planePoints, list):
        planePoint = planePoints[0]
    else:
        planePoint = planePoints

    ndotu = planeNormal.dot(ray.vec_norm)

    if ndotu != 0:
        w = ray.S - planePoint
        si = -planeNormal.dot(w) / ndotu
        intersect = w + si * ray.vec_norm + planePoint

        # print("planepoint", planePoint)
        # print("planeNormal", planeNormal)
        # print("rayPoint", ray.S)
        # print("rayDirection", ray.vec_norm)
        # print("intersect", intersect)

        distance_to_source = np.sqrt(
            sum(list(map(lambda x, y: (x - y)**2, ray.S, intersect))))
        intersects.append(
            Point(intersect, planeNormal, distance_to_source))
        # self.intersects_xes.append(intersect[0])
        # self.intersects_yes.append(intersect[1])
        # self.intersects_zes.append(intersect[2])

    # print(intersects)

    return intersects


# def _find_intersect_points_planes(ray, medium):
#         intersects = []
        
#         planePoints_x = np.zeros((medium.medium.shape[0] + 1, 3))
#         planePoints_x[0:, 0] = np.linspace(
#             medium.mesh_origin[0], medium.mesh_apex[0], num=medium.medium.shape[0] + 1)
        
#         planePoints_y = np.zeros((medium.medium.shape[1] + 1, 3))
#         planePoints_y[0:, 1] = np.linspace(
#             medium.mesh_origin[1], medium.mesh_apex[1], num=medium.medium.shape[1] + 1)

#         planePoints_z = np.zeros((medium.medium.shape[2] + 1, 3))
#         planePoints_z[0:, 2] = np.linspace(
#             medium.mesh_origin[2], medium.mesh_apex[2], num=medium.medium.shape[2] + 1)
        
#         xPlane = np.array([1, 0, 0])
#         yPlane = np.array([0, 1, 0])
#         zPlane = np.array([0, 0, 1])

#         z_intersects = _find_intersects_with_planes(
#             ray, np.asarray([0, 0, 1]), planePoints_z)
#         x_intersects = _find_intersects_with_planes(
#             ray, np.asarray([1, 0, 0]), planePoints_x)
#         y_intersects = _find_intersects_with_planes(
#             ray, np.asarray([0, 1, 0]), planePoints_y)

#         x_planes = np.tile(xPlane, (len(x_intersects), 1))
#         y_planes = np.tile(yPlane, (len(y_intersects), 1))
#         z_planes = np.tile(zPlane, (len(z_intersects), 1))

#         con = np.concatenate((x_intersects, y_intersects, z_intersects))
#         planes = np.concatenate((x_planes, y_planes, z_planes))

#         indices = np.where((medium.minX <= con[0:, 0]) & (medium.maxX >= con[0:, 0]) & (medium.minY <= con[0:, 1]) & (
#             medium.maxY >= con[0:, 1]) & (medium.minZ <= con[0:, 2]) & (medium.maxZ >= con[0:, 2]))
        
#         con = con[indices]
#         planes = planes[indices]

#         intersects, planes_sorted = _sort_by_distance_to_source(
#             con, planes, ray.S)
        
#         return intersects, planes_sorted



def _fill_depth_mesh(ray, intersect_points, intersect_planes, medium, checked_values, depth):
        
        traced_indices = []
        
        xPlane = np.array([1, 0, 0])
        yPlane = np.array([0, 1, 0])
        zPlane = np.array([0, 0, 1])
        
        a0 = intersect_points[0]


        idx = np.floor(
            (a0[0] - medium.mesh_origin[0])/medium.resolution[0]).astype(int)
        idy = np.floor(
            (a0[1] - medium.mesh_origin[1])/medium.resolution[1]).astype(int)
        idz = np.floor(
            (a0[2] - medium.mesh_origin[2])/medium.resolution[2]).astype(int)
        


        # Adresses the case of first intersect is at the mesh faces associated with the mesh apex point
        if np.fmod(a0[0], medium.resolution[0]) == 0 and ray.vec_norm[0] < 0:
            idx -= 1
        if np.fmod(a0[1], medium.resolution[1]) == 0 and ray.vec_norm[1] < 0:
            idy -= 1
        if np.fmod(a0[2], medium.resolution[2]) == 0 and ray.vec_norm[2] < 0:
            idz -= 1
        

        
        for intersect, plane in zip(intersect_points[1:], intersect_planes[1:]):


            a1 = intersect


            if 0 <= idx < medium.medium.shape[0] and 0 <= idy < medium.medium.shape[1] and 0 <= idz < medium.medium.shape[2]:
                

                
                traced_indices.append([idx, idy, idz])
                

                D = np.sqrt(((a1[0]-a0[0]))**2 +
                              ((a1[1]-a0[1]))**2 + ((a1[2]-a0[2]))**2)
                
                water_eq_density = medium.medium[idx, idy, idz]

                if not checked_values[idx, idy, idz]:
                    depth[idx, idy, idz] = (
                        ray.acc_depth + D/2) * water_eq_density
                    checked_values[idx, idy, idz] = 1
                elif checked_values[idx, idy, idz] and ray.acc_depth + D/2 < depth[idx, idy, idz]:
                    depth[idx, idy, idz] = (
                        ray.acc_depth + D/2) * water_eq_density

                ray.add_depth(D * water_eq_density)

            if (plane == xPlane).all() and idx <= depth.shape[0]:
                if ray.vec[0] > 0:
                    idx += 1
                elif ray.vec[0] < 0:
                    idx -= 1
            if (plane == yPlane).all() and idy <= depth.shape[1]:
                if ray.vec[1] > 0:
                    idy += 1
                elif ray.vec[1] < 0:
                    idy -= 1
            if (plane == zPlane).all() and idz <= depth.shape[2]:
                if ray.vec[2] > 0:
                    idz += 1
                elif ray.vec[2] < 0:
                    idz -= 1

            a0 = a1
