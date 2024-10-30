# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
import numpy as np
import math
from Ray import Ray


def dotproduct(v1, v2):
  return sum((a*b) for a, b in zip(v1, v2))

def length(v):
  return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
  return math.acos(dotproduct(v1, v2) / (length(v1) * length(v2)))



        

def between_0_and_1(I):
    
    return 0 <= I[0] <= 1 and  0<= I[1] <= 1 and 0<= I[2] <= 1

def between_0_and_minus1(I):
    
    return -1 <= I[0] <= 0 and  -1 <= I[1] <= 0 and -1 <= I[2] <= 0


def between_minus1_and_1(I):
    
    return -1 <= I[0] <= 1 and  -1 <= I[1] <= 1 and -1 <= I[2] <= 1


class Voxel():
    
    def __init__(self):
        
        self.xPlane = np.array([1, 0, 0])
        self.yPlane = np.array([0, 1, 0])
        self.zPlane = np.array([0, 0, 1])
        
        self.origin = np.array([0, 0, 0]) # Origin of voxel
        
        self.apex = np.array([1, 1, 1]) # Apex of voxel
        self.foot = np.array([-1, -1, -1]) # foot of voxel



    def calc_intersection(self, ray, a0):

        # if ray.vec_norm[0] >= 0 and ray.vec_norm[1] >= 0:
        #     apex = np.array([1, 1, 1])
        # elif ray.vec_norm[0] <= 0 and ray.vec_norm[1] < 0:
        #     apex = np.array([-1, -1, 0])
            
        # elif ray.vec_norm[0] > 0 and ray.vec_norm[1] < 0:
        #     apex = np.array([1, -1, 0])
        # elif ray.vec_norm[0] < 0 and ray.vec_norm[1] > 0:
        #     apex = np.array([-1, 1, 0])

        
        assert between_0_and_1(a0), "Entry point out of range"

        exit_candidates = []

        
        for i, val in enumerate(a0):
            if math.isclose(val, 0, abs_tol=1e-8) and val != 0:
                a0[i] = 0
                # print(val, i)
        
        
        # compensate = np.asarray([0,0,0])
        # if a0[0] == 0 and ray.vec_norm[0] < 0:
        #     compensate[0] = 1 
        
        
        # a0 += compensate
        
        print(a0)

        # P0: end point 1 of the segment P0P1
        # P1:  end point 2 of the segment P0P1
        P0 = ray.S
        P1 = ray.T
        
        #vec = ray.vec_norm
        
        # print(a0)
        # print(P0)
        # print(P1)

        # entrylowZ, entryhighZ = False, False
        # entrylowX, entryhighX = False, False
        # entrylowY, entryhighY = False, False        
                

        # if a0[2] == 0:
        #     entrylowZ = True
        # if a0[2] == 1:
        #     entryhighZ = True
        # if a0[0] == 0:
        #     entrylowX = True            
        # if a0[0] == 1:
        #     entryhighX = True                  
        # if a0[1] == 0:
        #     entrylowY = True            
        # if a0[1] == 1:
        #     entryhighY = True       


        #for side in range(1,7):
            
            # planeNormal: normal vector of the Plane 
            # V0: any point that belongs to the Plane             
            # if side == 1:
            #     planeNormal = self.xPlane
            #     V0 = self.origin
            # elif side == 2:
            #     planeNormal = self.yPlane
            #     V0 = self.origin
            # elif side == 3:
            #     planeNormal = self.zPlane
            #     V0 = self.origin         
            # elif side == 4:
            #     planeNormal = self.xPlane
            #     V0 = self.apex
            # elif side == 5:
            #     planeNormal = self.yPlane
            #     V0 = self.apex
            # elif side == 6:
            #     planeNormal = self.zPlane
            #     V0 = self.apex
            

            # w = P0 - V0;
            # u = P1-P0;
            # N = -np.dot(n,w);
            # D = np.dot(n,u)
            # if D != 0:
            #     sI = N / D
            #     I = P0+ sI*u
                
            #     if between_0_and_1(I):
            #         exit_candidates.append(I)
            #         print(I, side)
                    
                
        rayDirection = ray.vec_norm
        rayPoint = a0
        #planePoint = V0
        
        
        for planeNormal in [self.xPlane, self.yPlane, self.zPlane]:
            for planePoint in [self.origin, self.apex, self.foot]:
            
            
            # epsilon=1e-6
            
                ndotu = planeNormal.dot(rayDirection) 
                
                # if abs(ndotu) < epsilon:
                #     print ("no intersection or line is within plane")
                
                if ndotu != 0:
                    w = rayPoint - planePoint
                    si = -planeNormal.dot(w) / ndotu
                    Psi = w + si * rayDirection + planePoint
                    
                    
                    for i,val in enumerate(Psi):
                        if math.isclose(val, 0, abs_tol=1e-8):
                            Psi[i] = 0
                        if math.isclose(val, 1, abs_tol=1e-8):
                            Psi[i] = 1             
                    
                    if between_minus1_and_1(Psi):
                        exit_candidates.append(Psi)
                        print(Psi)
                    
                    # if between_0_and_1(Psi) and (planePoint == V0).any():
                    #     exit_candidates.append(Psi)
                    #     print(Psi, side)
                        
                    # if between_0_and_minus1(Psi) and (planePoint == self.foot).any():
                    #    exit_candidates.append(Psi)
                    #    print(Psi, side)                       
                        
                                

        print(exit_candidates)

        # Removes duplicate entries
        dict_tuple = {tuple(item): index for index, item in enumerate(exit_candidates)}
        exit_candidates = [list(itm) for itm in dict_tuple.keys()]
        
        if a0 in exit_candidates:
            exit_candidates.remove(a0)
            
            
        exit_candidates.sort()
            
            
        print("Before")
            
        print(exit_candidates)        
        #print(len(exit_candidates))
        
        
        
        #im_point = a0 + 1000000*ray.vec_norm
        
        
        remove_candidates = []
        
        for candidate in exit_candidates:
            
            # print("Starting", candidate, " Left: " , exit_candidates)
            
            remove = False
            
            if ray.vec_norm[0] > 0 and candidate[0] < a0[0]:
                remove = True
                #print("x",0)
            elif ray.vec_norm[1] > 0 and candidate[1] < a0[1]:
                remove = True
                #print("y",1)
            elif ray.vec_norm[2] > 0 and candidate[2] < a0[2]:
                remove = True            
                #print("z",2)
            
            
            
            # if ray.vec_norm[0] > 0 and candidate[0] < a0[0]:
            #     remove = True
            # elif ray.vec_norm[1] > 0 and candidate[1] < a0[1]:
            #     remove = True
            # elif ray.vec_norm[2] > 0 and candidate[2] < a0[2]:
            #     remove = True
                

                
                
            if ray.vec_norm[0] < 0 and candidate[0] > a0[0]:
                remove = True  
                #exit_candidates.remove(candidate)
            elif ray.vec_norm[1] < 0 and candidate[1] > a0[1]:
                remove = True  
                #exit_candidates.remove(candidate)                
            elif ray.vec_norm[2] < 0 and candidate[2] > a0[2]:
                remove = True  
                #exit_candidates.remove(candidate)            
                
                                
            if remove:
                # print("removing ", candidate)
                # exit_candidates.remove(candidate)
                remove_candidates.append(candidate)
                
        #     if candidate[2] == 0 and entrylowZ:
        #         exit_candidates.remove(candidate)
        #     elif candidate[2] == 1 and entryhighZ:
        #         exit_candidates.remove(candidate)
        #     elif candidate[0] == 0 and entrylowX:
        #         exit_candidates.remove(candidate)
        #     elif candidate[0] == 1 and entryhighX:
        #         exit_candidates.remove(candidate)
        #     elif candidate[1] == 0 and entrylowY:
        #         exit_candidates.remove(candidate)
        #     elif candidate[1] == 1 and entryhighY:
        #         exit_candidates.remove(candidate)
        
        
        
        
        print("Removes: ",remove_candidates)
        
        for remove in remove_candidates:
            if remove in exit_candidates:
                exit_candidates.remove(remove)
        
        
        # # pick shortest in case multiple intersections
        # if len(exit_candidates) > 1:
            
        #     if a0[0]
        
        
        
        print("After")                
        print(exit_candidates)        
        #print(len(exit_candidates))
        
        assert len(exit_candidates) != 0, "No exit candidates"
        assert len(exit_candidates) == 1, "More than 1 exit candidate position"
        
        
        return exit_candidates[0] # - compensate


    def calc_intersection2(self, ray, a0):
        
        planes = []
        


#Test 


if 1:
    voxel = Voxel()
    
    
    S = [0.5,10,0.5]
    T = [0.5,-3,0.5]
    ray = Ray(S,  T)
    a0 = [0.5 , 0 , 0.5]
    
    assert voxel.calc_intersection(ray, a0) == [0.5, -1. , 0.5]
    
    
    S = [0.5,0.5,-2]
    T = [0.5,0.5,-3]
    ray = Ray(S,  T)
    a0 = [0.5 , 0.5 , 0]
    
    assert voxel.calc_intersection(ray, a0) == [0.5, 0.5 , -1]
    
    
    
    a0 = [0.5 , 0 , 0.5]
    S = a0
    T = [0.5,2,2]
    ray = Ray(S,  T)
    
    assert voxel.calc_intersection(ray, a0) == [0.5, 0.6666666666666667, 1. ]
    
    
    
    a0 = [0 , 0 , 0]
    S = a0
    T = [2,2,2]
    ray = Ray(S,  T)
    
    assert voxel.calc_intersection(ray, a0) == [1, 1, 1. ]
    
    
    
    S = [0,-20,0]
    T = [20,70,0]
    a0 = [0.5, 0.5, 0.0]
    ray = Ray(S,  T)
    assert voxel.calc_intersection(ray, a0) == [0.6111111111111112, 1.0, 0.0]
    
    
    
    S = [7,2,5]
    
    T = [0,7,5]  
    a0 = [1.0, 0.0, 0.0]
    ray = Ray(S,  T)    
    voxel.calc_intersection(ray, a0)
    
    
    S = [70.5,20,50]
    T = [0,70,50]    
    a0 = [0.11999999999998323, 0.0, 0.0]
    ray = Ray(S,  T)    
    assert voxel.calc_intersection(ray, a0) == [0.0, 0.08510638297871151, 0.0]
    
    
    
    S = [70,70,50]
    T = [20,20,50]    
    a0 = [0, 0, 0]
    ray = Ray(S,  T)      
    voxel.calc_intersection(ray, a0)
    
    
    def func(a0, vec, planeNormal, planePoint):
        rayDirection = vec
        rayPoint = a0
        
        epsilon=1e-6
        
        ndotu = planeNormal.dot(rayDirection) 
        
        if abs(ndotu) < epsilon:
            print ("no intersection or line is within plane")
        
        w = rayPoint - planePoint
        si = -planeNormal.dot(w) / ndotu
        Psi = w + si * rayDirection + planePoint
        
        return np.asarray(Psi)
        #print ("intersection at", Psi)
    
    
    
    
    a0 = np.asarray([0.5, 0.5 , 0])
    vec = np.asarray([1,1,0])
    planeNormal = voxel.xPlane
    planePoint = voxel.apex
    
    assert (func(a0, vec, planeNormal, planePoint) == [1., 1., 0.]).any()
    
    
    
    S = [0.5,0.5,-2]
    T = [0.5,0.5,-3]
    ray = Ray(S,  T)
    a0 = [0.5 , 0.5 , 0]
    vec = ray.vec_norm
    planeNormal = voxel.zPlane
    planePoint = voxel.apex
    
    assert (func(a0, vec, planeNormal, planePoint) == [0.5, 0.5, 1. ]).any()
    
    
    
    S = [0,-20,0]
    
    T = [20,70,0]
    a0 = [0.5, 0.5, 0.0]
    ray = Ray(S,  T)
    vec = ray.vec_norm
    planeNormal = voxel.yPlane
    planePoint = voxel.apex
    
    assert (func(a0, vec, planeNormal, planePoint) == [0.61111111, 1.        , 0.        ]).any()
    
    
    
    
    
    S = [0,20,50]
    T = [70,70,50]
    a0 = [0.5, 0.5, 0.0]
    ray = Ray(S,  T)
    vec = ray.vec_norm
    planeNormal = voxel.xPlane
    planePoint = voxel.apex
    func(a0, vec, planeNormal, planePoint)
    
    
    
    
    S = [0.5,0.5,-2]
    T = [0.5,0.5,-3]
    ray = Ray(S,  T)
    a0 = [0.5 , 0.5 , 0]
    vec = ray.vec
    planeNormal = voxel.zPlane
    planePoint = voxel.apex
    func(a0, vec, planeNormal, planePoint)
    
    assert voxel.calc_intersection(ray, a0) == [0.5, 0.5 , 1]
    
    
    












# vec = np.asarray([1,1,0])
# vec_norm = vec/length(vec)



# fig = plt.figure()
# ax = fig.gca(projection='3d')
# ax.set_aspect("auto")

# # draw cube
# r = [0, 1]
# for s, e in combinations(np.array(list(product(r, r, r))), 2):
#     if np.sum(np.abs(s-e)) == r[1]-r[0]:
#         ax.plot3D(*zip(s, e), color="b")

# # # draw sphere
# # u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
# # x = np.cos(u)*np.sin(v)
# # y = np.sin(u)*np.sin(v)
# # z = np.cos(v)
# # ax.plot_wireframe(x, y, z, color="r")

# # draw a point
# ax.scatter(a0[0], a0[1], a0[2], color="g", s=100)




# # draw a vector
# a = Arrow3D([a0[0], 1], [a0[1], 1], [a0[2], 1], mutation_scale=20,
#             lw=1, arrowstyle="-|>", color="k")
# ax.add_artist(a)


# plt.xlabel("x")
# plt.ylabel("y")
# ax.set_zlabel("z")

# plt.show()








# #Based on http://geomalgorithms.com/a05-_intersect-1.html
# #from __future__ import print_function
# import numpy as np

# import matplotlib.pyplot as plt

# epsilon=1e-6

# #Define plane
# planeNormal = np.array([0, 0, 1])
# planePoint = np.array([0, 0, 5]) #Any point on the plane

# #Define ray
# rayDirection = np.array([0, -1, -1])
# rayPoint = np.array([0, 0, 10]) #Any point along the ray

# ndotu = planeNormal.dot(rayDirection) 

# if abs(ndotu) < epsilon:
#     print ("no intersection or line is within plane")

# w = rayPoint - planePoint
# si = -planeNormal.dot(w) / ndotu
# Psi = w + si * rayDirection + planePoint

# print ("intersection at", Psi)

# # create the figure
# fig = plt.figure()

# # add axes
# ax = fig.add_subplot(111,projection='3d')

# xx, yy = np.meshgrid(range(10), range(10))
# z = (9 - xx - yy) / 2 

# # plot the plane
# ax.plot_surface(xx, yy, z, alpha=0.5)

# plt.show()






# #
# # n: normal vector of the Plane 
# # V0: any point that belongs to the Plane 
# # P0: end point 1 of the segment P0P1
# # P1:  end point 2 of the segment P0P1
# n = np.array([1., 1., 1.])
# V0 = np.array([1., 1., -5.])
# P0 = np.array([-5., 1., -1.])
# P1 = np.array([1., 2., 3.])

# w = P0 - V0;
# u = P1-P0;
# N = -np.dot(n,w);
# D = np.dot(n,u)
# sI = N / D
# I = P0+ sI*u
# print(I)

# T = [-3.90909091,  1.18181818, -0.27272727]
# #assert I == np.asarray([-3.90909091,  1.18181818, -0.27272727])

# # for i in range(3):
# #     print(I[i] / T[i])
    
    

















