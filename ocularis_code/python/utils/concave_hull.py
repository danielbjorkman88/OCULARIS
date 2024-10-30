# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import bisect
from collections import OrderedDict
import math
#import numpy as np
import matplotlib.tri as tri
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.ops import linemerge


# import geometry


class ConcaveHull:
    """
    
    hull_in_plane = ConcaveHull(algo.config)
    hull_in_plane.loadpoints(points[:, :2])
    hull_in_plane.calculatehull(3.5)
    hull_in_plane.configure(np.asarray([0, 0, 0]))
    
    
    https://gist.github.com/AndreLester/589ea1eddd3a28d00f3d7e47bd9f28fb
    """
    
    
    
    def __init__(self, config):
        
        self.config = config
        
        self.triangles = {}
        self.crs = {}
        self.boundary = []
        
        self.boundary_xes = []
        self.boundary_yes = []
        self.boundary_zes = []
        
        self.boundary_in_plane_xes = []
        self.boundary_in_plane_yes = []
        
        self.boundary_points = []
        self.boundary_edges = []
        
    
    def configure(self, p_on_axis):
        
        from reference_frames import Reference_Frame
        
        # geo = Geometry.OptisGeometry(self.config)
        bev = Reference_Frame(self.config, 'Beams Eye View')
        
        
        self.boundary_in_plane_xes = self.boundary.boundary.xy[0]
        self.boundary_in_plane_yes = self.boundary.boundary.xy[1] 
        
        # xes = self.boundary.boundary.xy[0]
        # yes = self.boundary.boundary.xy[1]  
        
        for idx in range(len(self.boundary_in_plane_xes)):
            
            x = self.boundary_in_plane_xes[idx]
            y = self.boundary_in_plane_yes[idx]
            p_global = p_on_axis - bev.xvec_in_ref*x + bev.yvec_in_ref*y
            
            self.boundary_points.append(p_global)
            self.boundary_xes.append(p_global[0])
            self.boundary_yes.append(p_global[1])
            self.boundary_zes.append(p_global[2])
            
            
        for idx in range(len(self.boundary_points)):
            
            
            
            next_idx = idx + 1
            #print(next_idx , len(self.boundary_points))
            if next_idx >= len(self.boundary_points):
                next_idx -= len(self.boundary_points)
                
                
            xes = [self.boundary_points[idx][0], self.boundary_points[next_idx][0]]
            yes = [self.boundary_points[idx][1], self.boundary_points[next_idx][1]]
            zes = [self.boundary_points[idx][2], self.boundary_points[next_idx][2]]
                
            self.boundary_edges.append([xes, yes, zes])
            
    
    
    def loadpoints(self, points):
        #self.points = np.array(points)
        self.points = points
        
        
    def edge(self, key, triangle):
        '''Calculate the length of the triangle's outside edge
        and returns the [length, key]'''
        pos = triangle[1].index(-1)
        if pos==0:
            x1, y1 = self.points[triangle[0][0]]
            x2, y2 = self.points[triangle[0][1]]
        elif pos==1:
            x1, y1 = self.points[triangle[0][1]]
            x2, y2 = self.points[triangle[0][2]]
        elif pos==2:
            x1, y1 = self.points[triangle[0][0]]
            x2, y2 = self.points[triangle[0][2]]
        length = ((x1-x2)**2+(y1-y2)**2)**0.5
        rec = [length, key]
        return rec
        
    
    def triangulate(self):
        
        if len(self.points) < 2:
            raise Exception('CountError: You need at least 3 points to Triangulate')
        
        temp = list(zip(*self.points))
        x, y = list(temp[0]), list(temp[1])
        del(temp)
        
        triang = tri.Triangulation(x, y)
        
        self.triangles = {}
        
        for i, triangle in enumerate(triang.triangles):
            self.triangles[i] = [list(triangle), list(triang.neighbors[i])]
        

    def calculatehull(self, tol=50):
        
        self.tol = tol
        
        if len(self.triangles) == 0:
            self.triangulate()
        
        # All triangles with one boundary longer than the tolerance (self.tol)
        # is added to a sorted deletion list.
        # The list is kept sorted from according to the boundary edge's length
        # using bisect        
        deletion = []    
        self.boundary_vertices = set()
        for i, triangle in self.triangles.items():
            if -1 in triangle[1]:
                for pos, neigh in enumerate(triangle[1]):
                    if neigh == -1:
                        if pos == 0:
                            self.boundary_vertices.add(triangle[0][0])
                            self.boundary_vertices.add(triangle[0][1])
                        elif pos == 1:
                            self.boundary_vertices.add(triangle[0][1])
                            self.boundary_vertices.add(triangle[0][2])
                        elif pos == 2:
                            self.boundary_vertices.add(triangle[0][0])
                            self.boundary_vertices.add(triangle[0][2])
            if -1 in triangle[1] and triangle[1].count(-1) == 1:
                rec = self.edge(i, triangle)
                if rec[0] > self.tol and triangle[1].count(-1) == 1:
                    bisect.insort(deletion, rec)
                    
        while len(deletion) != 0:
            # The triangles with the longest boundary edges will be 
            # deleted first
            item = deletion.pop()
            ref = item[1]
            flag = 0
            
            # Triangle will not be deleted if it already has two boundary edges            
            if self.triangles[ref][1].count(-1) > 1:
                continue
                
            # Triangle will not be deleted if the inside node which is not
            # on this triangle's boundary is already on the boundary of 
            # another triangle
            adjust = {0: 2, 1: 0, 2: 1}            
            for i, neigh in enumerate(self.triangles[ref][1]):
                j = adjust[i]
                if neigh == -1 and self.triangles[ref][0][j] in self.boundary_vertices:
                    flag = 1
                    break
            if flag == 1:
                continue
           
            for i, neigh in enumerate(self.triangles[ref][1]):
                if neigh == -1:
                    continue
                pos = self.triangles[neigh][1].index(ref)
                self.triangles[neigh][1][pos] = -1
                rec = self.edge(neigh, self.triangles[neigh])
                if rec[0] > self.tol and self.triangles[rec[1]][1].count(-1) == 1:
                    bisect.insort(deletion, rec)
                    
            for pt in self.triangles[ref][0]:
                self.boundary_vertices.add(pt)
                                        
            del self.triangles[ref]
            
        self.polygon()
            
                    

    def polygon(self):
        
        edgelines = []
        for i, triangle in self.triangles.items():
            if -1 in triangle[1]:
                for pos, value in enumerate(triangle[1]):
                    if value == -1:
                        if pos==0:
                            x1, y1 = self.points[triangle[0][0]]
                            x2, y2 = self.points[triangle[0][1]]
                        elif pos==1:
                            x1, y1 = self.points[triangle[0][1]]
                            x2, y2 = self.points[triangle[0][2]]
                        elif pos==2:
                            x1, y1 = self.points[triangle[0][0]]
                            x2, y2 = self.points[triangle[0][2]]
                        line = LineString([(x1, y1), (x2, y2)])
                        edgelines.append(line)

        bound = linemerge(edgelines)
    
        self.boundary = Polygon(bound.coords)

        
        
        
        
        

        
        