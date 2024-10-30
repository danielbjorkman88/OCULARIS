# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import bisect
import matplotlib.tri as tri
from shapely.geometry import LineString
from shapely.geometry import Polygon
from shapely.ops import linemerge


class ConcaveHull:
    """
    https://gist.github.com/AndreLester/589ea1eddd3a28d00f3d7e47bd9f28fb
    """

    def __init__(self):
        self.triangles = {}
        self.crs = {}
        self.boundary = []
        
        self.triangulated = False

    def loadpoints(self, points):
        self.points = points

    def edge(self, key, triangle):
        '''Calculate the length of the triangle's outside edge
        and returns the [length, key]'''
        pos = triangle[1].index(-1)
        if pos == 0:
            x1, y1 = self.points[triangle[0][0]]
            x2, y2 = self.points[triangle[0][1]]
        elif pos == 1:
            x1, y1 = self.points[triangle[0][1]]
            x2, y2 = self.points[triangle[0][2]]
        elif pos == 2:
            x1, y1 = self.points[triangle[0][0]]
            x2, y2 = self.points[triangle[0][2]]
        length = ((x1-x2)**2+(y1-y2)**2)**0.5
        rec = [length, key]
        return rec

    def triangulate(self):

        if len(self.points) < 2:
            raise Exception(
                'CountError: You need at least 3 points to Triangulate')

        temp = list(zip(*self.points))
        x, y = list(temp[0]), list(temp[1])
        del(temp)

        # Fix for case if points are on straight line, not allowing the algorithm to triangulate
        try:
            triang = tri.Triangulation(x, y)
        except:
            x[0] = x[0]*(1 + 0.001)
            triang = tri.Triangulation(x, y)

        self.triangles = {}

        for i, triangle in enumerate(triang.triangles):
            self.triangles[i] = [list(triangle), list(triang.neighbors[i])]
        
        self.triangulated = True
        
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
                        if pos == 0:
                            x1, y1 = self.points[triangle[0][0]]
                            x2, y2 = self.points[triangle[0][1]]
                        elif pos == 1:
                            x1, y1 = self.points[triangle[0][1]]
                            x2, y2 = self.points[triangle[0][2]]
                        elif pos == 2:
                            x1, y1 = self.points[triangle[0][0]]
                            x2, y2 = self.points[triangle[0][2]]
                        line = LineString([(x1, y1), (x2, y2)])
                        edgelines.append(line)

        bound = linemerge(edgelines)

        self.boundary = Polygon(bound.coords)
