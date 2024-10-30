# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 16:24:34 2021

@author: bjoerk_c

This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import numpy as np
from abc import ABC
import logging
import os, sys
from pathlib import Path


class Medium(ABC):
    def __init__(self, config):
        self.config = config
        self.mesh_origin = np.asarray(
            [-config["Mesh dimensions"][0]/2, -config["Mesh dimensions"][1]/2,  -config["SID"]]) #-config["Mesh dimensions"][2]/2
        self.mesh_apex = np.asarray(
            [config["Mesh dimensions"][0]/2, config["Mesh dimensions"][1]/2, -config["SID"] + config["Mesh dimensions"][2]])
        self.resolution = np.asarray([self.config["Mesh dimensions"][0]/self.config["Nvoxels"][0], self.config["Mesh dimensions"]
                                     [1]/self.config["Nvoxels"][1], self.config["Mesh dimensions"][2]/self.config["Nvoxels"][2]])
        # self.isocentre = np.asarray([0, 0, 0])
        # self.vector = np.asarray([0, 0, 1])
        # self.azimuth = []
        # self.polar = []
        # self.torsion = []
        # self.medium = []
        # self.structures = []

        self.update_grid()

    def update_grid(self):

        self.minZ = self.mesh_origin[2]
        self.maxZ = self.mesh_apex[2]
        self.widthZ = self.resolution[2]

        self.minX = self.mesh_origin[0]
        self.maxX = self.mesh_apex[0]
        self.widthX = self.resolution[0]

        self.minY = self.mesh_origin[1]
        self.maxY = self.mesh_apex[1]
        self.widthY = self.resolution[1]

        Z, X = np.meshgrid(np.arange(self.minZ, self.maxZ + self.widthZ,
                                     self.widthZ), np.arange(self.minX, self.maxX + self.widthX, self.widthX))
        # # Correction for rounding error
        if Z.shape[0] == self.config["Nvoxels"][2]+2:
            Z = Z[0:-1,0:]
        if Z.shape[1] == self.config["Nvoxels"][2]+2:
            Z = Z[0:,0:-1]
        if X.shape[0] == self.config["Nvoxels"][0]+2:
            X = X[0:-1, 0:]
        if X.shape[1] == self.config["Nvoxels"][0]+2:
            X = X[0:, 0:-1]


        self.Z, self.X = Z, X 
        
        
        
        
        _, Y = np.meshgrid(np.arange(self.minZ, self.maxZ + self.widthZ,
                                     self.widthZ), np.arange(self.minY, self.maxY + self.widthY, self.widthY))
        
        if Y.shape[0] == self.config["Nvoxels"][1]+2:
            Y = Y[0:-1,0:]
        
        if Y.shape[1] == self.config["Nvoxels"][1]+2:
            Y = Y[0:,0:-1]
        
        self.Y = Y
        
        self.trans_X, self.trans_Y = np.meshgrid(np.arange(
            self.minX, self.maxX + self.widthX, self.widthX), np.arange(self.minY, self.maxY + self.widthY, self.widthY))
        
        # # Correction for rounding error
        if self.trans_X.shape[0] == self.config["Nvoxels"][0]+2:
            self.trans_X = self.trans_X[0:-1,0:]
        if self.trans_X.shape[1] == self.config["Nvoxels"][0]+2:
            self.trans_X = self.trans_X[0:,0:-1]
            
        if self.trans_Y.shape[0] == self.config["Nvoxels"][1]+2:
            self.trans_Y = self.trans_Y[0:-1,0:]
        if self.trans_Y.shape[1] == self.config["Nvoxels"][1]+2:
            self.trans_Y = self.trans_Y[0:,0:-1]            

        self.z_voxel_pos = np.linspace(
            self.minZ + self.widthZ/2, self.maxZ - self.widthZ/2, num=self.config["Nvoxels"][2])
        self.x_voxel_pos = np.linspace(
            self.minX + self.widthX/2, self.maxX - self.widthX/2, num=self.config["Nvoxels"][0])
        self.y_voxel_pos = np.linspace(
            self.minY + self.widthY/2, self.maxY - self.widthY/2, num=self.config["Nvoxels"][1])

    def rotate_azimuth(self):
        pass

    def rotate_polar(self):
        pass

    def rotate_torsion(self):
        pass

    def translate(self, vector):
        self.mesh_origin = self.mesh_origin + vector
        self.mesh_apex = self.mesh_apex + vector
        
        self.update_grid()
        
        self.logger = logging.getLogger("Algo.log")
        if not self.logger.handlers:
            self.logger.setLevel(logging.INFO)
            self.logger.addHandler(logging.StreamHandler())
            
        self.logger.info("Medium updated")
        


class Phantom(Medium):
    def __init__(self, config):
        super().__init__(config)
        self.medium = np.ones(
            [self.config["Nvoxels"][0], self.config["Nvoxels"][1], self.config["Nvoxels"][2]])
        
        if config["anterior_aligned"]:
            n_voxels = config["n_voxels_front_of_skinplane"]
            self.medium[0:, 0:,-n_voxels:] = 0
            
          

         
class Dicom(Medium):
    def __init__(self, config):
        super().__init__(config)
        self.medium = np.zeros(
            [self.config["Nvoxels"][0], self.config["Nvoxels"][1], self.config["Nvoxels"][2]])
        self.medium[self.config["Sclera"]] = 1
        self.medium[0:, 0:, 0:int(self.config["Nvoxels"][2]/1.5)] = 1
        

class GeoModel(Medium):
    pass
