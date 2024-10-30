# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:39:47 2021

@author: bjoerk_c
"""

# import copy
import os, sys
from pathlib import Path
import copy


currentdir = Path(os.getcwd())
parent_dir = currentdir.parent
newdir = currentdir.parent.parent / "pyproton"
if newdir not in sys.path:
    sys.path.insert(1, str(newdir))
if parent_dir not in sys.path:
    sys.path.insert(1, str(parent_dir))

from configuration import data
#from optis_tps.treatment_plan import TreatmentPlan
from generators.dose_engine_generators import generate_reference
from pyproton.volume.grid import Grid

import ocularis_tps
import pathlib
import numpy as np


class Patient:
    def __init__(self, patient_config):
        self.attributes = []
        self.plans = []
        self.patient_config = patient_config
        # self.BeamSet = optis_tps.BeamSet(self.patient_config )
        # self.Structures = [optis_tps.Structure(self.patient_config )]
        # self.TreatmentPlan = optis_tps.TreatmentPlan(self.patient_config)
        self._patient_model = None
        self._dose_engine_config = None
        self._grid = None
        self.n_voxels_front_of_skinplane = 7

        if type(patient_config['patient_path']) == pathlib.WindowsPath:
            self._patient_identifier = str(
                self.patient_config['patient_path']).split("\\")[-1]
        else:
            self._patient_identifier = self.patient_config['patient_path'].split(
                "\\")[-1]

    def dose_engine_ref_from_eyeplan_model(self, simple=False):

        if self._dose_engine_config == None:
            raise ValueError

        assert type(simple) == bool

        algo = generate_reference(
            data, self._dose_engine_config, self.patient_model.eyeplan_model, simple)

      

        assert self.grid.origin[0] == algo.medium.minX
        assert self.grid.origin[1] == algo.medium.minY
        assert self.grid.origin[2] == algo.medium.minZ

        assert (self.grid.spacing == algo.medium.resolution).all()
        assert (self.grid.size == np.asarray(algo.dose.shape)).all()

        # assert - 0.0001 < (algo.medium.maxZ -
        #                    algo.config["skin_plane"])/algo.config["skin_plane"] < 0.0001

        return algo
    
    

    def dose_engine_from_config(self, config, simple=False):

        if self._dose_engine_config == None:
            raise ValueError

        assert type(simple) == bool

        algo = generate_reference(
            data, config, self.patient_model.eyeplan_model, simple)

        assert self.grid.origin[0] == algo.medium.minX
        assert self.grid.origin[1] == algo.medium.minY
        assert self.grid.origin[2] == algo.medium.minZ

        # assert (self.grid.spacing == algo.medium.resolution).all()
        # assert (self.grid.size == np.asarray(algo.dose.shape)).all()

        # assert - 0.0001 < (algo.medium.maxZ -
        #                    algo.config["skin_plane"])/algo.config["skin_plane"] < 0.0001

        return algo

    @property
    def generate_config_copy(self):
        exclude_keys = ['surfaces']
        new_config = {k: self.dose_engine_config[k] for k in set(list(self.dose_engine_config.keys())) - set(exclude_keys)}
        config2 = copy.deepcopy(new_config)
        config2['surfaces'] = self.dose_engine_config['surfaces']
        return config2

    @property
    def grid(self):

        # if self._grid == None:
        #     self.define_grid()

        return self._grid

    @grid.setter
    def grid(self, grid):
        self._grid = grid

    def define_grid(self):

        if self._dose_engine_config == None:
            raise ValueError

        mesh_origin = np.asarray(
            [-self._dose_engine_config["Mesh dimensions"][0]/2, -self._dose_engine_config["Mesh dimensions"][1]/2,  -self._dose_engine_config["SID"]])  # -self._dose_engine_config["Mesh dimensions"][2]/2
        mesh_apex = np.asarray(
            [self._dose_engine_config["Mesh dimensions"][0]/2, self._dose_engine_config["Mesh dimensions"][1]/2, -self._dose_engine_config["SID"] + self._dose_engine_config["Mesh dimensions"][2]])
        resolution = np.asarray([self._dose_engine_config["Mesh dimensions"][0]/self._dose_engine_config["Nvoxels"][0], self._dose_engine_config["Mesh dimensions"]
                                 [1]/self._dose_engine_config["Nvoxels"][1], self._dose_engine_config["Mesh dimensions"][2]/self._dose_engine_config["Nvoxels"][2]])

        ep_model = self.patient_model.eyeplan_model

        # translation_x = 0
        # translation_y = 0

        
        # if self._dose_engine_config["anterior_aligned"]:
        # # Aligns mesh with anterior segment of eyeglobe + 0.5 mm
        # z_alignment = np.max(ep_model.structure_set_clips_registered["eyeglobe"].contour_coordinates[0:,2]) + 0.5
        
        n_voxels_front_of_skinplane = self.n_voxels_front_of_skinplane #int(np.ceil(3/ep_model.doseplane_h.resolution[2]))
        z_alignment = self._dose_engine_config['skin_plane_point'][2] + n_voxels_front_of_skinplane*ep_model.doseplane_h.resolution[2]
        
        
        print(f"Optis_tps.patient: Mesh aligned to {n_voxels_front_of_skinplane*ep_model.doseplane_h.resolution[2]} mm in front of skin plane")
        # else:
        #     # Aligns mesh with skin plane
        #     z_alignment = self._dose_engine_config["skin_plane_point"][2]
        #     print("Optis_tps: Mesh aligned to skin plane")
            
            
        translation_vector = np.zeros(3)
        
        translation_vector[2] = z_alignment - mesh_apex[2]
        translation_vector[1] = ep_model.doseplane_h.resolution[1]/2
        translation_vector[0] = ep_model.doseplane_h.resolution[0]/2
        
        
        
        # translation_z = self.dose_engine_config['skin_plane_point'][2] - mesh_apex[2]
        # translation = np.asarray([translation_x, translation_y, translation_z])

        # print(mesh_origin[2],  translation_z)
        mesh_origin += translation_vector

        grid = Grid(mesh_origin, resolution,
                    self._dose_engine_config["Nvoxels"], "mm")

        self._grid = grid

        self.patient_model.eyeplan_model.grid = grid

        print("Grid defined for patient and Eyeplan model")

    @property
    def dose_engine_config(self):
        """ Reference config from database"""
        return self._dose_engine_config

    @dose_engine_config.setter
    def dose_engine_config(self, config):

        # mesh_origin = np.asarray(
        #     [-config["Mesh dimensions"][0]/2, -config["Mesh dimensions"][1]/2,  -config["SID"]]) #-config["Mesh dimensions"][2]/2
        # # mesh_apex = np.asarray(
        # #     [config["Mesh dimensions"][0]/2, config["Mesh dimensions"][1]/2, -config["SID"] + config["Mesh dimensions"][2]])
        # resolution = np.asarray([config["Mesh dimensions"][0]/config["Nvoxels"][0], config["Mesh dimensions"]
        #                              [1]/config["Nvoxels"][1], config["Mesh dimensions"][2]/config["Nvoxels"][2]])

        # self.grid = Grid(mesh_origin, resolution, config["Nvoxels"], "mm" )

        self._dose_engine_config = config

        self.define_grid()

    @property
    def patient_identifier(self):
        return self._patient_identifier

    def __repr__(self):
        return "Patient instance " + self.patient_identifier

    @property
    def patient_model(self):

        if self._patient_model == None:
           self._patient_model = ocularis_tps.PatientModel(self.patient_config)
        return self._patient_model

    def define_treatment_plan(self):
        pass
        # patient_config = copy.deepcopy(self.patient_config)
        # patient_config["config"]["aperture_expansion"] = 2.5

        # self.plans.append(TreatmentPlan(patient_config, self.patient_model))

        # patient_config = copy.deepcopy(self.patient_config)
        # patient_config["config"]["aperture_expansion"] = 1

        # self.plans.append(TreatmentPlan(patient_config, self.patient_model))

        # patient_config = copy.deepcopy(self.patient_config)
        # patient_config["config"]["aperture_expansion"] = 3

        # self.plans.append(TreatmentPlan(patient_config, self.patient_model))

        # patient_config = copy.deepcopy(self.patient_config)
        # patient_config["config"]["aperture_expansion"] = 5

        # self.plans.append(TreatmentPlan(patient_config, self.patient_model))

        # patient_config = copy.deepcopy(self.patient_config)
        # patient_config["config"]["aperture_expansion"] = 8

        # self.plans.append(TreatmentPlan(patient_config, self.patient_model))

    def evaluate_treatment_goals(self):
        pass
