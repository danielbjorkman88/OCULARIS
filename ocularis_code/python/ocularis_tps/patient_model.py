# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:56:49 2021

@author: bjoerk_c
"""



import os
import sys
from pathlib import Path
from configuration import config

from pyproton.volume.volume import Volume
# from pyproton.handler.dicom_volume_handler import DicomVolumeHandler
from pyproton.handler.dicom_eye_structure_handler import DicomEyeStructureHandler



#import optis_tps

currentdir = Path(os.getcwd())
newdir = currentdir.parent #.parent / "pyproton"
if newdir not in sys.path:
    sys.path.insert(1, str(newdir))

import models


class PatientModel:
    def __init__(self, patient_config):
        # self.RegionsOfInterest = []
        self.OrganData = []
        self.patient_config = patient_config
        self._dcm_volume = None
        self._structure_set = None
        self._eyeplan_model = None
        self._grid = None
        
        
    @property
    def grid(self):
        return self._grid

    @property
    def structure_set(self):
        if self._structure_set == None:
            structure_handler = DicomEyeStructureHandler(
                self.patient_config["dicom_structures_path"])
            self._structure_set = structure_handler.structure_set

            volume = Volume(self.dcm_volume.pixel_array,
                            self.dcm_volume.voxel_to_world_matrix)

            self._structure_set.grid = volume.grid
        return self._structure_set

    @property
    def dcm_volume(self):
        if self._dcm_volume == None:
            currdir = os.getcwd()
            os.chdir(self.patient_config["dicom_volume_path"])
            list_dir = os.listdir(self.patient_config["dicom_volume_path"])
            list_dir = list(filter(lambda x: x[-4:] == ".dcm", list_dir))
            self._dcm_volume = DicomVolumeHandler(list_dir)
            os.chdir(currdir)
            print("Dicom model loaded")
        return self._dcm_volume

    @property
    def eyeplan_model(self):
        
        if self._eyeplan_model == None:
            
            if not type(self.patient_config) == dict:
                raise ValueError
            
            self._eyeplan_model = models.EyePlanModel(self.patient_config)
            self._eyeplan_model.config = config
            self._eyeplan_model.grid = self._grid
            print("Eyeplan model loaded")
            
        return self._eyeplan_model
        
        
        
        
        
        
        
        
        
        
        
        
        
