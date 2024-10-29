# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:06:45 2021

@author: bjoerk_c
"""

from patient_db.parser_db import parser_database

import ocularis_tps

import os
import copy
import logging
import numpy as np
import pandas as pd



class OCULARIS:
    def __init__(self, tps_config):
        self.tps_config = tps_config
        # self.patient = optis_tps.Patient(self.tps_config)
        self.patients = []
        self.patient_identifiers = []
        self.terminal = None
        self._dataframe = None
        
        
        self._ml_filename = "ml_data_sheet.csv"

        self.logger = logging.getLogger("TPS.log")
        if not self.logger.handlers:
            self.logger.setLevel(logging.INFO)
            self.logger.addHandler(logging.StreamHandler())
        
    
    def load_patients(self, patient_data_dir, fraction = 1):
        patient_folders = list(filter(lambda x: (x[0] == "P"), os.listdir(patient_data_dir))) 
        for folder in patient_folders[0:int(fraction*len(patient_folders))]:
            
            patient_config2 = copy.deepcopy(self.tps_config)
            
    
            
            patient_config2["patient_path"] = patient_data_dir + "\\" + folder
            
            patient_config2["patient_model_path"] = patient_config2["patient_path"] + "\\" + "model"
            
            patient_config2["eyeplan ptd filename"] = list(filter(lambda x: (x[-4:] == ".ptd"), os.listdir(patient_config2["patient_path"])))[0]
            
            patient = optis_tps.Patient(patient_config2)
            self.patients.append( patient)
            self.patient_identifiers.append(patient.patient_identifier)
            
    @property
    def dataframe(self, reload = 0):
        
        if 1:
            self._dataframe = _load_frame(self)
            self._dataframe.to_csv(self._ml_filename)
            return self._dataframe
        
        if self._ml_filename not in os.listdir():
            self._dataframe = _load_frame(self)
            self._dataframe.to_csv(self._ml_filename)
            return self._dataframe
        

        self._dataframe = pd.read_csv(self._ml_filename)
        self._dataframe.to_csv(self._ml_filename)
        return self._dataframe
            
    
    def load_patient(basedir, patient_identifier, anterior_aligned = False, wider = False):
        
        if type(patient_identifier) == int:
            patient_identifier = "P" + str(patient_identifier)
        
        return parser_database(basedir, patient_identifier, anterior_aligned, wider)
        #return database(patient_identifier)



def _load_frame(self):
            
    data_features = np.zeros([len(self.patients), len(self.patients[0].patient_model.eyeplan_model.feature_vector)])
    
    data_labels = np.zeros([len(self.patients), len(self.patients[0].patient_model.eyeplan_model.label_vector)])
    
    for i, patient in enumerate( self.patients):
        features = patient.patient_model.eyeplan_model.feature_vector
        labels = patient.patient_model.eyeplan_model.label_vector
        
        data_features[i, 0:] = features
        data_labels[i, 0:] = labels
        

    my_data = pd.DataFrame()
    
    my_data["patient"] = self.patient_identifiers
    
    # Add features
    for i, feature_name in enumerate(self.patients[0].patient_model.eyeplan_model.feature_labels):
        my_data[feature_name] = data_features[0:,i]
    
    
    # Add labels
    for i, label_name in enumerate(self.patients[0].patient_model.eyeplan_model.label_labels):
        my_data[label_name] = data_labels[0:, i]
    
    return my_data
