# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 08:54:38 2023

@author: bjoerk_c


This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.


"""


import os
import sys
from pathlib import Path

currentdir = Path(os.getcwd())
parent_dir = currentdir.parent.parent
newdir = currentdir.parent.parent.parent.parent / "pyproton"
if newdir not in sys.path:
    sys.path.insert(1, str(newdir))
if parent_dir not in sys.path:
    sys.path.insert(1, str(parent_dir))



from ocularis_tps.tps import OCULARIS
from plot_utils.eyeplan_comparison_utils import one_dimensional_eyeplan_comparison, compare_data, compare_lateral_profile_at_iso, compare_depth_profile_at_iso, compare_horizontal_plane, doseplane_structures, collimator_comparison, plot_wedge_use, compare_dvh






if __name__ == "__main__":

    basedir = Path(os.environ.get(
        "BASEDIR", "E:\OCULARIS_publish\patient_data"))
    
    cohort = [1,2,3,4,5]
    
    cohort = [5]
    
    
    for patient_number in cohort:
    
    # patient_number = cohort[0]
        patient_string = "P" + str(patient_number)
    
        pat1 = OCULARIS.load_patient(basedir, patient_string)
        tps_config = pat1.patient_config
        ep_model = pat1.patient_model.eyeplan_model
        config2 = pat1.dose_engine_config
        skin_plane = pat1.dose_engine_config["skin_plane_point"][2]
    
        algo = pat1.dose_engine_ref_from_eyeplan_model()
    
        collimator_comparison(algo, ep_model)
        one_dimensional_eyeplan_comparison(algo, ep_model)
        
        doseplane_structures(algo, ep_model)
        
        
        
        
        
