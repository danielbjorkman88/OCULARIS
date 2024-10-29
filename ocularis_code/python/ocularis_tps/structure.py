# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 10:57:10 2021

@author: bjoerk_c
"""

import Optis_TPS

class Structure:
    def __init__(self,tps_config ):
        self.tps_config = tps_config 
        self.OrganType = Optis_TPS.OrganData() #[] #Target, OAR etc
        self.dcm = []
        
        
    def set_OrganData(self):
        self.OrganData = Optis_TPS.OrganData()

