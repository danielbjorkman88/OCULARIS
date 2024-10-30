# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.
"""

import Optis_TPS

class Structure:
    def __init__(self,tps_config ):
        self.tps_config = tps_config 
        self.OrganType = Optis_TPS.OrganData() #[] #Target, OAR etc
        self.dcm = []
        
        
    def set_OrganData(self):
        self.OrganData = Optis_TPS.OrganData()

