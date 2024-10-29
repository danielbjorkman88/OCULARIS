# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 09:44:14 2021

@author: bjoerk_c
"""

# import os
# import sys
# import inspect


# currentdir = os.path.dirname(os.path.abspath(
#     inspect.getfile(inspect.currentframe())))
# parentdir = os.path.dirname(currentdir)
# if parentdir not in sys.path:
#     sys.path.insert(1, str(parentdir))

import dose_engine

class Beam:
    def __init__(self, beam_config):
        self.beam_config = beam_config
        self._theta_angle = beam_config["Gantry rot theta"]
        self._phi_angle = beam_config["Gantry rot phi"]
        self.norm_factor = 1
        self.title = f"Beam from theta {self._theta_angle} and phi {self._phi_angle} degrees"
        # self.target_setup = []
        # self.nozzle_setup = []

        self.dose_engine = dose_engine.BroadBeam(
            beam_config["data"], beam_config)

    def __repr__(self):
        return self.title

    def evaluate(self):
        pass

    def configure_dose_engine(self):
        pass

    def calc(self):
        pass
