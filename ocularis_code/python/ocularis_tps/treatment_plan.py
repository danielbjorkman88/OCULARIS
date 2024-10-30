# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.
"""


import optis_tps

import os
import sys
from pathlib import Path
import numpy as np

currentdir = Path(os.getcwd())
newdir = currentdir.parent.parent / "pyproton"
if newdir not in sys.path:
    sys.path.insert(1, str(newdir))


from pyproton.metrics.dvh import DVH

class TreatmentPlan:
    def __init__(self, tps_config, patient_model):
        self.tps_config = tps_config
        self._patient_model = patient_model
        #self._total_dose = None
        self._beam_set = None
        # self.BeamSet = optis_tps.BeamSet(self.tps_config)
        # self.Optimizer = optis_tps.Optimizer(tps_config)

        self._dvhs = {}

        self.type = []
        #assert self.type in ["Clinical", "Research"]

    @property
    def beam_set(self):
        if self._beam_set == None:
            self._beam_set = optis_tps.BeamSet(
                self.tps_config, self._patient_model)
        return self._beam_set

    @property
    def dvhs(self):
        if not self._dvhs:
            total_dose = self.total_dose
            for name in self._patient_model.structure_set.name_list:
                self._dvhs[name] = DVH(
                    total_dose, self._patient_model.structure_set[name].binary_mask*self.beam_set.beams[0].dose_engine.binary_mask)
        return self._dvhs

    @property
    def total_dose(self):
        total_dose = np.zeros(self._patient_model.structure_set.grid.size)
        for beam in self.beam_set:
            if beam.dose_engine.dose_calculated:
                total_dose += beam.dose_engine.dose
                print(f"Dose from {beam} added to total dose")
            else:
                print(f"Warning: Dose not added from: {beam}")

        return total_dose

    # def Evaluate_Goals(self):
    #     pass

    def evaluate(self):
        pass
