# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

from plot_utils.Plotter import Plotter

from scipy.special import erf
from scipy.interpolate import interp1d
from pathlib import Path
from configuration import config, data
from matplotlib.pyplot import cm
import numpy as np

import matplotlib.pyplot as plt


def format_func(value, tick_number):

    if value == 0:
        return 0
    else:
        return -value


class PlotDoseOnDicom(Plotter):

   def __init__(self, algo, patient):
       self.algo = algo
       self.patient = patient

   def plot(self, x_bin, y_bin):
       
       
       
       image_structures = np.zeros(self.patient.patient_model.structure_set["Sclera"].binary_mask[0:,y_bin, 0:].shape)
       image_structures[self.patient.patient_model.structure_set["Sclera"].binary_mask[0:,y_bin, 0:]] = 1
       image_structures[self.patient.patient_model.structure_set["GTV"].binary_mask[0:,y_bin, 0:]] = 2
       
       image_dose = self.algo.dose[0:, y_bin, 0:]

       dose_vector = self.algo.dose[x_bin, y_bin, 0:]

       fig = plt.figure()

       ax = plt.subplot(111)

       # self.algo.medium.trans_X , self.algo.medium.trans_Y,
       plt.pcolor(self.algo.medium.Z, self.algo.medium.X,  image_structures)

       masked_data = np.ma.masked_where(image_dose < 0.001, image_dose)
       plt.pcolor(self.algo.medium.Z, self.algo.medium.X,  masked_data, alpha=0.6,
                  cmap="jet", vmin=0, vmax=1)  # self.algo.medium.trans_X , self.algo.medium.trans_Y,
       # plt.pcolor(self.algo.medium.Z, self.algo.medium.X,  [], alpha = 1, cmap = "jet", vmin = 0, vmax = 1)
       cbar = plt.colorbar()
       cbar.ax.set_ylabel('Dose', rotation=270, size='xx-large')

       plt.plot(self.algo.central_axis.zes, self.algo.central_axis.xes,
                color="k", label="Central axis")

       plt.ylim(-20, 20)
       plt.xlim(-16, 20)

       plt.ylabel("x [mm]", fontsize=14)
       plt.xlabel("z [mm]", fontsize=14)

       plt.title(
           f'Aperture expanded by {self.algo.config["aperture_expansion"]}', fontsize=14)

       ax2 = ax.twinx()

       plt.plot(self.algo.medium.z_voxel_pos, dose_vector, color="w")

       # plt.xlim(-10, 80)
       # plt.ylabel("Depth [mm]", fontsize = 14)
       plt.xlabel("z [mm]", fontsize=14)

       plt.xlim(self.algo.medium.minZ, self.algo.medium.maxZ)

       plt.suptitle(f"{self.algo.nozzle.title},  y slice {y_bin}")

       xlength = 12
       fig.set_size_inches(xlength, xlength/1.61803398875)
       plt.show()
       try:
            plt.savefig(self.plot_path / "plot_dose_on_dicom.pdf",
                        bbox_inches='tight', pad_inches=0.1)
       except:
            pass


