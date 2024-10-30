# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


from scipy.special import erf
from scipy.interpolate import interp1d
from pathlib import Path
from configuration import config, data
from matplotlib.pyplot import cm
import numpy as np
from scipy.ndimage import zoom

import matplotlib.pyplot as plt
import reference_frames
import geometry
from plot_utils.Plotter import Plotter


class PlotPointCloudAndGaze(Plotter):

    def __init__(self, config, patient):
        self.config = config
        self.patient = patient

    def plot(self):

        ref = reference_frames.Reference_Frame(self.config, 'TreatmentRoom')

        fig = plt.figure()

        ax = fig.add_subplot(121, projection='3d')
        ax.view_init(-70, 90)

        ref.add_frame_to_plot(ax, ref, 10)

        for name in ['GTV', 'Lens', 'Sclera', "Lens_L", "Sklera", "clip4", "clip1", "clip2", "clip3"]:
            if name in self.patient.patient_model.structure_set.name_list:
                struct = self.patient.patient_model.structure_set[name]
            else:
                continue
            plt.plot([], [], [], color=struct.color, label=name)
            for i, j, k in struct.contour_coordinates[::2]:

                alpha = 1

                if name in ['Sclera', "Sklera"]:
                    alpha = 0.03
                elif name == 'Vitreous':
                    alpha = 0.01

                if name in ["clip4", "clip1", "clip2", "clip3"]:
                    color = "C1"
                else:
                    color = self.patient.patient_model.structure_set[name].color

                plt.plot(i, j, k, marker="*", color=color, alpha=alpha)

        gaze_vector = self.patient.patient_model.structure_set.gaze_vector

        plt.quiver(
            0, 0, 0, gaze_vector[0], gaze_vector[1], gaze_vector[2], color="k", length=15)

        plt.title(
            "com " + str(self.patient.patient_model.structure_set["GTV"].com))

        ax.set_xlabel('X', fontsize=12)
        # ax.set_xlim(0, 120)
        ax.set_ylabel('Y', fontsize=12)
        # ax.set_ylim(150, 200)
        ax.set_zlabel('Z', fontsize=12)

        plt.legend()

        ax = fig.add_subplot(122, projection='3d')
        ax.view_init(-70, 90)

        ref.add_frame_to_plot(ax, ref, 10)

        gaze_vector = self.patient.patient_model.structure_set.gaze_vector

        gaze_vector_voxel = self.patient.patient_model.structure_set.gaze_vector_voxel

        plt.quiver(0, 0, 0, gaze_vector[0], gaze_vector[1], gaze_vector[2],
                   color="C0", length=5, label=f"Point cloud, {gaze_vector}")

        plt.quiver(0, 0, 0, gaze_vector_voxel[0], gaze_vector_voxel[1], gaze_vector_voxel[2],
                   color="C1", length=5, label=f"Voxel distribution, {gaze_vector_voxel}")

        plt.title("Gaze orientation by model", fontsize=12)

        plt.legend()

        ax.set_xlabel('X', fontsize=12)
        ax.set_xlim(-10, 10)
        ax.set_ylabel('Y', fontsize=12)
        ax.set_ylim(-10, 10)
        ax.set_zlabel('Z', fontsize=12)
        ax.set_zlim(-10, 10)

        xlength = 16
        fig.set_size_inches(xlength, xlength/1.61803398875)
        plt.show()

        try:
            plt.savefig(self.algo.plot_path / "Plot_point_cloud_and_gaze.pdf",
                        bbox_inches='tight', pad_inches=0.1)
        except:
            pass
