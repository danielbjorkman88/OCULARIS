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
import reference_frames
from matplotlib.colors import LogNorm

import matplotlib.pyplot as plt

from pyproton.metrics.dvh import DVH

from utils.color_utils import my_colors


class CompareTransformedEye(Plotter):

   def __init__(self, algo, algo2, eyeplan_model, structure_set, slice_idx, path,  plot_idx):
       self.algo = algo
       self.algo2 = algo2
       self.eyeplan_model = eyeplan_model
       self.structure_set = structure_set
       self.slice_idx = slice_idx
       self.path = path
       self.plot_idx = plot_idx

   def plot(self):

       
        
        fig = plt.figure()
        ax = plt.subplot(221)
        
        # out = fraction1.points_at_idx(0)
        target_pts = self.algo.config["target_points"]  # out["target"]
        eyeglobe_pts = self.algo.config["eyeglobe_point"]  # out["eyeglobe"]
        
        image = self.algo.raytracer.depth[0:, self.slice_idx[1],
                                          0:] - self.algo2.raytracer.depth[0:, self.slice_idx[1], 0:]
        
        limit = np.max([abs(np.min(image)), np.max(image)])
        plt.pcolor(self.algo.medium.Z, self.algo.medium.X, image,
                   cmap="seismic", vmin=-limit, vmax=limit)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Range difference', rotation=270, size='large')
        
        
        ax.scatter([], [],
                   alpha=1, color="k", label="Treatment position")
        ax.scatter([], [], alpha=1,
                   color=my_colors[0], label="Transformed eye position")
        
        ax.scatter(eyeglobe_pts[0:, 2], eyeglobe_pts[0:, 0],
                   alpha=0.1, color="k")
        ax.scatter(target_pts[0:, 2], target_pts[0:, 0], alpha=0.1, color="k")
        
        # out = fraction1.points_at_idx(idx_of_interest)
        target_pts = self.algo2.config["target_points"]  # out["target"]
        eyeglobe_pts = self.algo2.config["eyeglobe_point"]  # out["eyeglobe"]
        
        ax.scatter(eyeglobe_pts[0:, 2], eyeglobe_pts[0:, 0], alpha=0.1,
                   color=my_colors[0])
        ax.scatter(target_pts[0:, 2], target_pts[0:, 0],
                   alpha=0.1, color=my_colors[0])
        
        plt.axvline(x=self.algo.config['skin_plane_point'][2], color="g", linewidth=4,
                    label=f"Skin plane at z = {round(self.algo.config['skin_plane_point'][2], 3)}")
        
        plt.xlabel("z [mm]", fontsize=12)
        plt.ylabel("x [mm]", fontsize=12)
        
        plt.title("Range", fontsize=14)
        plt.legend()
        
        ax = plt.subplot(222)
        
        # out = fraction1.points_at_idx(0)
        target_pts = self.algo.config["target_points"]  # out["target"]
        eyeglobe_pts = self.algo.config["eyeglobe_point"]  # out["eyeglobe"]
        
        image = self.algo.dose[0:, self.slice_idx[1],
                               0:] - self.algo2.dose[0:, self.slice_idx[1], 0:]
        
        limit = np.max([abs(np.min(image)), np.max(image)])
        plt.pcolor(self.algo.medium.Z, self.algo.medium.X, image,
                   cmap="seismic", vmin=-limit, vmax=limit)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel('Dose difference', rotation=270, size='large')
        
        ax.scatter(eyeglobe_pts[0:, 2], eyeglobe_pts[0:, 0],
                   alpha=0.1, color="k", label="Treatment position")
        ax.scatter(target_pts[0:, 2], target_pts[0:, 0], alpha=0.1, color="k")
        
        # out = fraction1.points_at_idx(idx_of_interest)
        target_pts = self.algo2.config["target_points"]  # out["target"]
        eyeglobe_pts = self.algo2.config["eyeglobe_point"]  # out["eyeglobe"]
        
        ax.scatter(eyeglobe_pts[0:, 2], eyeglobe_pts[0:, 0],
                   alpha=0.1, color=my_colors[0], label="Extreme")
        ax.scatter(target_pts[0:, 2], target_pts[0:, 0],
                   alpha=0.1, color=my_colors[0])
        
        plt.axvline(x=self.algo.config['skin_plane_point'][2], color="g", linewidth=4,
                    label=f"Skin plane at z = {round(self.algo.config['skin_plane_point'][2], 3)}")
        
        plt.xlabel("z [mm]", fontsize=12)
        plt.ylabel("x [mm]", fontsize=12)
        
        plt.title("Dose", fontsize=14)
        
        ax = plt.subplot(223)
        
        # fig = plt.figure()
        
        x = self.algo.medium.z_voxel_pos
        y = self.algo.medium.x_voxel_pos
        X, Y = np.meshgrid(x, y)
        
        struct = self.eyeplan_model.structure_set_clips_registered["eyeglobe"]
        struct.resample_contour_points_in_grid(1)
        
        ax.contour(X, Y, struct.binary_mask[0:, self.slice_idx[1], 0:], [
            1], colors='black')
        ax.contour(X, Y, self.eyeplan_model.structure_set_clips_registered["target"].binary_mask[0:, self.slice_idx[1], 0:], [
            1], colors='k')
        
        globe_struct = self.structure_set["eyeglobe"]
        globe_struct.resample_contour_points_in_grid(1)
        ax.contour(X, Y, globe_struct.binary_mask[0:, self.slice_idx[1], 0:], [
            1], colors='C0')
        
        target_struct = self.structure_set["target"]
        target_struct.resample_contour_points_in_grid(1)
        ax.contour(X, Y, target_struct.binary_mask[0:, self.slice_idx[1], 0:], [
            1], colors='C0')
        
        plt.plot([], [], color="k", label="Binary mask treatment position")
        plt.plot([], [], color="C0",
                 label="Binary mask transformed eye position")
        
        # target_pts = self.algo.config["target_points"] #out["target"]
        # eyeglobe_pts = self.algo.config["eyeglobe_point"] #out["eyeglobe"]
        
        # target_pts = np.asarray(list(filter(lambda x : abs(x[1]) < 0.5, target_pts )))
        # eyeglobe_pts = np.asarray(list(filter(lambda x : abs(x[1]) < 0.5, eyeglobe_pts )))
        
        # plt.scatter(eyeglobe_pts[0:,2], eyeglobe_pts[0:,0], color = "C0")
        # plt.scatter(target_pts[0:,2], target_pts[0:,0], color = "green")
        
        plt.xlabel("z [mm]", fontsize=12)
        plt.ylabel("x [mm]", fontsize=12)
        plt.grid(linewidth=0.3)
        
        plt.legend()
        
        plt.title("Binary masks", fontsize=14)
        
        ax = plt.subplot(224)
        
        # "lens", "cornea", "optdisk",
        for structure_name, color in zip(["target", "eyeglobe", "lens", "macula", "optdisk"], my_colors):
              struct = self.eyeplan_model.structure_set_clips_registered[structure_name]
              struct2 = self.structure_set[structure_name]
        
              dvh = DVH(self.algo.dose, struct.binary_mask)
              dvh2 = DVH(self.algo2.dose, struct2.binary_mask)
        
              # volume_fractions = np.linspace(0, 1, 100)
              dose_fractions = np.linspace(0, 1, 1000)
        
              volumes = dvh.V(dose_fractions)
              volumes2 = dvh2.V(dose_fractions)
        
              plt.plot(dose_fractions*100, volumes*100,
                      label=structure_name, color=color)
              plt.plot(dose_fractions*100, volumes2 *
                      100, color=color, linestyle="--")

        plt.plot([], [], color="k", linewidth=2, label="Reference")
        plt.plot([], [], color="k", linestyle="--",
                 label="Transformed position")
        
        plt.legend()
        
        plt.xlabel("% Dose", fontsize=12)
        plt.ylabel("% Volume or surface", fontsize=12)
        
        plt.grid(linewidth=0.3)
        
        plt.suptitle(
            f"Ref - transformed eye position, idx = {self.plot_idx}", fontsize=14)
        
        xlength = 15
        fig.set_size_inches(xlength, xlength/1.61803398875)
        plt.show()
        try:
            plt.savefig(self.path / f"plot_idx_{self.plot_idx}_",
                        bbox_inches='tight', pad_inches=0.1)
        except:
            pass
