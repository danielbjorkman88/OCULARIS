# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import os
import sys
import pandas as pd
import math

import numpy as np
from pathlib import Path

from scipy.interpolate import interp1d

from scipy.ndimage import uniform_filter1d

from utils.transformations import mapping, rotation_matrix_from_vectors
from utils.com import find_com


currentdir = Path(os.getcwd())
newdir = currentdir.parent.parent / "pyproton"
if newdir not in sys.path:
    sys.path.insert(1, str(newdir))


class Fraction():

    def __init__(self, eyeplan_model, fraction_number, avoid_filling_table = False):
        
        try:
            self.path = Path(
                eyeplan_model.patient_config["motion"] + "\\fraction" + str(fraction_number))
        except:
            self.path = Path(eyeplan_model.patient_config["motion"] / ("fraction" + str(fraction_number) ))
            
        

        self._patient_number = eyeplan_model.patient_config["patient_number"]
        self._fraction_number = fraction_number
        self.M_translation = None

        self._pps = pd.read_csv(
            self.path / f'Patient{self._patient_number}_fraction{self._fraction_number}_pp.csv', names=['x', 'y', 'z'])
        self._cps = pd.read_csv(
            self.path / f'Patient{self._patient_number}_fraction{self._fraction_number}_cp.csv', names=['x', 'y', 'z'])
        self.time = pd.read_csv(
            self.path / f'Patient{self._patient_number}_fraction{self._fraction_number}_time.csv', names=['t'])

        if not len(self._pps) == len(self._cps):
            raise ValueError
        
        nan_indices = []
        for df in [self._pps, self._cps]:
            for key in df.keys():
                    index = df[key].index[df[key].apply(np.isnan)]
                    for item in index:
                        if item not in nan_indices:
                            nan_indices.append(item)
                                                    
        if len(nan_indices) > 0:
            self._pps = self._pps.drop(nan_indices)
            self._cps = self._cps.drop(nan_indices)
            print("Nan values dropped from dataframe")
            

        self.N = len(self._pps)

        self.treatment_time = self.time["t"][0]

        self.interval = 1000 * self.treatment_time / self.N  # ms

        t = np.linspace(self.interval/2, self.interval *
                        self.N, self.N)

        self._df = pd.DataFrame({'t': t,
                                 'pps.x': uniform_filter1d(self._pps.x, size=5),
                                 'pps.y': uniform_filter1d(self._pps.y, size=5),
                                 'pps.z': uniform_filter1d(self._pps.z, size=5),
                                 'cps.x': uniform_filter1d(self._cps.x, size=5),
                                 'cps.y': uniform_filter1d(self._cps.y, size=5),
                                 'cps.z': uniform_filter1d(self._cps.z, size=5)})
        

        new_t = self._df['t'][::int(len(self._df)*1000/max(t))]
        self._df_interpolated = pd.DataFrame({'t': new_t,
                                              'pps.x': interp1d(self._df['t'], self._df['pps.x'])(new_t),
                                              'pps.y': interp1d(self._df['t'], self._df['pps.y'])(new_t),
                                              'pps.z': interp1d(self._df['t'], self._df['pps.z'])(new_t),
                                              'cps.x': interp1d(self._df['t'], self._df['cps.x'])(new_t),
                                              'cps.y': interp1d(self._df['t'], self._df['cps.y'])(new_t),
                                              'cps.z': interp1d(self._df['t'], self._df['cps.z'])(new_t)}).reset_index().drop(columns=['index'])
        
        self.N_resampled = len(self._df_interpolated)

        self.eyeplan_model = eyeplan_model

        self.eyeglobe_points_start = self.eyeplan_model.structure_set_clips_registered[
            "eyeglobe"].structure_points
        self.target_points_start = self.eyeplan_model.structure_set_clips_registered[
            "target"].structure_points
        self.lens_points_start = self.eyeplan_model.structure_set_clips_registered[
            "lens"].structure_points
        self.optdisk_points_start = self.eyeplan_model.structure_set_clips_registered[
            "optdisk"].structure_points
        self.macula_points_start = self.eyeplan_model.structure_set_clips_registered[
            "macula"].structure_points        
        
        self.target_com_start = find_com(self.target_points_start)
        
        if avoid_filling_table:
            pass
        else:
            self.fill_table()
            self.fill_table_resampled()
            
        self.define_transformation_matrices()

    def __repr__(self):
        return f"Patient {self._patient_number}, fraction {self._fraction_number}"

    def define_transformation_matrices(self):

        x, y, z, xp, yp, zp = self.get_arrow(0)
        # pos1 = np.asarray([x, y, z])
        vec = np.asarray([xp, yp, zp])
        vec /= np.linalg.norm(vec)

        orthogonal1 = np.asarray([1, 0, 0])
        orthogonal1 = orthogonal1 - orthogonal1.dot(vec) * vec
        orthogonal1 /= np.linalg.norm(orthogonal1)
        orthogonal2 = np.cross(vec, orthogonal1)

        M = np.eye(4)
        M[0:3, 0] = orthogonal1
        M[0:3, 1] = orthogonal2
        M[0:3, 2] = vec
        self.M = M

        M_translation = np.eye(4)
        M_translation[0:3, 3] = np.asarray([-x, -y, -z])
        self.M_translation = M_translation

    @property
    def df(self):

        p0 = self.get_point_at_idx(0)

        # if self._df == None:
        step_distances = np.zeros(self.N)
        distance_to_origin = np.zeros(self.N)
        for i in range(self.N):
            p1 = self.get_point_at_idx(i)
            distance_to_origin[i] = math.dist(p0, p1)
            if i + 1 < self.N:
                p2 = self.get_point_at_idx(i + 1)
                dist = math.dist(p1, p2)
                step_distances[i] = dist

        self._df['step'] = step_distances
        self._df['dist_origin'] = distance_to_origin

        return self._df

    @property
    def df_resampled(self):

        p0 = self.get_point_at_idx(0, True)

        #if self._df == None:
        step_distances = np.zeros(self.N_resampled)
        distance_to_origin = np.zeros(self.N_resampled)
        for i in range(self.N_resampled):
            p1 = self.get_point_at_idx(i, True)
            distance_to_origin[i] = math.dist(p0, p1)
            if i + 1 < self.N_resampled:
                p2 = self.get_point_at_idx(i + 1, True)
                dist = math.dist(p1, p2)
                step_distances[i] = dist

        self._df_interpolated['step'] = step_distances
        self._df_interpolated['dist_origin'] = distance_to_origin

        return self._df_interpolated

    def export_df(self, filename):
        self._df.to_csv(filename)

    def get_arrow(self, idx, downsampled = False):
        
        if downsampled:
            p_pp = np.asarray(
                [self._df["pps.x"][idx], self._df["pps.y"][idx], self._df["pps.z"][idx]])
            p_cp = np.asarray(
                [self._df["cps.x"][idx], self._df["cps.y"][idx], self._df["cps.z"][idx]])
        else:
            p_pp = np.asarray(
                [self._df["pps.x"][idx], self._df["pps.y"][idx], self._df["pps.z"][idx]])
            p_cp = np.asarray(
                [self._df["cps.x"][idx], self._df["cps.y"][idx], self._df["cps.z"][idx]])

        x = p_pp[0]
        y = p_pp[1]
        z = p_pp[2]

        vec = p_pp - p_cp

        u = vec[0]
        v = vec[1]
        w = vec[2]

        return x, y, z, u, v, w

    def get_point_at_idx(self, idx, downsampled=False):

        if downsampled:
            p_cp = np.asarray(
                [self._df_interpolated["cps.x"][idx], self._df_interpolated["cps.y"][idx], self._df_interpolated["cps.z"][idx]])
        else:
            p_cp = np.asarray(
                [self._df["cps.x"][idx], self._df["cps.y"][idx], self._df["cps.z"][idx]])
        return p_cp

    def get_vector_at_idx(self, idx, downsampled = False):
        
        if downsampled:
            p_pp = np.asarray(
                [self._df_interpolated["pps.x"][idx], self._df_interpolated["pps.y"][idx], self._df_interpolated["pps.z"][idx]])
            p_cp = np.asarray(
                [self._df_interpolated["cps.x"][idx], self._df_interpolated["cps.y"][idx], self._df_interpolated["cps.z"][idx]])
        else:
            p_pp = np.asarray(
                [self._df["pps.x"][idx], self._df["pps.y"][idx], self._df["pps.z"][idx]])
            p_cp = np.asarray(
                [self._df["cps.x"][idx], self._df["cps.y"][idx], self._df["cps.z"][idx]])

        vec = p_pp - p_cp
        return vec

    def points_at_idx(self, idx, downsampled=False):

        out = {}

        rotation_point = self.translate_points_according_to_idx(
            self.eyeplan_model.centre_of_model, idx, downsampled)

        target_pts = self.translate_points_according_to_idx(
            self.target_points_start, idx, downsampled)
        target_pts = self.rotate_points_according_to_idx(
            rotation_point, target_pts, idx, downsampled)

        eyeglobe_pts = self.translate_points_according_to_idx(
            self.eyeglobe_points_start, idx, downsampled)
        eyeglobe_pts = self.rotate_points_according_to_idx(
            rotation_point, eyeglobe_pts, idx, downsampled)

        lens_pts = self.translate_points_according_to_idx(
            self.lens_points_start, idx, downsampled)
        lens_pts = self.rotate_points_according_to_idx(
            rotation_point, lens_pts, idx, downsampled)

        optdisk_pts = self.translate_points_according_to_idx(
            self.optdisk_points_start, idx, downsampled)
        optdisk_pts = self.rotate_points_according_to_idx(
            rotation_point, optdisk_pts, idx, downsampled)

        macula_pts = self.translate_points_according_to_idx(
            self.macula_points_start, idx, downsampled)
        macula_pts = self.rotate_points_according_to_idx(
            rotation_point, macula_pts, idx, downsampled)

        out["target"] = target_pts
        out["eyeglobe"] = eyeglobe_pts
        out["lens"] = lens_pts
        out["optdisk"] = optdisk_pts
        out["macula"] = macula_pts

        return out

    def transform_points_by_idx(self, points, idx, downsampled=False):

        rotation_point = self.translate_points_according_to_idx(
            self.eyeplan_model.centre_of_model, idx, downsampled)

        my_points = self.translate_points_according_to_idx(
            points, idx, downsampled)
        my_points = self.rotate_points_according_to_idx(
            rotation_point, my_points, idx, downsampled)

        return my_points

    def fill_table(self):

        centers_of_target = np.zeros([self.N, 3])
        distances_target_to_start = np.zeros(self.N)

        for idx in list(range(0, self.N)):

            # ps = points_at_idx(idx)
            ps = self.points_at_idx(idx)
            target_pts = ps["target"]
            eyeglobe_pts = ps["eyeglobe"]

            com = find_com(target_pts)
            centers_of_target[idx, 0:] = np.asarray([com[0], com[1], com[2]])
            distances_target_to_start[idx] = math.dist(
                com, self.target_com_start)

        self._df["target.com.x"] = centers_of_target[0:, 0]
        self._df["target.com.y"] = centers_of_target[0:, 1]
        self._df["target.com.z"] = centers_of_target[0:, 2]

        self._df['distances target to start'] = distances_target_to_start
        
        
    
    def fill_table_resampled(self):

        centers_of_target = np.zeros([self.N_resampled, 3])
        distances_target_to_start = np.zeros(self.N_resampled)

        for idx in list(range(0, self.N_resampled)):
            # ps = points_at_idx(idx)
            ps = self.points_at_idx(idx, True)
            target_pts = ps["target"]
            eyeglobe_pts = ps["eyeglobe"]

            com = find_com(target_pts)
            centers_of_target[idx, 0:] = np.asarray([com[0], com[1], com[2]])
            distances_target_to_start[idx] = math.dist(
                com, self.target_com_start)

        self._df_interpolated["target.com.x"] = centers_of_target[0:, 0]
        self._df_interpolated["target.com.y"] = centers_of_target[0:, 1]
        self._df_interpolated["target.com.z"] = centers_of_target[0:, 2]

        self._df_interpolated['distances target to start'] = distances_target_to_start

    def get_relative_rotation_matrix(self, idx, downsampled = False):
        
        vec1 = self.get_vector_at_idx(0, downsampled)
        vec2 = self.get_vector_at_idx(idx, downsampled)
        
        if (vec1 == vec2).all():
            return np.eye(3)

        mat = rotation_matrix_from_vectors(
            vec1, vec2)
        return mat

    def get_relative_translation_matrix(self, idx, downsampled = False):
        
        p1 = self.get_point_at_idx(0, downsampled)
        p2 = self.get_point_at_idx(idx, downsampled)

        M_translation = np.eye(4)
        M_translation[0:3, 3] = np.asarray(
            [-(p1[0] - p2[0]), - (p1[1] - p2[1]), - (p1[2] - p2[2])])
        return M_translation

    def translate_points_according_to_idx(self, points, idx, downsampled = False):
        M = self.get_relative_translation_matrix(idx, downsampled)
        points = mapping(points, M)
        return points

    def rotate_points_according_to_idx(self, rot_point, points, idx, downsampled = False):

        # cntEYE = ep_model.structure_set_clips_registered["cor"].contour_coordinates[0]
        #t_c = ep_model.centre_of_model - cntEYE

        # T_t = np.eye(4)
        # T_t[0:3, 3] = -rot_point
        
        T_ = np.eye(4)
        T_[0:3, 3] = -rot_point

        T_r = np.eye(4)
        T_r[0:3, 0:3] = self.get_relative_rotation_matrix(idx, downsampled)

        T = np.eye(4)
        T[0:3, 3] = rot_point

        # points2 = mapping(points, T_t)

        points = mapping(points, T_)

        points = mapping(points, T_r)

        points = mapping(points, T)
        return points

    def get_arrow_nominal(self, idx):
        p_pp = np.asarray(
            [self._df["pps.x"][idx], self._df["pps.y"][idx], self._df["pps.z"][idx]])
        p_cp = np.asarray(
            [self._df["cps.x"][idx], self._df["cps.y"][idx], self._df["cps.z"][idx]])

        x = p_pp[0]
        y = p_pp[1]
        z = p_pp[2]
        pos = np.asarray([x, y, z])

        nominal_pos = _transform_point(self.M_translation, pos)

        vec = p_pp - p_cp
        vec /= np.linalg.norm(vec)
        vec = _transform_vec(self.M, vec)

        u = vec[0]
        v = vec[1]
        w = vec[2]

        return nominal_pos[0], nominal_pos[1], nominal_pos[2], u, v, w
    
    
    def get_arrow_relative(self, idx):
        p_pp = np.asarray(
            [self._df["pps.x"][idx], self._df["pps.y"][idx], self._df["pps.z"][idx]])
        p_cp = np.asarray(
            [self._df["cps.x"][idx], self._df["cps.y"][idx], self._df["cps.z"][idx]])

        x = p_pp[0]
        y = p_pp[1]
        z = p_pp[2]
        pos = np.asarray([x, y, z])

        nominal_pos = _transform_point(self.M_translation, pos)

        vec = p_pp - p_cp
        norm = np.linalg.norm(vec)
        vec /= norm
        vec = _transform_vec(self.M, vec)
        vec = norm *vec

        u = vec[0]
        v = vec[1]
        w = vec[2]

        return nominal_pos[0], nominal_pos[1], nominal_pos[2], u, v, w


def _transform_vec(M, my_vec):

    vec_extended = np.asarray([my_vec[0], my_vec[1], my_vec[2], 1])

    new_vec = M.T.dot(vec_extended)

    return np.asarray([new_vec[0], new_vec[1], new_vec[2]])


def _transform_point(M_translation, my_point):

    point_extended = np.asarray([my_point[0], my_point[1], my_point[2], 1])

    new_point = M_translation.dot(point_extended)

    return np.asarray([new_point[0], new_point[1], new_point[2]])
