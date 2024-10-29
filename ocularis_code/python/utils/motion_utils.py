# -*- coding: utf-8 -*-


import os
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
from scipy.spatial.transform import Rotation
import copy
import pyvista as pv
from scipy.interpolate import interp1d
from scipy.stats import linregress
from matplotlib.colors import LogNorm
from matplotlib import cm
import math
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
import pandas as pd
from utils.transformations import mapping, rotation_matrix_from_vectors, rotation_matrix_from_vectors_scipy
    
import matplotlib.cm as cmx
import matplotlib.colors as colors
import reference_frames
import plot_utils
from optis_tps.patient import Patient
from motion.fraction import Fraction
from time import process_time

# import pyvista as pv
# from pyproton.metrics.dvh import DVH
# from pyproton.volume.grid import Grid
# import numpy as np
# from matplotlib import pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.animation
# import pandas as pd
# import imageio



from configuration import config, data
import dose_engine
from pyproton.metrics.dvh import DVH
from pyproton.volume.grid import Grid

# import reference_frames
from utils.vector_utils import angle


from math import isclose



from utils.normalization_utils import find_80_20_R, find_80_20_L
# from matplotlib.pyplot import cm

# from patient_db.database import database
from generators.dose_engine_generators import  generate_reference


from shapely.geometry import Polygon



def rotate_points(points, rot_point, vec1, vec2):
    
    if (vec1 == vec2).all():
        mat = np.eye(3)
    else:
        mat = rotation_matrix_from_vectors(
            vec1, vec2)
        mat_scipy = rotation_matrix_from_vectors_scipy(
            vec1, vec2)        
    
    # print("mat", mat)
    # print("mat_scipy", mat_scipy)

    T_ = np.eye(4)
    T_[0:3, 3] = -rot_point
    
    T_r = np.eye(4)
    T_r[0:3, 0:3] = mat
    
    T = np.eye(4)
    T[0:3, 3] = rot_point
    
    points = mapping(points, T_)
    
    points = mapping(points, T_r)
    
    points = mapping(points, T)
    
    return points


def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))



def areas_of_80_20_contour(contour):
    
    # contour.allsegs[1][0]
    
    assert len(contour.allsegs) == 2
    
    iso_20 = contour.allsegs[0]
    iso_20_x = iso_20[0][0:,0]
    iso_20_y = iso_20[0][0:,1]
    
    iso_80 = contour.allsegs[1]
    iso_80_x = iso_80[0][0:,0]
    iso_80_y = iso_80[0][0:,1]
    
    # print(PolyArea(iso_20_x, iso_20_y))
    # print(Polygon(zip(iso_20_x, iso_20_y)).area)
    
    # print(PolyArea(iso_80_x, iso_80_y))
    # print(Polygon(zip(iso_80_x, iso_80_y)).area)
    
    assert isclose(PolyArea(iso_20_x, iso_20_y) , Polygon(zip(iso_20_x, iso_20_y)).area,  abs_tol=1e-6 )
    assert isclose(PolyArea(iso_80_x, iso_80_y) , Polygon(zip(iso_80_x, iso_80_y)).area,  abs_tol=1e-6 )
    
    
    return PolyArea(iso_20_x, iso_20_y), PolyArea(iso_80_x, iso_80_y)

# a1, a2 = areas_of_80_20_contour(contour2)




def generate_area_dics(alg, tot_dose, ax):
    
    x = alg.medium.z_voxel_pos
    y = alg.medium.x_voxel_pos
    X, Y = np.meshgrid(x, y)

    area_ratio_dict = {}
    areas_dict = {}
    
    for z_bin in range(alg.dose.shape[2]):
        
        # print(f"---------------{z_bin}-----------------------------")
        
        
        z_pos = alg.medium.resolution[2] *(z_bin + 0.5) + alg.medium.mesh_origin[2]
        
        # diffs_local = []
        
        image1 = alg.dose[0:, 0:, z_bin].T
        image2 = tot_dose[0:, 0:, z_bin].T
        
        if np.max(image1) < 0.8 and np.max(image2) < 0.8:
            continue
        
        # fig = plt.figure()
        
        # ax = plt.subplot(111)
        
        contour1 = ax.contour(X, Y, image1, [0.2, 0.8], colors= "C0", linestyles = [":", "--"])
        contour2 = ax.contour(X, Y, image2, [0.2, 0.8], colors= "C1", linestyles = [":", "--"])
        
        try:
            a1, a2 = areas_of_80_20_contour(contour1)
            tot_area_ref = a1 - a2
            
            
            a1, a2 = areas_of_80_20_contour(contour2)
            tot_area_acc = a1 - a2
            
            area_ratio_dict[z_pos] = tot_area_acc/tot_area_ref
            areas_dict[z_pos] = (tot_area_acc, tot_area_ref)
        except:
            pass
        
    return area_ratio_dict, areas_dict





   
def filenames_at_path(my_path):
        return sorted(list(filter(lambda filename: filename[0:len("subfraction")] == "subfraction", os.listdir(my_path))), key = lambda filename: int(filename[12:].split(".")[0]))
    


def accumulate_dose(basepath, algo, tps_config, N_subfractions):
    
    
    
    fractions = tps_config["fractions"]
    slice_idx = algo.config["Slice"]
    
    my_subfractions = []
    tot_acc_dose = np.zeros(algo.dose.shape)
    counter = 0
    for fraction_number in fractions:
        path = Path(basepath + str(fraction_number)) 
        filenames = filenames_at_path(path)
        for filename in filenames:
            dose = np.load(path / filename)
            my_subfractions.append(dose) #[0:,slice_idx[1], 0:]
            tot_acc_dose += dose/N_subfractions
            counter +=1

    assert counter == len(my_subfractions)
    
    return tot_acc_dose, my_subfractions



