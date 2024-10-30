# -*- coding: utf-8 -*-
"""

This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import numpy as np
import pandas as pd
import random
from scipy import interpolate
import matplotlib.pyplot as plt

# from patient_db.database import database

def get_displacements(path):

    cohort = [17027, 17029, 17034, 17111, 17114, 17115, 17116, 17117, 17186,
              17197, 17212, 17213, 17216, 17218, 17247, 17301]  # 17111 # 17186 #17115

    all_xes = []

    for patient_number in cohort:

        patient_string =  "P" + str(patient_number)
        pat1 = database(patient_string)
        tps_config = pat1.tps_config
        ep_model = pat1.patient_model.eyeplan_model
        config2 = pat1.dose_engine_config
        skin_plane = pat1.dose_engine_config["skin_plane_point"][2]

        df = pd.read_csv(
            path / f"Patient_{tps_config['patient_number']}_motion_dataframe_full.csv")

        centre = df["xs"][0]
        xes = df["xs"] - centre
        all_xes.extend(xes)

    all_yes = []

    for patient_number in cohort:

        patient_string =  "P" + str(patient_number)
        pat1 = database(patient_string)
        tps_config = pat1.tps_config
        ep_model = pat1.patient_model.eyeplan_model
        config2 = pat1.dose_engine_config
        skin_plane = pat1.dose_engine_config["skin_plane_point"][2]

        df = pd.read_csv(
            path / f"Patient_{tps_config['patient_number']}_motion_dataframe_full.csv")

        centre = df["ys"][0]
        xes = df["ys"] - centre
        all_yes.extend(xes)

    all_zes = []

    for patient_number in cohort:

        patient_string =  "P" + str(patient_number)
        pat1 = database(patient_string)
        tps_config = pat1.tps_config
        ep_model = pat1.patient_model.eyeplan_model
        config2 = pat1.dose_engine_config
        skin_plane = pat1.dose_engine_config["skin_plane_point"][2]

        df = pd.read_csv(
            path / f"Patient_{tps_config['patient_number']}_motion_dataframe_full.csv")

        centre = df["zs"][0]
        xes = df["zs"] - centre
        all_zes.extend(xes)

    return all_xes, all_yes, all_zes



def generate_mc_from_cumulative(f_rev_x, f_rev_y, f_rev_z, N_random_numbers):
    random_vals_x = np.zeros(N_random_numbers)
    random_vals_y = np.zeros(N_random_numbers)
    random_vals_z = np.zeros(N_random_numbers)
    for i in range(N_random_numbers):
        rand_variable_x = random.random()
        rand_variable_y = random.random()
        rand_variable_z = random.random()
        #random_sequence[i] = rand_variable
        random_vals_x[i] = f_rev_x(rand_variable_x)
        random_vals_y[i] = f_rev_y(rand_variable_y)
        random_vals_z[i] = f_rev_z(rand_variable_z)
        
    my_data = {"random_vals_x": random_vals_x, "random_vals_y": random_vals_y, "random_vals_z": random_vals_z}
    # centers = { "centers_x": centers_x, "centers_y": centers_y, "centers_z": centers_z}
    
    
       
    # df = pd.DataFrame(my_data)
    # df.to_csv(out_path / f"Patient_{tps_config['patient_number']}_dataframe.csv")  
        
        
    return pd.DataFrame(my_data) 








def get_intrafractional_displacements(path, N_bins,  N_random_numbers):
    
    """
    Intra-fractional motion
    
    """

    ax = plt.subplot(111)    

    all_xes, all_yes, all_zes = get_displacements(path)
    
    n_x, bins_x, patches = ax.hist(all_xes, N_bins, density=True, histtype='step',
                               cumulative=True)
    xes_centers = bins_x[0:-1] + (bins_x[1] - bins_x[0])/2
    f_rev_x = interpolate.interp1d(n_x, xes_centers, fill_value='extrapolate')



    n_y, bins_y, patches = ax.hist(all_yes, N_bins, density=True, histtype='step',
                               cumulative=True)
    yes_centers = bins_y[0:-1] + (bins_y[1] - bins_y[0])/2
    f_rev_y = interpolate.interp1d(n_y, yes_centers, fill_value='extrapolate')


    n_z, bins_z, patches = ax.hist(all_zes, N_bins, density=True, histtype='step',
                               cumulative=True)
    zes_centers = bins_z[0:-1] + (bins_z[1] - bins_z[0])/2
    f_rev_z = interpolate.interp1d(n_z, zes_centers, fill_value='extrapolate' )
    
    random_vals_df = generate_mc_from_cumulative(f_rev_x, f_rev_y, f_rev_z, N_random_numbers)
    
    return random_vals_df




def get_setup_displacements(N_iterations):
    """
    Setup uncertainty
    
    """
    
    mu = 0
    # sigma_x = 0.25 
    # sigma_y = 0.44
    # sigma_z = 0.37
    
    sigma_x = 0.17 
    sigma_y = 0.15
    sigma_z = 0.21    


    x_random = np.random.normal(mu, sigma_x, N_iterations)
    y_random = np.random.normal(mu, sigma_y, N_iterations)
    z_random = np.random.normal(mu, sigma_z, N_iterations)
    
    return x_random, y_random, z_random


