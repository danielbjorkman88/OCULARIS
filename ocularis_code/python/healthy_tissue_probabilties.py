# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import numpy as np
import matplotlib.pyplot as plt

def probability_visual_acuity(D2_macula):
    """
    NTCP Espensen et al
    
    probability_visual_acuity deteoration , ie maculopathy
    
    Linear interpolation function between the points

    Parameters:
    - D2_macula (float): Near maxium dose to macula. Input value between 0 and 1.

    Returns:
    - float: Interpolated value.
    """
    
    # Espensen treating with 52 Gy scaled to eGy and accounting for different fractionation scheme with factor 0.96
    end_espensen = round(52*0.96*1.1/60,2)
    
    # Define the points
    x0, y0 = 0, 0.42
    x1, y1 = end_espensen, 0.64
    
    slope = (y1 - y0) / (x1 - x0)

    # Perform linear interpolation/extrapolation
    interpolated_value = y0 + slope * (D2_macula - x0)
    
    y_max = y0 + slope * (1 - x0)

    
    if interpolated_value < y0:
        interpolated_value = y0
    elif interpolated_value > y_max:
        interpolated_value = y_max
        
    if interpolated_value > 1.0:
        interpolated_value = 1.0
    if interpolated_value < 0:
        interpolated_value = 0        

    return interpolated_value



def probability_neovascular_glaucoma_Espensen(Cornea_D20):
    """
    NTCP Espensen et al

    Parameters
    ----------
    Cornea_D20 : TYPE
        DESCRIPTION.

    Returns
    -------
    interpolated_value : TYPE
        DESCRIPTION.

    """
    
    # Espensen treating with 52 Gy scaled to eGy and accounting for different fractionation scheme with factor 0.96
    end_espensen = round(52*0.96*1.1/60,2)
    

    # Define the points
    x0, y0 = 0, 0.05
    x1, y1 = end_espensen, 0.16
    
    slope = (y1 - y0) / (x1 - x0)

    # Perform linear interpolation/extrapolation
    interpolated_value = y0 + slope * (Cornea_D20 - x0)
    
    y_max = y0 + slope * (1 - x0)

    # Perform linear interpolation
    #interpolated_value = y0 + (y1 - y0) * Cornea_D20
    
    if interpolated_value < y0:
        interpolated_value = y0
    elif interpolated_value > y_max:
        interpolated_value = y_max

    if interpolated_value > 1.0:
        interpolated_value = 1.0
    if interpolated_value < 0:
        interpolated_value = 0          

    return interpolated_value



def probability_ocular_hypertension(Ciliary_body_D20):
    """
    NTCP Espensen et al

    Parameters
    ----------
    Ciliary_body_D20 : TYPE
        DESCRIPTION.

    Returns
    -------
    interpolated_value : TYPE
        DESCRIPTION.

    """

    # Define the points
    x0, y0 = 0, 0.02
    x1, y1 = 1, 0.1

    # Perform linear interpolation
    interpolated_value = y0 + (y1 - y0) * Ciliary_body_D20
    
    if interpolated_value > 1.0:
        interpolated_value = 1.0
    if interpolated_value < 0:
        interpolated_value = 0          

    return interpolated_value
    
    
def probability_optic_neuropathy(Optic_disc_D20):
    """
    NTCP Espensen et al

    Parameters
    ----------
    Optic_disc_D20 : TYPE
        DESCRIPTION.

    Returns
    -------
    interpolated_value : TYPE
        DESCRIPTION.

    """
    # Espensen treating with 52 Gy scaled to eGy and accounting for different fractionation scheme with factor 0.96
    end_espensen = round(52*0.96*1.1/60,2)

    x0, y0 = 0, 0.08
    x1, y1 = end_espensen, 0.16
    
    slope = (y1 - y0) / (x1 - x0)

    # Perform linear interpolation/extrapolation
    interpolated_value = y0 + slope * (Optic_disc_D20 - x0)
    
    y_max = y0 + slope * (1 - x0)

    # Perform linear interpolation
    #interpolated_value = y0 + (y1 - y0) * Optic_disc_D20
    
    if interpolated_value < y0:
        interpolated_value = y0
    elif interpolated_value > y_max:
        interpolated_value = y_max
        
    if interpolated_value > 1.0:
        interpolated_value = 1.0
    if interpolated_value < 0:
        interpolated_value = 0          

    return interpolated_value
    
    
    


def probability_cataract(ciliary_body_V26):
    """
    NTCP Espensen et al
    

    Parameters
    ----------
    Ciliary_body_V26 : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    # Define the points
    x_values = [0, 0.25, 0.5, 0.75, 1]
    y_values = [0.19, 0.32, 0.55, 0.78, 0.95]

    # Perform linear interpolation
    interpolated_value = np.interp(ciliary_body_V26, x_values, y_values)

    return interpolated_value





def probability_retinal_detachment(retina_V95):
    """
    NTCP Espensen et al
    

    Parameters
    ----------
    retina_V95 : float
        Corresponds to Espensen 52 Gy.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    x_values = [0, 0.1,0.25, 0.5, 0.6, 0.75, 1]
    y_values = [0.12, 0.21 , 0.39, 0.75, 0.83 ,0.92, 0.98]


    # Perform linear interpolation
    interpolated_value = np.interp(retina_V95, x_values, y_values)

    return interpolated_value


def probability_cataract(ciliary_body_V28):
    """
    NTCP Espensen et al
    

    Parameters
    ----------
    ciliary_body_V28 : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """

    # Define the points
    def sigmoid(x):
        return 1 / (1 + np.exp(-x))

    # Define the points
    x_values = [0, 0.25, 0.5, 0.75, 1]
    y_values = [0.18, 0.32, 0.5, 0.75, 0.95]

    # Perform linear interpolation
    interpolated_value = np.interp(ciliary_body_V28, x_values, y_values)

    return interpolated_value




def increased_risk_neovascular_glaucoma_Mishra(macula_max, optic_disc_max):
    """
    NTCP Mishra et al

    Parameters
    ----------
    Optic_disc_D20 : TYPE
        DESCRIPTION.

    Returns
    -------
    interpolated_value : TYPE
        DESCRIPTION.

    """
    
    threshold = 0.28/60
    
    if macula_max > threshold and optic_disc_max > threshold:
        return 1

    return 0

    
    


def increased_risk_neovascular_glaucoma_Mishra(macula_V28, optic_disc_V28):
    """
    NTCP Mishra et al
    https://www.sciencedirect.com/science/article/pii/S0360301613006494

    Parameters
    ----------
    Optic_disc_D20 : TYPE
        DESCRIPTION.

    Returns
    -------
    interpolated_value : TYPE
        DESCRIPTION.

    """
    
    threshold = 0.28/60 # fraction of 28 Gy
    
    if macula_V28 > threshold and optic_disc_V28 > threshold:
        return 1

    return 0



def increased_risk_cataract_Thariat(lens_v5):
    """
    NTCP Thariat et al 2017

    Parameters
    ----------
    Optic_disc_D20 : TYPE
        DESCRIPTION.

    Returns
    -------
    interpolated_value : TYPE
        DESCRIPTION.

    """
    
    threshold = 10/60 # fraction of 10 Gy
    
    if lens_v5 > threshold:
        return 1

    return 0











# x_value = 0.5  # Replace with your desired value between 0 and 1
# result = probability_visual_acuity(x_value)
# print(f"For x = {x_value}, the interpolated value is: {result}")





# fig = plt.figure()

# # Generate a range of x values
# x_values = np.linspace(0, 1, 1000)

# # Calculate corresponding y values using the probability_cataract function
# y_values = [probability_cataract(x) for x in x_values]

# # Plot the results
# plt.plot(x_values, y_values, label='Sigmoid Interpolation')
# plt.xlabel('')
# plt.ylabel('Interpolated Value')
# plt.title('Sigmoid Interpolation Function')
# plt.legend()
# plt.show()

