# -*- coding: utf-8 -*-
"""
Created on Fri Dec  2 16:43:55 2022

@author: bjoerk_c
"""
import numpy as np

def find_com(my_points):
    return np.asarray([sum(row[i] for row in my_points)/len(my_points) for i in range(len(my_points[0]))])
