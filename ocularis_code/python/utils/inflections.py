# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 16:03:24 2022

@author: bjoerk_c
"""

from scipy.interpolate import interp1d

def find_inflexions(image, x_pos, z_pos):
    
    out = {}
    
    for idx, z_pos in enumerate(z_pos):
    
        xes, yes = x_pos, image[0:,idx]
        
        if (yes == 0).all():
            continue
        
        f = interp1d(xes, yes)
        f_rev = interp1d(yes, xes)
        
        subx_pre = xes[xes < 0.99*f_rev(max(yes))]
        suby_pre = f(subx_pre)
        
        subx_post = xes[xes > 0.99*f_rev(max(yes))]
        suby_post = f(subx_post)
        
        f_sub_pre = interp1d(suby_pre, subx_pre)
        f_sub_post = interp1d(suby_post, subx_post)    
        
        R50_pre = f_sub_pre(max(suby_pre)*0.5).max()
        R50_post = f_sub_post(max(suby_post)*0.5).max()
        
        out[z_pos] = (R50_pre, R50_post)
    
    return out

