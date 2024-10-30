# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

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

