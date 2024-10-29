# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 13:38:25 2021

@author: bjoerk_c
"""

# Treatment room reference system

config = {}
config["Target_range"] = 0 #mm # 12.0 #mm
config["Modulation_range"] = 0 #mm
config["Aperture"] = 35 #mm
config["VPS"] = 1750 #mm Virtual Point Source distance to isocentre
config["Nvoxels"] = [160,160,100] 
config["Mesh dimensions"] = [80,80,50] #mm
config["SID"] = 50 #mm Distance of phantom back surface before Isocentre 
config["AID"] = 70 #mm Aperture to Isocentre Distance
config["FID"] = 132.5 #mm Fixation light point to Isocentre Distance
config["Medium"] = "Phantom" # Phantom, Dicom or GeoModel
config["proximal_margin"] = 0 #mm
config["distal_margin"] = 0 #mm

config["Image"] = "3D" # 3D, xslice, yslice, zslice
config["Slice"] = [80,80,50]

config["Gantry rot theta"] = 0 # rotation in xz plane [Degrees]
config["Gantry rot phi"] = 0   # rotation in yz plane [Degrees]

config["wedge_insert_angle"] = 0 # azimuth angle in Treatment Room in which wedge is inserted
config["wedge_cover"] = 0 # 0 = no cover, 35 = full aperture covered [mm]
config["wedge_angle"] = 0 #degrees

config["skin_plane"] = 0 # z position of skin plane [mm]
config["collimator_ruler_length"] = 3.5
config["collimator_expansion"] = 2.5


config["eyeglobe_mesh_triangle_size"] = 12
config["eyelid_mesh_triangle_size"] = 12

# Needed to set part of medium to air
config["n_voxels_front_of_skinplane"] = 7 
# Note __init__ in optis_tps.patient

constants = {}
constants["SID"] = 50 #mm Distance of phantom back surface before Isocentre 
constants["AID"] = 70 #mm Aperture to Isocentre Distance
constants["FID"] = 132.5 #mm Fixation light point to Isocentre Distance
constants["VPS"] = 1750 #mm Virtual Point Source distance to isocentre


data = {}
data["tune_foil"] = "tune_foil"
data["DepthProfile_ESS16"] = "ESS16_z_31_05_23"
data["DepthProfile_ESS18"] = ""
data["DepthProfile_ESS20"] = ""

#data["Tunes"] = ["priscine_ISS20_z", "priscine_ISS18_z", "priscine_ISS16_z", "priscine_ISS12_z"]
#data["DepthDose"] = ["ESS12_combined_xes.npy", "ESS12_combined_yes.npy"]
data["DepthProfile"] = "R35v1_SBP_COLL35_RS045_SC2.3pos250mm.dat"
data["MW"] =["MWsR11v2","MWsR14v3", "MWsR17v2" ,"MWsR20v3", "MWsR26v2" ,"MWsR23v3", "MWsR29v2", "MWsR32v3" , "MWsR35v1"]
data["MWwatereq"] = "MWwatereq"
data["MWweights"] = "MWweights"
