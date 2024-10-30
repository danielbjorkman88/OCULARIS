
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""


import os
import pydicom
import numpy as np
from pyproton.utils.ioUtils import is_dicom, list_dir_absolute
from pyproton.handler.volume_handler import BaseVolumeHandler
from pyproton.handler.dicom_tags import *


class DicomVolumeHandler(BaseVolumeHandler):

    def _load_file(self, filename):

        if isinstance(filename, list):
            assert is_dicom(filename[0])
            return self._load_sliced_dicom_volume(filename)
        elif isinstance(filename, str):
            if os.path.isfile(filename):
                assert is_dicom(filename)
                return self._load_single_dicom_volume(filename)
            elif os.path.isdir(filename):
                return self._load_sliced_dicom_volume([file_path for file_path in list_dir_absolute(filename) if is_dicom(file_path)])
            else:
                raise ValueError("File {} does not exist".format(filename))
        else:
            raise ValueError("Unknown dicom format for {}".format(filename))

    def _load_single_dicom_volume(self):
        dcm = pydicom.dcmread(self.filename)
        self.pixel_array = self._read_pixel_array(dcm)
        self.voxel_to_world_matrix = self._read_voxel_to_world_matrix(dcm)

    def _read_pixel_array(self, dcm):
        array = np.transpose(dcm.pixel_array, axes=[2, 1, 0]).astype(np.float)
        array = convert_stored_values_to_array(array, dcm)
        return array

    def _read_voxel_to_world_matrix(self,dcm):
        voxel_spacing = get_voxel_spacing(dcm)
        image_position = get_image_position_patient(dcm)
        image_orientation = get_image_orientation_patient(dcm)
        return create_voxel_to_world_matrix(
            voxel_spacing, image_orientation, image_position
        )
        
    def _load_sliced_dicom_volume(self, file_list):
        pixel_array_list = []
        slice_coordinates = []
        base_dicom = pydicom.dcmread(file_list[0])
        
        for path in file_list:
            dcm_slice = pydicom.dcmread(path)
            if is_valid_slice(dcm_slice, base_dicom):
                pixel_array_list.append(
                    convert_stored_values_to_array(
                        dcm_slice.pixel_array.astype(np.float), dcm_slice)
                )
                slice_coordinates.append(dcm_slice.SliceLocation)
            else:
                raise ValueError("Invalid slice found at {}".format(path))

        assert is_valid_slice_list(
            sorted(slice_coordinates), base_dicom.SliceThickness)
        pixel_array_list = sort_slice_list(pixel_array_list, slice_coordinates)
        new_pixel_array = np.stack(pixel_array_list, axis=0)

        self.pixel_array = np.transpose(new_pixel_array, axes=[2, 1, 0]).astype(
            np.float
        )

        self.base_dicom = base_dicom
        voxel_spacing = get_voxel_spacing(self.base_dicom)
        image_position = get_image_position(self.base_dicom,slice_coordinates)
        image_orientation = get_image_orientation_patient(self.base_dicom)  
        self.voxel_to_world_matrix = create_voxel_to_world_matrix(
            voxel_spacing, image_orientation, image_position
        )


def convert_stored_values_to_array(array, base_dicom):
    if hasattr(base_dicom, RESCALE_SLOPE_NAME):
        array *= float(base_dicom.RescaleSlope)
    if hasattr(base_dicom, RESCALE_INTERCEPT_NAME):
        array = array + float(base_dicom.RescaleIntercept)
    if hasattr(base_dicom, DOSE_GRID_SCALING_NAME):
        array *= float(base_dicom.DoseGridScaling)
    return array


def get_voxel_spacing(dcm):
    voxel_spacing = list(dcm.PixelSpacing) + [dcm.SliceThickness]
    return np.array([float(spacing) for spacing in voxel_spacing])

def get_image_position_patient(dcm):
    return np.array(dcm.ImagePositionPatient)

def get_image_orientation_patient(dcm):
    # TODO: when going to non identity orientations, we have to check whether the vectors should be stacked horizontally or vertically
    orientation = dcm.ImageOrientationPatient
    xVec = np.array(orientation[:3])
    yVec = np.array(orientation[3:])
    zVec = np.cross(xVec, yVec)
    return np.stack([xVec, yVec, zVec])


def create_voxel_to_world_matrix(voxel_spacing, image_orientation, image_position):
    voxel_to_world_matrix = np.identity(4)
    voxel_to_world_matrix[:3, :3] = image_orientation * voxel_spacing
    voxel_to_world_matrix[:3, 3] = image_position
    return voxel_to_world_matrix


def is_valid_slice(slice, base_dcm):
    if not has_same_slice_thickness(slice, base_dcm):
        return False
    if not has_same_pixel_spacing(slice, base_dcm):
        return False
    if not has_same_slice_orientation(slice, base_dcm):
        return False
    return True


def has_same_slice_thickness(dcm_slice, base_dcm):
    return dcm_slice.SliceThickness == base_dcm.SliceThickness


def has_same_pixel_spacing(dcm_slice, base_dcm):
    return dcm_slice.PixelSpacing == base_dcm.PixelSpacing


def has_same_slice_orientation(dcm_slice, base_dcm):
    return dcm_slice.ImageOrientationPatient == base_dcm.ImageOrientationPatient


def is_valid_slice_list(sorted_slice_list, slice_thickness):
    sorted_slice_list = np.array(sorted_slice_list)
    diff = sorted_slice_list[1:] - sorted_slice_list[:-1]
    return np.allclose(diff, slice_thickness)


def sort_slice_list(pixel_array_list, slice_coordinates):
    return [
        slice_array
        for _, slice_array in sorted(zip(slice_coordinates, pixel_array_list))
    ]


def get_image_position(base_dicom, slice_coordinates):
    pixel_array_position = get_image_position_patient(base_dicom)
    pixel_array_position[2] = min(slice_coordinates)
    return pixel_array_position
