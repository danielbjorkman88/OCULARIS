
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

from pyproton.utils.ioUtils import (
    is_mha,
    is_nifti,
    is_nrrd,
    is_pickle,
    read_pickle,
)
from pyproton.handler.base_handler import BaseHandler
from pyproton.volume.volume import Volume
import numpy as np
import medpy.io
import nrrd
import nibabel

NRRD_ORIENTATION = "space directions"
NRRD_ORIGIN = "space origin"


class BaseVolumeHandler(BaseHandler):
    def __init__(self, filename: str):
        self.filename = filename
        self._load_file(filename)

    @property
    def volume(self) -> Volume:
        return Volume(self.pixel_array, self.voxel_to_world_matrix)

    def _load_file(self, file):
        raise NotImplementedError


class MhaVolumeHandler(BaseVolumeHandler):
    def _load_file(self, file):
        assert is_mha(file)
        mha = medpy.io.load(file)
        self.pixel_array = self._read_pixel_array(mha)
        self.voxel_to_world_matrix = self._read_voxel_to_world_matrix(mha)

    def _read_pixel_array(self, mha):
        array = mha[0]
        dim = len(array.shape)
        if dim == 4:
            array = np.transpose(array, axes=[3, 0, 1, 2])
        return array

    def _read_voxel_to_world_matrix(self, mha):
        header = mha[1]
        voxel_spacing = np.array(header.spacing)
        image_position = np.array(header.offset)
        image_orientation = header.direction
        voxel_to_world_matrix = np.identity(4)
        voxel_to_world_matrix[:3, :3] = image_orientation * voxel_spacing
        voxel_to_world_matrix[:3, 3] = image_position
        return voxel_to_world_matrix


class NrrdVolumeHandler(BaseVolumeHandler):
    def _load_file(self, file):
        assert is_nrrd(file)
        nrrdFile = nrrd.read(file)
        self.pixel_array = self._read_pixel_array(nrrdFile)
        self.voxel_to_world_matrix = self._read_voxel_to_world_matrix(nrrdFile)

    def _read_pixel_array(self, nrrd):
        return nrrd[0]

    def _read_voxel_to_world_matrix(self, nrrd):
        voxel_to_world_matrix = np.identity(4)
        voxel_to_world_matrix[:3, :3] = nrrd[1][NRRD_ORIENTATION]
        voxel_to_world_matrix[:3, 3] = nrrd[1][NRRD_ORIGIN]
        return voxel_to_world_matrix


class NiftyVolumeHandler(BaseVolumeHandler):
    def _load_file(self, file):
        assert is_nifti(file)
        nifty = nibabel.load(file)
        self.pixel_array = self._read_pixel_array(nifty)
        self.voxel_to_world_matrix = self._read_voxel_to_world_matrix(nifty)

    def _read_pixel_array(self, nifti):
        array = nifti.get_fdata()
        array = np.squeeze(array)  # get rid of time dimension
        if len(array.shape) == 4:
            array = np.transpose(array, [3, 0, 1, 2])
        return array

    def _read_voxel_to_world_matrix(self, nifti):
        voxel_to_world_matrix = nifti.get_affine()
        voxel_to_world_matrix[0, :] = -voxel_to_world_matrix[
            0, :
        ]  # inverted convention?
        voxel_to_world_matrix[1, :] = -voxel_to_world_matrix[1, :]
        return voxel_to_world_matrix


class PickleVolumeHandler(BaseVolumeHandler):
    """
    Used to load saved volume objects in a similar way as other volume files.
    Likely to change because a daughter class of volume might contain many more
    attributes which would be lost in this way
    """

    def _load_file(self, file):
        assert is_pickle(file)
        volume_object = read_pickle(file)
        self.pixel_array = volume_object.pixel_array
        self.voxel_to_world_matrix = volume_object.voxel_to_world_matrix

