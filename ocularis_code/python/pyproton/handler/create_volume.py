from multiprocessing.sharedctypes import Value
from pyproton.utils.ioUtils import (
    get_file_extension,
    is_mha,
    is_nifti,
    is_nrrd,
    is_pickle,
    is_dicom
)
from pyproton.handler.volume_handler import *
from pyproton.handler.dicom_volume_handler import DicomVolumeHandler


def create_volume(file) -> Volume:
    if is_dicom(file):
        return DicomVolumeHandler(file).volume
    if is_mha(file):
        return MhaVolumeHandler(file).volume
    if is_nrrd(file):
        return NrrdVolumeHandler(file).volume
    if is_nifti(file):
        return NiftyVolumeHandler(file).volume
    if is_pickle(file):
        return PickleVolumeHandler(file).volume
    raise ValueError(f'Filetype {get_file_extension(file)} not recognized for creating a Volume from')

