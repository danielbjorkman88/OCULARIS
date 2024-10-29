import os, pickle, json
from pyproton.utils.constants import *

"""
Contains input-output functionality with the filesystem. Many functions are small and used for improved code readability
in other parts of the codebase
"""


def list_dir_absolute(dirPath):
    return [os.path.join(dirPath, filename) for filename in list_dir_relative(dirPath)]


def list_dir_relative(dirPath):
    return os.listdir(dirPath)


def is_mha(file_path):
    return is_file_type(file_path, MHA_FILE_EXTENSION)


def is_numpy(file_path):
    return is_file_type(file_path, NUMPY_ARRAY_EXTENSION)


def is_gz(file_path):
    return is_file_type(file_path, GZ_FILE_EXTENSION)


def is_pickle(file_path):
    return is_file_type(file_path, PKL_FILE_EXTENSION, allow_zipped=True)


def is_nrrd(file_path):
    return is_file_type(file_path, NRRD_FILE_EXTENSION)


def is_nifti(file_path):
    return is_file_type(file_path, NIFTI_FILE_EXTENSION)


def is_bdf(file_path):
    return is_file_type(file_path, BDF_FILE_EXTENSION)


def is_json(file_path):
    return is_file_type(file_path, JSON_FILE_EXTENSTION)


def is_dicom(file_path):
    if isinstance(file_path, list):
        return all(is_dicom(path) for path in file_path)
    if is_file_type(file_path, DICOM_FILE_EXTENSION):
        return True
    if os.path.isdir(file_path):
        return all(is_dicom(path) for path in list_dir_absolute(file_path))
    return False


def is_file_type(file_path, type_extension, allow_zipped=False):
    ext = get_file_extension(file_path)
    if ext == type_extension:
        return True
    if allow_zipped:
        if get_zipped_file_extension(file_path) == type_extension:
            return True
    return False


def get_file_extension(file_path):
    return os.path.splitext(file_path)[1][1:]


def get_zipped_file_extension(file_path):
    return os.path.splitext(get_file_name(file_path))[1][1:]


def get_file_name(filepath):
    return os.path.split(os.path.splitext(filepath)[0])[1]


def read_json(file_path):
    with open(file_path) as json_file:
        json_object = json.load(json_file)
        json_file.close()
    return json_object


def write_json(file_path, data):
    with open(file_path, "w") as json_file:
        json.dump(data, json_file, indent=4)
        json_file.close()
    return


def read_pickle(file_path):
    with open(file_path, "rb") as pkl:
        pickle_object = pickle.load(pkl)
    return pickle_object


def write_pickle(file_path, data):
    with open(file_path, "wb") as pkl:
        pickle.dump(data, pkl)
    return
