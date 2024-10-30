
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import pydicom


def dicom_to_dict(dataset, extract_pixel_array=False, recursive=False):
    """
    Parameters
    ----------
    dataset: pydicom.Dataset
        The dataset to convert into a Python dict.
    """
    keys = dataset.keys()
    names = [dataset[k].name for k in keys]

    entry = {}

    for name, key in zip(names, keys):
        if name == "Pixel Data":
            continue
        val = dataset[key].value

        if recursive and type(val) == pydicom.sequence.Sequence:
            val_converted = []
            for v in val:
                v_conv = dicom_to_dict(v, False, recursive)
                val_converted.append(v_conv)
            val = val_converted

        entry[name] = val

    if extract_pixel_array:
        entry["Pixel Data"] = dataset.pixel_array

    return entry
