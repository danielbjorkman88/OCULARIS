# -*- coding: utf-8 -*-
"""
This file is part of the OCULARIS Ocular Proton Therapy Treatment Planning System. 
It is subject to the license terms in the LICENSE file located in the top-level directory of this distribution.

This program is not certified for clinical use and is provided WITHOUT ANY WARRANTY or implied warranty.
For accuracy, users should validate OCULARIS independently before drawing any conclusions.

"""

import pandas as pd

def save_orientation_data_row(new_row , path, filename):


    # Check if the file already exists
    try:
        df = pd.read_csv(path / filename)
    except FileNotFoundError:
        # If the file doesn't exist, create a new DataFrame
        df = pd.DataFrame()

    # Add new information to the DataFrame
    new_data = new_row

    # Use all fields from new_data DataFrame as columns
    if not df.empty:
        columns_to_keep = new_data.columns
        df = df[columns_to_keep]

    df = df.append(new_data, ignore_index=True)

    # Save the DataFrame to a CSV file
    df.to_csv(path / filename, index=False)
