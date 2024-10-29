# -*- coding: utf-8 -*-
"""
Created on Fri Jul 12 09:00:20 2024

@author: bjoerk_c
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
