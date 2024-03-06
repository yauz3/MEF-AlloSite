#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/02/2021
# Author: Sadettin Y. Ugurlu & David McDonald

import requests
import os

# -Training proteins:
training=['1AO0', '1RX2', '1DD7', '3IAD', '2CLH', '3AO1', '3IYD', '1V4S', '1I7S', '2YC3', '3MK6', '3IDB', '1Z8D', '3I0R', '2D5Z', '3F3V', '3EPS', '1EFA', '3MKS', '2EWN', '1XJE', '2BU2', '1COZ', '1HAK', '3GR4', '3RZ3', '1EGY', '2ZMF', '1PJ3', '3PTZ', '2XO8', '1SHJ', '1DB1', '1CE8', '1S9J', '1QTI', '2Q5O', '2OI2', '1ESM', '2POC', '2X1L', '1XTU', '2BND', '2I80', '3GCP', '2AL4', '1X88', '3O2M', '3CQD', '3FIG', '3HO8', '1LTH', '1FAP', '3HV8', '3GVU', '3PJG', '3H30', '1T49', '1RD4', '2V92', '2C2B', '3FZY', '3NJQ', '3UO9', '1W96', '2GS7', '3IJG', '1ZDS', '3F6G', '2PUC', '2R1R', '2VGI', '1KP8', '3OS8', '1W25', '3PEE', '3QEL', '1LDN', '1XLS', '1PFK', '3IRH', '1FTA', '2QF7', '3BEO', '3ZLK', '4AVB', '1QW7', '4B9Q', '1TUG', '1PEQ']

# -Filtered protein based on TM-scores:
# There is no protein pairs having higher than 0.5 TM-scores accross bencmarks and training set.
# Pairs with a TM-Score greater than 0.5 share the same fold, perhaps leading to a similarity in their three-dimensional structure.
test_1=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_2=['1J07', '1I72', '1YP2', '3BRK', '1ECB', '4CFH', '4EAG', '3KH5', '3L76', '1WQW', '4OP0', '1UWH', '1CKK', '3J41', '3I54', '2OZ6', '2FSZ', '1CSM', '4UUU', '4DQW', '1KMP', '1NV7', '2VD4', '1NE7', '4LZ5', '4U5B', '4PKN', '3E2A', '1VEA', '3TUV', '1L5G', '3BLW', '4MQT', '1O0S', '1S9I', '3D2P', '1FCJ', '1TBF', '2K31', '2PA3', '2H06', '2QMX', '3OF1', '2VK1', '1A3W', '4IP7', '1XMS', '3CMU', '1HK8', '3RSR', '2ONB', '2NW8', '1I6K', '4NES', '2JJX', '3NWY', '1PZO', '4OR2', '4RQZ', '1Q3E', '2PTM', '2VVT', '4Q0A', '4B6E', '3PMA', '3HWS', '1UM8', '2BTY', '1AZX', '1XXA', '3KJN', '1MC0', '2C18', '4LEG', '4TPW', '4DLR', '4QSK', '3UVV', '3HL8', '3LMH', '1OJ9', '3RHW', '3P4W', '2Q6H', '4JKT', '2QXL', '3QKU', '3AUX', '3AV0', '3THO', '4GQQ', '4M0Y', '1T4G', '2Y39', '3DBA', '3K8S', '3ATH', '4H39', '1BM7', '3FUD', '3JPY', '3F1O', '3ZM9', '1BJ4', '4TPT', '3CEV', '3ZFZ', '4FXY', '4M1P', '4NIL', '4OHF', '4CLL', '4PPU', '4Q9M', '3QAK', '4PHU', '3UT3', '2O8B', '2RDE', '4BXC', '4QPL', '3PNA']

# - All proteins used in the study to download all of them. They WILL NOT USE TOGETHER! That's why they are separate lists.
all_proteins=training+test_1+test_2


# - Get current_dir
current_dir = os.getcwd()

# Define the directory path by joining current_dir, "data", and "pdbs"
data_pdb_dir = os.path.join(current_dir, "data", "pdbs")

# - Check and Clean Data
# Check if the directory specified by 'data_pdb_dir' exists
if not os.path.exists(data_pdb_dir):
    # If it doesn't exist, create the directory
    os.makedirs(data_pdb_dir)
else:
    # If it exists, iterate over the files in the directory
    for filename in os.listdir(data_pdb_dir):
        # Get the full path of the file
        file_path = os.path.join(data_pdb_dir, filename)
        try:
            # Check if the file exists and is a regular file
            if os.path.isfile(file_path):
                # If so, unlink (delete) the file
                os.unlink(file_path)
        except Exception as e:
            # If an error occurs during file deletion, print the error message
            print(f"Error: {e}")


# - Download pdb
# Loop through each PDB ID in the list 'all_proteins'
for pdbID in all_proteins:
    # Construct the URL for downloading the PDB file using the RCSB URL and the current PDB ID
    rcsb = "https://files.rcsb.org/download/"
    url = rcsb + pdbID + ".pdb"
    # Send a GET request to the URL to download the PDB file
    r = requests.get(url)
    # Save the content of the downloaded PDB file to a local file with the same name as the PDB ID
    with open(os.path.join(data_pdb_dir, f"{pdbID}.pdb"), "wb") as f:
        f.write(r.content)
