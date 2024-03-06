#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/02/2021
# Author: Sadettin Y. Ugurlu & David McDonald
import requests
import os
from uti import protein_uti
from uti import fix_fpocket_pockets_working
import glob
import glob
import os
import pickle
import warnings
from uti.atomCount import atomCount
from uti.extract_fpocket_feature import parse_pocket_values,extract_pocket_number
import csv

# -Training proteins:
training=['1AO0', '1RX2', '1DD7', '3IAD', '2CLH', '3AO1', '3IYD', '1V4S', '1I7S', '2YC3', '3MK6', '3IDB', '1Z8D', '3I0R', '2D5Z', '3F3V', '3EPS', '1EFA', '3MKS', '2EWN', '1XJE', '2BU2', '1COZ', '1HAK', '3GR4', '3RZ3', '1EGY', '2ZMF', '1PJ3', '3PTZ', '2XO8', '1SHJ', '1DB1', '1CE8', '1S9J', '1QTI', '2Q5O', '2OI2', '1ESM', '2POC', '2X1L', '1XTU', '2BND', '2I80', '3GCP', '2AL4', '1X88', '3O2M', '3CQD', '3FIG', '3HO8', '1LTH', '1FAP', '3HV8', '3GVU', '3PJG', '3H30', '1T49', '1RD4', '2V92', '2C2B', '3FZY', '3NJQ', '3UO9', '1W96', '2GS7', '3IJG', '1ZDS', '3F6G', '2PUC', '2R1R', '2VGI', '1KP8', '3OS8', '1W25', '3PEE', '3QEL', '1LDN', '1XLS', '1PFK', '3IRH', '1FTA', '2QF7', '3BEO', '3ZLK', '4AVB', '1QW7', '4B9Q', '1TUG', '1PEQ']

# -Filtered protein based on TM-scores:
# There is no protein pairs having higher than 0.5 TM-scores accross bencmarks and training set
test_1=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_2=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_3=['1J07', '1I72', '1YP2', '3BRK', '1ECB', '4CFH', '4EAG', '3KH5', '3L76', '1WQW', '4OP0', '1UWH', '1CKK', '3J41', '3I54', '2OZ6', '2FSZ', '1CSM', '4UUU', '4DQW', '1KMP', '1NV7', '2VD4', '1NE7', '4LZ5', '4U5B', '4PKN', '3E2A', '1VEA', '3TUV', '1L5G', '3BLW', '4MQT', '1O0S', '1S9I', '3D2P', '1FCJ', '1TBF', '2K31', '2PA3', '2H06', '2QMX', '3OF1', '2VK1', '1A3W', '4IP7', '1XMS', '3CMU', '1HK8', '3RSR', '2ONB', '2NW8', '1I6K', '4NES', '2JJX', '3NWY', '1PZO', '4OR2', '4RQZ', '1Q3E', '2PTM', '2VVT', '4Q0A', '4B6E', '3PMA', '3HWS', '1UM8', '2BTY', '1AZX', '1XXA', '3KJN', '1MC0', '2C18', '4LEG', '4TPW', '4DLR', '4QSK', '3UVV', '3HL8', '3LMH', '1OJ9', '3RHW', '3P4W', '2Q6H', '4JKT', '2QXL', '3QKU', '3AUX', '3AV0', '3THO', '4GQQ', '4M0Y', '1T4G', '2Y39', '3DBA', '3K8S', '3ATH', '4H39', '1BM7', '3FUD', '3JPY', '3F1O', '3ZM9', '1BJ4', '4TPT', '3CEV', '3ZFZ', '4FXY', '4M1P', '4NIL', '4OHF', '4CLL', '4PPU', '4Q9M', '3QAK', '4PHU', '3UT3', '2O8B', '2RDE', '4BXC', '4QPL', '3PNA']

current_dir = os.path.dirname(os.path.abspath(__file__))



def prepare_fpocket_features(protein_list, input_path, output_file_name):
    # label pockets with 1 or 0
    # size = num of proteins
    # size of each element = num of pockets in each protein
    labels = []  # List[List[int]]

    # FPocket features
    # size = num of proteins
    # size of each element = num of pockets in each protein * 19 features
    features = []  # List[List[List[float]]]

    # at least need 9 heavy atoms to be labeled as 1
    # ref: mean value of heavy atoms in 20 amino acids = 9
    n_atoms = 9

    # the features will be storaged in the output dictionary
    output_dictionary = {}
    for pdb in protein_list:
        # take ground true values from ADS
        info = protein_uti.get_info_from_ads(pdb)[-1]
        # get the directory containing pocket files
        direction_pockets = f"{input_path}/{pdb}_cleaned_out/pockets/"
        # get the all pockets found by Fpocket
        pocket_files = glob.glob(direction_pockets + "*.pdb")

        # Report if no pockets
        if len(pocket_files) < 1:
            warnings.warn("no pocket detected in %s" % pdb)
            continue
        # while extracting fpocket features, we will label the pockets as well.
        for pocket in pocket_files:
            pocket_number = extract_pocket_number(pocket)
            counts = (atomCount(pocket, info))

            # Assign labels
            if counts >= n_atoms:
                cur_label = 1
            else:
                cur_label = 0

            features_collections = f"{input_path}/{pdb}_cleaned_out/{pdb}_cleaned_info.txt"
            cur_feature = parse_pocket_values(features_collections)

            pocket_features = cur_feature.get(str(pocket_number))
            output_dictionary[pdb, pocket_number] = pocket_features, cur_label
    #parent_directory = os.path.dirname(current_path)
    os.chdir(current_dir)
    with open(f'{output_file_name}.csv', 'w', newline='') as csvfile:
        fieldnames = ['protein_pocket', 'Score', 'Druggability Score', 'Number of Alpha Spheres', 'Total SASA', 'Polar SASA',
                      'Apolar SASA', 'Volume', 'Mean local hydrophobic density', 'Mean alpha sphere radius',
                      'Mean alp. sph. solvent access', 'Apolar alpha sphere proportion', 'Hydrophobicity score',
                      'Volume score', 'Polarity score', 'Charge score', 'Proportion of polar atoms',
                      'Alpha sphere density', 'Cent. of mass - Alpha Sphere max dist', 'Flexibility', 'Label']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for key, value in output_dictionary.items():
            pdb_key = f"{key[0]}_{key[1]}"
            row_data = {'protein_pocket': pdb_key, **value[0], 'Label': value[1]}
            writer.writerow(row_data)


# prepare training set
prepare_fpocket_features(protein_list=training,
                input_path=f"{current_dir}/data/training",
                output_file_name="training")

# prepare test_1
prepare_fpocket_features(protein_list=test_1,
                input_path=f"{current_dir}/data/test_1",
                output_file_name="test_1")
# prepare test_2
prepare_fpocket_features(protein_list=test_1,
                input_path=f"{current_dir}/data/test_2",
                output_file_name="test_2")

# prepare test_3
prepare_fpocket_features(protein_list=test_3,
                input_path=f"{current_dir}/data/test_3",
                output_file_name="test_3")