#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/02/2021
# Author: Sadettin Y. Ugurlu & David McDonald
import time
import os
from uti import prepare_biopython
from uti import prepare_pybiomed
import pandas as pd
from uti import prepare_mathfeature
# -Training proteins:
training=['1AO0', '1RX2', '1DD7', '3IAD', '2CLH', '3AO1', '3IYD', '1V4S', '1I7S', '2YC3', '3MK6', '3IDB', '1Z8D', '3I0R', '2D5Z', '3F3V', '3EPS', '1EFA', '3MKS', '2EWN', '1XJE', '2BU2', '1COZ', '1HAK', '3GR4', '3RZ3', '1EGY', '2ZMF', '1PJ3', '3PTZ', '2XO8', '1SHJ', '1DB1', '1CE8', '1S9J', '1QTI', '2Q5O', '2OI2', '1ESM', '2POC', '2X1L', '1XTU', '2BND', '2I80', '3GCP', '2AL4', '1X88', '3O2M', '3CQD', '3FIG', '3HO8', '1LTH', '1FAP', '3HV8', '3GVU', '3PJG', '3H30', '1T49', '1RD4', '2V92', '2C2B', '3FZY', '3NJQ', '3UO9', '1W96', '2GS7', '3IJG', '1ZDS', '3F6G', '2PUC', '2R1R', '2VGI', '1KP8', '3OS8', '1W25', '3PEE', '3QEL', '1LDN', '1XLS', '1PFK', '3IRH', '1FTA', '2QF7', '3BEO', '3ZLK', '4AVB', '1QW7', '4B9Q', '1TUG', '1PEQ']

# -Filtered protein based on TM-scores:
# There is no protein pairs having higher than 0.5 TM-scores accross bencmarks and training set
test_1=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_2=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_3=['1J07', '1I72', '1YP2', '3BRK', '1ECB', '4CFH', '4EAG', '3KH5', '3L76', '1WQW', '4OP0', '1UWH', '1CKK', '3J41', '3I54', '2OZ6', '2FSZ', '1CSM', '4UUU', '4DQW', '1KMP', '1NV7', '2VD4', '1NE7', '4LZ5', '4U5B', '4PKN', '3E2A', '1VEA', '3TUV', '1L5G', '3BLW', '4MQT', '1O0S', '1S9I', '3D2P', '1FCJ', '1TBF', '2K31', '2PA3', '2H06', '2QMX', '3OF1', '2VK1', '1A3W', '4IP7', '1XMS', '3CMU', '1HK8', '3RSR', '2ONB', '2NW8', '1I6K', '4NES', '2JJX', '3NWY', '1PZO', '4OR2', '4RQZ', '1Q3E', '2PTM', '2VVT', '4Q0A', '4B6E', '3PMA', '3HWS', '1UM8', '2BTY', '1AZX', '1XXA', '3KJN', '1MC0', '2C18', '4LEG', '4TPW', '4DLR', '4QSK', '3UVV', '3HL8', '3LMH', '1OJ9', '3RHW', '3P4W', '2Q6H', '4JKT', '2QXL', '3QKU', '3AUX', '3AV0', '3THO', '4GQQ', '4M0Y', '1T4G', '2Y39', '3DBA', '3K8S', '3ATH', '4H39', '1BM7', '3FUD', '3JPY', '3F1O', '3ZM9', '1BJ4', '4TPT', '3CEV', '3ZFZ', '4FXY', '4M1P', '4NIL', '4OHF', '4CLL', '4PPU', '4Q9M', '3QAK', '4PHU', '3UT3', '2O8B', '2RDE', '4BXC', '4QPL', '3PNA']

current_dir = os.path.dirname(os.path.abspath(__file__))

mathfeature_path=f"{current_dir}/bin/MathFeature/methods"

def prepare_all_features(protein_list,input_file_name,fpocket_filename):
    # prepare biopython related features
    biopython_data=prepare_biopython.biopython_descriptor(protein_list=protein_list,
                         fixed_pocket_path=f"{current_dir}/data/fix_pockets/{input_file_name}",
                         )

    ####################################################################################################################
    # prepare pybiomed related features
    qso_data = prepare_pybiomed.get_pybiomed_features_QSOSW_QSOgrant(
        fixed_pocket_path=f"{current_dir}/data/fix_pockets/{input_file_name}",
        protein_list=protein_list)
    cdt_data = prepare_pybiomed.get_cdt(fixed_pocket_path=f"{current_dir}/data/fix_pockets/{input_file_name}",
                                        protein_list=protein_list)
    aacomposition_data = prepare_pybiomed.get_aacomposition(
        fixed_pocket_path=f"{current_dir}/data/fix_pockets/{input_file_name}",
        protein_list=protein_list)
    conjoin_data = prepare_pybiomed.get_conjointtriad(fixed_pocket_path=f"{current_dir}/data/fix_pockets/{input_file_name}",
                                                      protein_list=protein_list)
    ####################################################################################################################
    os.chdir(current_dir)
    math_feature_fast = f"{input_file_name}.fasta"
    if os.path.exists(math_feature_fast):
        os.remove(math_feature_fast)
    _math_feature_csv = f"{input_file_name}_math_feature.csv"
    if os.path.exists(_math_feature_csv):
        os.remove(_math_feature_csv)


    prepare_mathfeature.math_feature(protein_list=protein_list,
             input_path=f"{current_dir}/data/fix_pockets/{input_file_name}",
             output_filename=input_file_name)
    time.sleep(10)
    os.chdir(current_dir)

    os.system(
        f"python3 {mathfeature_path}/Kgap.py -i {input_file_name}.fasta -l protein -k 1 -bef 1 -aft 1 -seq 3 -o {input_file_name}_math_feature.csv")
    math_feature_data =pd.read_csv(f"{input_file_name}_math_feature.csv")
    math_feature_data.drop(columns=["label"], inplace=True) # label here is PROTEIN (in up os.system code: -l protein)
    math_feature_data.rename(columns={"nameseq": "protein_pocket"}, inplace=True)

    ####################################################################################################################
    # Load dataframes"""
    os.chdir(current_dir)
    fpocket_data = pd.read_csv(fpocket_filename)

    # Merge the dataframes sequentially
    frames = [biopython_data,
              qso_data,
              cdt_data,
              aacomposition_data,
              conjoin_data,
              math_feature_data,
              ]

    # Merge dataframes sequentially, accumulating the result
    merged_df = fpocket_data  # Start with the first dataframe
    for df in frames:
        merged_df = pd.merge(merged_df, df, on="protein_pocket", how="outer")


    # Print the first few rows of merged dataframe to check if data is present
    print("First few rows of merged dataframe:")
    print(merged_df.head())
    print(len(merged_df))
    # Save merged dataframe to CSV
    merged_df.to_csv(f"{input_file_name}_ready.csv", index=False)

"""prepare_all_features(protein_list=training,
                     input_file_name="training",
                     fpocket_filename="training.csv")"""

prepare_all_features(protein_list=test_1,
                     input_file_name="test_1",
                     fpocket_filename="test_1.csv")

"""prepare_all_features(protein_list=test_2,
                     input_file_name="test_2",
                     fpocket_filename="test_2.csv")

prepare_all_features(protein_list=test_3,
                     input_file_name="test_3",
                     fpocket_filename="test_3.csv")"""