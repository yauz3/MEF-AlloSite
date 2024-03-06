#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/02/2021
# Author: Sadettin Y. Ugurlu & David McDonald

import os
from uti import prepare_biopython
import glob
import os
import pandas as pd
from PyBioMed.PyInteraction.PyInteraction import CalculateInteraction2
from PyBioMed.PyInteraction.PyInteraction import CalculateInteraction1
from PyBioMed.PyInteraction import PyInteraction
from PyBioMed.PyDNA.PyDNApsenac import GetPseDNC
from PyBioMed.PyProtein import CTD
import pandas as pd
import numpy as np
import csv
from uti import protein_similarity_search
from PyBioMed.PyProtein import AAComposition,Autocorrelation,ConjointTriad,PseudoAAC,PyProteinAAComposition,QuasiSequenceOrder
from PyBioMed import Pyprotein
from uti import fpocket_to_sequence
from Bio.PDB import PDBParser
from uti.extract_fpocket_feature import parse_pocket_values,extract_pocket_number


# -Training proteins:
training=['1AO0', '1RX2', '1DD7', '3IAD', '2CLH', '3AO1', '3IYD', '1V4S', '1I7S', '2YC3', '3MK6', '3IDB', '1Z8D', '3I0R', '2D5Z', '3F3V', '3EPS', '1EFA', '3MKS', '2EWN', '1XJE', '2BU2', '1COZ', '1HAK', '3GR4', '3RZ3', '1EGY', '2ZMF', '1PJ3', '3PTZ', '2XO8', '1SHJ', '1DB1', '1CE8', '1S9J', '1QTI', '2Q5O', '2OI2', '1ESM', '2POC', '2X1L', '1XTU', '2BND', '2I80', '3GCP', '2AL4', '1X88', '3O2M', '3CQD', '3FIG', '3HO8', '1LTH', '1FAP', '3HV8', '3GVU', '3PJG', '3H30', '1T49', '1RD4', '2V92', '2C2B', '3FZY', '3NJQ', '3UO9', '1W96', '2GS7', '3IJG', '1ZDS', '3F6G', '2PUC', '2R1R', '2VGI', '1KP8', '3OS8', '1W25', '3PEE', '3QEL', '1LDN', '1XLS', '1PFK', '3IRH', '1FTA', '2QF7', '3BEO', '3ZLK', '4AVB', '1QW7', '4B9Q', '1TUG', '1PEQ']

# -Filtered protein based on TM-scores:
# There is no protein pairs having higher than 0.5 TM-scores accross bencmarks and training set
test_1=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_2=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_3=['1J07', '1I72', '1YP2', '3BRK', '1ECB', '4CFH', '4EAG', '3KH5', '3L76', '1WQW', '4OP0', '1UWH', '1CKK', '3J41', '3I54', '2OZ6', '2FSZ', '1CSM', '4UUU', '4DQW', '1KMP', '1NV7', '2VD4', '1NE7', '4LZ5', '4U5B', '4PKN', '3E2A', '1VEA', '3TUV', '1L5G', '3BLW', '4MQT', '1O0S', '1S9I', '3D2P', '1FCJ', '1TBF', '2K31', '2PA3', '2H06', '2QMX', '3OF1', '2VK1', '1A3W', '4IP7', '1XMS', '3CMU', '1HK8', '3RSR', '2ONB', '2NW8', '1I6K', '4NES', '2JJX', '3NWY', '1PZO', '4OR2', '4RQZ', '1Q3E', '2PTM', '2VVT', '4Q0A', '4B6E', '3PMA', '3HWS', '1UM8', '2BTY', '1AZX', '1XXA', '3KJN', '1MC0', '2C18', '4LEG', '4TPW', '4DLR', '4QSK', '3UVV', '3HL8', '3LMH', '1OJ9', '3RHW', '3P4W', '2Q6H', '4JKT', '2QXL', '3QKU', '3AUX', '3AV0', '3THO', '4GQQ', '4M0Y', '1T4G', '2Y39', '3DBA', '3K8S', '3ATH', '4H39', '1BM7', '3FUD', '3JPY', '3F1O', '3ZM9', '1BJ4', '4TPT', '3CEV', '3ZFZ', '4FXY', '4M1P', '4NIL', '4OHF', '4CLL', '4PPU', '4Q9M', '3QAK', '4PHU', '3UT3', '2O8B', '2RDE', '4BXC', '4QPL', '3PNA']

# Get the absolute path of the current directory
current_dir = os.path.dirname(os.path.abspath(__file__))

# Get the absolute path of the parent directory
parent_dir = os.path.dirname(current_dir)


def sequnce_from_PDBParser(pdbid):
    # You can use a dict to convert three letter code to one letter code
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


    # run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdbid)

    # iterate each model, chain, and residue
    # printing out the sequence for each chain

    sequence = ""
    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                try: # sometimes DNA or RNA may exist in the structure with ATOM label. To valid Amino Acid is provided above
                    seq.append(d3to1[residue.resname])
                except:
                    continue
            sequence = sequence + (''.join(seq))
    return sequence


def get_pybiomed_features_QSOSW_QSOgrant(fixed_pocket_path, protein_list):
    """
    The code needs fixed_pocket_path to find input files. Also, feature will be prepared for protein in protein list
    :param fixed_pocket_path: fixed_pocket_path after completing residues for pocket found by Fpocket
    :param protein_list: A protein list
    :return: CSV file or Dataframe
    """
    # the varibales will be stored in the dictionary
    output_dictionary={}
    for protein in protein_list:
        # find the fixed pockets by completing missing atoms on residues
        os.chdir(f"{fixed_pocket_path}/{protein}_fixed")
        # take all pocket file in the directory
        pocket_list=glob.glob("*.pdb")
        # for each pocket calculates desriptors
        for pocket in pocket_list:
            # Rarely (< 1%) sequence of pocket is too small to procude feature or Fpocket finds nucleic acid as a pocket,
            # they are not an allosteric binding site, so to overcome this kind of disturbance, we used try.
            try:
                # take sequence for pocket
                sequence_of_pocket=sequnce_from_PDBParser(f"{pocket}")
                # calculate protein/pocket descriptors
                protein_desriptors= QuasiSequenceOrder.GetQuasiSequenceOrder(sequence_of_pocket)
                # store them
                pocket_number=extract_pocket_number(pocket)
                output_dictionary[f"{protein}_{pocket_number}"]=(protein_desriptors.values())
                # take column names
                colums = list(protein_desriptors.keys())

            except:
                # if pockets are too small or pockets created by DNA or RNA, return NONE.
                # Autogluon can handle these during training and prediction
                pocket_number=extract_pocket_number(pocket)
                output_dictionary[f"{protein}_{pocket_number}"]= [np.nan] * len(colums)

    out_data=pd.DataFrame.from_dict(output_dictionary, orient='index',columns=colums)
    out_data.index.name = 'protein_pocket'
    #os.chdir(parent_dir)
    #out_data.to_csv(f"{output_name}_QuasiSequenceOrder.GetQuasiSequenceOrder.csv")
    return out_data

def get_cdt(fixed_pocket_path, protein_list):
    """
    The code needs fixed_pocket_path to find input files. Also, feature will be prepared for protein in protein list
    :param fixed_pocket_path: fixed_pocket_path after completing residues for pocket found by Fpocket
    :param protein_list: A protein list
    :return: CSV file or Dataframe
    """
    # variable will be stored in the dictionary
    output_dictionary={}
    for protein in protein_list:
        # go to the directory of the fixed protein
        os.chdir(f"{fixed_pocket_path}/{protein}_fixed")
        # take all pocket file in the directory
        pocket_list=glob.glob("*.pdb")
        # for each pocket calculates desriptors
        for pocket in pocket_list:
            # Rarely (< 1%) sequence of pocket is too small to procude feature or Fpocket finds nucleic acid as a pocket,
            # they are not an allosteric binding site, so to overcome this kind of disturbance, we used try.
            try:
                # take sequence for pocket
                sequence_of_pocket=sequnce_from_PDBParser(f"{pocket}")
                # calculate protein/pocket descriptors
                protein_desriptors= CTD.CalculateCTD(sequence_of_pocket)
                # store them
                pocket_number=extract_pocket_number(pocket)
                output_dictionary[f"{protein}_{pocket_number}"]=(protein_desriptors.values())
                # take column names
                colums = list(protein_desriptors.keys())
            except:
                # if pockets are too small or pockets created by DNA or RNA, return NONE.
                # Autogluon can handle these during training and prediction
                pocket_number = extract_pocket_number(pocket)
                output_dictionary[f"{protein}_{pocket_number}"] = [np.nan] * len(colums)

    out_data=pd.DataFrame.from_dict(output_dictionary, orient='index',columns=colums)
    out_data.index.name = 'protein_pocket'
    #os.chdir(parent_dir)
    #out_data.to_csv("QuasiSequenceOrder.GetQuasiSequenceOrder_test_1.csv")
    return out_data

def get_aacomposition(fixed_pocket_path, protein_list):
    """
    The code needs fixed_pocket_path to find input files. Also, feature will be prepared for protein in protein list
    :param fixed_pocket_path: fixed_pocket_path after completing residues for pocket found by Fpocket
    :param protein_list: A protein list
    :return: CSV file or Dataframe
    """
    # variable will be stored in the dictionary
    output_dictionary={}
    for protein in protein_list:
        # go to the directory of the fixed protein
        os.chdir(f"{fixed_pocket_path}/{protein}_fixed")
        # take all pocket file in the directory
        pocket_list=glob.glob("*.pdb")
        # for each pocket calculates desriptors
        for pocket in pocket_list:
            # Rarely (< 1%) sequence of pocket is too small to procude feature or Fpocket finds nucleic acid as a pocket,
            # they are not an allosteric binding site, so to overcome this kind of disturbance, we used try.
            try:
                # take sequence for pocket
                sequence_of_pocket=sequnce_from_PDBParser(f"{pocket}")
                # calculate protein/pocket descriptors
                protein_desriptors= AAComposition.CalculateAADipeptideComposition(sequence_of_pocket)
                # store them
                pocket_number=extract_pocket_number(pocket)
                output_dictionary[f"{protein}_{pocket_number}"]=(protein_desriptors.values())
                # take column names
                colums = list(protein_desriptors.keys())

            except:
                # if pockets are too small or pockets created by DNA or RNA, return NONE.
                # Autogluon can handle these during training and prediction
                pocket_number = extract_pocket_number(pocket)
                output_dictionary[f"{protein}_{pocket_number}"] = [np.nan] * len(colums)

    out_data=pd.DataFrame.from_dict(output_dictionary, orient='index',columns=colums)
    out_data.index.name = 'protein_pocket'
    #os.chdir(parent_dir)
    #out_data.to_csv("QuasiSequenceOrder.GetQuasiSequenceOrder_test_1_deneme.csv")
    return out_data

def get_conjointtriad(fixed_pocket_path, protein_list):
    """
    The code needs fixed_pocket_path to find input files. Also, feature will be prepared for protein in protein list
    :param fixed_pocket_path: fixed_pocket_path after completing residues for pocket found by Fpocket
    :param protein_list: A protein list
    :return: CSV file or Dataframe
    """
    # variable will be stored in the dictionary
    output_dictionary={}
    for protein in protein_list:
        # go to the directory of the fixed protein
        os.chdir(f"{fixed_pocket_path}/{protein}_fixed")
        # take all pocket file in the directory
        pocket_list=glob.glob("*.pdb")
        # for each pocket calculates desriptors
        for pocket in pocket_list:
            # Rarely (< 1%) sequence of pocket is too small to procude feature or Fpocket finds nucleic acid as a pocket,
            # they are not an allosteric binding site, so to overcome this kind of disturbance, we used try.
            try:
                # take sequence for pocket
                sequence_of_pocket=sequnce_from_PDBParser(f"{pocket}")
                # calculate protein/pocket descriptors
                protein_desriptors= ConjointTriad.CalculateConjointTriad(sequence_of_pocket)
                # store them
                pocket_number=extract_pocket_number(pocket)
                output_dictionary[f"{protein}_{pocket_number}"]=(protein_desriptors.values())
                # take column names
                colums = list(protein_desriptors.keys())

            except:
                # if pockets are too small or pockets created by DNA or RNA, return NONE.
                # Autogluon can handle these during training and prediction
                pocket_number = extract_pocket_number(pocket)
                output_dictionary[f"{protein}_{pocket_number}"] = [np.nan] * len(colums)

    out_data=pd.DataFrame.from_dict(output_dictionary, orient='index',columns=colums)
    out_data.index.name = 'protein_pocket'
    #os.chdir(parent_dir)
    #out_data.to_csv("QuasiSequenceOrder.GetQuasiSequenceOrder_test_1_deneme.csv")
    return out_data

