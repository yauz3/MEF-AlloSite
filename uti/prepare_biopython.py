#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/02/2021
# Author: Sadettin Y. Ugurlu & David McDonald

import requests
import os
from uti import fix_fpocket_pockets_working
import glob
import glob
import os
import pickle
import warnings
from uti.atomCount import atomCount
from uti.extract_fpocket_feature import parse_pocket_values,extract_pocket_number
import csv
import os
import sys
import glob
import os
from uti import psvina
from uti import vina_uti
from uti import smina_feature
import pandas as pd
import numpy as np
import csv
import subprocess
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB import PDBParser

# it will off the warning
warnings.filterwarnings('ignore')

# -Training proteins:
training=['1AO0', '1RX2', '1DD7', '3IAD', '2CLH', '3AO1', '3IYD', '1V4S', '1I7S', '2YC3', '3MK6', '3IDB', '1Z8D', '3I0R', '2D5Z', '3F3V', '3EPS', '1EFA', '3MKS', '2EWN', '1XJE', '2BU2', '1COZ', '1HAK', '3GR4', '3RZ3', '1EGY', '2ZMF', '1PJ3', '3PTZ', '2XO8', '1SHJ', '1DB1', '1CE8', '1S9J', '1QTI', '2Q5O', '2OI2', '1ESM', '2POC', '2X1L', '1XTU', '2BND', '2I80', '3GCP', '2AL4', '1X88', '3O2M', '3CQD', '3FIG', '3HO8', '1LTH', '1FAP', '3HV8', '3GVU', '3PJG', '3H30', '1T49', '1RD4', '2V92', '2C2B', '3FZY', '3NJQ', '3UO9', '1W96', '2GS7', '3IJG', '1ZDS', '3F6G', '2PUC', '2R1R', '2VGI', '1KP8', '3OS8', '1W25', '3PEE', '3QEL', '1LDN', '1XLS', '1PFK', '3IRH', '1FTA', '2QF7', '3BEO', '3ZLK', '4AVB', '1QW7', '4B9Q', '1TUG', '1PEQ']

# -Filtered protein based on TM-scores:
# There is no protein pairs having higher than 0.5 TM-scores accross bencmarks and training set
test_1=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_2=['3QOP', '11BG', '2XJC', '1H9G', '3QH0', '4BO2', '1UXV', '4I1R', '4AW0', '2Q8M', '1NJJ', '3F9N', '2HIM', '1DKU', '1W0F', '1OF6', '3GCD', '2I7N', '2BE9', '3N1V', '3LAJ', '2HVW', '4ETZ', '4HSG', '3O96', '4OO9', '3HNC', '4NBN', '1JLR', '1FX2', '3E5U', '4EBW', '3HO6', '1FIY', '2JFN', '3PYY', '4BZB', '4MBS', '4B1F', '3KGF', '4C7B', '2VPR', '1HKB', '2A69', '4BQH', '2Y0P', '3LW0', '3LU6', '3KF0', '3PXF', '1M8P', '2RD5', '4JAF', '4EO6', '2YLO', '3KCG']
test_3=['1J07', '1I72', '1YP2', '3BRK', '1ECB', '4CFH', '4EAG', '3KH5', '3L76', '1WQW', '4OP0', '1UWH', '1CKK', '3J41', '3I54', '2OZ6', '2FSZ', '1CSM', '4UUU', '4DQW', '1KMP', '1NV7', '2VD4', '1NE7', '4LZ5', '4U5B', '4PKN', '3E2A', '1VEA', '3TUV', '1L5G', '3BLW', '4MQT', '1O0S', '1S9I', '3D2P', '1FCJ', '1TBF', '2K31', '2PA3', '2H06', '2QMX', '3OF1', '2VK1', '1A3W', '4IP7', '1XMS', '3CMU', '1HK8', '3RSR', '2ONB', '2NW8', '1I6K', '4NES', '2JJX', '3NWY', '1PZO', '4OR2', '4RQZ', '1Q3E', '2PTM', '2VVT', '4Q0A', '4B6E', '3PMA', '3HWS', '1UM8', '2BTY', '1AZX', '1XXA', '3KJN', '1MC0', '2C18', '4LEG', '4TPW', '4DLR', '4QSK', '3UVV', '3HL8', '3LMH', '1OJ9', '3RHW', '3P4W', '2Q6H', '4JKT', '2QXL', '3QKU', '3AUX', '3AV0', '3THO', '4GQQ', '4M0Y', '1T4G', '2Y39', '3DBA', '3K8S', '3ATH', '4H39', '1BM7', '3FUD', '3JPY', '3F1O', '3ZM9', '1BJ4', '4TPT', '3CEV', '3ZFZ', '4FXY', '4M1P', '4NIL', '4OHF', '4CLL', '4PPU', '4Q9M', '3QAK', '4PHU', '3UT3', '2O8B', '2RDE', '4BXC', '4QPL', '3PNA']

current_dir = os.path.dirname(os.path.abspath(__file__))



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
    #print(sequence)
    return sequence


dna_rna_pockets={}
d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}
def biopython_descriptor(protein_list,fixed_pocket_path,):
    """
    The code will take protein in a list to produce sequence, then produce biopython related features. The features are
    mainly proportion of amino acids.
    :param protein_list: proteins in a protein list
    :param fixed_pocket_path: to find where input protein are
    :param output_path: the place where the code will write CSV output file
    :param output_name: the name of output CSV file
    :return: It will save CSV output file.
    """
    # the features will be store in the dictionary
    output_dictionary={}
    for protein in protein_list:
        # go to dictionary where pockets have been stored
        os.chdir(f"{fixed_pocket_path}/{protein}_fixed")
        # get all .pdb files into a list
        pocket_list=glob.glob("*.pdb")
        for pocket in pocket_list:
            # Fpocket can find DNA, RNA without any protein fragment as a pocket, or pocket can be too much small to produce
            # seqeunce based features, which breaks the code.
            try:
                # take sequence for pocket.pdb file
                sequence = sequnce_from_PDBParser(pocket)
                X = ProteinAnalysis(sequence)
                percent_list=[]
                # d3to1 is a dictionary for residues vs their label
                for key, value in d3to1.items():
                        percent=X.get_amino_acids_percent()[value] # percentage of residue
                        percent_list.append(percent)
                # calculate molecular weight for the pocket
                mw = X.molecular_weight()
                # calculate aromaticy for the pocket
                aromaticty = X.aromaticity()
                # calculate instability for the pocket
                instability= X.instability_index()
                # calculate isoelectiric point for the pocket
                isoelectric_point=X.isoelectric_point()
                # get secondary structure features for the pocket
                # [helix, turn, sheet]
                sec_structure=X.secondary_structure_fraction()
                helix=sec_structure[0] # with reduced cysteines
                turn=sec_structure[1] # with disulfid bridges
                sheet=sec_structure[2]
                # calculate gravy for the pocket
                gravy=X.gravy()
                # calculate charge at pH 7
                charge_at_pH=X.charge_at_pH(7)
                # calculate molar extinction coefficient for the protein
                molar_extinction_coefficient=X.molar_extinction_coefficient()
                molar_extinction_coefficient=str(molar_extinction_coefficient).replace(")","").replace("(","").split(",")
                # each varible have been stored in the dictionary
                pocket_number=extract_pocket_number(pocket)
                output_dictionary[f"{protein}_{pocket_number}"]=(
                                                   percent_list[0],percent_list[1],percent_list[2],percent_list[3],
                                                   percent_list[4],percent_list[5],percent_list[6],percent_list[7],
                                                   percent_list[8],percent_list[9],percent_list[10],percent_list[11],
                                                   percent_list[12],percent_list[13],percent_list[14],percent_list[15],
                                                   percent_list[16],percent_list[17],percent_list[18],percent_list[19],
                                                   mw,aromaticty,instability,isoelectric_point,
                                                   helix,turn,sheet,gravy,charge_at_pH,molar_extinction_coefficient[0],
                                                   molar_extinction_coefficient[1],
                                                  )
            except:
                # over than 98% Fpocket find a protein as a pocket or find large enough to produce biopython features.
                # if not, they will return None, which AutoGluon can handle during training a model or prediction.
                output_dictionary[f"{protein}_{pocket_number}"] = (None,None,None,None,None,None,None,None,None,None,None,
                                                                 None,None,None,None,None,None,None,None,None,None,None,
                                                                 None,None,None,None,None,None,None,None,None,)
    columns_names_of_csv=['biopython.csvCYS', 'biopython.csvASP', 'biopython.csvSER', 'biopython.csvGLN', 'biopython.csvLYS',
                          'biopython.csvILE', 'biopython.csvPRO', 'biopython.csvTHR', 'biopython.csvPHE', 'biopython.csvASN',
                          'biopython.csvGLY', 'biopython.csvHIS', 'biopython.csvLEU', 'biopython.csvARG', 'biopython.csvTRP',
                          'biopython.csvALA', 'biopython.csvVAL', 'biopython.csvGLU', 'biopython.csvTYR', 'biopython.csvMET',
                          'biopython.csvmw', 'biopython.csvaromaticty', 'biopython.csvinstability', 'biopython.csvisoelectric_point',
                          'biopython.csvhelix', 'biopython.csvturn', 'biopython.csvsheet', 'biopython.csvgravy',
                          'biopython.csvcharge_at_pH', 'biopython.csvmolar_extinction_coefficient_0', 'biopython.csvmolar_extinction_coefficient_1']

    out_data=pd.DataFrame.from_dict(output_dictionary, orient='index',columns=columns_names_of_csv)
    out_data.index.name = 'protein_pocket'
    #os.chdir(output_path)
    #out_data.to_csv(f"{output_name}.csv")
    return out_data
# exampmle usage
"""biopython_descriptor(protein_list=training,
                     fixed_pocket_path="/home/yavuz/yavuz_proje/allosteric_feature_selected/data/training",
                     output_path="/home/yavuz/yavuz_proje/allosteric_feature_selected",
                     output_name="training_biopython")"""

