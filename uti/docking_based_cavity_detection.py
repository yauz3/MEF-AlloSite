#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/02/2021
# Author: Sadettin Y. Ugurlu & David McDonald


import glob
import os
import shutil
from script import protein_similarity_search
from script import psvina
import pickle
import sys
from Bio import SeqIO
from Bio.PDB import PDBParser
import subprocess
from pymol import cmd




ALLOSTERIC_PATH="/home/yavuz/yavuz_proje/allosteric_binding_site"
PDBind_PATH="/home/yavuz/yavuz_proje/allosteric_binding_site/database/PDBbind_v2020_other_PL/v2020-other-PL"
Fpocket_out="/home/yavuz/yavuz_proje/allosteric_binding_site/passer_200/PASSer2.0-main/data/pockets"
ENLARGED_POCKETS_PATH="/home/yavuz/yavuz_proje/allosteric_binding_site/passer_200/PASSer2.0-main/data/two_residue"

# take input
def docking_based_cavity_detection(input_pdbid):
    # go to Fpocket outputs, since passer selecting chains
    fpocet_out_path=f"{Fpocket_out}/{input_pdbid}_out"
    os.chdir(fpocet_out_path)
    check_clean_pdb=fpocet_out_path+f"{input_pdbid}_clean.pdb"
    if os.path.isfile(check_clean_pdb) is not True:
        cmd.reinitialize()
        cmd.load(f"{input_pdbid}_out.pdb")
        cmd.remove("het")
        cmd.save(f"{input_pdbid}_clean.pdb")
        cmd.remove("all")

    # get input sequnce:
    input_sequence=protein_similarity_search.sequnce_from_PDBParser(f"{input_pdbid}_clean.pdb")
    #print(input_sequence)
    # read all seqeunce for PDBind
    os.chdir(ALLOSTERIC_PATH)
    a_file = open("../allosteric_binding_site/bin/sequence_dictionary.pkl", "rb")
    sequnce_dictionary = pickle.load(a_file)
    #print("len",len(sequnce_dictionary)) # 14126
    # we are going to check whether data bases has the protein or not
    if str(input_pdbid).lower() in list(sequnce_dictionary):
        selected_pdbid=input_pdbid.lower()
    # search the most similar protein in PDBind database
    else:
        similarty_output={}
        for key,value in sequnce_dictionary.items():
            # calculate similarty for each pair
            similarty_score=protein_similarity_search.sequence_similarty(input_sequence,value)
            #print(similarty_score)
            similarty_output[key]=similarty_score
        similarty_output=dict(sorted(similarty_output.items(), key=lambda item: item[1]))
        #print(similarty_output)
        selected_pdbid=list(similarty_output)[-1]
    # pocket da var burada RMSD hesaplayabilirsin
    try:
        os.mkdir(f"/home/yavuz/yavuz_proje/allosteric_binding_site/docking_200/{input_pdbid}")
    except:
        print("there is already file")
    shutil.copy(f"{PDBind_PATH}/{selected_pdbid}/{selected_pdbid}_ligand.mol2",
                f"/home/yavuz/yavuz_proje/allosteric_binding_site/docking_200/{input_pdbid}/{input_pdbid}_ligand.mol2")
    shutil.copy(f"{PDBind_PATH}/{selected_pdbid}/{selected_pdbid}_pocket.pdb",
                f"/home/yavuz/yavuz_proje/allosteric_binding_site/docking_200/{input_pdbid}/{input_pdbid}_pocket_from_library.pdb")
    os.chdir(f"{PDBind_PATH}/{selected_pdbid}")

    # we will use the pockets structure instead of whole protein
    os.chdir(f"{ENLARGED_POCKETS_PATH}/{input_pdbid}")
    pocket_list=glob.glob("*.pdb")

    for pocket in pocket_list:
        shutil.copy(f"{ENLARGED_POCKETS_PATH}/{input_pdbid}/{pocket}",
                    f"/home/yavuz/yavuz_proje/allosteric_binding_site/docking_200/{input_pdbid}/{pocket}")


if __name__ == '__main__':
    """input_pdbid="1AO0"""
    docking_based_cavity_detection(
        input_pdbid=input_pdbid)
