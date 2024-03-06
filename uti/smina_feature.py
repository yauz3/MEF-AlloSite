import datetime
import math
import os
import statistics
import subprocess
import sys
import textwrap
import timeit
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from pymol.cgo import *
from scipy.optimize import minimize
import re
import glob
import shutil


"""
os.system("./smina.static -r pocket17_atm_fixed.pdb -l pocket17_pose_1.pdb --score_only --scoring")
Options are:
Lin_F7
Lin_F9  #makale bunu kullanmış, iki üç tanesini kullanabilirsin
ad4_scoring
default
dkoes_fast
dkoes_scoring
dkoes_scoring_old
vina
vinardo
"""
### eğer sonuçlar istediğin gibi çıkmazsa msms den sasa feature çekebilirisin
# RunXGB.py den bakabilirsin nasıl yaptıklarına

# hatta olmazsa 6. water_features i da kullanabilirsin
# RunXGB.py den bakabilirsin nasıl yaptıklarına

BIN_PATH="/home/yavuz/yavuz_proje/allosteric_binding_site/bin"

# iki tane output alacağız.
"""
os.system("./smina_feature -r pocket17_atm_fixed.pdb -l pocket17_pose_1.pdb --score_only --custom_scoring sf_vina.txt")
"""


def smina_feature(pocket_pose, ligand_pose):
    os.chdir(BIN_PATH)
    output_feature=str(subprocess.check_output(f"./smina_feature -r {pocket_pose} -l {ligand_pose} --score_only --custom_scoring sf_vina.txt",
                                                          shell=True))
    features=(output_feature.replace("\\n","").split(" "))
    clean_features=[]
    for i in features:
        if "'" in i:
            i=i.replace("'","")
        try:
            clean_features.append(float(i))
            b=b+1
        except:
            continue
    #print(clean_features)
    #print(len(clean_features))
    # features:
    """"m_bond_0", "m_bond_1", "m_bond_2", "m_bond_3", "m_bond_4", "m_bond_5", "polar_polar_0",
    "polar_polar_1", "polar_polar_2", "polar_polar_3", "polar_polar_4", "polar_polar_5", "polar_polar_6", "polar_nonpolar_0",
    "polar_nonpolar_1", "polar_nonpolar_2", "polar_nonpolar_3", "polar_nonpolar_4", "polar_nonpolar_5", "polar_nonpolar_6",
    "nonpolar_0", "nonpolar_1", "nonpolar_2", "nonpolar_3", "nonpolar_4", "nonpolar_5", "nonpolar_6", "h_bond_0",
    "h_bond_1", "h_bond_2", "h_bond_3", "h_bond_4", "anti_h_bond_0", "anti_h_bond_1", "anti_h_bond_2", "anti_h_bond_3",
    "anti_h_bond_4", "repulsion", "ad4_solvation_0", "ad4_solvation_1", "electrostatic_0", "electrostatic_1", "num_tors_add_0",
    "num_rotors_add_1", "num_heavy_atoms", "num_hydrophobic_atoms", "ligand_max_num_h_bonds", "ligand_length"""
    return clean_features

def smina_static(pocket_pose, ligand_pose):
    os.chdir(BIN_PATH)
    # lin_f9
    static_lin_f9=str(subprocess.check_output(f"./smina.static -r {pocket_pose} -l {ligand_pose} --score_only --scoring Lin_F9",
                                                          shell=True))
    lin_f9_out=static_lin_f9.split(r"\n##  ")[1].split(r"\n")[0].split(" ") # 9 feature
    # vina
    vina=str(subprocess.check_output(f"./smina.static -r {pocket_pose} -l {ligand_pose} --score_only --scoring vina",
                                                          shell=True))
    vina_out=vina.split(r"\n##  ")[1].split(r"\n")[0].split(" ") # 7 feature
    # ad4_scoring
    ad4_scoring=str(subprocess.check_output(f"./smina.static -r {pocket_pose} -l {ligand_pose} --score_only --scoring ad4_scoring",
                                                          shell=True))
    ad4_scoring_out=ad4_scoring.split(r"\n##  ")[1].split(r"\n")[0].split(" ") # 5 features
    # dkoes_scoring
    dkoes_scoring=str(subprocess.check_output(f"./smina.static -r {pocket_pose} -l {ligand_pose} --score_only --scoring dkoes_scoring",
                                                          shell=True))
    dkoes_scoring_out=dkoes_scoring.split(r"\n##  ")[1].split(r"\n")[0].split(" ") # 5 features
    # vinardo it gives error
    """vinardo=str(subprocess.check_output(f"./smina.static -r {pocket_pose} -l {ligand_pose} --score_only --scoring vinardo",
                                                          shell=True))
    vinardo_out=vinardo.split(r"\n##  ")[1].split(r"\n")[0].split(" ") # 5 features"""
    """print(lin_f9_out)
    print(vina_out)
    print(dkoes_scoring_out)
    print(vinardo_out)"""
    return lin_f9_out,vina_out,ad4_scoring_out,dkoes_scoring_out


if __name__ == "__main__":
    """pocket_pose="/home/yavuz/yavuz_proje/allosteric_binding_site/pocket17_atm_fixed.pdb"
    ligand_pose="/home/yavuz/yavuz_proje/allosteric_binding_site/pocket17_pose_1.pdb"""
    smina_feature(
        pocket_pose=pocket_pose,
        ligand_pose=ligand_pose)
    smina_static(
        pocket_pose=pocket_pose,
        ligand_pose=ligand_pose)



