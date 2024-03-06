

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
import shutil

def psoVina_local_docking(receptor,
                          ligand,
                          x,
                          y,
                          z,
                          size_x,
                          size_y,
                          size_z,
                          MGLTOOLS_PATH,
                          PSO_VINA_PATH):
    prepare_receptor_path=f"{MGLTOOLS_PATH}/MGLToolsPckgs/AutoDockTools/Utilities24"
    ### Paths ####
    file_path_0 = os.path.dirname(sys.argv[0])
    print('\nfile path =', os.path.abspath(file_path_0))
    file_path = os.path.abspath(file_path_0).split("/")
    first = file_path[0]
    last = file_path[-1]
    file_path.remove(last)
    file_path.remove(first)
    output_file = "/"
    only_receptor_name = receptor.split("/")[-1].split(".")[0]
    only_ligand_name = ligand.split("/")[-1].split(".")[0]
    receptor_path=receptor.replace(f"{only_receptor_name}.pdb","")
    ouput_file_path = receptor_path
    os.chdir("/home/yavuz/yavuz_proje/allosteric_binding_site/script")
    code=f"./pythonsh prepare_receptor4.py -r {receptor} -o {receptor_path}{only_receptor_name}.pdbqt"
    code_1=f"./pythonsh prepare_ligand4.py -l {ligand} -o {receptor_path}{only_ligand_name}.pdbqt"
    subprocess.call(code, shell=True)
    subprocess.call(code_1, shell=True)
    os.chdir(ouput_file_path)
    vina_code = f"{PSO_VINA_PATH}/psovina2 --receptor {receptor_path}{only_receptor_name}.pdbqt --ligand " \
                f"{receptor_path}{only_ligand_name}.pdbqt --out {only_receptor_name}_{only_ligand_name}.out " \
                f"--log {only_receptor_name}_{only_ligand_name}.log --center_x {x} --center_y {y} --center_z {z} " \
                f"--size_x {size_x} --size_y {size_y} --size_z {size_z} "
    print("vina_code",vina_code)
    subprocess.call(vina_code, shell=True)
    #shutil.copy(receptor, f"{ouput_file_path}/{only_receptor_name}.pdb")

if __name__ == '__main__':
    """MGLTOOLS_PATH="/home/yavuz/yavuz_proje/protac/bin/MGLTools-1.5.6"
    receptor = "/home/yavuz/yavuz_proje/protac/test/receptor.pdb"  # protac structure should be here
    ligand = "/home/yavuz/yavuz_proje/protac/test/protac.pdb"
    PSO_VINA_PATH="/home/yavuz/yavuz_proje/protac/bin"
    x=12
    y=12
    z=12
    size_x =12
    size_y =12
    size_z =12"""
    psoVina_local_docking(
        receptor=receptor,
        ligand=ligand,
        x=x,
        y=y,
        z=z,
        size_x=size_x,
        size_y=size_y,
        size_z=size_z,
        MGLTOOLS_PATH=MGLTOOLS_PATH,
        PSO_VINA_PATH=PSO_VINA_PATH
    )
