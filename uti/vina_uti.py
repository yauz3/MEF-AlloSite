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



def read_vina_log_file(
    vina_log_file_path,
    verbose: bool = True,
    ):
    vina_log_file = vina_log_file_path.split("/")[-1].split(".")[0]
    print("vina_log_file",vina_log_file)
    vina_log_path =vina_log_file_path.replace(f"{vina_log_file}.log","")
    if verbose:
        print ("Reading Vina energies from Vina log file located at", vina_log_file)
    lines = []
    print(vina_log_path)
    os.chdir(vina_log_path)
    with open(f"{vina_log_file}.log", "r") as f:
        for line in map(str.rstrip, f.readlines()):
            if re.match(r"^ {2,3}[1-9]+", line):
                lines.append(line.split())

    modes = {}
    for mode, energy, rmsd_lb, rmsd_ub in lines:
        mode = int(mode)
        energy = float(energy)
        rmsd_lb = float(rmsd_lb)
        rmsd_ub = float(rmsd_ub)
        modes[mode] = {
            "mode": mode,
            "energy": energy,
            "rmsd_lb": rmsd_lb,
            "rmsd_ub": rmsd_ub,
        }
    print(modes)
    return modes

def convert_and_prepare_vina_complex(
    vina_output_filename,
    ):
    vina_out_file = vina_output_filename.split("/")[-1].split(".")[0]
    receptor_name= f'{vina_out_file.split("_")[0]}_atm_fixed.pdb'
    pocket_number=vina_out_file.split("_")[0]
    vina_out_path =vina_output_filename.replace(f"{vina_out_file}.out","")
    os.chdir(vina_out_path)
    result_list_docking=glob.glob(f"{pocket_number}*.out")
    print(f"vina_out_path, {vina_out_path}") # /home/yavuz/yavuz_proje/allosteric_binding_site/docking/1AO0/
    print("receptor_name",receptor_name) # pocket10_atm_fixed
    print("pocket_number",pocket_number)
    print("result_list_docking",result_list_docking)
    for i in range(1, 2):
        try:
            results = (open(result_list_docking[0])).read()
            start = results.find("MODEL {}".format(i))
            i = i + 1
            end = results.find("MODEL {}".format(i))
            atoms = results[start:end]
            gece = open(f"{pocket_number}_pose_{i - 1}.pdbqt", "w")
            for at in atoms:
                gece.writelines("%s" % at)
            gece.close()
            # convertion
            os.system(f"obabel -ipdbqt {pocket_number}_pose_{i - 1}.pdbqt -opdb -O {pocket_number}_pose_{i - 1}.pdb")
            # complex preparation
            cmd.reinitialize()
            cmd.load(receptor_name)
            cmd.load(f"{pocket_number}_pose_{i - 1}.pdb")
            cmd.save(f"{pocket_number}_complex_{i - 1}.pdb")
            # cmd.remove(f"{pocket_number}_pose_{i - 1}.pdb")
            cmd.remove("all")
        except:
            continue

if __name__ == "__main__":
    """vina_log_file_path='/home/yavuz/yavuz_proje/protac/outputs/Vina_OUT/receptor_clean/receptor_clean_protac.log'
    input_path="/home/yavuz/yavuz_proje/protac/test"
    vina_output_filename='/home/yavuz/yavuz_proje/protac/outputs/Vina_OUT/receptor_clean/receptor_clean_protac.out'
    """
    read_vina_log_file(
        vina_log_file_path=vina_log_file_path)
    convert_and_prepare_vina_complex(
        vina_output_filename=vina_output_filename,)
