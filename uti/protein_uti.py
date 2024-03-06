#!/usr/bin/env python
# -*- coding: utf-8 -*-
# 22/02/2021
# Author: Sadettin Y. Ugurlu & David McDonald


import os
import datetime
import math
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
import os
from collections import OrderedDict
import math
import numpy as np
from pymol.cgo import *
from scipy.optimize import minimize
from pymol import cmd
import pandas as pd
import sys
import argparse
import timeit
import statistics
import textwrap
import datetime
import subprocess
import csv
import time
import scoria
import math
import numpy as np
from pymol.cgo import *
from pathlib import Path
from Bio.PDB import PDBParser
from scipy.optimize import minimize
from pymol.cgo import *
import pandas as pd
import sys
import argparse
import timeit
import statistics
import textwrap
import datetime
import os.path
import time
import csv
import os
import shutil
import numpy
import scoria
import numpy as np
import mdtraj as md
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBParser, PDBIO, Select
import glob
import re


def get_coordinate_of_protein_using_mdtraj(pdbfile):
    import mdtraj as md
    traj = md.load_pdb(pdbfile)
    df = pd.DataFrame(traj.xyz[0])
    return df
    #############################################################3
def get_coordinate_of_CA(pdbid):
    from Bio.PDB.PDBParser import PDBParser
    parser = PDBParser(PERMISSIVE=1)
    only_pdbid = pdbid.split("/")[-1].split(".")[0]
    structure_id = only_pdbid
    filename = f"{only_pdbid}.pdb"
    input_file_path=pdbid.replace(filename,"")
    print(input_file_path)
    os.chdir(input_file_path)
    structure = parser.get_structure(structure_id, filename)
    model = structure[0]
    try:
        output_dictionary={}
        b=0
        for chain in model.get_list():
            for residue in chain.get_list():
                ca = residue["CA"]
                #print(ca.get_coord())
                output_dictionary[b]=ca.get_coord()
                b=b+1
        #print(output_dictionary)
    except:
        print("Protein structure should be cleaned")
    output_dataframe=pd.DataFrame.from_dict(output_dictionary, orient='index')
    print(output_dataframe)
    return output_dataframe
def get_coordinate_of_protein(pdbid):
    #############################################################3
    only_pdbid_name = pdbid.split("/")[-1].split(".")[0]
    structure_id = only_pdbid_name
    filename = f"{only_pdbid_name}.pdb"
    input_file_path = pdbid.replace(filename, "")
    os.chdir(input_file_path)
    p = PDBParser()
    s = p.get_structure(structure_id, filename)
    result = []
    for chains in s:
        for chain in chains:
            for residue in chain:
                for atom in residue:
                    a = atom.get_coord()
                    result.append(atom.get_coord())  # for dönsügünüsü [] yazdırmak lazım

    coordinate_dataframe = pd.DataFrame(result)
    print(coordinate_dataframe)
    return coordinate_dataframe

def load_scoria_molecule(
        mol_path: str,
        verbose: bool = False,
):
    """Use the `scoria` library to build a molecule model from the 3D structure file.
    Parameters
    ----------
    mol_path : str
        The path to the 3D file. Will be converted to PDB format, if it is not already.
    verbose : bool, optional
        Flag to print updates to a console, by default False
    Returns
    -------
    scoria.Molecule
        A `scoria` molecule object constructed from the given file
    """

    if verbose:
        print("Building scoria molecule from", mol_path)

    # convert file if necessary
    if not mol_path.endswith(".pdb"):
        stem, ext = os.path.splitext(mol_path)
        ext = ext.replace(".", "")  # remove .
        mol_path = obabel_convert(
            input_format=ext,
            input_filename=mol_path,
            output_format="pdb",
            output_filename=stem + ".pdb")

    return scoria.Molecule(mol_path)

def identify_centre_of_mass(
        mol_path: str = None,
        mol: scoria.Molecule = None,
        precision: int = 3,
        geometric: bool = True,
        verbose: bool = False,
):
    """Use scoria to compute the center of mass / geometric center for a PDB file or pre-constructed scoria Molecule.
    Parameters
    ----------
    mol_path : str, optional
        Path of 3D molecule file, by default None
    mol : scoria.Molecule, optional
        Pre-constructed scoria Molecule, by default None
    precision : int, optional
        Rounding precision, by default 3
    geometric : bool, optional
        Flag to compute geometric mean, rather than center of mass, by default True
    verbose : bool, optional
        Flag to print updates to the console, by default False
    Returns
    -------
    tuple
        Computed x, y, z co-ordinates of the center
    """

    if verbose:
        print("Using scoria to compute center of molecule")

    if mol is None:  # construct molecule
        assert mol_path is not None and isinstance(mol_path, str)
        mol = load_scoria_molecule(mol_path)

    if not geometric:
        if verbose:
            print("Calculating mass center")
        try:
            center_x, center_y, center_z = mol.get_center_of_mass().data
        except Exception as e:
            print("Mass center failed")
            print("Exception was", e)
            geometric = True

    if geometric:
        if verbose:
            print("Calculating geometric center", )
        try:
            center_x, center_y, center_z = mol.get_geometric_center()
        except Exception as e:
            print("geometric center calculation error", e)
            return None

    # round using precision
    if isinstance(precision, int):
        center_x = round(center_x, precision)
        center_y = round(center_y, precision)
        center_z = round(center_z, precision)
    #print(center_x, center_y, center_z)
    return center_x, center_y, center_z
def get_bounding_box_size(
        mol: scoria.Molecule,
        scale: float = 1.,
        allowance: float = 0.,
        verbose: bool = False,
):
    """Use `scoria` to compute the size of the bounding box required to contain the molecule.
    Optionally, scale the box by `scale`, and add `allowance` to all cartesian co-ordinates.
    Parameters
    ----------
    mol : scoria.Molecule
        A pre-contructed scoria Molecule, or a path to a 3D structure file.
    scale : float, optional
        Scalar to apply to all dimensions, by default 1.
    allowance : float, optional
        Extra Angstroms to add to all dimensions, by default 3.
    verbose : bool, optional
        Flag to print updates to the console, by default False
    Returns
    -------
    np.ndarray
        A numpy array of shape (3,) containing the x, y, z bounding box sizes
    """
    if not isinstance(mol, scoria.Molecule):  # pdb file input
        mol = load_scoria_molecule(mol, verbose=verbose)
    if verbose:
        print("determining bounding box of mol")
    bounding_box = mol.get_bounding_box()
    #print(np.ceil(np.abs(bounding_box[0] - bounding_box[1])) * scale + allowance)
    return np.ceil(np.abs(bounding_box[0] - bounding_box[1])) * scale + allowance

def get_closes_coordinates (
    pdb_1,
    pdb_2
):
    protein_1_coordinates = get_coordinate_of_protein_using_mdtraj(pdb_1)
    protein_2_coordinates = get_coordinate_of_protein_using_mdtraj(pdb_2)
    #print(protein_2_coordinates)
    #print(protein_1_coordinates)
    print("calculating distance of atom requires amount of time")
    distance_dictionary={}
    dist_list=[]
    for index, value in protein_1_coordinates.iterrows():
        for index_2, value_2 in protein_2_coordinates.iterrows():
            point_1=(value[0],value[1],value[2])
            point_2=(value_2[0],value_2[1],value_2[2])
            dist = math.sqrt(sum([(a - b) ** 2 for a, b in zip(point_1, point_2)]))
            distance_dictionary[dist]=point_1,point_2
            dist_list.append(dist)
    #print(dist_list)
    #print(distance_dictionary[min(dist_list)])
    #print(min(dist_list))
    return distance_dictionary[min(dist_list)]


def checkIfDuplicates_1(listOfElems):
    ''' Check if given list contains any duplicates '''
    if len(listOfElems) == len(set(listOfElems)):
        return False
    else:
        return True

def CA_distance(target_pdb,ligand_pdb):
    from Bio.PDB.ResidueDepth import ResidueDepth
    from Bio.PDB.PDBParser import PDBParser
    from Bio.PDB.PDBParser import PDBParser
    import numpy
    from Bio.PDB.PDBParser import PDBParser
    only_target_name = target_pdb.split("/")[-1].split(".")[0]
    only_ligand_name = ligand_pdb.split("/")[-1].split(".")[0]
    input_path=target_pdb.replace(f"{only_target_name}.pdb","")
    os.chdir(input_path)
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{only_target_name}", f"{only_target_name}.pdb")
    structure_2 = parser.get_structure(f"{only_ligand_name}", f"{only_ligand_name}.pdb")
    residues_1 = [r for r in structure_1.get_residues()]
    residues_2 =[r for r in structure_2.get_residues()]
    out_put_dictionary={}
    out_put_dictionary_0={}
    distance_list=[]
    for number_1 in range(len(residues_1)):
        for number_2 in range(len(residues_2)):
            one  = residues_1[number_1]["CA"].get_coord()
            two = residues_2[number_2]["CA"].get_coord()
            #print('{} - {} = {}'.format(one,two, numpy.linalg.norm(one-two)))
            carbon_distance=numpy.linalg.norm(one-two)
            distance_list.append(carbon_distance)
            out_put_dictionary[one[0], one[1], one[2], two[0], two[1], two[2], carbon_distance] = \
                (one[0], one[1], one[2], two[0], two[1], two[2], carbon_distance)
            out_put_dictionary_0[carbon_distance] = one[0], one[1], one[2], two[0], two[1], two[2]
            # print(distance)
    lowest_count = distance_list.count(min(distance_list))
    if lowest_count > 1:
        print("more than one point detected")
        for key, value in out_put_dictionary.items():
            if min(distance) == value[6]:
                print(key, min(distance_list), value[6])
        exit()
    # print(out_put_dictionary_0)
    print(out_put_dictionary_0.get(min(distance_list)))
    print(min(distance_list))


    # it will return target_closest_coordinates,
    print("(target_x, target_y, target_z, ligand_x, ligand_y, ligand_z), min_distance, mean_distance, median_distance")
    print(out_put_dictionary_0.get(min(distance_list)), min(distance_list), statistics.mean(distance_list), statistics.median(distance_list))
    return out_put_dictionary_0.get(min(distance_list)), min(distance_list), statistics.mean(distance_list), statistics.median(distance_list)

def chain_distance(target_pdb,ligand_pdb):
    only_target_name = target_pdb.split("/")[-1].split(".")[0]
    only_ligand_name = ligand_pdb.split("/")[-1].split(".")[0]
    input_path = target_pdb.replace(f"{only_target_name}.pdb", "")
    os.chdir(input_path)
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{only_target_name}", f"{only_target_name}.pdb")
    structure_2 = parser.get_structure(f"{only_ligand_name}", f"{only_ligand_name}.pdb")

    # For each chain
    b=0
    output_dictionary={}
    distance_list=[]
    for chain_1 in structure_1.get_chains():
        b=b+1
        c=0
        for chain_2 in structure_2.get_chains():
            one=(chain_1.center_of_mass())
            two=(chain_2.center_of_mass())
            chain_distance=numpy.linalg.norm(one - two)
            c=c+1
            #print(one)
            #print("chain",chain_2)
            #print(f"chain_mass_distance {chain_distance} between {chain_1} and {chain_2}")
            output_dictionary[chain_distance]=chain_1,chain_2
            distance_list.append(chain_distance)
    print(output_dictionary.get(min(distance_list)),min(distance_list))
    return output_dictionary.get(min(distance_list)),min(distance_list)

def prepare_cloest_chains(target_pdb,center_x,center_y,center_z):
    only_target_name = target_pdb.split("/")[-1].split(".")[0]
    target_path=target_pdb.replace(f"{only_target_name}.pdb","")
    os.chdir(target_path)
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{only_target_name}", f"{only_target_name}.pdb")
    # For each chain
    chain_out_dict={}
    for chain_1 in structure_1.get_chains():
        point1 = (chain_1.center_of_mass())
        print(point1)
        point2 = center_x,center_y,center_z
        print(point2)
        chain_distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(point1, point2)]))
        print(chain_distance)
        chain_1=chain_1.id
        chain_out_dict[chain_1]=chain_distance
    chain_out_dict={k: v for k, v in sorted(chain_out_dict.items(), key=lambda item: item[1])}
    print("chain_out")
    print(chain_out_dict)
    #print(chain_dis_list)
    cmd.load(f"{only_target_name}.pdb")
    cmd.remove("resn HOH")
    # cmd.h_add("all")
    # cmd.h_add(selection="acceptor or donors")
    closest_chain = list(chain_out_dict.keys())[0]
    second_closest_chain = list(chain_out_dict.keys())[1]
    print("the cloest chain: ", closest_chain)
    print("the second closest chain: ", second_closest_chain)
    b=0
    for key,value in chain_out_dict.items():
        if b >1:
            print(f"chain {key} is removed")
            cmd.remove(f"chain {key}")
            cmd.remove("het")
            cmd.save("protein_protein_complex.pdb")
        b=b+1
    cmd.remove(f"chain {closest_chain}")
    cmd.save("protein_1.pdb")
    cmd.remove("all")
    cmd.load("protein_protein_complex.pdb")
    cmd.remove(f"chain {second_closest_chain}")
    cmd.save("protein_2.pdb")
    data_1=pd.read_csv("protein_1.pdb")
    data_2=pd.read_csv("protein_2.pdb")
    print(len(data_2))
    if len(data_1) > len(data_2):
        os.rename("protein_1.pdb", "target.pdb")
        os.rename("protein_2.pdb", "ligand.pdb")
    else:
        os.rename("protein_2.pdb", "target.pdb")
        os.rename("protein_1.pdb", "ligand.pdb")
    return print("target and ligand pdb files are ready!")

def residue_distance(target_pdb,ligand_pdb):
    from Bio.PDB import PDBParser
    only_target_name = target_pdb.split("/")[-1].split(".")[0]
    only_ligand_name = ligand_pdb.split("/")[-1].split(".")[0]
    input_path = target_pdb.replace(f"{only_target_name}.pdb", "")
    os.chdir(input_path)
    out_put_dictionary={}
    out_put_dictionary_0= {}
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{only_target_name}", f"{only_target_name}.pdb")
    structure_2 = parser.get_structure(f"{only_ligand_name}", f"{only_ligand_name}.pdb")
    #com = structure_1.center_of_mass()
    distance=[]
    # For each chain
    b=0
    for residue_1 in structure_1.get_residues():
        b=b+1
        c=0
        for residue_2 in structure_2.get_residues():
            one=(residue_1.center_of_mass())
            two=(residue_2.center_of_mass())
            residue_distance=numpy.linalg.norm(one - two)
            #print(residue_distance)
            c=c+1
            distance.append(float(residue_distance))
            out_put_dictionary[one[0],one[1],one[2],two[0],two[1],two[2],residue_distance]=\
                (one[0],one[1],one[2],two[0],two[1],two[2],residue_distance)
            out_put_dictionary_0[residue_distance] =one[0], one[1], one[2], two[0], two[1], two[2]
    #print(distance)
    lowest_count=distance.count(min(distance))
    if lowest_count > 1:
        print("more than one point detected")
        for key,value in out_put_dictionary.items():
            if min(distance) == value[6]:
                print(key,min(distance),value[6])
        exit()
    #print(out_put_dictionary_0)
    print(out_put_dictionary_0.get(min(distance)))
    print(min(distance))

    print("structure 1 has ", b," residue")
    print("structure 2 has ",c," residue")
    # it will return target_closest_coordinates,
    print("(target_x, target_y, target_z, ligand_x, ligand_y, ligand_z), min_distance, mean_distance, median_distance")
    print(out_put_dictionary_0.get(min(distance)), min(distance), statistics.mean(distance),statistics.median(distance))
    return out_put_dictionary_0.get(min(distance)), min(distance), statistics.mean(distance),statistics.median(distance)
#merge_protein_files("complex_3","1nwl","1pne") # it shold be complex, target_protein and ligand_protein.
# this function is for megadock output
def merge_protein_files(complex_pdb,target_pdb,ligand_pdb,):
    cmd.load(f"{target_pdb}.pdb") # it is fixed
    cmd.load(f"{ligand_pdb}.pdb")
    cmd.load(f"{complex_pdb}.pdb")
    align_1=cmd.align(f"{ligand_pdb}",f"{complex_pdb}_0002")
    print("align_1 :",align_1)
    align_2=cmd.align(f"{target_pdb}",f"{complex_pdb}_0001")
    print("align_2 :",align_2)
    if align_1[0] and align_2[0] < 1:
        cmd.remove(f"{complex_pdb}_0002")
        cmd.remove(f"{complex_pdb}_0001")
        cmd.save("deneme.mol")
        os.system("obabel -imol deneme.mol -opdb -O deneme.pdb")
    else:
        align_1 = cmd.align(f"{target_pdb}", f"{complex_pdb}_0002")
        print("align_1 :", align_1)
        align_2 = cmd.align(f"{ligand_pdb}", f"{complex_pdb}_0001")
        print("align_2 :", align_2)
        cmd.remove(f"{complex_pdb}_0002")
        cmd.remove(f"{complex_pdb}_0001")
        cmd.save("deneme.mol")
        os.system("obabel -imol deneme.mol -opdb -O deneme.pdb")
    # TODO:  cheking whether formed pdb is the same or not using aliging is not possible. Current solution is to use
    # pymol to see exactly the same.

def is_het(residue):
    res = residue.id[0]
    return res != " " and res != "W"
class ResidueSelect(Select):
    def __init__(self, chain, residue):
        self.chain = chain
        self.residue = residue

    def accept_chain(self, chain):
        return chain.id == self.chain.id

    def accept_residue(self, residue):
        """ Recognition of heteroatoms - Remove water molecules """
        return residue == self.residue and is_het(residue)
def extract_ligands(pdb_code,output_path):
    """ Extraction of the heteroatoms of .pdb files """
    only_target_name = pdb_code.split("/")[-1].split(".")[0]
    input_path = pdb_code.replace(f"{only_target_name}.pdb", "")
    os.chdir(input_path)
    pdb = PDBParser().get_structure(only_target_name,f"{only_target_name}.pdb")
    io = PDBIO()
    io.set_structure(pdb)
    for model in pdb:
        for chain in model:
            for residue in chain:
                if not is_het(residue):
                    continue
                ligand_symbol = residue.get_resname()
                print(f"saving {chain} {residue} {ligand_symbol}")
                if len(ligand_symbol) > 2:
                    os.chdir(output_path)
                    io.save(f"{only_target_name}_{ligand_symbol}.pdb", ResidueSelect(chain, residue))
                    print(f"{only_target_name}_{ligand_symbol}.pdb is saved")

def calculate_RMSD_pymol(
        reference_filename,
        model_filename,
        return_if_fail=float("inf"),
):
    if not reference_filename.endswith(".pdb"):
        reference_basename, ext = os.path.splitext(reference_filename)
        ext = ext.replace(".", "")
        reference_filename = obabel_convert(
            input_format=ext,
            input_filename=reference_filename,
            output_format="pdb",
            output_filename=reference_basename + ".pdb"
        )

    if not model_filename.endswith(".pdb"):
        model_basename, ext = os.path.splitext(model_filename)
        ext = ext.replace(".", "")
        model_filename = obabel_convert(
            input_format=ext,
            input_filename=model_filename,
            output_format="pdb",
            output_filename=model_basename + ".pdb"
        )

    print("Using PyMol to compute RMSD for:", reference_filename, model_filename)

    # refresh program completely to allow re-reading previously read files
    cmd.reinitialize()

    # load the files into PyMol
    cmd.load(reference_filename, object="_reference",
             # state=0,
             )
    cmd.load(model_filename, object="_model",
             # state=0,
             )

    # set residue ID the same
    cmd.alter("all", "resn='LIG'")
    cmd.alter("all", "resi=''")
    cmd.alter("all", "resv=1")

    # set name to element
    cmd.alter("all", "name=elem")
    # remove chain
    cmd.alter("all", "chain=''")
    # remove sequence id
    cmd.alter("all", "segi=''")

    # remove hydrogen
    cmd.remove("hydro")
    # add hydrogen
    # cmd.h_add("all")

    cmd.sort("all")

    # cmd.save("test.pdb")

    try:
        # ret = cmd.align(
        #     "_reference", "_model",
        #     cutoff=100, # only relevant for outlier rejection
        #     cycles=0, # consider all atoms (no outlier rejection)
        #     transform=0,
        # )
        ret = cmd.rms_cur("_reference", "_model")
        print("obtained", ret)

        cmd.remove("all")

        # return ret[0] # RMSD
        return ret  # RMSD
    except Exception as e:
        print("PYMOL RMSD exception", e)
        # raise e
        return return_if_fail

def prepare_chain(target_pdb,chain_list):
    only_target_name = target_pdb.split("/")[-1].split(".")[0]
    target_path = target_pdb.replace(f"{only_target_name}.pdb", "")
    os.chdir(target_path)
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{only_target_name}", f"{only_target_name}.pdb")
    # For each chain
    chain_drop_list=[]
    for chain_1 in structure_1.get_chains():
        chain_1=chain_1.id
        if chain_1 not in chain_list:
            chain_drop_list.append(chain_1)
    print("chain_drop",chain_drop_list)
    cmd.load(f"{only_target_name}.pdb")
    cmd.remove("resn HOH")
    # cmd.h_add("all")
    # cmd.h_add(selection="acceptor or donors")
    for key in chain_drop_list:
        print(f"chain {key} is removed")
        cmd.remove(f"chain {key}")
        cmd.remove("het")
        cmd.save("protein_protein_complex.pdb")
    cmd.remove(f"chain {chain_list[0]}")
    cmd.save("protein_1.pdb")
    cmd.remove("all")
    cmd.load("protein_protein_complex.pdb")
    cmd.remove(f"chain {chain_list[1]}")
    cmd.save("protein_2.pdb")
    data_1=pd.read_csv("protein_1.pdb")
    data_2=pd.read_csv("protein_2.pdb")
    if len(data_1) > len(data_2):
        os.rename("protein_1.pdb", "target.pdb")
        os.rename("protein_2.pdb", "ligand.pdb")
    else:
        os.rename("protein_2.pdb", "target.pdb")
        os.rename("protein_1.pdb", "ligand.pdb")
    return print("target and ligand pdb files are ready!")

def download_pdb(pdbid):
    # run terminal to downloand pdb files
    try:
        os.system("pdb_fetch -biounit %s > %s.pdb" % (pdbid, pdbid))
        # os.system("pdb_fetch %s > %s.pdb" % (pdbid, pdbid))
    except:
        print("Please use pdbid without extention")
    return print("PDBid downloaded")

def clean_pdb(input_path, pdbid):
    os.chdir(input_path)  # Change the current working directory to the specified input path
    cmd.reinitialize()
    cmd.load(pdbid)  # Load the PDB file with the given PDB ID into PyMOL
    cmd.remove("resn HOH")  # Remove all residues with the name 'HOH' (water molecules) from the loaded structure
    # cmd.h_add("all")  # (This line is commented out) Add hydrogen atoms to all atoms in the structure (not currently active)
    cmd.remove("HETATM")  # Remove all heteroatoms (non-standard residues) from the loaded structure
    # cmd.h_add(selection="acceptor or donors")  # (This line is commented out) Add hydrogen atoms to atoms identified as acceptors or donors (not currently active)
    cmd.save(
        f"{pdbid.split('.')[0]}_cleaned.pdb")  # Save the cleaned structure with the original PDB ID plus '_cleaned' as the file name
    cmd.remove("all")

def get_ligand_chain_ids(example):
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{example}", f"{example}.pdb")
    b = 0
    chain_size={}
    chain_length_size=[]
    for chain_1 in structure_1.get_chains():
        size=[]
        for residue in chain_1:
            size.append(residue)
            chain_length_size.append(residue)
        chain_size[chain_1.id]=len(size)
    try:
        for key, value in chain_size.items():
            if value == 1:
                selected_chain=key
    except:
        for key, value in chain_size.items():
            if value == min(chain_length_size):
                selected_chain=key
    return selected_chain

def prepare_desired_chain(target_pdb,chain_list):
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{target_pdb}", f"{target_pdb}.pdb")
    # For each chain
    chain_drop_list=[]
    for chain_1 in structure_1.get_chains():
        chain_1=chain_1.id
        if chain_1 not in chain_list:
            chain_drop_list.append(chain_1)
    print("chain_drop",chain_drop_list)
    cmd.load(f"{target_pdb}.pdb")
    cmd.remove("resn HOH")
    # cmd.h_add("all")
    # cmd.h_add(selection="acceptor or donors")
    for key in chain_drop_list:
        print(f"chain {key} is removed")
        cmd.remove(f"chain {key}")
        cmd.remove("het")
    cmd.save(f"{target_pdb}_ligands.pdb")
    cmd.remove("all")
    print("desired chains are ready!")


def get_coordinate_pdb(pdbid):
    #############################################################3
    structure_id = pdbid.split(".")[0]
    filename = pdbid

    p = PDBParser()
    s = p.get_structure(structure_id, filename)
    result = []
    for chains in s:
        for chain in chains:
            for residue in chain:
                for atom in residue:
                    a = atom.get_coord()
                    result.append(atom.get_coord())  # for dönsügünüsü [] yazdırmak lazım

    cordinate_dataframe = pd.DataFrame(result)
    #print(c)
    return cordinate_dataframe

def get_coordinate_pdb_without_C(pdbid):
    #############################################################3
    structure_id = pdbid.split(".")[0]
    filename = pdbid

    p = PDBParser()
    s = p.get_structure(structure_id, filename)
    result = []
    for chains in s:
        for chain in chains:
            for residue in chain:
                for atom in residue:
                    if "C" not in str(atom):
                        a = atom.get_coord()
                        result.append(atom.get_coord())  # for dönsügünüsü [] yazdırmak lazım

    cordinate_dataframe = pd.DataFrame(result)
    #print(c)
    return cordinate_dataframe


def residue_mass_center(pdbid):
    from Bio.PDB import PDBParser
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{pdbid}", f"{pdbid}.pdb")
    output_dictionary={}
    b=1
    for residue_1 in structure_1.get_residues():
        one=(residue_1.center_of_mass())
        output_dictionary[b]=one
        b += 1
    output_dataframe=pd.DataFrame.from_dict(output_dictionary, orient='index')
    #print(output_dataframe)
    return output_dataframe

def get_chain_list(target_pdb):
    parser = PDBParser()
    structure_1 = parser.get_structure(f"{target_pdb}", f"{target_pdb}.pdb")
    # For each chain
    chain_list=[]
    for chain_1 in structure_1.get_chains():
        chain_list.append (chain_1.id)
    print("chain_list",chain_list)
    return chain_list

def promissing_pose_selection_receptor(input_file):
    test_path="/home/yavuz/yavuz_proje/protac/test"
    mega_out_path="/home/yavuz/yavuz_proje/protac/outputs/Megadock_OUT"

    os.chdir(mega_out_path)
    megadock_output_file = glob.glob(f"{input_file}_receptor_complex_*")
    output_dictionary={}
    for megadock_out in megadock_output_file:
        os.chdir(f"{mega_out_path}/{megadock_out}")
        output_file = glob.glob("*.out")
        with open(str(output_file[0])) as f:
            lines = f.readlines()
            home_line_number = 0
            pose_number=0
            for line in lines:
                if home_line_number == 2:
                    MegadockScore = (line.split()[-1])
                    pose_number +=1
                if "home" in str(line):
                    home_line_number += 1
                if pose_number > 0:
                    output_dictionary[megadock_out,pose_number]=MegadockScore
    ordered_dictionary={k: v for k, v in sorted(output_dictionary.items(), key=lambda item: item[1],reverse=True)}
    ordered_list = list(ordered_dictionary.keys())
    if os.path.exists(f'../{input_file}_selected_pose') == False:
        os.makedirs(f'../{input_file}_selected_pose')
    distance_dictionary={}
    for promissing_pose in range(5):
        result_path=ordered_list[promissing_pose][0]
        promissing_pose_number=ordered_list[promissing_pose][1]
        #print(result_path) # 6BOY_receptor_complex_2
        #print(promissing_pose_number) # 1
        coordinate_protac=get_coordinate_pdb_without_C(f"{mega_out_path}/{input_file}_receptor_protac/ligand_{promissing_pose_number}.pdb")
        os.chdir(f"{mega_out_path}/{result_path}")
        predicted_coordinate =get_coordinate_pdb_without_C(f"ligand_{promissing_pose_number}.pdb")
        """print(coordinate_protac)
        print(predicted_coordinate)"""
        dist = 0
        for key, value in predicted_coordinate.iterrows():
            point_1 = value[0], value[1], value[2]
            for key_2, value_2 in coordinate_protac.iterrows():
                point_2 = value_2[0], value_2[1], value_2[2]
                distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(point_1, point_2)]))
                dist = dist + distance
        distance_dictionary[f"{result_path}_{promissing_pose_number}"]=dist
    ordered_distance_dictionary={k: v for k, v in sorted(distance_dictionary.items(), key=lambda item: item[1])}
    ordered_distance_list = list(ordered_distance_dictionary.keys())
    print("ordered_distance")
    print(ordered_distance_dictionary)
    #half_list_distance=ordered_distance_list[0:(int(len(ordered_distance_list) / 2))]
    half_list_distance=ordered_distance_list
    final_dictionary={}
    b=1
    for promissing_pose in range(5):
        result_path=ordered_list[promissing_pose][0]
        promissing_pose_number=ordered_list[promissing_pose][1]
        if f"{result_path}_{promissing_pose_number}" in half_list_distance:
            os.chdir(f"{mega_out_path}/{result_path}")
            cmd.reinitialize()
            cmd.load(f"{result_path}.pdb")
            cmd.load(f"ligand_{promissing_pose_number}.pdb")
            os.chdir(f'../{input_file}_selected_pose')
            cmd.save(f"{result_path}_{promissing_pose_number}.pdb")
            #cmd.save(f"{input_file}_ligand_{b}.pdb")
            cmd.remove("all")
            final_dictionary[f"{result_path}_{promissing_pose_number}"]=b
            b+=1
    print(final_dictionary)
    return final_dictionary


def promissing_pose_selection_target(input_file):
    test_path="/home/yavuz/yavuz_proje/protac/test"
    mega_out_path="/home/yavuz/yavuz_proje/protac/outputs/Megadock_OUT"

    os.chdir(mega_out_path)
    megadock_output_file = glob.glob(f"{input_file}_target_complex_*")
    output_dictionary={}
    for megadock_out in megadock_output_file:
        os.chdir(f"{mega_out_path}/{megadock_out}")
        output_file = glob.glob("*.out")
        with open(str(output_file[0])) as f:
            lines = f.readlines()
            home_line_number = 0
            pose_number=0
            for line in lines:
                if home_line_number == 2:
                    MegadockScore = (line.split()[-1])
                    pose_number +=1
                if "home" in str(line):
                    home_line_number += 1
                if pose_number > 0:
                    output_dictionary[megadock_out,pose_number]=MegadockScore
    ordered_dictionary={k: v for k, v in sorted(output_dictionary.items(), key=lambda item: item[1],reverse=True)}
    ordered_list = list(ordered_dictionary.keys())
    if os.path.exists(f'../{input_file}_selected_pose') == False:
        os.makedirs(f'../{input_file}_selected_pose')
    distance_dictionary={}
    for promissing_pose in range(10):
        result_path=ordered_list[promissing_pose][0]
        promissing_pose_number=ordered_list[promissing_pose][1]
        #print(result_path) # 6BOY_receptor_complex_2
        #print(promissing_pose_number) # 1
        coordinate_protac=get_coordinate_pdb_without_C(f"{mega_out_path}/{input_file}_target_protac/ligand_{promissing_pose_number}.pdb")
        os.chdir(f"{mega_out_path}/{result_path}")
        predicted_coordinate = get_coordinate_pdb_without_C(f"ligand_{promissing_pose_number}.pdb")
        """print(coordinate_protac)
        print(predicted_coordinate)"""
        dist = 0
        for key, value in predicted_coordinate.iterrows():
            point_1 = value[0], value[1], value[2]
            for key_2, value_2 in coordinate_protac.iterrows():
                point_2 = value_2[0], value_2[1], value_2[2]
                distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(point_1, point_2)]))
                dist = dist + distance
        distance_dictionary[f"{result_path}_{promissing_pose_number}"]=dist
    ordered_distance_dictionary={k: v for k, v in sorted(distance_dictionary.items(), key=lambda item: item[1])}
    ordered_distance_list = list(ordered_distance_dictionary.keys())
    #half_list_distance=ordered_distance_list[0:(int(len(ordered_distance_list) / 2))]
    half_list_distance=ordered_distance_list
    final_dictionary = {}
    b = 1
    for promissing_pose in range(10):
        result_path = ordered_list[promissing_pose][0]
        promissing_pose_number = ordered_list[promissing_pose][1]
        if f"{result_path}_{promissing_pose_number}" in half_list_distance:
            os.chdir(f"{mega_out_path}/{result_path}")
            cmd.reinitialize()
            cmd.load(f"{result_path}.pdb")
            cmd.load(f"ligand_{promissing_pose_number}.pdb")
            os.chdir(f'../{input_file}_selected_pose')
            cmd.save(f"{result_path}_{promissing_pose_number}.pdb")
            #cmd.save(f"{input_file}_ligand_{b}.pdb")
            cmd.remove("all")
            final_dictionary[f"{result_path}_{promissing_pose_number}"] = b
            b += 1
    return final_dictionary

def get_atom_size(pdbid):
    os.system(f"grep 'ATOM' {pdbid} > temp.txt")
    with open("temp.txt") as tem:
        lines=tem.readlines()
        print(len(lines))
    if len(lines) > 2:
        os.system("rm temp.txt")
        return  len(lines)
    else:
        os.system(f"grep 'HETATM' {pdbid} > temp.txt")
        with open("temp.txt") as tem:
            lines = tem.readlines()
            print(len(lines))
        os.system("rm temp.txt")
        return len(lines)


def get_chain_id(protein):

    # Get the absolute path to the file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Find ASD_Release_201909_AS in the directory
    file_path = os.path.join(script_dir, "..", "sourceData", "ASD_Release_201909_AS.txt")
    # Check if the file exists
    if not os.path.exists(file_path):
        print("Error: File not found.")
        return None

    # Read the file and extract the chain IDs
    with open(file_path, 'r') as fh:
        alloInfoOri = fh.readlines()
        # take database as a list for each protein (element)
        alloInfo = [line.split("\t") for line in alloInfoOri]
        # search database to find chain_id for a given protein
        for protein_info in alloInfo:
            if protein in protein_info:
                for item in protein_info:
                    if "Chain" in item:
                        # Zincir isimlerini bulmak için düzenli ifade kullanma
                        chain_list = re.findall(r'Chain ([A-Z]):', item)

    return chain_list


def get_info_from_ads(protein):
    # Get the absolute path to the file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, "..", "sourceData", "ASD_Release_201909_AS.txt")
    # Check if the file exists
    if not os.path.exists(file_path):
        print("Error: File not found.")
        return None

    # Read the file and extract the chain ID
    with open(file_path, 'r') as fh:
        alloInfoOri = fh.readlines()
        alloInfo = [line.split("\t") for line in alloInfoOri]
        for i in alloInfo:
            if protein in i:
                return i

    return None



def run_fpocket(protein_list, output_name ,chain_selection=True):
    # Find current directory
    current_dir = os.path.dirname(os.path.abspath(__file__))

    # Find parent directory
    parent_dir = os.path.dirname(current_dir)

    # Find data directory
    data_dir = os.path.join(parent_dir, "data")

    # run Fpocket for each protein in protein_list
    for protein in protein_list:
        # clean proteins before run Fpocket, since P2rank mentioning Fpocket have a potential data leak when proteins are not cleaned!
        clean_pdb(input_path=f"{data_dir}/pdbs",
                  pdbid=f"{protein}.pdb")
        if chain_selection:
            # get chain_id for protein from ADS
            chain_id = get_chain_id(protein)
            formatted_chain_ids = ','.join(str(item) for item in chain_id)
            print("chain_id",formatted_chain_ids)
            # run FPocket with chain selection option: -k
            os.system(f"fpocket -f {data_dir}/pdbs/{protein}_cleaned.pdb -k {formatted_chain_ids}")
        else:
            # run FPocket without chain selection option: -k
            os.system(f"fpocket -f {data_dir}/pdbs/{protein}_cleaned.pdb")
        # save outputs
        training_pockets_dir=f"{data_dir}/{output_name}/"
        # check whether output_dic exist. If not, create one.
        if not os.path.exists(training_pockets_dir):
            os.makedirs(training_pockets_dir)
        # move outputs to labelled directory
        os.system(f"mv {data_dir}/pdbs/{protein}_cleaned_out {training_pockets_dir}")



