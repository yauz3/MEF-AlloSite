import os
import shutil
import re
import subprocess

from functools import partial

from concurrent.futures import ProcessPoolExecutor as Pool

import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(__file__), 
    os.pardir,
    os.pardir,
    ))

if __name__ == "__main__":

    import sys
    import os.path

    sys.path.insert(1, PROJECT_ROOT)

from utils.sys_utils import execute_system_command
from utils.io.io_utils import copy_file, dataframe_to_dict, delete_directory, delete_file, load_json, write_json

from utils.molecules.pdb_utils import identify_centre_of_mass, get_bounding_box_size, download_pdb_structure_using_pdb_fetch

FPOCKET_OUT = os.path.join(
    "cavities", 
    "FPOCKET_out",
    )
POCKET_ANALYSIS_N_PROC = int(os.environ.get("POCKET_ANALYSIS_N_PROC", default=1))

def call_fpocket(
    input_filename: str, 
    main_job_id: int = None,
    verbose: bool = False,
    ):
    if verbose:
        print ("Running fpocket on file", input_filename)
    
    cmd = f'''
    fpocket -f {input_filename}
    '''
    if os.path.exists(input_filename):
        try:
            execute_system_command(cmd, main_job_id=main_job_id, verbose=verbose)
        except Exception as e:
            print ("Exception calling fpocket", e)
            return None

    return input_filename

def show_pocket_scores(
    pocket_score_filename: str,
    ):
    return dataframe_to_dict(pocket_score_filename)

def select_best_scoring_pocket(
    pocket_score_filename: str, 
    by: str = "drug_score",
    verbose: bool = False,
    ):

    pocket_scores = pd.read_csv(pocket_score_filename, sep=" ", index_col=0)

    if verbose:
        print("Using column", by, "to select best pocket")

    return pocket_scores[by].idxmax()

def return_pocket(
    best_scoring_pocket: int, 
    pocket_dir: str,
    ):
    pocket_filename = os.path.join(
        pocket_dir, 
        f"pocket{best_scoring_pocket}_atm.pdb")
    assert os.path.exists(pocket_filename)
    return pocket_filename

def read_fpocket_pdb_file(
    pdb_filename: str,
    verbose: bool = False,
    ):

    if verbose:
        print ("Reading fpocket output PDB file", pdb_filename)
    all_header_values = {}
    try:
        lines = subprocess.check_output(
            f"awk '/^HEADER [0-9]+/' {pdb_filename}", shell=True).decode("utf-8").split("\n")
        for line in lines:
            if line == "":
                continue
            line = "-".join(line.split("-")[1:])
            key, value = line.split(":")
            # leading whitespace
            key = re.sub(r"^\s+", "", key)
            # trailing whitespace
            key = re.sub(r"\s+$", "", key)
            # convert to lowercase and remove spaces
            key = key.lower().replace(" ", "_")
            # remove .
            key = key.replace(".", "")
            value = float(value)
            all_header_values[key] = value
    except:
        all_header_values = {}
    return all_header_values

def run_fpocket_and_collate_single_target(
    target_identifier: str, # pdb_id / target_name
    output_dir: str,
    min_pocket_score: float = 0,
    existing_pdb_filename: str = None,
    delete_output_dir: bool = True,
    verbose: bool = False,
    ):

    if verbose:
        print ("Running fpocket on target", target_identifier, "and collecting pocket data into a dictionary")

    output_dir = os.path.join(
        output_dir,
        FPOCKET_OUT)
    os.makedirs(output_dir, exist_ok=True)

    target_output_dir = os.path.join(
        output_dir,
        f"{target_identifier}_out")

    if not os.path.isdir(target_output_dir):
        pdb_filename = os.path.join(
            output_dir,
            f"{target_identifier}.pdb")

        if not os.path.exists(pdb_filename):
            if existing_pdb_filename is not None and os.path.exists(existing_pdb_filename):
                if existing_pdb_filename != pdb_filename:
                    copy_file(existing_pdb_filename, pdb_filename) # because the file is deleted after fpocket is called
            else:
                pdb_filename = download_pdb_structure_using_pdb_fetch(
                    pdb_id=target_identifier,
                    pdb_filename=pdb_filename,
                    verbose=verbose,
                )

        # call fpocket            
        call_fpocket(pdb_filename, verbose=verbose)
        # delete the pdb file
        delete_file(pdb_filename, verbose=verbose)

    target_pocket_dir = os.path.join(
        target_output_dir,
        "pockets",
        )

    fpocket_target_data = {}

    if os.path.isdir(target_pocket_dir):

        pocket_pdb_filenames = [file 
        for file in os.listdir(target_pocket_dir) 
        if file.endswith("atm.pdb") and "env" not in file]

        for pocket_pdb_filename in pocket_pdb_filenames:
            pocket_number = pocket_pdb_filename.split("_")[0]
            # remove `pocket`
            pocket_number = pocket_number.split("pocket")[1]

            pocket_pdb_filename = os.path.join(
                target_pocket_dir, 
                pocket_pdb_filename)
            
            # read pdb file as text file
            fpocket_pdb_data = read_fpocket_pdb_file(pocket_pdb_filename, verbose=verbose)
            if "pocket_score" not in fpocket_pdb_data:
                continue # skip pocket 
            pocket_score = fpocket_pdb_data["pocket_score"] 

            if min_pocket_score is not None and pocket_score < min_pocket_score:
                continue
        
            center_x, center_y, center_z = identify_centre_of_mass(pocket_pdb_filename, verbose=verbose)
            size_x, size_y, size_z = get_bounding_box_size(pocket_pdb_filename, verbose=verbose)

            fpocket_target_data[pocket_number] = {
                "center_x": center_x,
                "center_y": center_y,
                "center_z": center_z,
                "size_x": size_x,
                "size_y": size_y,
                "size_z": size_z,
                **fpocket_pdb_data,
                "pocket_pdb_filename": pocket_pdb_filename,
            }

    if len(fpocket_target_data) == 0:
        # NO POCKETS
        print ("Fpocket failed to identify targets for", target_identifier)
        fpocket_target_data[None] = {
            "center_x": None,
            "center_y":None,
            "center_z": None,
            "pocket_score": None,
            "drug_score": None,
            "number_of_alpha_spheres": None,
            "mean_alpha-sphere_radius": None,
            "mean_alpha-sphere_solvent_acc.": None,
            "mean_b-factor_of_pocket_residues": None,
            "hydrophobicity_score": None,
            "polarity_score": None,
            "amino_acid_based_volume_score": None,
            "pocket_volume_(monte_carlo)": None,
            "pocket_volume_(convex_hull)": None,
            "charge_score": None,
            "local_hydrophobic_density_score": None,
            "number_of_apolar_alpha_sphere": None,
            "proportion_of_apolar_alpha_sphere": None,
            "pocket_pdb_filename": None,
        }

    if delete_output_dir:
        delete_directory(target_output_dir, verbose=verbose)

    # return {target_identifier: fpocket_target_data}
    return fpocket_target_data

# def collate_fpocket_data(
#     targets,
#     output_filename,
#     output_dir,
#     min_pocket_score=0,
#     ):

#     if os.path.exists(output_filename):
#         print (output_filename, "ALREADY EXISTS -- LOADING IT")
#         fpocket_data = load_json(output_filename)
#     else:

#         print ("RUNNNING FPOCKET FOR TARGETS", targets, "USING MIN_SCORE", min_pocket_score)

#         # os.makedirs(FPOCKET_OUT, exist_ok=True)

#         fpocket_data = {}

#         with Pool(processes=POCKET_ANALYSIS_N_PROC) as p:
            
#             fpocket_data_all_targets = p.map(
#                 partial(
#                     run_fpocket_and_collate_single_target,
#                     output_dir=output_dir,
#                     min_pocket_score=min_pocket_score,
#                     existing_pdb_filename=None,   
#                 ),
#                 targets,
#             )

#             for fpocket_data_target in fpocket_data_all_targets:
#                 fpocket_data.update(fpocket_data_target)
            
#         print ("WRITING FPOCKET DATA TO", output_filename)
#         write_json(fpocket_data, output_filename)
    
#     return fpocket_data

if __name__ == "__main__":


    # input_filename = "pdbs/1uk0/1uk0.pdb"

    # output_dir = "test"
    # os.makedirs(output_dir, exist_ok=True)
    # pdb_id = "2OQK"
    # from utils.pdb_utils import download_pdb_file
    # pdb_filename = download_pdb_file(pdb_id, output_dir)

    # pocket_score_filename, pocket_dir = call_fpocket(pdb_filename)

    # best_scoring_pocket = select_best_scoring_pocket(pocket_score_filename)


    # print (best_scoring_pocket)

    # print (return_pocket(best_scoring_pocket, pocket_dir))

    # targets = ["2OQK", "6NNA"]

    # docking_directory = os.path.abspath("Docking")

    # fpocket_data = collate_fpocket_data(
    #     targets=targets,
    #     output_filename="FPOCKET_OUT.json",
    #     )


    # write_json(fpocket_data, os.path.join(PROJECT_ROOT, "FPOCKET_OUT.json")) 

    # d = read_fpocket_pdb_file("/home/david/aiengine/aiengine/cavities/FPOCKET_out/5ESM_out/pockets/pocket1_atm.pdb")

    # for k, v in d.items():
    #     print (k, v)

    print(run_fpocket_and_collate_single_target(
        target_identifier="2SHK",
        output_dir="fpocket_out_test",
        min_pocket_score=None,
        existing_pdb_filename="2SHK.pdb",

    ))
