import os
import shutil
from functools import partial

from multiprocessing.pool import Pool

import numpy as np
import pandas as pd

PROJECT_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(__file__), 
    os.pardir,
    os.pardir))

if __name__ == "__main__":

    import sys
    import os.path

    sys.path.insert(1, PROJECT_ROOT)


from utils.sys_utils import execute_system_command
from utils.io.io_utils import copy_file, delete_directory, delete_file, write_json, load_json
from utils.molecules.pdb_utils import download_pdb_structure_using_pdb_fetch


P2RANK_OUT = os.path.join(
    # PROJECT_ROOT,
    "cavities", 
    "PRANK_out",
    )
    
P2RANK_DIR = os.path.join(
    PROJECT_ROOT,
    "bin",
    "location_analysis",
    "p2rank_2.3.1",
    )
P2RANK_LOCATION = os.path.join(P2RANK_DIR, "prank")

POCKET_ANALYSIS_N_PROC = int(os.environ.get("POCKET_ANALYSIS_N_PROC", default=1))

def call_p2rank(
    input_filename, 
    output_dir,
    generate_visualisations=0,
    n_threads=1,
    main_job_id=None,
    verbose: bool = False,
    ):

    if output_dir is None:
        stem, ext = os.path.splitext(input_filename)
        output_dir = os.path.basename(stem) + "_prank_out"

    if verbose:
        print ("CALLING p2rank ON FILE", input_filename)
        print ("Outputting to directory", output_dir)

    output_filename = os.path.join(output_dir, os.path.basename(input_filename) + "_predictions.csv")

    cmd = f"{P2RANK_LOCATION} predict -f {input_filename} -o {output_dir} -visualizations {generate_visualisations} -threads {n_threads}"
    try:
        return_code = execute_system_command(cmd, main_job_id=main_job_id, verbose=verbose)
    except Exception as e:
        print ("Error with p2rank", e)
        pass

    return output_filename

def run_p2rank_and_collate_single_target(
    target_identifier,
    output_dir,
    min_pocket_score,
    existing_pdb_filename,
    generate_visualisations=0,
    delete_output_dir: bool = True,
    verbose: bool = False,
    ):
   

    prank_target_output_dir = os.path.join(
        output_dir,
        P2RANK_OUT, 
        target_identifier)
    os.makedirs(prank_target_output_dir, exist_ok=True)

    p2rank_prediction_filename = os.path.join(
        prank_target_output_dir, 
        f"{target_identifier}_predictions.csv",
        )

    if not os.path.isfile(p2rank_prediction_filename):
        # run p2rank
        pdb_filename = os.path.join(prank_target_output_dir, f"{target_identifier}.pdb")
        
        if not os.path.exists(pdb_filename):
            if existing_pdb_filename is not None:
                if existing_pdb_filename != pdb_filename:
                    copy_file(existing_pdb_filename, pdb_filename)
            else:
                pdb_filename = download_pdb_structure_using_pdb_fetch(
                    pdb_id=target_identifier,
                    pdb_filename=pdb_filename,
                    verbose=verbose,
                )

        prank_output_filename = call_p2rank(
            pdb_filename,
            output_dir=prank_target_output_dir,
            generate_visualisations=generate_visualisations,
            verbose=verbose,
            )
        if os.path.exists(prank_output_filename):
            os.rename(
                prank_output_filename, 
                p2rank_prediction_filename,
                )

        delete_file(pdb_filename, verbose=verbose)

    p2rank_target_data = {}

    p2rank_prediction_df = None
    
    if os.path.exists(p2rank_prediction_filename):
        if verbose:
            print("Reading P2Rank predictions from", p2rank_prediction_filename)
        p2rank_prediction_df = pd.read_csv(p2rank_prediction_filename, index_col=0)
        # remove whitespace from columns TODO: WHERE DID THAT COME FROM?
        p2rank_prediction_df.columns = map(lambda col: col.strip(), p2rank_prediction_df.columns)

    if p2rank_prediction_df is not None:
        
        if verbose:
            print ("P2Rank identfied", p2rank_prediction_df.shape[0], "pockets")

        for pocket_number, row in p2rank_prediction_df.iterrows():

            # remove whitespace
            pocket_number = pocket_number.strip()
            # remove `pocket`
            pocket_number = pocket_number.split("pocket")[1]

            pocket_score = row["score"]

            # pocket_score must be bigger than min_pocket_score 
            if min_pocket_score is not None and pocket_score < min_pocket_score:
                continue

            center_x = row["center_x"]
            center_y = row["center_y"]
            center_z = row["center_z"]

            size_x = None # TODO
            size_y = None 
            size_z = None 

            sas_points = row["sas_points"]
            surf_atoms = row["surf_atoms"]

            p2rank_target_data[pocket_number] = {
                "center_x": center_x,
                "center_y": center_y,
                "center_z": center_z,
                "size_x": size_x,
                "size_y": size_y,
                "size_z": size_z,
                "pocket_score": pocket_score,
                "sas_points": sas_points,
                "surf_atoms": surf_atoms,
            }

    if len(p2rank_target_data) == 0:
        # no identified pockets
        if verbose:
            print ("P2Rank failed to identify pockets for target", target_identifier)

        p2rank_target_data[None] = {
            "center_x": None,
            "center_y":None,
            "center_z": None,
            "size_x": None,
            "size_y":None,
            "size_z": None,
            "pocket_score": None,
            "sas_points": None,
            "surf_atoms": None,
        }

    if delete_output_dir:
        delete_directory(prank_target_output_dir, verbose=verbose)


    # return {target_identifier: p2rank_target_data}
    return p2rank_target_data


# def collate_p2rank_data(
#     targets,
#     output_filename,
#     output_dir,
#     min_pocket_score=0,
#     verbose: bool = True,
#     ):

#     if os.path.exists(output_filename):

#         print (output_filename, "already exists -- loading it")
#         p2rank_data = load_json(output_filename)

#     else:

#         print ("Running P2Rank for targets", targets, "using minimum score threshold:", min_pocket_score)
#         p2rank_data = {}

#         with Pool(processes=POCKET_ANALYSIS_N_PROC) as p:
#             p2rank_data_all_targets = p.map(
#                 partial(
#                     run_p2rank_and_collate_single_target,
#                     output_dir=output_dir,
#                     min_pocket_score=min_pocket_score,
#                     existing_pdb_filename=None,
#                     verbose=verbose,
#                 ),
#                 targets,
#             )

#             for p2rank_data_target in p2rank_data_all_targets:
#                 p2rank_data.update(p2rank_data_target)

#         print ("Writing P2rank data to", output_filename)
#         write_json(p2rank_data, output_filename)

#     return p2rank_data

if __name__ == "__main__":



    # targets = [
    #     # "6NNA", 
    #     # "2OQK", 
    #     "1HTP",
    # ]

    # target = "test_target"
    # min_pocket_score = 0
    # existing_pdb_filename = f"{DATA_ROOT_DIR}/edock_benchmark/target_structures/coach/3erkA_BS01_SB4.pdb"

    output_filename = "P2RANK_OUT.json"


    # target = "5LXI"
    # target = "2WZS"
    target = "5HJB"
    # existing_pdb_filename = "5LXI.pdb"
    # existing_pdb_filename = "2WZS.pdb"
    existing_pdb_filename = "5HJB.pdb"

    # from utils.pdb_utils import select_chain

    # existing_pdb_filename = select_chain(
    #     pdb_id=target,
    #     pdb_filename=existing_pdb_filename,
    #     chain_id={
    #         "C", 
    #         "E",
    #     }
    # )

    output_dir = "p2rank_test"


    # from ai_blind_docking.blind_docking_utils.foldx_utils import run_foldx

    # existing_pdb_filename = run_foldx(
    #     existing_pdb_filename,
    #     output_dir=output_dir,
    # )



    p2rank_data = run_p2rank_and_collate_single_target(
        target,
        output_dir=output_dir,
        min_pocket_score=0,
        existing_pdb_filename=existing_pdb_filename,
    )

    write_json(p2rank_data, output_filename)
