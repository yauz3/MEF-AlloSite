"""
https://github.com/akiyamalab/MEGADOCK


 Available options:
  -o filename    : set the output filename (default to "$R-$L.out")
  -O             : output docking detail files
  -N integer     : set the number of output predictions (default to 2000)
  -t integer     : set the number of predictions per each rotation (default to 1)
  -F integer     : set the number of FFT point (default to none)
  -v float       : set the voxel size (default to 1.2)
  -D             : set the 6 deg. (54000 angles) of rotational sampling
                   (default to none, 15 deg. (3600 angles) of rotational sampling)
  -r integer     : set the number of rotational sampling angles
                   (54000: 54000 angles, 1: 1 angles, 24: 24 angles, default to 3600 angles)
  -e float       : set the electrostatics term ratio (default to 1.0)
  -d float       : set the hydrophobic term ratio (default to 1.0)
  -a float       : set the rPSC receptor core penalty (default to -45.0)
  -b float       : set the rPSC ligand core penalty (default to 1.0)
  -f 1/2/3       : set function
                   (1: rPSC, 2: rPSC+Elec, 3: rPSC+Elec+RDE, default to 3)
  -h             : show this message


"""

import os
import subprocess
import sys
import shutil
import time
from pymol import cmd

def megadock(
        target, # should be paht/target
        ligand, # should be paht/target
        MEGA_DOCK_PATH,
        output_number,
        ):
    file_path_0 = os.path.dirname(sys.argv[0])  # To store files in the exact dictionaries.
    only_target_name=target.split("/")[-1].split(".")[0]
    only_ligand_name = ligand.split("/")[-1].split(".")[0]
    output_file="/home/yavuz/yavuz_proje/protac/outputs/Megadock_OUT"
    print("output_filee", output_file)
    try:
        """os.chdir(output_file)
        if os.path.isfile("Megadock_OUT") == False:
            os.mkdir(f"{output_file}/Megadock_OUT")"""
        os.chdir(f"{output_file}")
        os.mkdir(f"{output_file}/{only_target_name}")
    except:
        print("There is already output file in the dictionary.")
    ouput_file_path=f"{output_file}/{only_target_name}"

    os.chdir(MEGA_DOCK_PATH)
    os.system(f"./megadock-gpu -R {target} -L {ligand} -N {output_number} -o {ouput_file_path}/{only_target_name}-{only_ligand_name}.out")
    """os.system(
        f"./megadock-gpu -R {target} -L {ligand} -N {output_number}")
    print(f"{MEGA_DOCK_PATH}/{only_target_name}-{only_ligand_name}.out", f"{ouput_file_path}/{only_target_name}-{only_ligand_name}.out")
    shutil.move(f"{MEGA_DOCK_PATH}/{only_target_name}-{only_ligand_name}.out", f"{ouput_file_path}/{only_target_name}-{only_ligand_name}.out")
    shutil.move(f"{MEGA_DOCK_PATH}/{only_target_name}-{only_ligand_name}.detail", f"{ouput_file_path}/{only_target_name}-{only_ligand_name}.detail")
    shutil.move(f"{MEGA_DOCK_PATH}/{only_target_name}-{only_ligand_name}.csv", f"{ouput_file_path}/{only_target_name}-{only_ligand_name}.csv")"""
    shutil.copy(f"{ligand}", f"{ouput_file_path}/{only_ligand_name}.pdb")
    shutil.copy(f"{target}", f"{ouput_file_path}/{only_target_name}.pdb")

    for i in range(1,output_number):
        os.chdir(MEGA_DOCK_PATH)
        os.system(f"./decoygen {ouput_file_path}/ligand_{i}.pdb {ouput_file_path}/{only_ligand_name}.pdb {ouput_file_path}/{only_target_name}-{only_ligand_name}.out {i}")
        """os.chdir(ouput_file_path)
        cmd.load(f"ligand_{i}.pdb")
        cmd.load(f"{only_target_name}.pdb")
        cmd.save(f"complex_{i}.pdb")
        #cmd.remove(f"ligand_{i}") # to save as complex
        cmd.remove("all")"""
        # this can breake align
        #os.system(f"cat {only_target_name}.pdb ligand_{i}.pdb > complex_{i}.pdb")
    #print(f"./decoygen {ouput_file_path}/ligand.pdb {ouput_file_path}/{only_ligand_name}.pdb {ouput_file_path}/{only_target_name}-{only_ligand_name}.out 2")
if __name__ == "__main__":
    """MEGA_DOCK_PATH = "/home/yavuz/yavuz_proje/protac/bin/MEGADOCK-master"
    target = "/home/yavuz/yavuz_proje/protac/test/receptor_protac.pdb"  # protac structure should be here
    ligand = "/home/yavuz/yavuz_proje/protac/test/target.pdb"""
    megadock(
        target=target,
        ligand=ligand,
        MEGA_DOCK_PATH=MEGA_DOCK_PATH,
        output_number=output_number,
    )
