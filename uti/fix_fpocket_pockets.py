import os
from Bio.PDB import PDBParser
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.PDBIO import Select
from pymol import cmd

class ChainResidueSelect(Select):
    def __init__(
        self,
        chain_residues,
        ):
        """Initialiser

        Parameters
        ----------
        chain_residues : set
            Set of chain_id, residue pairs to keep.
        """
        super(ChainResidueSelect, self).__init__()
        if isinstance(chain_residues, list):
            chain_residues = set(chain_residues)
        self.chain_residues = chain_residues

    def accept_residue(
        self,
        residue,
        ):
        """Custom select residue function.

        Parameters
        ----------
        residue : Residue
            Biopython residue object to check

        Returns
        -------
        bool
            Flag to keep residue
        """

        chain_id = residue.get_parent().id
        _, residue_id, _ = residue.id

        return (chain_id, residue_id) in self.chain_residues


def neigh_residue(pocket_name,protein_input,input_path_pocket,input_path_protein,output_path):
    # sometimes pockets has HETATOM; therefore, we clean it
    os.chdir(input_path_pocket)
    cmd.reinitialize()
    cmd.load(f"{pocket_name}.pdb")
    cmd.remove("het")
    cmd.save(f"{pocket_name}.pdb")
    cmd.remove("all")

    # we can add couple more residue to complete pocket. However, we need to know that residue exist on the structure
    # it will return all residue
    os.chdir(input_path_protein)
    # set target structure
    parser = PDBParser()
    residue_check = parser.get_structure(
        protein_input,  # can be any string
        f"{protein_input}.pdb",  # pdb filename of full target structure
    )
    all_residue_list=[]
    for chains in residue_check:
        for chain in chains:
            for residue in chain:
                residue_info=residue.get_full_id()
                all_residue_list.append(residue_info[3][1])
    os.chdir(input_path_pocket)
    # initialise
    parser = PDBParser()
    # load pocket structure
    pocket_structure = parser.get_structure(pocket_name, pocket_name + ".pdb")
    # build list of residue IDs
    residue_ids = []
    for pocket_residue in pocket_structure.get_residues():
            chain_id = pocket_residue.get_parent().id
            _, residue_id, _ = pocket_residue.id
            residue_ids.append((chain_id, int(residue_id)))
            # The pocket size can be adjustable to have more informative sequence based features
            # However, we did not use the approach to keep focus on model improvement.
            """if (int(residue_id)-1) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id)-1))
            if (int(residue_id) +1) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id) + 1))"""
            """if (int(residue_id) - 2) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id) - 2))
            if (int(residue_id) + 2) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id) + 2))
            if (int(residue_id) + 3) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id) + 3))
            if (int(residue_id) - 3) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id) - 3))
            if (int(residue_id) - 4) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id) - 4))
            if (int(residue_id) + 4) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id) + 4))
            if (int(residue_id) - 5) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id) - 5))
            if (int(residue_id) + 5) in all_residue_list:
                residue_ids.append((chain_id, int(residue_id) + 5))
    """
    os.chdir(input_path_protein)
    # set target structure
    parser = PDBParser()
    target_structure = parser.get_structure(
         protein_input, # can be any string
         f"{protein_input}.pdb", # pdb filename of full target structure
    )
    # write pocket
    io = PDBIO()
    io.set_structure(target_structure)
    # save only residues listed in selected_pocket_residue_is
    # target pocket_pdb_filename is the output filename for the pocket
    os.chdir(output_path)
    io.save(f"{pocket_name}_fixed.pdb",ChainResidueSelect(residue_ids))

if __name__ == "__main__":
    """pocket_name="pocket1_atm"
    protein_input="1A3W"
    input_path_pocket="/home/yavuz/yavuz_proje/allosteric_binding_site/passer_test/PASSer2.0-main/data/pockets/1A3W_out/pockets"
    input_path_protein="/home/yavuz/yavuz_proje/allosteric_binding_site/passer_test/PASSer2.0-main/data/pdbs"
    output_path="/home/yavuz/yavuz_proje/allosteric_binding_site/passer_test/PASSer2.0-main/data"""
    neigh_residue(
    pocket_name=pocket_name,
    protein_input=protein_input,
    input_path_pocket=input_path_pocket,
    input_path_protein=input_path_protein,
    output_path=output_path)