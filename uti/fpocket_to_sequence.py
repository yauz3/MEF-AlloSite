
from Bio.PDB import PDBParser


# map three letter AA codes to one letter AA code
# taken from PyBioMed.PyGetMol.GetProtein
AA3_TO_AA1 = [
    ("ALA", "A"),
    ("CYS", "C"),
    ("ASP", "D"),
    ("GLU", "E"),
    ("PHE", "F"),
    ("GLY", "G"),
    ("HIS", "H"),
    ("HSE", "H"),
    ("HSD", "H"),
    ("ILE", "I"),
    ("LYS", "K"),
    ("LEU", "L"),
    ("MET", "M"),
    ("MSE", "M"),
    ("ASN", "N"),
    ("PRO", "P"),
    ("GLN", "Q"),
    ("ARG", "R"),
    ("SER", "S"),
    ("THR", "T"),
    ("VAL", "V"),
    ("TRP", "W"),
    ("TYR", "Y"),
]
# convert to dict
AA3_TO_AA1 = {
    k: v 
    for k, v in AA3_TO_AA1
}

def fpocket_pdb_file_to_sequence(
    pocket_pdb_filename: str,
    ):

    print ("Reading FPocket pocket PDB file", pocket_pdb_filename, "and converting to seqeunce")

    # initialise PDBParser
    parser = PDBParser()
    # load in pdb file
    pocket_structure = parser.get_structure(
        id="my-pocket",
        file=pocket_pdb_filename,
    )

    # build list of residue IDs as (chain_id, resname, residue_id) pairs 
    residue_ids = []
    for pocket_residue in pocket_structure.get_residues():
        chain_id = pocket_residue.get_parent().id
        _, residue_id, _ = pocket_residue.id
        resname = pocket_residue.resname

        # (chain_id, resname, residue_id) tuple
        residue_ids.append((chain_id, resname, residue_id))

    # use residue_ids to define sequence (should be ordered by residue_id)
    chain_to_residue = {}
    for chain_id, resname, residue_id in residue_ids:
        if chain_id not in chain_to_residue:
            chain_to_residue[chain_id] = []
        chain_to_residue[chain_id].append((resname, residue_id))

    # map chain ID to sequence string (sorted by residue_id)
    chain_to_sequence = []
    for chain_id, chain_residues in chain_to_residue.items():
        # sort by residue id
        chain_residues = sorted(chain_residues, key=lambda residue: residue[1])
        # add to chain_to_sequence
        chain_to_sequence = "".join(
            (
                # convert 3-character amino acid string to single character
                AA3_TO_AA1[resname] if resname in AA3_TO_AA1 else "X" # also taken from PyBioMed.PyGetMol.GetProtein
                for resname, residue_id in chain_residues
            )
        )

    return chain_to_sequence

if __name__ == "__main__":

    #pocket_pdb_filename = "fpocket_out/cavities/FPOCKET_out/3ETO_out/pockets/pocket1_atm.pdb"

    chain_to_sequence = fpocket_pdb_file_to_sequence(pocket_pdb_filename)

    # prints dict (chain_id -> sequence string)
    #print (chain_to_sequence)

    # convert to single string for all chains (sorted by chain_id)
    #sequence_all_chains = "".join((chain_to_sequence[chain] for chain in sorted(chain_to_sequence)))
    #print (sequence_all_chains)