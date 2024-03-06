import glob
import os
import subprocess
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Align
import pickle
from Bio import SeqIO
from Bio.PDB import PDBParser




def sequence_similarty(seq_1,seq_2):
    alignments = pairwise2.align.globalxx(seq_1, seq_2)
    for alignment in alignments:
        #print(format_alignment(*alignment).split("Score=")[1])
        similarty_score=int(format_alignment(*alignment).split("Score=")[1].replace("\n",""))
    return similarty_score


def get_sequence(pdbid):
    import requests
    data = requests.get(f'https://www.ebi.ac.uk/pdbe/api/pdb/entry/molecules/{pdbid}').json()[pdbid.lower()]
    #print(data[0]['sequence'])
    return data[0]['sequence']

PDBIND_PATH="/home/yavuz/yavuz_proje/allosteric_binding_site/database/PDBbind_v2020_other_PL/v2020-other-PL"
def sequence_dictionary_pdbind():
    pdbid_list=os.listdir(PDBIND_PATH)
    pdb_ligand_dictionary={}
    for pdbid in pdbid_list:
        if "index" and "readme" not in pdbid:
            sequence=get_sequence(pdbid)
            pdb_ligand_dictionary[pdbid]=str(sequence)
            print(pdbid)
    print(pdb_ligand_dictionary)
    a_file = open("sequence_dictionary.pkl", "wb")
    pickle.dump(pdb_ligand_dictionary, a_file)
    a_file.close()
#sequence_dictionary_pdbind()

def sequnce_from_PDBParser(pdbid):
    # You can use a dict to convert three letter code to one letter code
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}


    # run parser
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdbid)

    # iterate each model, chain, and residue
    # printing out the sequence for each chain

    sequence = ""
    for model in structure:
        for chain in model:
            seq = []
            for residue in chain:
                try: # sometimes DNA or RNA may exist in the structure with ATOM label. To valid Amino Acid is provided above
                    seq.append(d3to1[residue.resname])
                except:
                    continue
            sequence = sequence + (''.join(seq))
    #print(sequence)
    return sequence

def get_sequnce_from_pdb2fast(pdbid):
    os.chdir("/home/yavuz/yavuz_proje/allosteric_binding_site/bin")
    output=(str(subprocess.check_output(f"./pdb2fasta {pdbid}",shell=True)).replace("'","")).split(r"\n")
    sequnce=""
    if len(output) > 2:
        for i in range(len(output)):
            if ">" not in output[i]:
                sequnce=sequnce+str(output[i])
    return sequnce

def sequnce_from_pdb_file(pdbid):
    for record in SeqIO.parse(pdbid, "pdb-atom"):
        print(record.seq)


if __name__ == "__main__":
    """ seq_1="MLMPKKERKVEGDEVIRVPLPEGNQLFGVVEQALGAGWMDVRCEDGKIRRCRIPGKLRRRVWIRVGDLVIVQPWPVQSDKRGDIVYRYTQTQVD"
        seq_2="MLMPKKERKVEGDEVIRVPLPEGNQLFGVVEQALGAGWMDVRCEDGKIRRCRIPGKLRRRVWIRVGDLVIVQPWPVQSDKRGDIVYRYTQTQVDWLLRKGKITQEFLTGGSLLVE"""
    sequence_similarty(
    seq_1=seq_1,
    seq_2=seq_2,
    )
    sequence_dictionary_pdbind(

    )
    get_sequence(
        pdbid=pdbid
    )
    get_sequnce_from_pdb2fast(
        pdbid=pdbid
    )
    sequnce_from_PDBParser(
        pdbid=pdbid
    )