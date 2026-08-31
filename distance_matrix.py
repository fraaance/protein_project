import numpy as np
from scipy.spatial.distance import cdist
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd
import sys
from pathlib import Path
import re
from RNA_Reader import FastaReader 

def print_matrix(seq, dist_matrix):
    print("", *seq, sep="\t")
    for _, i in enumerate(seq):
        print(i, dist_matrix[_])

# extracts the name, sequence and the x,y,z-coordinates 
def read_cif_file(file_path):
    cif = MMCIF2Dict(file_path)
    atom_df = pd.DataFrame({
            "atom": cif["_atom_site.label_atom_id"],
            "amino_acid": cif["_atom_site.label_comp_id"],
            "entity_id": cif["_atom_site.label_entity_id"],
            "seq_id": cif["_atom_site.label_seq_id"],
            "x": cif["_atom_site.Cartn_x"],
            "y": cif["_atom_site.Cartn_y"],
            "z": cif["_atom_site.Cartn_z"],
            "occupacy" : cif["_atom_site.occupancy"]
    })

    # remove alternative chiral centres and only focus on the atom with highest occupacy
    atom_df_chiral = (atom_df[atom_df["atom"] == "CA"]
                     .sort_values("occupacy", ascending=False)
                     .drop_duplicates(subset=["atom", "amino_acid", "entity_id", "seq_id"], keep="first")
                     #.sort_values(["chain", "position"])
                    ).sort_index()
    seq = cif["_entity_poly.pdbx_seq_one_letter_code"][0].replace("\n", "")
    seq_name = cif["_entry.id"][0].replace("\n", "")
    
    return seq_name, seq, atom_df_chiral

# calculates and returns the distant matrix based on the extracted chiral centres from the cif file
def calculate_dist(seq, atom_df):
    vectors = atom_df[["x", "y", "z"]].to_numpy(dtype=float)
    dist_matrix = cdist(vectors, vectors, metric='euclidean')
    #print_matrix(seq, dist_matrix)
    #print(dist_matrix)
    return dist_matrix

# creates the fasta file for the data set structure
def print_in_fasta(seq_name, seq, file_path):
    with open(file_path, "a") as f:
       f.write(f">{seq_name}\n") 
       f.write(f"{seq}\n")

def validate_paths(dictionary_path, output_dir):
    dictionary_path = Path(dictionary_path)
    output_dir = Path(output_dir)

    if not dictionary_path.is_dir():
        return False
    
    output_dir.mkdir(parents=True, exist_ok=True)
    return True

# only accepts sequences with length > 0 and capital letter
def validate_seq(seq):
    if len(seq) <= 0 or not re.match(r"[A-Z]", seq):
        return False
    return True


# opens the directory and extracts ever cif-file, calculates the seq dist matrix
# and stores the seq_name and seq in a common fasta-file
def build_dataset(dir_name, output_dir):
    directory = Path(dir_name)
    output_dir = Path(output_dir)

    distance_dir = output_dir / "distances"
    distance_dir.mkdir(parents=True, exist_ok=True)
    fasta_file = output_dir / "sequences.faa"

    duplic_sequences = set()
    reader = FastaReader(str(fasta_file))
    for prot in reader.prot_list:
        duplic_sequences.add(prot.header)
    # here: add on, that RNA_Reader() gets called, to read the Fasta and adds all seq to the duplic_seq set

    for file in directory.glob("*.cif"):
        seq_name, seq, atom_df = read_cif_file(file)
        if not validate_seq(seq=seq):
            continue
        
        if seq_name in duplic_sequences:
            continue
        duplic_sequences.add(seq_name)

        dist_matrix = calculate_dist(seq=seq, atom_df=atom_df)
        #print(dist_matrix)

        # Store seq_name and seq in overall FASTA-file
        print_in_fasta(seq_name, seq, fasta_file)

        # Store distant-matrix as numpy file
        dist_path = distance_dir / f"{seq_name}.npy"
        np.save(dist_path, dist_matrix) 

############################################################
####################### program start ######################
############################################################

# version one: create a new dataset
# distance_matrix.py -d dictionary_path -f fasta_file_path.faa 
dictionary_path = sys.argv[1]
output_dir = sys.argv[2]
if validate_paths(dictionary_path=dictionary_path, output_dir=output_dir):
    build_dataset(dictionary_path, output_dir)
print(f"FASTA-file stored at {output_dir}/sequences.faa")
print(f"Distance Matrices stored in {output_dir}/distances/")
        
# version two: data set already exists -> only calculate the distant matrices
# distance_matrix.py -d dictionary_path.faa