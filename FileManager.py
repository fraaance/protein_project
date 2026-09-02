from dataclasses import dataclass
import re
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist


# manager.read_file(file_path) -> returns a tuple list of all found sequences 
#                                 as a (seq_name, seq)
# manager.write_file(file_path, seq_list) -> checks if the file already exists and only 
#                                            adds new sequences from the list -> no duplicates 
class FastaManager:
    def __init__(self):
        num_stars = 5

    def read_file(self, file_path):
        if self.validate_path(file_path):
            path = Path(file_path)
            if not path.is_file():
                return None
            prot_list = []
            with open(file_path) as f:
                liste = re.split(">", f.read())[1:]
                
                print("Found sequences in FASTA File:")

                for e in liste: 
                    lines = e.splitlines()
                    header = lines[0]
                    seq = "".join(lines[1:])
                    print(f"Seq: {header}")
                    prot_list.append((header, seq))
            return prot_list
            

    def write_file(self, file_path, seq_list):
        path = file_path
        duplicate_seq = set()
        prot_list = self.read_file(file_path)
        if prot_list is not None:
            for p in prot_list:
                duplicate_seq.add(p[0])
        
        if self.validate_path(path) and self.validate_seq_list(seq_list):
            with open(path, "a") as f: 
                for name, seq in seq_list:
                    if name in duplicate_seq:
                        print(f"Sequence {name} already in FASTA File")
                        continue
                    if not self.validate_seq(512, seq):
                        continue

                    f.write(f">{name}\n") 
                    f.write(f"{seq}\n")

    def validate_path(self, file_path):
        path = Path(file_path)
        if path.is_file() or str(file_path).endswith(".faa"):
            return True
        return False

    def validate_seq_list(self, seq_list):
        if len(seq_list) == 0:
            return False
        return True
    
    def validate_seq(self, max_length, seq):
        if len(seq) <= 0 or len(seq) > max_length or not re.match(r"[A-Z]", seq):
            return False
        return True



# manager.cifManager(directory_path, output_dir_path) -> directory_path has to exist and all the files have to be .cif format
# manager.read_cif_file(file_path) -> returns list with all the amino acid seqences in the file, incl the information for distance matrix
#                                  -> stores the distance matrix of every sequence as a seq_name.npy file under the "dist_matrices" directory
    
class cifManager:
    def __init__(self, directory_path, output_dir_path):
        self.directory_path = directory_path
        self.output_dir_path = output_dir_path
        num_stars = 50

    def validate_paths(self, directory_path, output_dir_path):
        directory_path = Path(directory_path)
        output_dir = Path(output_dir_path)

        if not directory_path.is_dir():
            return False
        
        output_dir.mkdir(parents=True, exist_ok=True)
        return True

    # returns list with: seq_name, seq, atom_df
    def read_cif_file(self, directory_path):
        #if not self.validate_cif_path(file_path):
        #    return False

        directory_path = Path(directory_path)
        output_list = []

        for file in directory_path.glob("*.cif"):

            print(f" Reading File {str(file)[-8:]} ".center(50, "*"))

            cif = MMCIF2Dict(file)

            prot_df = pd.DataFrame(
                {
                    "entity_id": cif["_struct_ref.id"],
                    #"entity_id": cif["_entity_poly.entity_id"],
                    "seq_name": cif["_struct_ref.pdbx_db_accession"],
                    "seq": cif["_struct_ref.pdbx_seq_one_letter_code"]
                    #"seq": cif["_entity_poly.pdbx_seq_one_letter_code"],
                    #"chains": cif["_entity_poly.pdbx_strand_id"],
                    #"type": cif["_entity_poly.type"]
                }
            )

            print(prot_df)

            # remove non amino acid sequences
            prot_df = self.remove_non_aa_sequences(prot_df)
            
            atom_df = pd.DataFrame(
                {
                    "atom": cif["_atom_site.label_atom_id"],
                    "amino_acid": cif["_atom_site.label_comp_id"],
                    "entity_id": cif["_atom_site.label_entity_id"],
                    #"chain": cif["_atom_site.label_asym_id"],
                    "seq_id": cif["_atom_site.label_seq_id"],
                    "x": cif["_atom_site.Cartn_x"],
                    "y": cif["_atom_site.Cartn_y"],
                    "z": cif["_atom_site.Cartn_z"],
                    #"occupacy" : cif["_atom_site.occupancy"],
                    "alt_id": cif["_atom_site.label_alt_id"]
                }
            )
            
            atom_df = atom_df[(atom_df["atom"] == "CA")].copy()

            atom_df_chiral = atom_df[
                atom_df["alt_id"].isin([".", "A"])
            ].copy()
            atom_df_chiral = atom_df_chiral.drop_duplicates(
                subset=["entity_id", "seq_id"],
                keep="first"
            )

            for _, row in prot_df.iterrows():
                entity_id = row["entity_id"]
                seq = row["seq"]
                seq_name = row["seq_name"]

                entity_c_atoms = atom_df_chiral[atom_df_chiral["entity_id"] == str(entity_id)]

                # only keeps sequences with valid length 
                if len(seq) == len(entity_c_atoms):
                    output_list.append((seq_name, seq, entity_c_atoms))
                    #print(len(seq), entity_c_atoms.shape)
        #print(output_list)
        return output_list

    # returns seq_name, dist_matrix
    def calculate_dist_matrix(self, seq_name, seq, atom_df):
        id = atom_df["entity_id"]
        print(f"LENGHT OF SEQ: {len(seq)}".center(100, "#"))
        if len(seq) != len(atom_df) or len(seq) > 512:
            print(f" Skipped prot: {seq_name}; Length mismatch ".center(50, "*"))
            return None, None

        vectors = atom_df[["x", "y", "z"]].to_numpy(dtype=float)
        dist_matrix = cdist(vectors, vectors, metric='euclidean')

        # limit distance_matrix and allowed seq size to 512
        dummy = np.zeros((512, 512), dtype=float)
        n = len(seq)
        dummy[:n, :n] = dist_matrix
        
        return seq_name, dummy

    # saves dist_matrix as: seq_name.npy under output_dir/dist_matrices/seq_name.npy
    def save_dist_matrix(self, seq_name, dist_matrix, output_dir_path):
        dist_matrix_path = output_dir_path / f"{seq_name}.npy"
        np.save(dist_matrix_path, dist_matrix)

    # to add: .cif.gz files -> for unzipping
    def validate_cif_path(self, file_path):
        path = Path(file_path)
        return (file_path.endswith(".cif")) and path.is_file()
    
    def remove_non_aa_sequences(self, df):
        alph = set("ACDEFGHIKLMNPQRSTVWY")

        df["seq"] = df["seq"].str.replace("\n", "", regex=False)
        df = df[df["seq"].apply(lambda seq: set(seq).issubset(alph))].copy()

        return df


#class MMSeqManager:

#manager = FastaManager()

#lsite = manager.read_file("/Users/franzweisel/Downloads/project_output/sequences.faa")
#print(len(lsite))
#print(lsite)
#manager.write_file("/Users/franzweisel/Downloads/project_output/sequences2.faa", lsite)

#manager = cifManager("/Users/franzweisel/Downloads/testtt", "/Users/franzweisel/Downloads/testtt_output")
#path = Path("/Users/franzweisel/Downloads/testtt")
#x = manager.read_cif_file(path)
#print("end".center(50, "+"))
#print(x)
#for file in path.glob("*.cif"):
    #manager.read_cif_file(file)
