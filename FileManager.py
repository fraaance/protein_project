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
        num_stars = 70

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
        self.num_stars = 70

        aa_3to1 = {
        "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
        "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
        "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"
        }
        aa_1to3 = {v: k for k, v in aa_3to1.items()}  

    def validate_paths(self, directory_path, output_dir_path):
        directory_path = Path(directory_path)
        output_dir = Path(output_dir_path)

        if not directory_path.is_dir():
            return False
        
        output_dir.mkdir(parents=True, exist_ok=True)
        return True

    # returns list with: seq_name, seq, atom_df
    def read_cif_file(self, directory_path):
        directory_path = Path(directory_path)
        output_list = []

        for file in directory_path.glob("*.cif"):
            print(f" Reading File {str(file)[-8:]} ".center(self.num_stars, "*"))
            cif = MMCIF2Dict(file)

            prot_df = pd.DataFrame(
                {
                    "entity_id": cif["_struct_ref.id"],
                    "seq_name": cif["_struct_ref.pdbx_db_accession"],
                    "seq": cif["_struct_ref.pdbx_seq_one_letter_code"]
                }
            )

            #remove sequences that are not an amino acid
            prot_df = self.remove_non_aa_sequences(prot_df)

            # store seq_name and sequence for later call up
            seq_dict = dict(prot_df[["seq_name", "seq"]].itertuples(index=False, name=None))

            # transform this dataframe from entity_id, seq_name, seq into
            # list: [entity_id, seq_name, pos, aminoacid]

            # transform linear amino acid seq into dataframe with every amino acid being a own row
            prot_list = [
                [row["entity_id"], row["seq_name"], pos, aa]
                for _, row in prot_df.iterrows()
                for pos, aa in enumerate(row["seq"], start=1)
            ]

            prot_res_df = pd.DataFrame(prot_list, columns=["entity_id", "seq_name", "pos", "seq_amino_acid"])
            
            atom_df = pd.DataFrame(
                {
                    "atom": cif["_atom_site.label_atom_id"],
                    "amino_acid": cif["_atom_site.label_comp_id"],
                    "entity_id": cif["_atom_site.label_entity_id"],
                    "seq_id": cif["_atom_site.label_seq_id"],
                    "x": cif["_atom_site.Cartn_x"],
                    "y": cif["_atom_site.Cartn_y"],
                    "z": cif["_atom_site.Cartn_z"],
                    # alt_id to remove alternative chains -> only the first is kept
                    "alt_id": cif["_atom_site.label_alt_id"]
                }
            )

            # only extract chiral centres and remove alternative positions
            atom_df_chiral = (atom_df[
                (atom_df["atom"] == "CA") & (atom_df["alt_id"].isin([".", "A"]))
            ].drop_duplicates(subset=["entity_id", "seq_id"], keep="first").copy()
            )

            # merge for entity = entity, seq_id = pos
            inval_seq = atom_df_chiral.loc[
                pd.to_numeric(atom_df_chiral["seq_id"], errors="coerce").isna(), "entity_id"
            ].unique()

            atom_df_chiral = atom_df_chiral[~atom_df_chiral["entity_id"].isin(inval_seq)].copy()
            atom_df_chiral["seq_id"] = atom_df_chiral["seq_id"].astype(int)
            
            # to prevent bug with non numerical seq_id's:
            atom_df_chiral["seq_id"] = atom_df_chiral["seq_id"].astype(int)

            merged_df = prot_res_df.merge(
                atom_df_chiral[
                    ["entity_id", "seq_id", "amino_acid", "x", "y", "z"]
                ],
                left_on=["entity_id", "pos"],
                right_on=["entity_id", "seq_id"],
                how="left"
            ).drop(columns=["seq_id", "amino_acid"])

            for (seq_name, protein_df) in merged_df.groupby("seq_name"):
                missing_values = protein_df["x"].isna().sum()
                if (len(protein_df) - missing_values) < len(protein_df) * 0.7:
                    print(f" sequence {seq_name} has less than 70% of coordinates; skipped ".center(self.num_stars, "*"))
                    continue

                protein_df = protein_df.drop(columns=["seq_name"])
                output_list.append((seq_name, seq_dict.get(seq_name),protein_df))

        return output_list

    # returns seq_name, dist_matrix
    def calculate_dist_matrix(self, seq_name, seq, atom_df):
        if len(seq) > 512 or len(atom_df) > 512:
            print(f" skipped {seq_name} ".center(self.num_stars, "*"))
            return None
        
        vectors = atom_df[["x", "y", "z"]].to_numpy(dtype=float)
        valid_residues = ~np.isnan(vectors).any(axis=1)

        dist_matrix = np.zeros((512, 512), dtype=np.float32)
        mask = np.zeros((512, 512), dtype=np.float32)
        indeces = np.where(valid_residues)[0]

        valid_coords = vectors[valid_residues]
        dists = cdist(valid_coords, valid_coords)

        dist_matrix[np.ix_(indeces, indeces)] = dists
        mask[np.ix_(indeces, indeces)] = 1.0

        return seq_name, dist_matrix, mask

    # saves dist_matrix as: seq_name.npy under output_dir/dist_matrices/seq_name.npy
    def save_dist_matrix(self, seq_name, dist_matrix, mask, output_dir_path):
        dist_matrix_path = output_dir_path / f"{seq_name}.npz"
        np.savez(dist_matrix_path, dist_matrix=dist_matrix, mask=mask)

    # to add: .cif.gz files -> for unzipping
    def validate_cif_path(self, file_path):
        path = Path(file_path)
        return (file_path.endswith(".cif")) and path.is_file()
    
    def remove_non_aa_sequences(self, df):
        alph = set("ACDEFGHIKLMNPQRSTVWY")

        df["seq"] = df["seq"].str.replace("\n", "", regex=False)
        df = df[df["seq"].apply(lambda seq: set(seq).issubset(alph))].copy()
        return df

    def translate_aa(self, item):
        item = item.upper()
        return self.aa_3to1.get(item, self.aa_1to3.get(item))
