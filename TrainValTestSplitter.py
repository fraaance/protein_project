import subprocess
from pathlib import Path
import pandas as pd
from sklearn.model_selection import train_test_split
import numpy as np
from FileManager import FastaManager



class TrainValTestSplitter:
    def __init__(self, dir_path):
        self.dir_path = Path(dir_path)

    def train_val_test_split(self):
        directory = self.dir_path

        fasta_path = directory / "sequences.faa"
        output_dir = directory / "mmseq2_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        clusters = output_dir / "clusters"
        tmp = output_dir / "mmseq2_tmp"

        # Run MMSEQ2 to find clusters and prevent to similar train-val-test sets
        subprocess.run([
            "mmseqs",
            "easy-cluster",
            str(fasta_path),
            str(clusters),
            str(tmp),
            "--min-seq-id",
            "0.2",
            "-c",
            "0.8"
        ], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

        # Extract found clusters
        cluster_df = pd.read_csv(str(output_dir / "clusters_cluster.tsv"),
                            sep="\t",
                            names=["cluster", "prot_id"])
        #print(clusters)

        clusters = cluster_df["cluster"].unique()

        # train, val, test split of clusters
        train_cl, temp_cl = train_test_split(
            clusters, test_size=0.2, random_state=42
        )

        val_cl, test_cl = train_test_split(
            temp_cl, test_size=0.5, random_state=42
        ) 


        # Combine with Distance Matrices as tupel-list
        train_set = cluster_df[cluster_df["cluster"].isin(train_cl)]["prot_id"].tolist()
        test_set = cluster_df[cluster_df["cluster"].isin(test_cl)]["prot_id"].tolist()
        val_set = cluster_df[cluster_df["cluster"].isin(val_cl)]["prot_id"].tolist()

        # get the sequences, numerical encoded
        reader = FastaManager()
        prot_list = reader.read_file(fasta_path)

        prot_dict = {header: seq for header, seq in prot_list}

        distances_path = directory / "dist_matrices"
        train_list = [(self.aa_encoder(prot_dict[el]), np.load(f"{distances_path}/{el}.npy")) for el in train_set]
        val_list = [(self.aa_encoder(prot_dict[el]), np.load(f"{distances_path}/{el}.npy")) for el in val_set]
        test_list = [(self.aa_encoder(prot_dict[el]), np.load(f"{distances_path}/{el}.npy")) for el in test_set]
        
        return train_list, val_list, test_list
    
    def aa_encoder(self, seq):
        aa_to_int = {
        "A": 1, "C": 2, "D": 3,
        "E": 4, "F": 5, "G": 6,
        "H": 7, "I": 8, "K": 9,
        "L": 10, "M": 11, "N": 12,
        "P": 13, "Q": 14, "R": 15,
        "S": 16, "T": 17, "V": 18,
        "W": 19, "Y": 20
        }
        encod_seq = [aa_to_int[c] for c in seq]

        encod_seq += [0] * (512 - len(encod_seq))

        return np.array(encod_seq, dtype=np.int32)
 #train_set, val_set, test_set = train_val_test_split("/Users/franzweisel/Downloads/project_output")
