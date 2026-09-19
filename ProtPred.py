import sys
import pandas as pd 
from pathlib import Path
import numpy as np
from FileManager import FastaManager, cifManager
from TrainValTestSplitter import TrainValTestSplitter
from ML_Model import forward_net, run_model
import keras
import matplotlib.pyplot as plt

class ProtPredictor:
    def __init__(self, db_path, cluster_path, keras_path, output_dir_path, model_mode):
        self.db_path = db_path
        self.cluster_path = cluster_path
        self.keras_path = keras_path
        self.output_dir_path = output_dir_path
        self.model_mode = model_mode
        self.num_stars = 70

        self.build_directory(db_path, cluster_path, output_dir_path)
        train_list, val_list, test_list = self.build_train_dataset(cluster_path, output_dir_path)
        
        if keras_path is None:
            distances_path = output_dir_path / "dist_matrices"
            print(f" training model ".center(self.num_stars, "*"))
            model = forward_net()
            save_model, fig = run_model(model, train_list, val_list, test_list, distances_path)
            save_model.save(f"{output_dir_path}/model.keras")
            fig.savefig(f"{output_dir_path}/training.png")
            plt.show()
            
    # directory with:
            # mmseq2_output -> cluster information
            # dist_matrices -> distance matrices stored as .npy
            # sequences.faa -> all valid seq of provided database
    def build_directory(self, db_path, cluster_path, output_dir_path):
        #fasta_manager = FastaManager()
        cif_manager = cifManager(db_path, output_dir_path)

        prot_db = []
        if cif_manager.validate_paths(db_path, output_dir_path):
            prot_db = cif_manager.read_cif_file(db_path)
    
        dist_matrix_dir = output_dir_path / "dist_matrices"
        dist_matrix_dir.mkdir(parents=True, exist_ok=True)

        # test, if the provided database has a sufficient size
        if len(prot_db) < 10:
            raise ValueError(f" the provided database is too small [{len(prot_db)}]; min: 10")
        
        for seq_name, seq, atom_df in prot_db:
            print(f" Calculating distance matrix for {seq_name} ")
            output = cif_manager.calculate_dist_matrix(seq_name, seq, atom_df)
            if output is not None:
                name, matrix, mask = output
                cif_manager.save_dist_matrix(name, matrix, mask, dist_matrix_dir)

        ###### store the sequences of the cif-directory in the sequences.faa file
        self.store_sequences(prot_db, output_dir_path)

    # split of dataset into 80/10/10, based on cluster profile
    def build_train_dataset(self, cluster_path, output_dir_path):
        print(f"splitting dataset into 80/10/20")
        splitter = TrainValTestSplitter(output_dir_path)
        return splitter.train_val_test_split(cluster_path)
    
    # stores the seq found in the database in output_dir_path/sequences.faa
    def store_sequences(self, prot_db, output_dir_path):
        fasta_manager = FastaManager()
        seq_list = []
        for seq_name, seq, atom_df in prot_db:
            seq_list.append((seq_name, seq))

        fasta_file_path = Path(output_dir_path / "sequences.faa")
        if fasta_manager.validate_path(fasta_file_path):
            print(f" Writing valid sequences in sequences.faa ".center(self.num_stars, "*"))
            fasta_manager.write_file(fasta_file_path, seq_list)
