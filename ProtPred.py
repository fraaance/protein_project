import sys
import pandas as pd 
from pathlib import Path
import numpy as np
from FileManager import FastaManager, cifManager
from TrainValTestSplitter import TrainValTestSplitter
import ML_Model

class ProtPredictor:
    def __init__(self, db_path, cluster_path, keras_path, output_dir_path, model_mode):
        self.db_path = db_path
        self.cluster_path = cluster_path
        self.keras_path = keras_path
        self.output_dir_path = output_dir_path
        self.model_mode = model_mode
        self.num_stars = 50

        self.build_directory(db_path, cluster_path, output_dir_path)
        train_list, val_list, test_list = self.build_train_dataset(output_dir_path)
        self.train_model(train_list, val_list, test_list)
        
            
    def build_directory(self, db_path, cluster_path, output_dir_path):
        fasta_manager = FastaManager()
        cif_manager = cifManager(db_path, output_dir_path)

        prot_db = []
        if cif_manager.validate_paths(db_path, output_dir_path):
            prot_db = cif_manager.read_cif_file(db_path)
    
        dist_matrix_dir = output_dir_path / "dist_matrices"
        dist_matrix_dir.mkdir(parents=True, exist_ok=True)

        for seq_name, seq, atom_df in prot_db:
            print(f" Calculating distance matrix for {seq_name} ")
            name, matrix = cif_manager.calculate_dist_matrix(seq_name, seq, atom_df)
            if name is not None:
                cif_manager.save_dist_matrix(name, matrix, dist_matrix_dir)

        ###### store the sequences of the cif-directory in the sequences.faa file
        seq_list = []
        for seq_name, seq, atom_df in prot_db:
            seq_list.append((seq_name, seq))

        fasta_file_path = Path(output_dir_path / "sequences.faa")
        if fasta_manager.validate_path(fasta_file_path):
            print(f" Writing valid sequences in sequences.faa ".center(self.num_stars, "*"))
            fasta_manager.write_file(fasta_file_path, seq_list)
            
            print(f" NUM of Training Seq: {len(fasta_manager.read_file(fasta_file_path))}".center(100, "O"))

    def build_train_dataset(self, output_dir_path):
        print(f"splitting dataset into 80/10/20")
        splitter = TrainValTestSplitter(output_dir_path)
        return splitter.train_val_test_split()
    
    def train_model(self, train_list, val_list, test_list):
        model = ML_Model.forward_net()
        ML_Model.run_model(model, train_list, val_list, test_list)
        
        x, y = test_list[0]
        x = np.expand_dims(x, axis=0)

        print(x.shape)
        print(y)
        print(model.predict(x))
