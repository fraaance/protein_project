# ProteinPredictor

Keras based model predicting distance matrix based on a single amino acid sequence.

> **Status:** work in progress; the current implementations focuses on dataset preprocessing and model training; next step: ... positional connections

## Pipeline
- Read protein structures form PDB **mmCif** (`.cif`) files
- Extract amino-acid sequences and C_a coordinates 
- Filter sequences: 
    - canonical amino acids only
    - max. length: **512 AA**
    - at least **70% available C_a coordinates**
- Calculate C_a - C_a distance matrices
- Pad matrices to **512 x 512** 
- Store distance matrix with validity mask as `.npz` file
- Cluster sequences with **MMseqs2** for **80/10/10 - Train/Val/Test Split**
- Encoding AA as integers (`10-20`, padding `0`)
- Train and evaluate Keras model

## Current Model Design
- Input: integer-encoded sequence `(512, )`
- Embedding `21 -> 256`
- Two dense layers (ReLU)
- Output: predicted distance Matrix `(512, 512)`

## Requirements: 
- python 3
- Tensorflow / Keras
- Numpy
- pandas
- SciPy
- scikit-learn
- Matplotlib
- MMseqs2

## Usage
    ```bash 
    python ProtStructure_Predictor.py \
    -db <cif_directory> \
    -m learn \
    [-clusters <clusters.tsv>] \
    [-o <output_directory>]
    ```

## Arguments
- `-db` — directory containing `.cif` files
- `-m learn` — preprocess data and train a new model
- `-clusters` — optional existing MMseqs2 cluster TSV
- `-o` — optional output directory
- `-m no_learn` / `-m <model.keras>`(still under development)

## Structure of Project
- ProtStructure_Predictor.py -> argument handling
- ProtPred.py -> pipeline orchestration
- FileManager.py -> mmCIF/ FASTA parsing and distance matrices
- TrainValTestSplitter.py -> MMseqs2 clustering, splitting, sequence encoding
- ML_Model.py -> Keras model, streaming datasets, losses and training