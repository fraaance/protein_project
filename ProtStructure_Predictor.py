from pathlib import Path
import argparse
from ProtPred import ProtPredictor

num_stars = 70

parser = argparse.ArgumentParser(description="Protein distance matrix prediction")

parser.add_argument(
    "-db",
    required=True,
    help="Directory containing .cif files"
)

parser.add_argument(
    "-clusters",
    required=False,
    default=None,
    help="Optional existing MMseq2 cluster .tsv file"
)

parser.add_argument(
    "-m",
    required=True,
    help="'learn', 'no_learn' or path to .keras model"
)

parser.add_argument(
    "-o",
    required=False,
    default=None,
    help="Optional output-dir; default is dir within cif-DB directory"
)

args = parser.parse_args()

######### check Database
db_path = Path(args.db)

if not db_path.is_dir():
    raise ValueError("-db must be a directory")

cif_files = list(db_path.glob("*.cif"))

if not cif_files:
    raise ValueError("No .cif files found in database directory")

######## check Cluster
cluster_path = None
if args.clusters is None:
    print(f"No cluster file provided - MMseqs2 cluster")

else:
    cluster_path = Path(args.clusters)

    if not cluster_path.is_file():
        raise ValueError("Cluster .tsv does not exist")
    
    if cluster_path.suffix != ".tsv":
        raise ValueError("Cluster file must be a .tsv file")
    
######## check Model
keras_path = None
model_mode = args.m

if args.m == "learn":
    print(f" Training new model on .cif directory".center(num_stars, "*") )

elif args.m == "no_learn":
    print(f" Skipping training of new model, using pretrained model ".center(num_stars, "*"))
    keras_path = Path(__file__).resolve().parent / "pretrained_models/base_model.keras"

else:
    keras_path = Path(args.m)

    if not keras_path.is_file():
        raise ValueError("Model .keras file does not exist")
    if keras_path.suffix != ".keras":
        raise ValueError("Model must be .keras file")
    
    print(f" Loading existing model ".center(num_stars, "*"))

######## check Output Directory
if args.o is None:
    output_dir_path = db_path.parent / "output"
    output_dir_path.mkdir(parents=True, exist_ok=True)
    print(f"Output dir with path {output_dir_path} created")

else:
    output_dir_path = Path(args.o)

    if not output_dir_path.is_dir():
        raise ValueError(f"{output_dir_path} is not a directory")
    
print(db_path, cluster_path, keras_path, output_dir_path, model_mode)

predictor = ProtPredictor(
    db_path = db_path,
    cluster_path = cluster_path,
    keras_path = keras_path,
    output_dir_path = output_dir_path,
    model_mode = model_mode
)