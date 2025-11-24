import os
from ftplib import FTP
import csv
from FPSim2.io import create_db_file
from rdkit import Chem
from rdkit import __version__ as rdkit_version
from rdkit.Chem.Scaffolds import MurckoScaffold

print(rdkit_version)
assert rdkit_version == "2025.03.3", "Please use RDKit 2025.3.3"

root = os.path.dirname(os.path.abspath(__file__))

dest_dir = os.path.join(root, "..", "data")


print("Indexing ChEMBL...")
with open(os.path.join(root, "..", "data", "chembl_35_chemreps.txt"), "r") as f: #file to be downloaded online
    reader = csv.reader(f, delimiter="\t")
    next(reader)
    smiles_list = []
    for r in reader:
        smiles_list += [r[1]]

smiles_list = sorted(set(smiles_list))

mols = [[smiles, i] for i, smiles in enumerate(smiles_list)]

print("Creating a database file with Morgan fingerprints")

create_db_file(
    mols_source=mols,
    filename=os.path.join(dest_dir, "fpsim2_database_chembl.h5"),
    mol_format='smiles',
    fp_type='Morgan',
    fp_params={'radius': 2, 'fpSize': 1024}
)

with open(os.path.join(dest_dir, "fpsim2_database_chembl_smiles.csv"), "w") as f:
    writer = csv.writer(f)
    writer.writerow(["smiles", "index"])
    for smiles, i in mols:
        writer.writerow([smiles, i])