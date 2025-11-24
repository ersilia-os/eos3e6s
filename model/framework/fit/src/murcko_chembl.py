from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
import pandas as pd
import os
from tqdm import tqdm

root = os.path.dirname(os.path.abspath(__file__))
dest_dir = os.path.join(root, "..", "data")
df = pd.read_csv(os.path.join(dest_dir, "fpsim2_database_chembl_smiles.csv"))
scaff_data = []

def scaffold_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    scaff_mol = MurckoScaffold.GetScaffoldForMol(mol)
    if scaff_mol is None or scaff_mol.GetNumAtoms() == 0:
        return None
    scaff_smi = Chem.MolToSmiles(scaff_mol, isomericSmiles=False)
    return scaff_smi

scaffolds = []
for smi in tqdm(df["smiles"].tolist()):
    scaffold = scaffold_smiles(smi)
    scaffolds += [scaffold]
df["scaffold"] = scaffolds

df.to_csv(os.path.join(dest_dir, "chembl_scaffolds.csv"), index=False)


