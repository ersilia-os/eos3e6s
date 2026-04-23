import pandas as pd
import os
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors

root = os.path.dirname(os.path.abspath(__file__))

dest_dir = os.path.join(root, "..", "data")
df = pd.read_csv(os.path.join(dest_dir, "fpsim2_database_chembl_smiles.csv"))

prop_dict = {
    "smiles": [],
    "MW": [],
    "LogP": [],
    "Hba": [],
    "Hbd": [],
    "rb": [],
    "HeavyAtoms": []
}

def calc_physchem(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print("Mol None")
        prop_dict["smiles"] += [smiles]
        prop_dict["MW"]+= [None]
        prop_dict["LogP"]+= [None]
        prop_dict["Hba"]+= [None]
        prop_dict["Hbd"]+= [None]
        prop_dict["rb"]+= [None]
        prop_dict["HeavyAtoms"]+= [None]
    else:
        prop_dict["smiles"] += [smiles]
        prop_dict["MW"]+= [Descriptors.MolWt(mol)]
        prop_dict["LogP"]+= [Crippen.MolLogP(mol)]
        prop_dict["Hba"]+= [rdMolDescriptors.CalcNumLipinskiHBA(mol)]
        prop_dict["Hbd"]+= [rdMolDescriptors.CalcNumLipinskiHBD(mol)]
        prop_dict["rb"]+= [rdMolDescriptors.CalcNumRotatableBonds(mol)]
        prop_dict["HeavyAtoms"]+= [Descriptors.HeavyAtomCount(mol)]

for smi in tqdm(df["smiles"].tolist()):
    calc_physchem(smi)
df_props = pd.DataFrame(prop_dict)
df_props.to_csv("../data/chembl_props.csv", index=False)