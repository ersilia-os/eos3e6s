import os
from FPSim2 import FPSim2Engine
import h5py
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


root = os.path.dirname(os.path.abspath(__file__))

class ChemblLowSimilarity(object):

    def __init__(self, threshold=0.25):
        self.threshold = threshold
        self.fp_database = os.path.join(root, "..", "fit", "data", "fpsim2_database_chembl.h5")
        self.fp_smiles = os.path.join(root, "..", "fit", "data", "fpsim2_database_chembl_smiles.csv")
        self.fpe = FPSim2Engine(self.fp_database, in_memory_fps=False)

    def get_low_similarity(self, smiles):
        high_sim = self.fpe.on_disk_similarity(
            smiles,
            threshold=self.threshold,
            metric="tanimoto",
            n_workers=1
        )
        high_ids = {res[0] for res in high_sim}
        input_smiles = pd.read_csv(self.fp_smiles)
        low_smiles = input_smiles[~input_smiles["index"].isin(high_ids)]
        low_smiles = low_smiles["smiles"].tolist()
        return low_smiles


class ChemblProps(object):
    def __init__(self):
        df = pd.read_csv(os.path.join(root, "..", "fit", "data", "chembl_props.csv"))

        self.smiles_to_props = {}
        for row in df.itertuples(index=False):
            smi = str(row.smiles)
            if pd.isna(row.MW) or pd.isna(row.LogP) or pd.isna(row.Hba) \
               or pd.isna(row.Hbd) or pd.isna(row.rb):
                continue
            self.smiles_to_props[smi] = {
                "MW": float(row.MW),
                "LogP": float(row.LogP),
                "Hba": int(row.Hba),
                "Hbd": int(row.Hbd),
                "rb": int(row.rb),
            }

    def _calc_physchem(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        self.prop_dict = {}
        self.prop_dict["MW"] = Descriptors.MolWt(mol)
        self.prop_dict["LogP"] = Crippen.MolLogP(mol)
        self.prop_dict["Hba"] = rdMolDescriptors.CalcNumLipinskiHBA(mol)
        self.prop_dict["Hbd"] = rdMolDescriptors.CalcNumLipinskiHBD(mol)
        self.prop_dict["rb"] = rdMolDescriptors.CalcNumRotatableBonds(mol)
        return self.prop_dict

    def _passes_global_constraints(self, prop_dict):
        return (
            100 < prop_dict["MW"] < 1000 and
            -5 < prop_dict["LogP"] < 10 and
            prop_dict["rb"] < 20 and
            prop_dict["Hba"] < 20 and
            prop_dict["Hbd"] < 20
        )

    def _build_physchem_window(self, smiles):
        self._calc_physchem(smiles)
        mw = self.prop_dict["MW"]
        logp = self.prop_dict["LogP"]
        hba = self.prop_dict["Hba"]
        hbd = self.prop_dict["Hbd"]
        rb  = self.prop_dict["rb"]
        window = {
            "MW_min": mw - 50,
            "MW_max": mw + 50,
            "LogP_min": logp - 0.5,
            "LogP_max": logp + 0.5,
            "Hba_min": max(hba - 1, 0),
            "Hba_max": hba + 1,
            "Hbd_min": max(hbd - 1, 0),
            "Hbd_max": hbd + 1,
            "rb_min": max(rb - 1, 0),
            "rb_max": rb + 1,
        }
        return window
    
    def filter_by_physchem(self, candidate_smiles_list, smiles):
        window = self._build_physchem_window(smiles)
        filtered = []
        for smi in candidate_smiles_list:
            props = self.smiles_to_props.get(smi)
            if props is None:
                continue
            if not self._passes_global_constraints(props):
                continue
            if not (
                window["MW_min"]   <= props["MW"]   <= window["MW_max"]   and
                window["LogP_min"] <= props["LogP"] <= window["LogP_max"] and
                window["Hba_min"]  <= props["Hba"]  <= window["Hba_max"]  and
                window["Hbd_min"]  <= props["Hbd"]  <= window["Hbd_max"]  and
                window["rb_min"]   <= props["rb"]   <= window["rb_max"]
            ):
                continue
            filtered += [smi]
        return filtered
    
class ScaffoldDiversity(object):
    def __init__(self):
        df = pd.read_csv(os.path.join(root, "..", "fit", "data", "chembl_scaffolds.csv"))
        self.smiles_to_scaffold = dict(
            zip(df["smiles"].astype(str), df["scaffold"].astype(str))
        )

    def _get_precalculated_scaffold(self, smiles):
        return self.smiles_to_scaffold.get(smiles)

    def _calculate_scaffold(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        scaff_mol = MurckoScaffold.GetScaffoldForMol(mol)
        if scaff_mol is None:
            return None
        scaff_smi = Chem.MolToSmiles(scaff_mol, isomericSmiles=False)
        return scaff_smi

    def filter_by_scaffold(self, candidate_smiles_list, input_smiles):
        input_scaffold = self._calculate_scaffold(input_smiles)
        filtered = []
        for smi in candidate_smiles_list:
            scaff= self._get_precalculated_scaffold(smi)
            if input_scaffold is None or scaff is None or scaff != input_scaffold:
                filtered += [smi]
        return filtered
