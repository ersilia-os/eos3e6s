# imports
import os
import csv
import sys
import random
import numpy as np

# current file directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(root)
from similarity import ChemblLowSimilarity, ChemblProps, ScaffoldDiversity

# parse arguments
input_file = sys.argv[1]
output_file = sys.argv[2]

def filter_pipeline(smiles_list):
    R = []
    len_sim = []
    len_phys = []
    len_scaff = []
    sim_filter = ChemblLowSimilarity()
    physchem_filter = ChemblProps()
    scaffold_filter = ScaffoldDiversity()
    active_set = set(smiles_list)
    used_decoys = set()
    for smi in smiles_list:
        low_smiles = sim_filter.get_low_similarity(smi)
        print("Similarity filter: ", len(low_smiles))
        len_sim += [len(low_smiles)]
        low_smiles = [s for s in low_smiles if s not in used_decoys and s not in active_set]
        chembl_physchem_filtered = physchem_filter.filter_by_physchem(low_smiles,smi)
        print("Physchem filter: ", len(chembl_physchem_filtered))
        len_phys += [len(chembl_physchem_filtered)]
        chembl_scaffold_filtered = scaffold_filter.filter_by_scaffold(chembl_physchem_filtered, smi)
        print("Scaffold filter: ", len(chembl_scaffold_filtered))
        len_scaff += [len(chembl_scaffold_filtered)]
        chembl_active_filtered = [s for s in chembl_scaffold_filtered if s not in active_set]
        if len(chembl_active_filtered) < 100:
            n = 100-len(chembl_active_filtered)
            fillin_subset = [s for s in low_smiles if s not in chembl_active_filtered] #select from low similarity without other constraints
            if len(fillin_subset) >= n:
                extra = random.sample(fillin_subset, n)
            else: #make extra sure that even if we dont have enough low sim molecules code does not crash
                extra = fillin_subset + ([np.nan] * (n - len(fillin_subset)))
            r = chembl_active_filtered+extra
        if len(chembl_active_filtered)>=100:
            r = random.sample(chembl_active_filtered,100)
        used_decoys.update([s for s in r if isinstance(s, str)])
        R += [r]
    return R, len_sim, len_phys, len_scaff

# read SMILES from .csv file, assuming one column with header
with open(input_file, "r") as f:
    reader = csv.reader(f)
    next(reader)  # skip header
    smiles_list = [r[0] for r in reader]

# run model
outputs, len_sim, len_phys, len_scaff = filter_pipeline(smiles_list)

#check input and output have the same lenght
input_len = len(smiles_list)
output_len = len(outputs)
assert input_len == output_len

# write output in a .csv file
with open(output_file, "w") as f:
    writer = csv.writer(f)
    writer.writerow([f"smi_{i:02d}" for i in range(100)])
    for r in outputs:
        writer.writerow(r)