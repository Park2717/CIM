import os
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem
from rdkit.Chem.Fingerprints import FingerprintMols

def Tofingerprint(smiles):
    smi = Chem.CanonSmiles(smiles)
    m = Chem.MolFromSmiles(smi)
    ep1 = AllChem.GetMorganFingerprint(m, 1)
    ep2 = AllChem.GetMorganFingerprint(m, 2)
    ep3 = AllChem.GetMorganFingerprint(m, 3)
    fp1 = AllChem.GetMorganFingerprint(m, 1, useFeatures=True)
    fp2 = AllChem.GetMorganFingerprint(m, 2, useFeatures=True)
    fp3 = AllChem.GetMorganFingerprint(m, 3, useFeatures=True)
    return [ep1, ep2, ep3, fp1, fp2, fp3]

A_Smiles = ['OCCc1cn(-c2ccccc2)c(NCc2cccc(-c3cn[nH]c3)c2)n1']
B_Smiles = \
['O=C(O)CCCCCc1cn(-c2ccccc2)c(NC(=O)c2cccc(-c3cn[nH]c3)c2)n1',\
'O=C(O)CCCc1cn(-c2ccccc2)c(NC(=O)c2cccc(-c3cn[nH]c3)c2)n1',\
'O=C(O)CCc1cn(-c2ccccc2)c(NC(=O)c2cccc(-c3cn[nH]c3)c2)n1']

for A in A_Smiles:
    for B in B_Smiles:
        molA = Tofingerprint(A)
        molB = Tofingerprint(B)
        for mol in range(6):
            sim = DataStructs.DiceSimilarity(molA[mol], molB[mol])
            print(f'{A} : {B} = {sim}')

# from joblib import Parallel
# import multiprocessing
#
# num_cores = multiprocessing.cpu_count()
# # results = Parallel(n_jobs=num_cores)(delayed(file_ref)(i) for i in file_list_ref_py)
#
# for cutoff in [0.5, 0.6, 0.7, 0.8, 0.9]:
#     IRAK4_Pair = pd.DataFrame(columns=['molA', 'molB'])
#
#     for inactive in IRAK4_inactive.itertuples():
#         for active in IRAK4_active.itertuples():
#             molA = Tofingerprint(inactive.Smiles)
#             molB = Tofingerprint(active.Smiles)
#             for mol in range(6):
#                 sim = DataStructs.DiceSimilarity(molA[mol], molB[mol])
#                 if sim >= cutoff:
#                     pair = pd.DataFrame([{'molA' : inactive.Smiles, 'molB' : active.Smiles}])
#                     IRAK4_Pair = pd.concat([IRAK4_Pair, pair])
#                     continue
#
#     print(f'cutoff ({cutoff}) pair data : {len(IRAK4_Pair)}')
#     IRAK4_Pair.to_csv(f'./IRAK4_cutoff_{cutoff}.txt', header=False, index=False)

