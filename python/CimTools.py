import pandas as pd
import numpy as np
from tqdm import tqdm
from CIM_Module.Refining import convert_format as cf
from rdkit import Chem
from rdkit.Chem import PandasTools, AllChem
from CIM_Module.calculation import *

def MoleculeMinimizer(inputfile, forcefield='MMFF'):
    df = cf(inputfile)()
    df['ROMol'] = df['Smiles'].apply(lambda x: Chem.MolFromSmiles(x))
    df['ROMol'] = df['ROMol'].apply(lambda x: AllChem.AddHs(x))
    df['ROMol'].apply(lambda x: AllChem.EmbedMolecule(x))
    if forcefield.upper() == 'MMFF':
        df['ROMol'].apply(lambda x: AllChem.MMFFOptimizeMolecule(x))
    elif forcefield.upper() == 'UFF':
        df['ROMol'].apply(lambda x: AllChem.UFFOptimizeMolecule(x))
    else:
        print(f'Forcefield "{forcefield}" is not supported.\nAvailable forcefields: "MMFF", "UFF".')
    PandasTools.WriteSDF(df, inputfile.split('.')[0] + '_minimized.sdf', molColName='ROMol', idName='ID')
    print('successfully minimized')

def fragment_search(DB, fragment):
    DF = cf(DB)()
    dataframe = pd.DataFrame()
    match = DF['Smiles'].apply(lambda x: find_fragment(x, fragment))
    result = DF[match]
    if len(result) != 0:
        print(f'There are {len(result)} compounds matching {fragment}.')
    dataframe = pd.concat([dataframe, result], ignore_index=True)
    return dataframe

def superimpose_proteins(Reference, inputfile, All=True): # only pdb format available.
    # calculate focus of Reference protein
    f1, d1 = analyze_pdb(Reference)
    # calculate focus of target protein to superimpose
    f2, d2 = analyze_pdb(inputfile)
    # parallel movement
    pm = (f2[0] - f1[0], f2[1] - f1[1], f2[2] - f1[2])
    # cal abr dif
    a1, b1, r1 = cal_rot_angle(f1, list(d1['Coordinate'])[0], list(d1['Coordinate'])[3])
    a2, b2, r2 = cal_rot_angle(f2, list(d2['Coordinate'])[0], list(d2['Coordinate'])[3])
    angle = (a1, a2, b1, b2, r1, r2)
    # superimpose
    with open(inputfile, 'r') as f:
        lines = f.readlines()
        new_lines = []
        for line in tqdm(lines, desc='superimposing'):
            if line.startswith('ATOM'):
                Coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                pm_Coord = (Coord[0] - pm[0], Coord[1] - pm[1], Coord[2] - pm[2])
                m_Coord = rotation(f1, pm_Coord, angle)
                final_Coord = (round(m_Coord[0], 3), round(m_Coord[1], 3), round(m_Coord[2], 3))
                new_line = line[:30] + "{:>8}".format(final_Coord[0]) + "{:>8}".format(final_Coord[1]) + "{:>8}".format(final_Coord[2]) + line[54:]
                new_lines += [new_line]
            elif All and line.startswith('HETATM'):
                Coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                pm_Coord = (Coord[0] - pm[0], Coord[1] - pm[1], Coord[2] - pm[2])
                m_Coord = rotation(f1, pm_Coord, angle)
                final_Coord = (round(m_Coord[0], 3), round(m_Coord[1], 3), round(m_Coord[2], 3))
                new_line = line[:30] + "{:>8}".format(final_Coord[0]) + "{:>8}".format(final_Coord[1]) + "{:>8}".format(final_Coord[2]) + line[54:]
                new_lines += [new_line]
            else:
                new_lines += [line]
    fr = open(f'{inputfile.split(".")[0]}_superimposed.pdb', 'w')
    for new_line in tqdm(new_lines, desc='save'):
        fr.write(new_line)
    fr.close()
    print('successfully superimposed')

def remove_proteins(inputfile): # only pdb format available.
    file = inputfile
    with open(file, 'r') as f:
        lines = f.readlines()
        new_lines = []
        for line in lines:
            if line.startswith('ATOM'):
                if line[17:20] != 'SOL':
                    continue
                else:
                    new_lines += [line]
            else:
                new_lines += [line]
    fr = open(f'{file.split(".")[0]}_removed.pdb', 'w')
    for new_line in new_lines:
        fr.write(new_line)
    fr.close()
    print('successfully removed')

def compounds_search(DB, compounds): # only sdf format available.
    # input_file
    data = PandasTools.LoadSDF(compounds, smilesName='Smiles')
    # DB to search
    DB_data = PandasTools.LoadSDF(DB, smilesName='Smiles')
    # smi filter
    smis = data['Smiles']
    canon_smis = smis.apply(lambda x: Chem.CanonSmiles(x))
    DB_smis = DB_data['Smiles']
    DB_canon_smis = DB_smis.apply(lambda x: Chem.CanonSmiles(x))
    n = 0
    for smi in canon_smis:
        if smi in DB_canon_smis:
            n += 1
    print(f"Already exist smis : {n}")

def devide_sdf_mols(inputfile, cutoff):
    path = '/'.join(inputfile.split('/')[:-1])
    file = inputfile.split('/')[-1]
    with open(f'{path}/{file}', 'r') as f:
        lines = f.readlines()
        n = 1
        m = 1
        new_lines = []
        for idx, line in enumerate(lines):
            new_lines += [line]
            if line.startswith('$$$$'):
                n += 1
            if n == cutoff:
                fr = open(f'{path}/{file.split(".")[0]}_sep{m}.sdf', 'w')
                for new_line in new_lines:
                    fr.write(new_line)
                fr.close()
                n = 1
                new_lines = []
                m += 1
            if idx == len(lines) - 1:
                fr = open(f'{path}/{file.split(".")[0]}_sep{m}.sdf', 'w')
                for new_line in new_lines:
                    fr.write(new_line)
                fr.close()
    print('complete')