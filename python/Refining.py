import os
import pandas as pd
import numpy as np
from tqdm import tqdm
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import PandasTools
from rdkit.Chem import AllChem

class convert_format:
    def __init__(self, inputfile):
        self.inputfile = inputfile
        lst = inputfile.split('/')
        path = '/'.join(lst[:-1])
        name = lst[-1].split('.')[0]
        self.path = path
        self.name = name

    def xlsx_to_sdf(self):
        dataframe = pd.read_excel(self.inputfile)
        dataframe = dataframe.loc[:, ~dataframe.columns.str.contains('^Unnamed')]
        ID = input(f'Your columns are {list(dataframe.columns)}. \n'
                   'Which column do you want to select to ID?:')
        PandasTools.AddMoleculeColumnToFrame(dataframe, 'Smiles', 'Molecule')
        auto = input('Do you want to set properties Automatically? [y/n]:')
        if auto == 'y':
            columns = [x for x in dataframe.columns if not (x == ID or x == 'Smiles')]
        else:
            columns = input('Which columns do you want to use as properties?:').split(', ')
        PandasTools.WriteSDF(dataframe, f'{self.path}/{self.name}.sdf', molColName='Molecule',idName=str(ID), properties=columns)
        print(f'Your file "{self.name}" is successfully converted from excel to sdf file.')

    def sdf_to_xlsx(self):
        dataframe = PandasTools.LoadSDF(self.inputfile)
        smi_list = []
        for idx, mol in enumerate(list(dataframe['ROMol'])):
            try:
                smi = Chem.MolToSmiles(mol)
                smi_list.append(smi)
            except:
                print(f'Lost Smiles \n index: {idx}, ID: {dataframe["ID"][idx]}')
                continue
        dataframe['Smiles'] = smi_list
        result = pd.DataFrame(dataframe, columns=['ID'] + ['Smiles'] + list(dataframe.columns[:-3]))
        result.to_excel(f'{self.path}/{self.name}.xlsx')
        print(f'Your file "{self.name}" is successfully converted from sdf to excel file.')

    def sdf_to_dataframe(self):
        dataframe = PandasTools.LoadSDF(self.inputfile)
        smi_list = []
        for idx, mol in enumerate(list(dataframe['ROMol'])):
            try:
                smi = Chem.MolToSmiles(mol)
                smi_list.append(smi)
            except:
                print(f'Lost Smiles \n index: {idx}, ID: {dataframe["ID"][idx]}')
                continue
        dataframe['Smiles'] = smi_list
        result = pd.DataFrame(dataframe, columns=['ID'] + ['Smiles'] + list(dataframe.columns[:-3]))
        return result

if __name__ == "__main__":
    convert_format('D:/New_Target/ULK1/ULK1_purecompounds.xlsx').xlsx_to_sdf()