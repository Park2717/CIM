import os
import pandas as pd
import numpy as np
import math
from tqdm import tqdm
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Draw, PandasTools, MolFromSmiles

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

    def xlsx_to_png(self):
        dataframe = pd.read_excel(self.inputfile)
        dataframe = dataframe.loc[:, ~dataframe.columns.str.contains('^Unnamed')]
        ID = input(f'Your columns are {list(dataframe.columns)}. \n'
                   'Which column do you want to select as "ID"? [import column name]:\n'
                   'Or do you want to save only structure img? [import "only structure"]:')
        try:
            os.mkdir(f"{self.path}/{self.name}_draw")
        except:
            print('Already exist folder')
        if ID.upper() == 'ONLY STRUCTURE':
            if len(list(dataframe['Smiles'])) <= 8:
                if len(list(dataframe['Smiles'])) <= 4:
                    mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=len(list(dataframe['Smiles'])),
                                               subImgSize=(250, 250), legends=None)
                    img.save(f"{self.path}/{self.name}_draw/{self.name}.png")
                else:
                    mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                               subImgSize=(250, 250), legends=None)
                    img.save(f"{self.path}/{self.name}_draw/{self.name}.png")
            else:
                split = input('The file have too many structures.\n'
                              'How many structures do you want to include in one png file?:')
                lst = np.array_split(dataframe, math.ceil(len(dataframe)/int(split)))
                for idx, i in enumerate(lst):
                    mols = [MolFromSmiles(mol) for mol in list(i['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                               subImgSize=(250, 250), legends=None)
                    img.save(f"{self.path}/{self.name}_draw/{self.name}_{idx}.png")
        else:
            if len(list(dataframe['Smiles'])) <= 8:
                if len(list(dataframe['Smiles'])) <= 4:
                    mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=len(list(dataframe['Smiles'])),
                                               subImgSize=(250, 250), legends=list(dataframe[ID]))
                    img.save(f"{self.path}/{self.name}_draw/{self.name}.png")
                else:
                    mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                               subImgSize=(250, 250), legends=list(dataframe[ID]))
                    img.save(f"{self.path}/{self.name}_draw/{self.name}.png")
            else:
                split = input('The file have too many structures.\n'
                              'How many structures do you want to include in one png file?:')
                lst = np.array_split(dataframe, math.ceil(len(dataframe)/int(split)))
                for idx, i in enumerate(lst):
                    mols = [MolFromSmiles(mol) for mol in list(i['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                               subImgSize=(250, 250), legends=list(i[ID]))
                    img.save(f"{self.path}/{self.name}_draw/{self.name}_{idx+1}.png")
        print(f'Your file "{self.name}" is successfully converted from excel to png file.')

    def sdf_to_png(self):
        dataframe = PandasTools.LoadSDF(self.inputfile)
        smi_list = []
        for idx, mol in enumerate(list(dataframe['ROMol'])):
            try:
                smi = Chem.MolToSmiles(mol)
                smi_list.append(smi)
            except:
                print(f'Lost Smiles \n index: {idx}, ID: {dataframe["ID"][idx]}')
                dataframe.drop([idx])
                continue
        dataframe['Smiles'] = smi_list
        ID = input(f'Do you want to naming ID? [y/n]:')
        try:
            os.mkdir(f"{self.path}/{self.name}_draw")
        except:
            print('Already exist folder')
        if ID.upper() == 'N':
            if len(list(dataframe['Smiles'])) <= 8:
                if len(list(dataframe['Smiles'])) <= 4:
                    mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=len(list(dataframe['Smiles'])),
                                               subImgSize=(250, 250), legends=None)
                    img.save(f"{self.path}/{self.name}_draw/{self.name}.png")
                else:
                    mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                               subImgSize=(250, 250), legends=None)
                    img.save(f"{self.path}/{self.name}_draw/{self.name}.png")
            else:
                split = input('The file have too many structures.\n'
                              'How many structures do you want to include in one png file?:')
                lst = np.array_split(dataframe, math.ceil(len(dataframe) / int(split)))
                for idx, i in enumerate(lst):
                    mols = [MolFromSmiles(mol) for mol in list(i['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                               subImgSize=(250, 250), legends=None)
                    img.save(f"{self.path}/{self.name}_draw/{self.name}_{idx}.png")
        else:
            if len(list(dataframe['Smiles'])) <= 8:
                if len(list(dataframe['Smiles'])) <= 4:
                    mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=len(list(dataframe['Smiles'])),
                                               subImgSize=(250, 250), legends=list(dataframe['ID']))
                    img.save(f"{self.path}/{self.name}_draw/{self.name}.png")
                else:
                    mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                               subImgSize=(250, 250), legends=list(dataframe['ID']))
                    img.save(f"{self.path}/{self.name}_draw/{self.name}.png")
            else:
                split = input('The file have too many structures.\n'
                              'How many structures do you want to include in one png file?:')
                lst = np.array_split(dataframe, math.ceil(len(dataframe) / int(split)))
                for idx, i in enumerate(lst):
                    mols = [MolFromSmiles(mol) for mol in list(i['Smiles'])]
                    img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                               subImgSize=(250, 250), legends=list(i['ID']))
                    img.save(f"{self.path}/{self.name}_draw/{self.name}_{idx + 1}.png")
        print(f'Your file "{self.name}" is successfully converted from sdf to png file.')

def dataframe_to_sdf(dataframe):
    ID = input(f'Your columns are {list(dataframe.columns)}. \n'
               'Which column do you want to select to ID?:')
    PandasTools.AddMoleculeColumnToFrame(dataframe, 'Smiles', 'Molecule')
    columns = [x for x in dataframe.columns if not (x == ID or x == 'Smiles')]
    path = input('Indicate your path to save:')
    name = input('Indicate the file name:')
    PandasTools.WriteSDF(dataframe, f'{path}/{name}.sdf', molColName='Molecule', idName=str(ID),
                         properties=columns)
    print(f'Your file "{name}" is successfully converted from dataframe to sdf file.')

def dataframe_to_png(dataframe):
    ID = input(f'Your columns are {list(dataframe.columns)}. \n'
               'Which column do you want to select as "ID"? [import column name]:\n'
               'Or do you want to save only structure img? [import "only structure"]:')
    path = input('Indicate your path to save:')
    name = input('Indicate the file name:')
    try:
        os.mkdir(f"{path}/{name}_draw")
    except:
        print('Already exist folder')
    if ID.upper() == 'ONLY STRUCTURE':
        if len(list(dataframe['Smiles'])) <= 8:
            if len(list(dataframe['Smiles'])) <= 4:
                mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                img = Draw.MolsToGridImage(mols, molsPerRow=len(list(dataframe['Smiles'])),
                                           subImgSize=(250, 250), legends=None)
                img.save(f"{path}/{name}_draw/{name}.png")
            else:
                mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                           subImgSize=(250, 250), legends=None)
                img.save(f"{path}/{name}_draw/{name}.png")
        else:
            split = input('The file have too many structures.\n'
                          'How many structures do you want to include in one png file?:')
            lst = np.array_split(dataframe, math.ceil(len(dataframe) / int(split)))
            for idx, i in enumerate(lst):
                mols = [MolFromSmiles(mol) for mol in list(i['Smiles'])]
                img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                           subImgSize=(250, 250), legends=None)
                img.save(f"{path}/{name}_draw/{name}_{idx}.png")
    else:
        if len(list(dataframe['Smiles'])) <= 8:
            if len(list(dataframe['Smiles'])) <= 4:
                mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                img = Draw.MolsToGridImage(mols, molsPerRow=len(list(dataframe['Smiles'])),
                                           subImgSize=(250, 250), legends=list(dataframe[ID]))
                img.save(f"{path}/{name}_draw/{name}.png")
            else:
                mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
                img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                           subImgSize=(250, 250), legends=list(dataframe[ID]))
                img.save(f"{path}/{name}_draw/{name}.png")
        else:
            split = input('The file have too many structures.\n'
                          'How many structures do you want to include in one png file?:')
            lst = np.array_split(dataframe, math.ceil(len(dataframe) / int(split)))
            for idx, i in enumerate(lst):
                mols = [MolFromSmiles(mol) for mol in list(i['Smiles'])]
                img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                           subImgSize=(250, 250), legends=list(i[ID]))
                img.save(f"{path}/{name}_draw/{name}_{idx + 1}.png")
    print(f'Your file "{name}" is successfully converted from dataframe to png file.')

if __name__ == "__main__":
    convert_format('D:/New_Target/ULK1/ULK1_purecompounds.xlsx').xlsx_to_png()
    #dataframe = pd.read_excel('D:/New_Target/ULK1/ULK1_purecompounds.xlsx')
    #dataframe_to_sdf(dataframe)
    #dataframe = convert_format('D:/python_coding_files/OTAVA_90_IRAK4_ACD.sdf').sdf_to_dataframe()
    #print(dataframe)