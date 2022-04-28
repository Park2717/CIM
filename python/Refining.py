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
        if type(inputfile) == type(pd.DataFrame()): # input path and data when inputformat is dataframe.
            self.fm = 'dataframe'
            path = input('PATH:')
            name = input('NAME:')
            self.path = path
            self.name = name
        else:
            lst = inputfile.split('/') # else list the path, name and format.
            path = '/'.join(lst[:-1])
            name = lst[-1].split('.')[0]
            fm = lst[-1].split('.')[1]
            self.path = path
            self.name = name
            self.fm = fm

    def to_dataframe(self, sep='\t'): # convert your file to dataframe regardless of the file format.
        if type(self.inputfile) == type(pd.DataFrame()):
            dataframe = self.inputfile
            dataframe = dataframe.loc[:, ~dataframe.columns.str.contains('^Unnamed')]
            return dataframe
        elif self.fm == 'txt':
            dataframe = pd.read_csv(self.inputfile, sep=sep)
            return dataframe
        elif self.fm == 'xlsx':
            dataframe = pd.read_excel(self.inputfile)
            dataframe = dataframe.loc[:, ~dataframe.columns.str.contains('^Unnamed')]
            return dataframe
        elif self.fm == 'sdf':
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
            result = pd.DataFrame(dataframe, columns=['ID'] + ['Smiles'] + list(dataframe.columns[:-3]))
            return result
        else:
            print('This format is not available yet.')

    def to_txt(self, sep='\t'):
        dataframe = self.to_dataframe()
        dataframe.to_csv(f'{self.path}/{self.name}.txt', sep=sep)
        print(f'Your file "{self.name}" is successfully converted from {self.fm} to text file.')

    def to_xlsx(self):
        dataframe = self.to_dataframe()
        dataframe.to_excel(f'{self.path}/{self.name}.xlsx')
        print(f'Your file "{self.name}" is successfully converted from {self.fm} to excel file.')

    def to_sdf(self):
        dataframe = self.to_dataframe()
        ID = input(f'Your columns are {list(dataframe.columns)}. \n' # show your columns.
                   'Which column do you want to select to ID?:')     # select ID column.
        PandasTools.AddMoleculeColumnToFrame(dataframe, 'Smiles', 'Molecule')
        auto = input('Do you want to set properties Automatically? [Y/N]:')
        # if yes, input all of your propertiees except ID and Smiles.
        if auto.upper() == 'Y':
            columns = [x for x in dataframe.columns if not (x == ID or x == 'Smiles')]
        # Or you can choose your properties.
        else:
            columns = input('Which columns do you want to use as properties?:').split(', ')
        PandasTools.WriteSDF(dataframe, f'{self.path}/{self.name}.sdf', molColName='Molecule',idName=str(ID), properties=columns)
        print(f'Your file "{self.name}" is successfully converted from {self.fm} to sdf file.')

    def to_png(self):
        dataframe = self.to_dataframe()
        ID = input(f'Your columns are {list(dataframe.columns)}. \n'
                   'Which column do you want to select as "ID"? [import column name]:\n'
                   'Or do you want to save only structure img? [import "only structure"]:')
        # make new folder for images.
        try:
            os.mkdir(f"{self.path}/{self.name}_draw")
        except:
            print('Already exist folder')
        # indicate structure name.
        if ID.upper() == 'ONLY STRUCTURE':
            structure_name = None
        else:
            structure_name = list(dataframe[ID])
        # make splits if you want.
        if len(list(dataframe['Smiles'])) >= 5:
            split = input(f'The file have {len(dataframe)} structures.\n'
                          'How many structures do you want to include in one png file?:')
            lst = np.array_split(dataframe, math.ceil(len(dataframe) / int(split)))
            for idx, i in enumerate(lst):
                mols = [MolFromSmiles(mol) for mol in list(i['Smiles'])]
                img = Draw.MolsToGridImage(mols, molsPerRow=4,
                                           subImgSize=(250, 250), legends=structure_name)
                img.save(f"{self.path}/{self.name}_draw/{self.name}_{idx}.png")
        else:
            mols = [MolFromSmiles(mol) for mol in list(dataframe['Smiles'])]
            img = Draw.MolsToGridImage(mols, molsPerRow=len(list(dataframe['Smiles'])),
                                       subImgSize=(250, 250), legends=structure_name)
            img.save(f"{self.path}/{self.name}_draw/{self.name}.png")

    def abstract(self, start=0, finish=None):
        dataframe = self.to_dataframe()
        # select rows which you indicated from input.
        dataframe = dataframe[start:finish]
        # select columns you want.
        columns = input(f'Your columns are {list(dataframe.columns)}. \n'
                        'Which column do you want to abstract?:').split(', ')
        if len(columns) == 1:
            columns = columns[0]
            print(f'{list(dataframe[columns])}\n[type: list, length: {len(list(dataframe[columns]))}]')
            return list(dataframe[columns])
        else:
            result = pd.DataFrame(dataframe, columns=columns)
            print(f'{result}\n[type: dataframe, length: {len(result)}]')
            return result

if __name__ == "__main__":
    #dataframe = convert_format('D:/New_Target/ULK1/ULK1_purecompounds.xlsx').abstract()
    #dataframe = pd.read_excel('D:/New_Target/ULK1/ULK1_purecompounds.xlsx')
    #dataframe_to_sdf(dataframe)
    #dataframe = convert_format('D:/python_coding_files/OTAVA_90_IRAK4_ACD.sdf').sdf_to_dataframe()
    #print(dataframe)
    #df = pd.DataFrame()
    #type(dataframe) == type(pd.DataFrame())
    #type(1) == int
    #convert_format(dataframe).to_png()
    convert_format(dataframe).to_png()