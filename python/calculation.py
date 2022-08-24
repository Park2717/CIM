import numpy as np
import math
import pandas as pd
from tqdm import tqdm
from rdkit import Chem

def find_fragment(smiles, fragment):
    s = Chem.MolFromSmiles(Chem.CanonSmiles(smiles))
    fragment = Chem.CanonSmiles(fragment)
    f = Chem.MolFromSmiles(fragment)
    match = s.HasSubstructMatch(f)
    return match

def analyze_pdb(file):
    with open(file, 'r') as f:
        lines = f.readlines()
        n=0
        m=0
        k=0
        result = pd.DataFrame()
        for line in tqdm(lines, desc=f'Analyzing {file.split("/")[-1]}'):
            if line.startswith('TER'):
                n+=1
            if line.startswith('ATOM'):
                if k==0:
                    m = n
                    atm_type = line[11:17].strip()
                    rsd_type = line[17:20]
                    rsd_num = int(line[23:26])
                    Coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    tmp = {'ATOM Type': atm_type, 'RESIDUE Type': rsd_type, 'RESIDUE NUM': rsd_num, 'Coordinate': Coord}
                    if atm_type in ['N', 'C', 'O']:
                        result = pd.concat([result, pd.DataFrame([tmp])], ignore_index=True)
                    k += 1
                elif k!=0 and m!=n:
                    break
                else:
                    atm_type = line[11:17].strip()
                    rsd_type = line[17:20]
                    rsd_num = int(line[23:26])
                    Coord = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
                    tmp = {'ATOM Type': atm_type, 'RESIDUE Type': rsd_type, 'RESIDUE NUM': rsd_num, 'Coordinate': Coord}
                    if atm_type in ['N', 'C', 'O']:
                        result = pd.concat([result, pd.DataFrame([tmp])], ignore_index=True)
        Coordinates = list(result['Coordinate'])
        focus = tuple(np.average(Coordinates, axis=0))
        return focus, result

def rotation(focus, Coord, angle):
    (x, y, z) = (Coord[0] - focus[0], Coord[1] - focus[1], Coord[2] - focus[2])
    (a1, a2, b1, b2, r1, r2) = angle
    # alpha inverse rotation
    n = x * math.cos(-a2) - y * math.sin(-a2)
    m = x * math.sin(-a2) + y * math.cos(-a2)
    (x, y, z) = (n, m ,z)
    # gamma inverse rotation
    n = z * math.cos(-r2) - x * math.sin(-r2)
    m = z * math.sin(-r2) + x * math.cos(-r2)
    (x, y, z) = (m, y, n)
    # beta inverse rotation
    n = x * math.cos(-b2) - y * math.sin(-b2)
    m = x * math.sin(-b2) + y * math.cos(-b2)
    (x, y, z) = (n, m, z)
    # beta rotation
    n = x * math.cos(b1) - y * math.sin(b1)
    m = x * math.sin(b1) + y * math.cos(b1)
    (x, y, z) = (n, m, z)
    # gamma rotation
    n = z * math.cos(r1) - x * math.sin(r1)
    m = z * math.sin(r1) + x * math.cos(r1)
    (x, y, z) = (m, y, n)
    # alpha rotation
    n = x * math.cos(a1) - y * math.sin(a1)
    m = x * math.sin(a1) + y * math.cos(a1)
    (x, y, z) = (n + focus[0], m + focus[1], z + focus[2])
    return x, y, z

def cal_rot_angle(focus, Coord1, Coord2):
    (x1, y1, z1) = (Coord1[0] - focus[0], Coord1[1] - focus[1], Coord1[2] - focus[2])
    (x2, y2, z2) = (Coord2[0] - focus[0], Coord2[1] - focus[1], Coord2[2] - focus[2])
    alpha = math.atan2(y1, x1)
    n1 = x1 * math.cos(-alpha) - y1 * math.sin(-alpha)
    m1 = x1 * math.sin(-alpha) + y1 * math.cos(-alpha)
    (x1, y1, z1) = (n1, m1, z1)
    n2 = x2 * math.cos(-alpha) - y2 * math.sin(-alpha)
    m2 = x2 * math.sin(-alpha) + y2 * math.cos(-alpha)
    (x2, y2, z2) = (n2, m2, z2)
    gamma = math.atan2(x1, z1)
    n2 = z2 * math.cos(-gamma) - x2 * math.sin(-gamma)
    m2 = z2 * math.sin(-gamma) + x2 * math.cos(-gamma)
    (x2, y2, z2) = (m2, y2, n2)
    beta = math.atan2(y2, x2)
    return alpha, beta, gamma

