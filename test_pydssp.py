import sys
from Bio.PDB import PDBParser
import pydssp
import numpy as np

parser = PDBParser(QUIET=True)
struct = parser.get_structure('1AFW', '1AFW.pdb')
coords = []
res_list = []
for model in struct:
    for chain in model:
        for res in chain:
            if res.id[0].strip(): continue
            try:
                n = res['N'].coord
                ca = res['CA'].coord
                c = res['C'].coord
                o = res['O'].coord
                coords.append([n, ca, c, o])
                res_list.append(res)
            except KeyError:
                pass
                
coord_arr = np.array(coords)
print(f"Shape: {coord_arr.shape}")
ss = pydssp.assign(coord_arr, out_type='c3')
print(f"SS Shape: {ss.shape}")
print(f"First 10 SS: {ss[:10]}")
