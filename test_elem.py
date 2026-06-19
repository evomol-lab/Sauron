from Bio.PDB import PDBParser
p = PDBParser(QUIET=True)
s = p.get_structure('s', '4ake.pdb')
atoms = list(s.get_atoms())
print("Element of CA atom:", atoms[1].element)
print("Name of CA atom:", atoms[1].name)
