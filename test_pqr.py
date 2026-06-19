from Bio.PDB import PDBParser
import subprocess

subprocess.run([".venv/bin/pdb2pqr", "--ff=PARSE", "--keep-chain", "1ake.pdb", "test.pqr"], check=True)
p = PDBParser(QUIET=True)
s = p.get_structure('s', 'test.pqr')
a = list(s.get_atoms())[0]
print("occupancy (should be charge):", a.occupancy)
print("bfactor (should be radius):", a.bfactor)
