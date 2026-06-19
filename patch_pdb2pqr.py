with open('sauron.py', 'r') as f:
    content = f.read()

content = content.replace(
    'cmd = ["pdb2pqr", "--ff=PARSE", "--keep-chain", input_file, output_file]',
    'import os\n            pdb2pqr_path = os.path.join(os.path.dirname(sys.executable), "pdb2pqr")\n            if not os.path.exists(pdb2pqr_path): pdb2pqr_path = "pdb2pqr"\n            cmd = [pdb2pqr_path, "--ff=PARSE", "--keep-chain", input_file, output_file]'
)

with open('sauron.py', 'w') as f:
    f.write(content)

