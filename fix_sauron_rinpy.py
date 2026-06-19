import sys

with open('sauron.py', 'r') as f:
    content = f.read()

# First remove the bad insertion
import re
# We need to find the block that was inserted and move it up.
lines = content.split('\n')
start_idx = -1
end_idx = -1
for i, line in enumerate(lines):
    if "if args.calc_method == 'rinpy':" in line and "from rinpy import RINProcess" in content:
        start_idx = i
        break

if start_idx != -1:
    for i in range(start_idx, len(lines)):
        if "continue" in lines[i]:
            end_idx = i
            break

if start_idx != -1 and end_idx != -1:
    block = lines[start_idx:end_idx+1]
    # Remove the block
    lines = lines[:start_idx] + lines[end_idx+1:]
    
    # Unindent the block by 4 spaces
    unindented_block = []
    for line in block:
        if line.startswith('    '):
            unindented_block.append(line[4:])
        else:
            unindented_block.append(line)
            
    # Insert it before `if args.add_h:`
    insert_idx = -1
    for i, line in enumerate(lines):
        if "if args.add_h:" in line:
            insert_idx = i
            break
            
    if insert_idx != -1:
        # Also change pqr_file to temp_pdb in the block
        final_block = []
        for line in unindented_block:
            line = line.replace('pqr_file', 'temp_pdb')
            final_block.append(line)
            
        lines = lines[:insert_idx] + final_block + lines[insert_idx:]
        
        with open('sauron.py', 'w') as f:
            f.write('\n'.join(lines))
        print("Fixed!")
    else:
        print("Could not find insert index.")
else:
    print("Could not find block.")
