import re

with open('sauron.py', 'r') as f:
    lines = f.readlines()

def replace_energy_in_line(i, v1, v2):
    lines[i] = re.sub(r"'Energy': 'nan'", f"'Energy': f\"{{calc_energy_terms({v1}, {v2}):.4f}}\" if calc_energy else 'nan'", lines[i])

# SSBOND
for i, line in enumerate(lines):
    if "'Interaction': 'SSBOND:SC_SC'" in line:
        replace_energy_in_line(i+1, "r1['SG']", "r2['SG']")
    elif "'Interaction': 'IONIC:SC_SC'" in line:
        replace_energy_in_line(i+1, "a1", "a2")
    elif "'Interaction': f'METAL:" in line:
        replace_energy_in_line(i+1, "a1", "a2")
    elif "'Interaction': f'HALOGEN:" in line:
        replace_energy_in_line(i+1, "a1", "a2")
    elif "'Interaction': f'HBOND:" in line:
        replace_energy_in_line(i+1, "a1", "a2")
    elif "'Interaction': f'VDW:" in line:
        replace_energy_in_line(i+1, "a1", "a2")
    elif "'Interaction': f'WATERBRIDGE:" in line:
        replace_energy_in_line(i+1, "pa1", "pa2")

with open('sauron.py', 'w') as f:
    f.writelines(lines)

