import re

with open('sauron.py', 'r') as f:
    content = f.read()

# 1. Update def calculate_rin
content = content.replace(
    'def calculate_rin(structure, strict_angle=False, remove_multiples=False, model_num=1, calc_method="standard", pqr_file=None):',
    'def calculate_rin(structure, strict_angle=False, remove_multiples=False, model_num=1, calc_method="standard", pqr_file=None, calc_energy=False, include_vdw=False):'
)

# 2. Add Energy helper inside calculate_rin
helper = """    def calc_energy_terms(a1, a2):
        if not calc_energy: return 'nan'
        try:
            q1 = getattr(a1, 'occupancy', 0.0) or 0.0
            q2 = getattr(a2, 'occupancy', 0.0) or 0.0
            r1 = getattr(a1, 'bfactor', 0.0) or 0.0
            r2 = getattr(a2, 'bfactor', 0.0) or 0.0
            dist = np.linalg.norm(a1.getCoords() if hasattr(a1, 'getCoords') else a1.coord - (a2.getCoords() if hasattr(a2, 'getCoords') else a2.coord))
            if dist == 0: return 0.0
            e_elec = 332.0637 * q1 * q2 / (4.0 * dist * dist)
            e_vdw = 0.0
            if include_vdw and r1 > 0 and r2 > 0:
                ratio = (r1 + r2) / dist
                ratio6 = ratio**6
                e_vdw = 0.1 * (ratio6**2 - 2 * ratio6)
            return e_elec + e_vdw
        except:
            return 'nan'
"""
content = content.replace('    edges = []\n    nodes = set()', helper + '\n    edges = []\n    nodes = set()')

# 3. Add Energy to add_edge in insty
content = content.replace(
    "def add_edge(n1, n2, interaction, dist, angle, a1='', a2=''):",
    "def add_edge(n1, n2, interaction, dist, angle, a1='', a2='', energy='nan'):"
)
content = content.replace(
    "'Angle': f\"{angle:.3f}\" if angle != 'nan' else 'nan', \n                    'Atom1': a1",
    "'Angle': f\"{angle:.3f}\" if angle != 'nan' else 'nan', \n                    'Energy': f\"{energy:.4f}\" if energy != 'nan' else 'nan', 'Atom1': a1"
)

# 4. We need to add 'Energy': 'nan' to all other edges.append
# We can use regex to inject 'Energy': 'nan', inside the dictionary
import re
content = re.sub(
    r"'Angle':([^,]+), 'Atom1'",
    r"'Angle':\1, 'Energy': 'nan', 'Atom1'",
    content
)

# 5. We need to update the cols
content = content.replace(
    "cols = ['NodeId1', 'Interaction', 'NodeId2', 'Distance', 'Angle', 'Atom1', 'Atom2', 'Donor', 'Positive', 'Cation', 'Orientation', 'Model']",
    "cols = ['NodeId1', 'Interaction', 'NodeId2', 'Distance', 'Angle', 'Energy', 'Atom1', 'Atom2', 'Donor', 'Positive', 'Cation', 'Orientation', 'Model']"
)

with open('sauron.py', 'w') as f:
    f.write(content)

