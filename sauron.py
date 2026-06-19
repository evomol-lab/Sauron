import sys
import os
import argparse
try:
    import rinpy
except ImportError:
    rinpy = None

import subprocess
import numpy as np
import pandas as pd
import networkx as nx
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, Selection

# VDW Radii based on generic values
VDW_RADII = {
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'H': 1.20, 'P': 1.80
}

MC_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}

def get_vdw_radius(atom_name):
    element = atom_name[0]
    return VDW_RADII.get(element, 1.70)

def format_node_id(res):
    ins = res.id[2]
    if ins == ' ': ins = '_'
    return f"{res.parent.id}:{res.id[1]}:{ins}:{res.resname}"

def get_contact_type(a1, a2):
    t1 = 'LIG' if a1.parent.id[0].startswith('H') else ('MC' if a1.name in MC_ATOMS else 'SC')
    t2 = 'LIG' if a2.parent.id[0].startswith('H') else ('MC' if a2.name in MC_ATOMS else 'SC')
    return f"{t1}_{t2}"

def get_centroid(res):
    if res.resname in ['PHE', 'TYR']:
        atoms = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
    elif res.resname == 'TRP':
        atoms = ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
    elif res.resname == 'HIS':
        atoms = ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
    else:
        return None
    coords = [res[a].coord for a in atoms if a in res]
    if not coords: return None
    return np.mean(coords, axis=0)

def calc_angle(p1, p2, p3):
    v1 = p1 - p2
    v2 = p3 - p2
    cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
    angle = np.arccos(np.clip(cosine_angle, -1.0, 1.0))
    return np.degrees(angle)

def get_attached_hydrogens(atom):
    hydrogens = []
    for a in atom.parent:
        if a.name.startswith('H') or (len(a.name) > 1 and a.name[0].isdigit() and a.name[1] == 'H'):
            dist = np.linalg.norm(a.coord - atom.coord)
            if dist < 1.3:
                hydrogens.append(a)
    return hydrogens

def check_hbond_angle(a1, a2):
    # Check if a1 is donor
    h1 = get_attached_hydrogens(a1)
    for h in h1:
        angle = calc_angle(a1.coord, h.coord, a2.coord)
        if angle > 120: return True, angle
    
    # Check if a2 is donor
    h2 = get_attached_hydrogens(a2)
    for h in h2:
        angle = calc_angle(a2.coord, h.coord, a1.coord)
        if angle > 120: return True, angle
        
    return False, np.nan

def add_hydrogens(input_file):
    output_file = input_file + ".pqr"
    
    print(f"Adding hydrogens using pdb2pqr to {input_file}...")
    try:
        if getattr(sys, 'frozen', False):
            cmd = [sys.executable, "--run-pdb2pqr", "--ff=PARSE", "--keep-chain", input_file, output_file]
        else:
            import os
            pdb2pqr_path = os.path.join(os.path.dirname(sys.executable), "pdb2pqr")
            if not os.path.exists(pdb2pqr_path): pdb2pqr_path = "pdb2pqr"
            cmd = [pdb2pqr_path, "--ff=PARSE", "--keep-chain", input_file, output_file]
            
        subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        return output_file
    except Exception as e:
        print(f"pdb2pqr failed or not available: {e}. Falling back to original file.")
        return input_file

def parse_structure(filepath):
    if filepath.endswith('.cif'):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    return parser.get_structure('struct', filepath)

def calculate_rin(structure, strict_angle=False, remove_multiples=False, model_num=1, calc_method="standard", pqr_file=None, calc_energy=False, include_vdw=False):
    def calc_energy_terms(a1, a2):
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

    edges = []
    nodes = set()
    
    atoms = list(structure.get_atoms())
    
    if calc_method == "insty":
        try:
            import prody
        except ImportError:
            prody = None

        try:
            import rinpy
        except ImportError:
            rinpy = None

            import warnings
            import re
            warnings.filterwarnings('ignore')
            prody_atoms = prody.parsePQR(pqr_file) if pqr_file and pqr_file.endswith('.pqr') else prody.parsePDB(pqr_file)
            if not prody_atoms:
                raise ValueError("ProDy failed to parse structure.")
            res = prody.calcProteinInteractions(prody_atoms)
            
            def parse_prody_node(res_str, chid):
                m = re.match(r'^(.*?)(\-?\d+)$', res_str)
                if not m: return None, None, None
                resname = m.group(1)
                pos = m.group(2)
                return f"{chid}:{pos}:_:{resname}", resname, pos

            def add_edge(n1, n2, interaction, dist, angle, a1='', a2='', energy='nan'):
                if n1 > n2: n1, n2 = n2, n1
                edges.append({
                    'NodeId1': n1, 'Interaction': interaction, 'NodeId2': n2,
                    'Distance': f"{dist:.3f}", 'Angle': f"{angle:.3f}" if angle != 'nan' else 'nan', 
                    'Energy': f"{energy:.4f}" if energy != 'nan' else 'nan', 'Atom1': a1, 'Atom2': a2, 'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                })
            
            mc_atoms = ['N', 'C', 'CA', 'O']

            def get_prody_contact_type(a1_name, a2_name):
                if a1_name in mc_atoms and a2_name in mc_atoms: return 'MC_MC'
                elif a1_name in mc_atoms or a2_name in mc_atoms: return 'SC_MC'
                else: return 'SC_SC'

            if len(res) > 0 and res[0]:
                for hb in res[0]:
                    n1, _, _ = parse_prody_node(hb[0], hb[2])
                    n2, _, _ = parse_prody_node(hb[3], hb[5])
                    if n1 and n2:
                        a1 = hb[1].split('_')[0] if hb[1] else ''
                        a2 = hb[4].split('_')[0] if hb[4] else ''
                        ctype = get_prody_contact_type(a1, a2)
                        add_edge(n1, n2, f'HBOND:{ctype}', hb[6], hb[7], a1, a2)
            if len(res) > 1 and res[1]:
                for sb in res[1]:
                    n1, _, _ = parse_prody_node(sb[0], sb[2])
                    n2, _, _ = parse_prody_node(sb[3], sb[5])
                    if n1 and n2:
                        a1 = sb[1].split('_')[0] if sb[1] else ''
                        a2 = sb[4].split('_')[0] if sb[4] else ''
                        ctype = get_prody_contact_type(a1, a2)
                        add_edge(n1, n2, f'IONIC:{ctype}', sb[6], 'nan', a1, a2)
            if len(res) > 3 and res[3]:
                for pi in res[3]:
                    n1, _, _ = parse_prody_node(pi[0], pi[2])
                    n2, _, _ = parse_prody_node(pi[3], pi[5])
                    if n1 and n2:
                        add_edge(n1, n2, 'PIPISTACK:SC_SC', pi[6], 'nan')
            if len(res) > 4 and res[4]:
                for pc in res[4]:
                    n1, _, _ = parse_prody_node(pc[0], pc[2])
                    n2, _, _ = parse_prody_node(pc[3], pc[5])
                    if n1 and n2:
                        add_edge(n1, n2, 'PICATION:SC_SC', pc[6], 'nan')
            if len(res) > 5 and res[5]:
                for hy in res[5]:
                    n1, _, _ = parse_prody_node(hy[0], hy[2])
                    n2, _, _ = parse_prody_node(hy[3], hy[5])
                    if n1 and n2:
                        add_edge(n1, n2, 'VDW:SC_SC', hy[6], 'nan')
            if len(res) > 6 and res[6]:
                for ss in res[6]:
                    n1, _, _ = parse_prody_node(ss[0], ss[2])
                    n2, _, _ = parse_prody_node(ss[3], ss[5])
                    if n1 and n2:
                        add_edge(n1, n2, 'SSBOND:SC_SC', ss[6], 'nan', 'SG', 'SG')
            
            try:
                wb_res = prody.calcWaterBridges(prody_atoms)
                for bridge in wb_res:
                    if hasattr(bridge, 'proteins') and len(bridge.proteins) == 2:
                        a1, a2 = bridge.proteins
                        n1 = f"{a1.getChid()}:{a1.getResnum()}:_:{a1.getResname()}"
                        n2 = f"{a2.getChid()}:{a2.getResnum()}:_:{a2.getResname()}"
                        dist = np.linalg.norm(a1.getCoords() - a2.getCoords())
                        ctype = get_prody_contact_type(a1.getName(), a2.getName())
                        add_edge(n1, n2, f'WATERBRIDGE:{ctype}', dist, 'nan', a1.getName(), a2.getName())
                        # add nodes
                        nodes.add(orig_struct[0][a1.getChid()][(' ', a1.getResnum(), ' ')]) # Approximation, actual node adding is handled below
            except Exception as e:
                pass
            
            for a in atoms:
                if a.parent.id[0] != 'W':
                    nodes.add(a.parent)
                    

        except Exception as e:
            print(f"ProDy calculation failed: {e}. Falling back to standard method.")
            calc_method = "standard"

    pairs = []
    if calc_method == "voronoi":
        from scipy.spatial import Voronoi
        filtered_atoms = [a for a in atoms if a.parent.id[0] != 'W']
        coords = np.array([a.coord for a in filtered_atoms])
        if len(coords) >= 4:
            try:
                vor = Voronoi(coords)
                pairs_set = set()
                for ridge in vor.ridge_points:
                    pairs_set.add(tuple(sorted([ridge[0], ridge[1]])))
                for i, j in pairs_set:
                    a1, a2 = filtered_atoms[i], filtered_atoms[j]
                    if np.linalg.norm(a1.coord - a2.coord) <= 6.5:
                        pairs.append((a1, a2))
            except:
                calc_method = "standard"
                
    if calc_method == "standard":
        ns = NeighborSearch(atoms)
        pairs = ns.search_all(6.5)
        
    res_pairs = {}
    for a1, a2 in pairs:
        r1, r2 = a1.parent, a2.parent
        if r1 == r2: continue
        
        if r1.id[0] == 'W' or r2.id[0] == 'W': continue
        
        if abs(r1.id[1] - r2.id[1]) <= 1 and r1.parent == r2.parent:
            if a1.name in ['N', 'C', 'CA', 'O'] and a2.name in ['N', 'C', 'CA', 'O']:
                dist = np.linalg.norm(a1.coord - a2.coord)
                if dist < 2.0: continue
                
        key = tuple(sorted([r1, r2], key=lambda x: format_node_id(x)))
        if key not in res_pairs:
            res_pairs[key] = []
        res_pairs[key].append((a1, a2))
        
    for (r1, r2), atom_pairs in res_pairs.items():
        n1 = format_node_id(r1)
        n2 = format_node_id(r2)
        nodes.add(r1)
        nodes.add(r2)
        
        # Check specific interactions
        # 1. SSBOND
        if r1.resname == 'CYS' and r2.resname == 'CYS':
            if 'SG' in r1 and 'SG' in r2:
                dist = np.linalg.norm(r1['SG'].coord - r2['SG'].coord)
                if dist <= 2.5:
                    edges.append({
                        'NodeId1': n1, 'Interaction': 'SSBOND:SC_SC', 'NodeId2': n2,
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(r1['SG'], r2['SG']):.4f}" if calc_energy else 'nan', 'Atom1': 'SG', 'Atom2': 'SG',
                        'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })
                    continue
        
        # 2. IONIC (Salt bridge)
        pos = {'ARG': ['NH1', 'NH2'], 'LYS': ['NZ'], 'HIS': ['ND1', 'NE2']}
        neg = {'ASP': ['OD1', 'OD2'], 'GLU': ['OE1', 'OE2']}
        
        for ax, ay in atom_pairs:
            if format_node_id(ax.parent) > format_node_id(ay.parent):
                a1, a2 = ay, ax
            else:
                a1, a2 = ax, ay
                
            if r1.resname in pos and a1.name in pos[r1.resname] and r2.resname in neg and a2.name in neg[r2.resname]:
                dist = np.linalg.norm(a1.coord - a2.coord)
                if dist <= 4.0:
                    edges.append({
                        'NodeId1': n1, 'Interaction': 'IONIC:SC_SC', 'NodeId2': n2,
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                        'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })
                    break
            elif r2.resname in pos and a2.name in pos[r2.resname] and r1.resname in neg and a1.name in neg[r1.resname]:
                dist = np.linalg.norm(a1.coord - a2.coord)
                if dist <= 4.0:
                    edges.append({
                        'NodeId1': n1, 'Interaction': 'IONIC:SC_SC', 'NodeId2': n2,
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                        'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })
                    break

        # 3. PIPISTACK
        aromatics = ['PHE', 'TYR', 'TRP', 'HIS']
        if r1.resname in aromatics and r2.resname in aromatics:
            c1 = get_centroid(r1)
            c2 = get_centroid(r2)
            if c1 is not None and c2 is not None:
                dist = np.linalg.norm(c1 - c2)
                if dist <= 6.5:
                    edges.append({
                        'NodeId1': n1, 'Interaction': 'PIPISTACK:SC_SC', 'NodeId2': n2,
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': '', 'Atom2': '',
                        'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': 'T', 'Model': model_num
                    })

        # 4. PICATION
        if r1.resname in aromatics and r2.resname in pos:
            c1 = get_centroid(r1)
            c2_coords = [r2[a].coord for a in pos[r2.resname] if a in r2]
            if c1 is not None and c2_coords:
                c2_cm = np.mean(c2_coords, axis=0)
                dist = np.linalg.norm(c1 - c2_cm)
                if dist <= 5.0:
                    edges.append({
                        'NodeId1': n1, 'Interaction': 'PICATION:SC_SC', 'NodeId2': n2,
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': '', 'Atom2': '',
                        'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })
        elif r2.resname in aromatics and r1.resname in pos:
            c2 = get_centroid(r2)
            c1_coords = [r1[a].coord for a in pos[r1.resname] if a in r1]
            if c2 is not None and c1_coords:
                c1_cm = np.mean(c1_coords, axis=0)
                dist = np.linalg.norm(c2 - c1_cm)
                if dist <= 5.0:
                    edges.append({
                        'NodeId1': n1, 'Interaction': 'PICATION:SC_SC', 'NodeId2': n2,
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': '', 'Atom2': '',
                        'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })

        # 5. METAL, HALOGEN, PIHBOND, HBOND & VDW
        donors_acceptors = {'N', 'O', 'S'}
        halogens = {'F', 'CL', 'BR', 'I'}
        metals = {'ZN', 'MG', 'CA', 'FE', 'CU', 'MN', 'CO', 'NI', 'K', 'NA'}
        
        pihbond_found = False
        aromatics = ['PHE', 'TYR', 'TRP', 'HIS']
        if r1.resname in aromatics or r2.resname in aromatics:
            for r_arom, r_donor in [(r1, r2), (r2, r1)]:
                if r_arom.resname in aromatics:
                    c_arom = get_centroid(r_arom)
                    if c_arom is not None:
                        for a_donor in r_donor:
                            if getattr(a_donor, 'element', a_donor.name[0]).upper() in donors_acceptors:
                                dist = np.linalg.norm(c_arom - a_donor.coord)
                                if dist <= 4.0:
                                    edges.append({
                                        'NodeId1': n1, 'Interaction': f'PIHBOND:SC_{"MC" if a_donor.name in MC_ATOMS else "SC"}', 'NodeId2': n2,
                                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': '', 'Atom2': a_donor.name,
                                        'Donor': format_node_id(r_donor), 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                                    })
                                    pihbond_found = True
                                    break
                if pihbond_found: break

        hbond_found = False
        vdw_found = False
        min_vdw_dist = 999.0
        best_vdw_pair = None
        metal_found = False
        halogen_found = False
        
        for ax, ay in atom_pairs:
            if format_node_id(ax.parent) > format_node_id(ay.parent):
                a1, a2 = ay, ax
            else:
                a1, a2 = ax, ay
                
            dist = np.linalg.norm(a1.coord - a2.coord)
            el1 = getattr(a1, 'element', a1.name[0]).upper()
            el2 = getattr(a2, 'element', a2.name[0]).upper()
            
            # METAL Coordination
            if not metal_found:
                if (el1 in metals and el2 in donors_acceptors) or (el2 in metals and el1 in donors_acceptors):
                    if dist <= 3.0:
                        edges.append({
                            'NodeId1': format_node_id(a1.parent), 'Interaction': f'METAL:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                            'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                            'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                        })
                        metal_found = True
            
            # HALOGEN Bond
            if not halogen_found:
                if (el1 in halogens and el2 in donors_acceptors) or (el2 in halogens and el1 in donors_acceptors):
                    if dist <= 4.0:
                        edges.append({
                            'NodeId1': format_node_id(a1.parent), 'Interaction': f'HALOGEN:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                            'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                            'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                        })
                        halogen_found = True

            # Check HBOND
            if el1 in donors_acceptors and el2 in donors_acceptors:
                if dist <= 3.5:
                    angle_val = 150.000
                    if strict_angle:
                        valid, computed_angle = check_hbond_angle(a1, a2)
                        if not valid:
                            continue
                        angle_val = computed_angle
                        
                    edges.append({
                        'NodeId1': format_node_id(a1.parent), 'Interaction': f'HBOND:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                        'Distance': f"{dist:.3f}", 'Angle': f"{angle_val:.3f}", 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                        'Donor': format_node_id(a1.parent), 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })
                    hbond_found = True
                    
            # Check VDW
            vdw_dist = get_vdw_radius(a1.name) + get_vdw_radius(a2.name) + 0.5
            if dist <= vdw_dist:
                if dist < min_vdw_dist:
                    min_vdw_dist = dist
                    best_vdw_pair = (a1, a2)
                    
        if not hbond_found and not metal_found and not halogen_found and not pihbond_found and best_vdw_pair:
            a1, a2 = best_vdw_pair
            edges.append({
                'NodeId1': format_node_id(a1.parent), 'Interaction': f'VDW:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                'Distance': f"{min_vdw_dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
            })

    if calc_method != "insty":
        # Calculate Water Bridges manually
        waters = [a for a in atoms if a.parent.id[0] == 'W' and getattr(a, 'element', a.name[0]).upper() == 'O']
        if waters:
            donors_acceptors = {'N', 'O', 'S'}
            try:
                water_ns = NeighborSearch(atoms)
                for w_atom in waters:
                    close_atoms = water_ns.search(w_atom.coord, 3.5)
                    prot_atoms = [a for a in close_atoms if a.parent.id[0] != 'W' and getattr(a, 'element', a.name[0]).upper() in donors_acceptors]
                    if len(prot_atoms) >= 2:
                        for i in range(len(prot_atoms)):
                            for j in range(i+1, len(prot_atoms)):
                                pa1, pa2 = prot_atoms[i], prot_atoms[j]
                                r1, r2 = pa1.parent, pa2.parent
                                if r1 == r2: continue
                                n1 = format_node_id(r1)
                                n2 = format_node_id(r2)
                                dist = np.linalg.norm(pa1.coord - pa2.coord)
                                if format_node_id(r1) > format_node_id(r2):
                                    n1, n2 = n2, n1
                                    pa1, pa2 = pa2, pa1
                                edges.append({
                                    'NodeId1': n1, 'Interaction': f'WATERBRIDGE:{get_contact_type(pa1, pa2)}', 'NodeId2': n2,
                                    'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(pa1, pa2):.4f}" if calc_energy else 'nan', 'Atom1': pa1.name, 'Atom2': pa2.name,
                                    'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                                })
                                nodes.add(r1)
                                nodes.add(r2)
            except Exception as e:
                pass

    # Return edges df
    edges_df = pd.DataFrame(edges)
    # Reorder to match exact RING format
    cols = ['NodeId1', 'Interaction', 'NodeId2', 'Distance', 'Angle', 'Energy', 'Atom1', 'Atom2', 'Donor', 'Positive', 'Cation', 'Orientation', 'Model']
    if not edges_df.empty:
        edges_df = edges_df[cols]
        def extract_sort_key(node_id):
            parts = node_id.split(':')
            if len(parts) >= 2:
                chain = parts[0]
                try:
                    pos = int(parts[1])
                except ValueError:
                    pos = 0
                return (chain, pos)
            return (node_id, 0)
            
        if remove_multiples:
            edges_df['Distance_val'] = edges_df['Distance'].astype(float)
            edges_df = edges_df.sort_values('Distance_val')
            edges_df = edges_df.drop_duplicates(subset=['NodeId1', 'Interaction', 'NodeId2'], keep='first')
            edges_df = edges_df.drop(columns=['Distance_val'])
        else:
            edges_df = edges_df.drop_duplicates()
            
        # Ensure consistent order by Chain and numeric Position
        edges_df['sort_key'] = edges_df['NodeId1'].apply(extract_sort_key)
        edges_df = edges_df.sort_values(by=['sort_key'])
        edges_df = edges_df.drop(columns=['sort_key'])
        
    return edges_df, nodes

def build_nodes_df(nodes, edges_df, pdbFileName="", structure=None, model_num=1):
    G = nx.Graph()
    if not edges_df.empty:
        for idx, row in edges_df.iterrows():
            G.add_edge(row['NodeId1'], row['NodeId2'])
            
    dssp_dict = {}
    if structure is not None:
        coords = []
        res_list = []
        for res in structure.get_residues():
            if res.id[0] == 'W': continue
            try:
                n = res['N'].coord
                ca = res['CA'].coord
                c = res['C'].coord
                o = res['O'].coord
                coords.append([n, ca, c, o])
                res_list.append(res)
            except KeyError:
                pass
        if coords:
            try:
                import pydssp
                import numpy as np
                coord_arr = np.array(coords)
                ss = pydssp.assign(coord_arr, out_type='c3')
                for res, s in zip(res_list, ss):
                    dssp_dict[res] = s if s != '-' else 'C'
            except Exception as e:
                print(f"DSSP calculation failed: {e}")
            
    nodes_data = []
    for r in nodes:
        node_id = format_node_id(r)
        chain = r.parent.id
        pos = r.id[1]
        resname = r.resname
        ca_bfactor = 0.0
        x, y, z = 0.0, 0.0, 0.0
        if 'CA' in r:
            ca_bfactor = r['CA'].bfactor
            x, y, z = r['CA'].coord
        else:
            coords = [a.coord for a in r]
            if coords:
                c = np.mean(coords, axis=0)
                x, y, z = c[0], c[1], c[2]
            
        deg = G.degree(node_id) if G.has_node(node_id) else 0
        dssp_val = dssp_dict.get(r, ' ')
        node_type = 'LIG' if r.id[0].startswith('H') else 'RES'
        
        nodes_data.append({
            'NodeId': node_id, 'Chain': chain, 'Position': pos, 'Residue': resname,
            'Type': node_type, 'Dssp': dssp_val, 'Degree': deg, 'Bfactor_CA': f"{ca_bfactor:.3f}",
            'x': f"{x:.3f}", 'y': f"{y:.3f}", 'z': f"{z:.3f}", 'pdbFileName': f"{pdbFileName}#{pos}.{chain}", 'Model': model_num
        })
        
    df = pd.DataFrame(nodes_data)
    cols = ['NodeId', 'Chain', 'Position', 'Residue', 'Type', 'Dssp', 'Degree', 'Bfactor_CA', 'x', 'y', 'z', 'pdbFileName', 'Model']
    if not df.empty:
        # Sort by Chain, then by numeric Position
        df['Position_int'] = pd.to_numeric(df['Position'], errors='coerce').fillna(0).astype(int)
        df = df.sort_values(by=['Chain', 'Position_int'])
        df = df[cols]
    return df

def analyze_network(edges_df, nodes_df, prefix):
    if edges_df.empty:
        print("Empty network, skipping metrics.")
        return
        
    G = nx.Graph()
    for idx, row in edges_df.iterrows():
        G.add_edge(row['NodeId1'], row['NodeId2'])
        
    # Degree
    degree_dict = dict(G.degree())
    # Clustering Coefficient
    cc_dict = nx.clustering(G)
    # Betweenness
    bw_dict = nx.betweenness_centrality(G)
    # Eigenvector Centrality
    try:
        ec_dict = nx.eigenvector_centrality(G, max_iter=1000)
    except:
        ec_dict = {n: 0.0 for n in G.nodes()}
        
    metrics = []
    for node in G.nodes():
        metrics.append({
            'NodeId': node,
            'Degree': degree_dict.get(node, 0),
            'Clustering_Coefficient': cc_dict.get(node, 0.0),
            'Betweenness_Centrality': bw_dict.get(node, 0.0),
            'Eigenvector_Centrality': ec_dict.get(node, 0.0)
        })
        
    df = pd.DataFrame(metrics)
    
    # Output full table
    df.to_csv(f"{prefix}_network_metrics.tsv", sep='\t', index=False)
    
    # Top 25 tables
    top_deg = df.nlargest(25, 'Degree').reset_index(drop=True)
    top_cc = df.nlargest(25, 'Clustering_Coefficient').reset_index(drop=True)
    top_bw = df.nlargest(25, 'Betweenness_Centrality').reset_index(drop=True)
    top_ec = df.nlargest(25, 'Eigenvector_Centrality').reset_index(drop=True)
    
    N = len(top_deg)
    if N > 0:
        merged_top25 = pd.DataFrame({
            'Rank': range(1, N + 1),
            'Degree_NodeId': top_deg['NodeId'],
            'Degree': top_deg['Degree'],
            'Betweenness_NodeId': top_bw['NodeId'],
            'Betweenness': top_bw['Betweenness_Centrality'],
            'Clustering_NodeId': top_cc['NodeId'],
            'Clustering': top_cc['Clustering_Coefficient'],
            'Eigenvector_NodeId': top_ec['NodeId'],
            'Eigenvector': top_ec['Eigenvector_Centrality']
        })
        merged_top25.to_csv(f"{prefix}_top25_metrics.tsv", sep='\t', index=False)
        print(f"Network metrics saved to {prefix}_network_metrics.tsv and {prefix}_top25_metrics.tsv")
    return degree_dict

def main():
    parser = argparse.ArgumentParser(description="Calculate RIN using RING 4.0 thresholds")
    parser.add_argument("input", help="Input PDB or CIF file")
    parser.add_argument("--add-h", action="store_true", help="Add hydrogens using pdb2pqr")
    parser.add_argument("--strict-angle", action="store_true", help="Enforce strict angle constraints for Hydrogen Bonds (e.g. >120 deg)")
    parser.add_argument("--remove-multiples", action="store_true", help="Remove multiple interactions of the same type between the same residue pair")
    parser.add_argument("--calc-method", choices=["standard", "voronoi", "insty", "rinpy"], default="standard", help="Method for RIN calculation")
    parser.add_argument("--chains", type=str, help="Comma-separated list of chains to calculate (e.g. A,B,C)")
    parser.add_argument("--energy", action="store_true", help="Calculate Interaction Energy (PENs). Requires a PQR file or --add-h.")
    parser.add_argument("--vdw-energy", action="store_true", help="Include generic Lennard-Jones VDW term in energy calculation. If not set, only Electrostatics is computed.")
    args = parser.parse_args()
    
    input_file = args.input
    base_prefix = os.path.splitext(os.path.basename(input_file))[0]
    
    orig_struct = parse_structure(input_file)
    models = list(orig_struct.get_models())
    from Bio.PDB import PDBIO, Structure
    io = PDBIO()
    
    for i, orig_model in enumerate(models):
        model_num = orig_model.id if orig_model.id is not None and orig_model.id >= 0 else i + 1
        prefix = f"{base_prefix}_model_{model_num}" if len(models) > 1 else base_prefix
        
        temp_pdb = f"{prefix}_temp.pdb"
        temp_struct = Structure.Structure("temp")
        temp_struct.add(orig_model)
        io.set_structure(temp_struct)
        io.save(temp_pdb)
        
        if args.calc_method == 'rinpy':
            try:
                from rinpy import RINProcess
            except ImportError:
                RINProcess = None
            if RINProcess is None:
                print("RinPy is not installed. Please install it using 'pip install rinpy'.")
                sys.exit(1)
            
            import shutil
            
            rinpy_in = f"{prefix}_rinpy_in"
            rinpy_out = f"{prefix}_rinpy_out"
            
            os.makedirs(rinpy_in, exist_ok=True)
            shutil.copy(temp_pdb, os.path.join(rinpy_in, os.path.basename(temp_pdb)))
            
            proc = RINProcess(output_path=rinpy_out, input_path=rinpy_in, num_workers=1)
            try:
                proc.start_process()
            except Exception as e:
                print(f"RinPy calculation failed: {e}")
                sys.exit(1)
                
            print(f"RinPy calculation completed. Results are in {rinpy_out}")
            
            # Write Log
            import datetime
            with open(f"{prefix}_sauron.log", 'w') as logf:
                logf.write(f"Sauron RIN Calculation Log\n")
                logf.write(f"Timestamp: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                logf.write(f"Input File: {args.input}\n")
                logf.write(f"Method: RINPY\n")
                
            # Extract Top 10 and Hinges
            import json
            import glob
            
            rinpy_summary = {
                "degree": [],
                "betweenness": [],
                "closeness": [],
                "hinges": {}
            }
            
            def parse_rinpy_csv(filepath):
                if not os.path.exists(filepath): return []
                with open(filepath, 'r') as f:
                    lines = f.readlines()
                data = []
                for line in lines:
                    parts = line.strip().split()
                    if len(parts) >= 5:
                        val = float(parts[1])
                        res = parts[2]
                        chain = parts[3]
                        num = parts[4]
                        data.append({"Residue": f"{chain}:{num}:{res}", "Value": val})
                # sort descending
                data.sort(key=lambda x: x["Value"], reverse=True)
                return data[:10]
                
            rinpy_summary["degree"] = parse_rinpy_csv(os.path.join(rinpy_out, f"{prefix}_temp", f"{prefix}_temp_centrality_degree.csv"))
            rinpy_summary["betweenness"] = parse_rinpy_csv(os.path.join(rinpy_out, f"{prefix}_temp", f"{prefix}_temp_centrality_betweenness.csv"))
            rinpy_summary["closeness"] = parse_rinpy_csv(os.path.join(rinpy_out, f"{prefix}_temp", f"{prefix}_temp_centrality_closeness.csv"))
            
            hinge_files = glob.glob(os.path.join(rinpy_out, f"{prefix}_temp", "hinge_modes", "*_hinge_residues.txt"))
            for hf in hinge_files:
                mode_name = os.path.basename(hf).replace(f"{prefix}_temp_laplacian_", "").replace("_hinge_residues.txt", "")
                hinges = []
                with open(hf, 'r') as f:
                    for line in f:
                        parts = line.strip().split(';')
                        if len(parts) >= 3:
                            hinges.append(f"{parts[1]}:{parts[2]}:{parts[0]}")
                rinpy_summary["hinges"][mode_name] = hinges
                
            with open(f"{prefix}_rinpy_summary.json", "w") as f:
                json.dump(rinpy_summary, f)
                
            if os.path.exists(temp_pdb): os.remove(temp_pdb)
            if temp_pdb != temp_pdb and os.path.exists(temp_pdb): os.remove(temp_pdb)
            shutil.rmtree(rinpy_in)
            continue
        if args.add_h:
            pqr_file = add_hydrogens(temp_pdb)
            
            try:
                struct_to_calc = list(parse_structure(pqr_file).get_models())[0]
            except:
                struct_to_calc = orig_model
                
            if pqr_file != temp_pdb:
                for orig_chain in orig_model:
                    if orig_chain.id not in struct_to_calc: continue
                    chain = struct_to_calc[orig_chain.id]
                    for orig_res in list(orig_chain):
                        if orig_res.id[0] != ' ' and orig_res.id[0] != 'W':
                            if orig_res.id not in chain:
                                orig_chain.detach_child(orig_res.id)
                                chain.add(orig_res)
        else:
            struct_to_calc = orig_model
            pqr_file = temp_pdb
            
        if args.chains:
            valid_chains = [c.strip() for c in args.chains.split(',')]
            for chain in list(struct_to_calc):
                if chain.id not in valid_chains:
                    struct_to_calc.detach_child(chain.id)
            
        edges_df, nodes = calculate_rin(struct_to_calc, strict_angle=args.strict_angle, remove_multiples=args.remove_multiples, model_num=model_num, calc_method=args.calc_method, pqr_file=pqr_file, calc_energy=args.energy, include_vdw=args.vdw_energy)
        nodes_df = build_nodes_df(nodes, edges_df, os.path.basename(args.input), structure=struct_to_calc, model_num=model_num)
        
        edges_out = f"{prefix}.edges"
        nodes_out = f"{prefix}.nodes"
        
        if not edges_df.empty:
            edges_df.to_csv(edges_out, sep='\t', index=False)
        if not nodes_df.empty:
            nodes_df.to_csv(nodes_out, sep='\t', index=False)
            
        print(f"Generated {edges_out} and {nodes_out}")
        degree_dict = analyze_network(edges_df, nodes_df, prefix)
        
        # Write Log
        import datetime
        with open(f"{prefix}_sauron.log", 'w') as logf:
            logf.write(f"Sauron RIN Calculation Log\n")
            logf.write(f"Timestamp: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            logf.write(f"Input File: {args.input}\n")
            logf.write(f"Method: {args.calc_method}\n")
            logf.write(f"Strict Angle: {args.strict_angle}\n")
            logf.write(f"Remove Multiples: {args.remove_multiples}\n")
            logf.write(f"Energy (PEN): {args.energy}\n")
            logf.write(f"Include VDW: {args.vdw_energy}\n")
            logf.write(f"Chains: {args.chains if args.chains else 'ALL'}\n")
            logf.write(f"Total Nodes: {len(nodes)}\n")
            logf.write(f"Total Edges: {len(edges_df)}\n")

        # Write Degree PDB
        from Bio.PDB import PDBIO
        
        # Helper to format node ID for structure traversal
        def get_node_id_for_atom(res):
            chain_id = res.parent.id
            res_id = res.id[1]
            ins_code = res.id[2].strip() if res.id[2].strip() else '_'
            res_name = res.resname
            return f"{chain_id}:{res_id}:{ins_code}:{res_name}"

        # `struct_to_calc` contains the subset of chains or the whole structure
        # But `orig_struct` has the exact coordinates. We can annotate `orig_struct`.
        for model in orig_struct:
            for chain in model:
                for res in chain:
                    node_id = get_node_id_for_atom(res)
                    deg = degree_dict.get(node_id, 0)
                    for atom in res:
                        atom.set_bfactor(deg)
        io = PDBIO()
        io.set_structure(orig_struct)
        io.save(f"{prefix}_degree.pdb")
        
        if os.path.exists(temp_pdb): os.remove(temp_pdb)
        if pqr_file != temp_pdb and os.path.exists(pqr_file): os.remove(pqr_file)

if __name__ == "__main__":
    main()
