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

# VDW Radii (Bondi / CSD-derived generic values), keyed by *element*, not atom name.
VDW_RADII = {
    'C': 1.70, 'N': 1.55, 'O': 1.52, 'S': 1.80, 'H': 1.20, 'P': 1.80,
    'F': 1.47, 'CL': 1.75, 'BR': 1.85, 'I': 1.98,
    'ZN': 1.39, 'MG': 1.73, 'CA': 2.31, 'MN': 1.73, 'FE': 1.72, 'CU': 1.40,
    'NI': 1.63, 'CO': 1.53, 'NA': 2.27, 'K': 2.75, 'CD': 1.58, 'HG': 1.55,
    'LI': 1.82, 'CS': 3.43,
}

# Effective metal-ion coordination distance cutoffs (Angstrom), approximating typical
# M-O/N/S coordination bond lengths (Harding 2001; CheckMyMetal, Zheng et al. 2017).
# A flat 3.0 A cutoff (the previous behaviour) systematically misses coordination
# bonds to larger cations such as K+ or Ca2+, so these are widened individually.
METAL_COORD_DIST = {
    'NA': 3.0, 'K': 3.4, 'MG': 2.6, 'CA': 3.2, 'MN': 2.7, 'FE': 2.7,
    'CO': 2.6, 'NI': 2.5, 'CU': 2.6, 'ZN': 2.7, 'CD': 2.9, 'HG': 3.0,
    'LI': 2.6, 'CS': 3.6,
}
DEFAULT_METAL_DIST = 3.0

MC_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}

# Two-letter element symbols that commonly appear in PDB/PQR atom names (metals and
# halogens). Used to disambiguate atom names like "CA" (alpha-carbon vs. calcium ion)
# or "CL"/"NA" when the file's own element column is missing or was stripped (e.g.
# after pdb2pqr, which does not carry the standard PDB element field).
TWO_LETTER_ELEMENTS = set(VDW_RADII.keys()) & {
    'CL', 'BR', 'ZN', 'MG', 'CA', 'MN', 'FE', 'CU', 'NI', 'CO', 'NA', 'K', 'CD', 'HG', 'LI', 'CS'
}

def get_element(atom):
    """Best-effort chemical element for an atom, preferring the parsed PDB element
    field and only falling back to name-based guessing (correctly handling two-letter
    symbols) when that field is absent."""
    el = (getattr(atom, 'element', '') or '').strip().upper()
    if el and el != 'X':
        return el
    name = atom.name.strip().upper()
    if name[:2] in TWO_LETTER_ELEMENTS:
        return name[:2]
    return name[0] if name else 'C'

def get_vdw_radius(element):
    return VDW_RADII.get(element, 1.70)

def get_metal_coord_dist(element):
    return METAL_COORD_DIST.get(element, DEFAULT_METAL_DIST)

# Hydrogen-bond donor/acceptor atom roles for the 20 standard amino acids. A blanket
# "any N/O/S near any N/O/S" test (the previous behaviour) also matches chemically
# impossible pairs such as two backbone carbonyl oxygens or two amide nitrogens -
# neither of which can hydrogen-bond to its own kind - which was the single largest
# source of over-counted HBOND edges relative to RING.
STANDARD_AA = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL'
}
HBOND_DONOR_ATOMS = {
    'ARG': {'NE', 'NH1', 'NH2'}, 'ASN': {'ND2'}, 'GLN': {'NE2'},
    'HIS': {'ND1', 'NE2'}, 'LYS': {'NZ'}, 'SER': {'OG'}, 'THR': {'OG1'},
    'TYR': {'OH'}, 'TRP': {'NE1'}, 'CYS': {'SG'},
}
HBOND_ACCEPTOR_ATOMS = {
    'ASN': {'OD1'}, 'ASP': {'OD1', 'OD2'}, 'GLN': {'OE1'}, 'GLU': {'OE1', 'OE2'},
    'HIS': {'ND1', 'NE2'}, 'SER': {'OG'}, 'THR': {'OG1'}, 'TYR': {'OH'},
    'CYS': {'SG'}, 'MET': {'SD'},
}

def get_hbond_role(atom):
    """Returns (is_donor, is_acceptor) for an atom. Standard amino acids use their
    known backbone (N donor except Pro; O/OXT acceptor) and side-chain chemistry.
    Anything else (ligands, nucleic acids, unknown HETATMs) falls back to a
    permissive "any polar atom can play either role" rule, since we don't have a
    chemical-component dictionary to type them precisely."""
    res = atom.parent
    resname = res.resname
    name = atom.name
    if resname not in STANDARD_AA:
        el = get_element(atom)
        is_polar = el in ('N', 'O', 'S')
        return is_polar, is_polar
    is_donor = (name == 'N' and resname != 'PRO') or name in HBOND_DONOR_ATOMS.get(resname, ())
    is_acceptor = name in ('O', 'OXT') or name in HBOND_ACCEPTOR_ATOMS.get(resname, ())
    return is_donor, is_acceptor

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

def get_ring_normal(res):
    """Unit normal vector of an aromatic side-chain ring, from the plane best fit
    through its ring atoms (first three ring atoms present is enough for planar
    aromatic rings)."""
    if res.resname in ['PHE', 'TYR']:
        atoms = ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ']
    elif res.resname == 'TRP':
        atoms = ['CG', 'CD1', 'CD2', 'NE1', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2']
    elif res.resname == 'HIS':
        atoms = ['CG', 'ND1', 'CD2', 'CE1', 'NE2']
    else:
        return None
    coords = np.array([res[a].coord for a in atoms if a in res])
    if len(coords) < 3: return None
    centroid = coords.mean(axis=0)
    # SVD of the centered ring coordinates: the singular vector with the smallest
    # singular value is the plane normal (robust to which 3 atoms one might pick).
    _, _, vh = np.linalg.svd(coords - centroid)
    return vh[2]

def classify_pipistack_orientation(r1, r2):
    """Classify a pi-pi stacking contact as Parallel ('P', roughly face-to-face or
    offset-stacked, normals aligned) or T-shaped ('T', edge-to-face, normals ~perpendicular)."""
    n1 = get_ring_normal(r1)
    n2 = get_ring_normal(r2)
    if n1 is None or n2 is None:
        return 'T'
    cosang = np.clip(abs(np.dot(n1, n2)), -1.0, 1.0)
    angle = np.degrees(np.arccos(cosang))
    return 'P' if angle <= 45 else 'T'

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

def calculate_rin(structure, strict_angle=False, remove_multiples=False, model_num=1, calc_method="standard", pqr_file=None, calc_energy=False, include_vdw=False, water_bridges=False):
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

    # The +0.5 A VDW tolerance below is meant to compensate for a heavy-atom-only
    # structure's *missing* hydrogens (an implicit H's own reach isn't otherwise
    # accounted for). When the structure already carries explicit hydrogens (NMR
    # ensembles, or after --add-h), that compensation is no longer needed and just
    # widens the net on top of already-real atoms - verified against a 20-model NMR
    # ensemble (2N8N) with reference RING output: at the default 0.5 A tolerance,
    # Sauron called a VDW contact on 97 residue pairs in one representative model
    # vs RING's 42; dropping to a small residual (thermal-motion-only) tolerance
    # once hydrogens are present brings that back in line.
    has_explicit_h = any(get_element(a) == 'H' for a in atoms)
    # With explicit hydrogens, the classic "sum of vdW radii" touching-distance
    # definition (no added slack) already matches RING closely - empirically tuned
    # against the 2N8N ensemble (VDW ratio dropped from 1.73x at +0.5 A to ~1.1x at
    # +0.0 A; more negative isn't a physically meaningful "contact" anymore).
    VDW_TOLERANCE = 0.0 if has_explicit_h else 0.5

    if calc_method == "insty":
        try:
            import prody
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

        # Exclude "trivial" local backbone contacts. Any two residues within the
        # same chain closer than 3 apart in sequence are geometrically forced into
        # contact by the fixed peptide-bond geometry alone (e.g. O(i)...N(i+1) sits
        # at ~2.2 A purely from the C=O/C-N bond angle, well inside every distance
        # cutoff below) - this is not a meaningful non-covalent interaction. RING
        # never reports edges for |i-j|<3 (confirmed against reference RING output:
        # zero HBOND/VDW edges at sequence separation 1 or 2); without this filter
        # every peptide bond in the structure was being mis-reported as an HBOND.
        if r1.parent == r2.parent and abs(r1.id[1] - r2.id[1]) < 3:
            continue

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
        # Both the cation and the anion are treated as point groups rather than as
        # loose "any qualifying atom pair" contacts: the cation is represented by its
        # single charge-bearing atom (Arg CZ, Lys NZ - the guanidinium/ammonium
        # center) and the anion by the centroid of its two carboxylate oxygens. Using
        # individual atom pairs instead (the previous approach) let salt bridges fire
        # from whichever single terminal atom happened to be closest, over-counting
        # relative to this stricter group-center definition. His is excluded as a
        # cation here since its charge is tautomer/pH-dependent and ambiguous without
        # explicit protonation state.
        CATION_ATOM = {'ARG': 'CZ', 'LYS': 'NZ'}
        ANION_ATOMS = {'ASP': ('OD1', 'OD2'), 'GLU': ('OE1', 'OE2')}

        def ionic_group_center(res, spec):
            if isinstance(spec, str):
                return res[spec].coord if spec in res else None
            coords = [res[a].coord for a in spec if a in res]
            return np.mean(coords, axis=0) if coords else None

        cation_res, anion_res = None, None
        if r1.resname in CATION_ATOM and r2.resname in ANION_ATOMS:
            cation_res, anion_res = r1, r2
        elif r2.resname in CATION_ATOM and r1.resname in ANION_ATOMS:
            cation_res, anion_res = r2, r1
        if cation_res is not None:
            cation_atom = cation_res[CATION_ATOM[cation_res.resname]]
            cp = cation_atom.coord
            ap = ionic_group_center(anion_res, ANION_ATOMS[anion_res.resname])
            if ap is not None:
                dist = np.linalg.norm(cp - ap)
                if dist <= 4.0:
                    # For the energy term (if requested), use the actual closer anion
                    # oxygen rather than its centroid, since calc_energy_terms needs
                    # real atoms (occupancy/bfactor) to work with.
                    anion_atom = min(
                        (anion_res[a] for a in ANION_ATOMS[anion_res.resname] if a in anion_res),
                        key=lambda a: np.linalg.norm(a.coord - cp)
                    )
                    edges.append({
                        'NodeId1': n1, 'Interaction': 'IONIC:SC_SC', 'NodeId2': n2,
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(cation_atom, anion_atom):.4f}" if calc_energy else 'nan', 'Atom1': cation_atom.name, 'Atom2': anion_atom.name,
                        'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })

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
                        'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': classify_pipistack_orientation(r1, r2), 'Model': model_num
                    })

        # 4. PICATION (same single-atom cation convention as IONIC above)
        pication_pair = None
        if r1.resname in aromatics and r2.resname in CATION_ATOM:
            pication_pair = (r1, r2)
        elif r2.resname in aromatics and r1.resname in CATION_ATOM:
            pication_pair = (r2, r1)
        if pication_pair is not None:
            r_arom, r_cation = pication_pair
            c_arom = get_centroid(r_arom)
            c_cation = ionic_group_center(r_cation, CATION_ATOM[r_cation.resname])
            if c_arom is not None and c_cation is not None:
                dist = np.linalg.norm(c_arom - c_cation)
                if dist <= 5.0:
                    # A favorable cation-pi contact has the cation sitting roughly
                    # above the ring face (along its normal), not off to the side in
                    # the ring plane - a distance-only check (the previous behaviour)
                    # also fires on off-axis geometries a true cation-pi interaction
                    # wouldn't count. Angle is measured from the ring normal (0 deg =
                    # directly above the ring center).
                    normal = get_ring_normal(r_arom)
                    angle_ok = True
                    cp_angle = 'nan'
                    if normal is not None:
                        v = c_cation - c_arom
                        vnorm = np.linalg.norm(v)
                        if vnorm > 1e-6:
                            cosang = np.clip(abs(np.dot(v / vnorm, normal)), -1.0, 1.0)
                            cp_angle = np.degrees(np.arccos(cosang))
                            angle_ok = cp_angle <= 45.0
                    if angle_ok:
                        edges.append({
                            'NodeId1': n1, 'Interaction': 'PICATION:SC_SC', 'NodeId2': n2,
                            'Distance': f"{dist:.3f}", 'Angle': f"{cp_angle:.3f}" if cp_angle != 'nan' else 'nan', 'Energy': 'nan', 'Atom1': '', 'Atom2': '',
                            'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                        })

        # 5. METAL, HALOGEN, PIHBOND, HBOND & VDW
        donors_acceptors = {'N', 'O', 'S'}
        halogens = {'F', 'CL', 'BR', 'I'}
        metals = set(METAL_COORD_DIST.keys())

        pihbond_found = False
        aromatics = ['PHE', 'TYR', 'TRP', 'HIS']
        if r1.resname in aromatics or r2.resname in aromatics:
            for r_arom, r_donor in [(r1, r2), (r2, r1)]:
                if r_arom.resname in aromatics:
                    c_arom = get_centroid(r_arom)
                    if c_arom is not None:
                        for a_donor in r_donor:
                            # Require a genuine H-bond donor (an atom that actually
                            # carries a polar hydrogen), not just any N/O/S - the same
                            # "any polar atom qualifies" mistake that over-counted HBOND.
                            is_donor, _ = get_hbond_role(a_donor)
                            if is_donor:
                                dist = np.linalg.norm(c_arom - a_donor.coord)
                                # Heavy-atom-to-ring-centroid distance; ~1 A longer than
                                # the H...centroid geometric criterion since no explicit
                                # hydrogen position is assumed (Steiner & Koellner, 2001).
                                if dist <= 4.0:
                                    edges.append({
                                        'NodeId1': n1, 'Interaction': f'PIHBOND:SC_{"MC" if a_donor.name in MC_ATOMS else "SC"}', 'NodeId2': n2,
                                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': '', 'Atom2': a_donor.name,
                                        'Donor': format_node_id(r_donor), 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                                    })
                                    pihbond_found = True
                                    break
                if pihbond_found: break

        # METAL/HALOGEN/HBOND are each independent (at most one edge per category, per
        # residue pair). VDW is instead treated as the *fallback* contact: it is only
        # reported when the pair has no more specific typed interaction at all. RING
        # does occasionally report VDW alongside a typed interaction on the same pair
        # (usually because it places explicit hydrogens, exposing a separate close
        # H...atom contact we cannot see on a heavy-atom-only structure), but validating
        # against reference RING output showed the fallback rule reproduces RING's VDW
        # counts far more closely than treating VDW as fully independent, which -
        # without RING's explicit-hydrogen geometry to naturally limit it - roughly
        # doubled true VDW edges by re-reporting the same heavy-atom contact under a
        # second label for nearly every already-typed pair.
        hbond_found = False
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
            el1 = get_element(a1)
            el2 = get_element(a2)

            # METAL Coordination (per-metal distance cutoff, see METAL_COORD_DIST)
            if not metal_found:
                if el1 in metals and el2 in donors_acceptors:
                    metal_el, other_el = el1, el2
                elif el2 in metals and el1 in donors_acceptors:
                    metal_el, other_el = el2, el1
                else:
                    metal_el = None
                if metal_el is not None and dist <= get_metal_coord_dist(metal_el):
                    edges.append({
                        'NodeId1': format_node_id(a1.parent), 'Interaction': f'METAL:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                        'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })
                    metal_found = True

            # HALOGEN Bond (sum of vdW radii + small tolerance, rather than a flat cutoff,
            # so it scales correctly across F/Cl/Br/I instead of under-counting the larger halogens)
            if not halogen_found:
                if (el1 in halogens and el2 in donors_acceptors) or (el2 in halogens and el1 in donors_acceptors):
                    halogen_dist = get_vdw_radius(el1) + get_vdw_radius(el2) + 0.3
                    if dist <= halogen_dist:
                        edges.append({
                            'NodeId1': format_node_id(a1.parent), 'Interaction': f'HALOGEN:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                            'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                            'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                        })
                        halogen_found = True

            # Check HBOND: require a genuine donor/acceptor pairing (not just "both
            # polar"), matching the atom's actual role in its residue (see get_hbond_role).
            if dist <= 3.5:
                d1, acc1 = get_hbond_role(a1)
                d2, acc2 = get_hbond_role(a2)
                donor, acceptor = (None, None)
                if d1 and acc2:
                    donor, acceptor = a1, a2
                elif d2 and acc1:
                    donor, acceptor = a2, a1
                if donor is not None:
                    angle_val = 150.000
                    if strict_angle:
                        valid, computed_angle = check_hbond_angle(donor, acceptor)
                        if not valid:
                            continue
                        angle_val = computed_angle

                    edges.append({
                        'NodeId1': format_node_id(a1.parent), 'Interaction': f'HBOND:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                        'Distance': f"{dist:.3f}", 'Angle': f"{angle_val:.3f}", 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                        'Donor': format_node_id(donor.parent), 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })
                    hbond_found = True

            # Track the closest VDW-qualifying atom pair (element-correct radii; see get_element)
            vdw_dist = get_vdw_radius(el1) + get_vdw_radius(el2) + VDW_TOLERANCE
            if dist <= vdw_dist and dist < min_vdw_dist:
                min_vdw_dist = dist
                best_vdw_pair = (a1, a2)

        # VDW is the fallback contact type (see note above): only reported when this
        # residue pair has no more specific interaction already.
        if not hbond_found and not metal_found and not halogen_found and not pihbond_found and best_vdw_pair:
            a1, a2 = best_vdw_pair
            edges.append({
                'NodeId1': format_node_id(a1.parent), 'Interaction': f'VDW:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                'Distance': f"{min_vdw_dist:.3f}", 'Angle': 'nan', 'Energy': f"{calc_energy_terms(a1, a2):.4f}" if calc_energy else 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
            })

    # Water-mediated bridges are a Sauron-specific extension - RING has no equivalent
    # interaction type - so they're opt-in (--water-bridges) rather than computed by
    # default. A water molecule has at most ~4 usable H-bond positions (2 donor
    # hydrogens, 2 acceptor lone pairs), so each water is limited to its WATER_BRIDGE_MAX_PARTNERS
    # closest genuine donor/acceptor atoms (get_hbond_role, not "any N/O/S") before
    # pairing them up - reporting every combination among *all* nearby polar atoms
    # (the previous behaviour) let a single crowded water spawn many edges between
    # residues that may not be simultaneously bridged by it at all.
    WATER_BRIDGE_MAX_PARTNERS = 4
    if water_bridges and calc_method != "insty":
        waters = [a for a in atoms if a.parent.id[0] == 'W' and get_element(a) == 'O']
        if waters:
            try:
                water_ns = NeighborSearch(atoms)
                for w_atom in waters:
                    close_atoms = water_ns.search(w_atom.coord, 3.5)
                    prot_atoms = [
                        a for a in close_atoms
                        if a.parent.id[0] != 'W' and any(get_hbond_role(a))
                    ]
                    if len(prot_atoms) < 2:
                        continue
                    prot_atoms.sort(key=lambda a: np.linalg.norm(a.coord - w_atom.coord))
                    prot_atoms = prot_atoms[:WATER_BRIDGE_MAX_PARTNERS]
                    for i in range(len(prot_atoms)):
                        for j in range(i + 1, len(prot_atoms)):
                            pa1, pa2 = prot_atoms[i], prot_atoms[j]
                            r1, r2 = pa1.parent, pa2.parent
                            if r1 == r2:
                                continue
                            n1 = format_node_id(r1)
                            n2 = format_node_id(r2)
                            # "Distance" is the total water-mediated path length
                            # (atom1-water + water-atom2), not the direct atom1-atom2
                            # distance - the two bridged atoms can otherwise sit on
                            # opposite sides of the water, up to ~7 A apart, which
                            # read as a nonsensically long "contact".
                            leg1 = np.linalg.norm(pa1.coord - w_atom.coord)
                            leg2 = np.linalg.norm(pa2.coord - w_atom.coord)
                            dist = leg1 + leg2
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
    parser.add_argument("--water-bridges", action="store_true", help="Detect water-mediated bridges (WATERBRIDGE edges) between residues sharing a crystallographic water. This is a Sauron-specific extension with no equivalent RING interaction type, off by default.")
    args = parser.parse_args()
    
    input_file = args.input
    base_prefix = os.path.splitext(os.path.basename(input_file))[0]
    
    orig_struct = parse_structure(input_file)
    models = list(orig_struct.get_models())
    from Bio.PDB import PDBIO, Structure
    io = PDBIO()
    
    for i, orig_model in enumerate(models):
        # Biopython always assigns Model.id starting at 0, regardless of the actual
        # model numbers in the source file (e.g. mmCIF's own pdbx_PDB_model_num,
        # which is 1-indexed) - using it directly mislabeled every model off-by-one
        # against the file itself and against tools like RING (verified against a
        # 20-model NMR ensemble: Sauron's "model_0..19" lined up exactly with RING's
        # "Model 1..20"). The loading order (1-indexed) is what actually matches.
        model_num = i + 1
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
            
        edges_df, nodes = calculate_rin(struct_to_calc, strict_angle=args.strict_angle, remove_multiples=args.remove_multiples, model_num=model_num, calc_method=args.calc_method, pqr_file=pqr_file, calc_energy=args.energy, include_vdw=args.vdw_energy, water_bridges=args.water_bridges)
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
            logf.write(f"Water Bridges: {args.water_bridges}\n")
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
