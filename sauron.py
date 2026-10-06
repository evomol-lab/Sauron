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

# Van der Waals contact classes follow the all-atom contact analysis of Word et
# al. (1999, J. Mol. Biol. 285:1735), used by MolProbity/Probe: by the gap
# between the two vdW surfaces, "close contact" (0 < gap <= 0.25 A), "wide
# contact" (0.25 < gap <= 0.5 A), "small overlap" (overlap < 0.4 A) and "bad
# overlap" (overlap >= 0.4 A). A VDW edge needs a close contact or a small
# overlap (--vdw-contacts close, the default) or, optionally, also accepts wide
# contacts (--vdw-contacts wide). A bad overlap is a steric clash (Chen et al.
# 2010, Acta Cryst. D66:12) - a structural defect, not an interaction - so it
# is reported separately and never counted as a VDW edge. Radii are Bondi
# (1964), see VDW_RADII.
VDW_CONTACT_MAX_GAP = {'close': 0.25, 'wide': 0.5}
DEFAULT_VDW_CONTACTS = 'close'
CLASH_MIN_OVERLAP = 0.4
# Heavy atoms of different residues closer than this are covalently linked
# (glycosylation, covalent inhibitors, cross-links), not in non-bonded contact.
COVALENT_LINK_MAX = 2.2

# Halogen bonds are directional: the C-X...A angle is close to linear (~165 deg)
# and X...A is shorter than the sum of vdW radii (Auffinger et al. 2004, PNAS
# 101:16789; IUPAC definition, Desiraju et al. 2013, Pure Appl. Chem. 85:1711).
HALOGEN_MIN_ANGLE = 140.0

# X-H...pi hydrogen bonds: the donor's hydrogen points at the aromatic ring
# (Steiner & Koellner 2001, J. Mol. Biol. 305:535).
PIHBOND_MAX_H_CENTROID = 3.0
PIHBOND_MIN_ANGLE = 120.0
# Legacy heavy-atom-only criterion, used only when hydrogen geometry is not checked.
PIHBOND_MAX_DONOR_CENTROID = 4.0

# Residue-residue interaction energy (--energy): Coulomb with a distance-dependent
# dielectric eps(r) = 4r and PQR partial charges, plus (--vdw-energy) a generic
# 12-6 Lennard-Jones term with a uniform well depth and r_min = sum of Bondi radii.
COULOMB_CONSTANT = 332.0637  # kcal*A/(mol*e^2)
LJ_WELL_DEPTH = 0.1          # kcal/mol, generic

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
# neither of which can hydrogen-bond to its own kind. Roles follow standard
# amino-acid chemistry (e.g. Baker & Hubbard 1984; McDonald & Thornton 1994).
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

def check_hbond_angle(donor, acceptor):
    """Best donor-H...acceptor angle over the donor's own hydrogens; valid if > 120 deg.
    Only the donor's hydrogens are tried: testing the acceptor's hydrogens as well
    (the previous behaviour) accepted bonds in the reverse direction even when the
    other atom cannot act as an acceptor (e.g. a backbone N)."""
    best = np.nan
    for h in get_attached_hydrogens(donor):
        angle = calc_angle(donor.coord, h.coord, acceptor.coord)
        if np.isnan(best) or angle > best:
            best = angle
    if not np.isnan(best) and best > 120:
        return True, best
    return False, best

def polar_hydrogen_parent(atom):
    """The N/O/S atom a hydrogen is bonded to, or None (non-hydrogen, or C-H)."""
    if get_element(atom) != 'H':
        return None
    for a in atom.parent:
        if a is not atom and get_element(a) in ('N', 'O', 'S') and np.linalg.norm(a.coord - atom.coord) < 1.3:
            return a
    return None

def load_pqr_parameters(model, pqr_file):
    """Attach PQR partial charge and radius to each atom (atom.xtra['charge'],
    atom.xtra['pqr_radius']) and return how many atoms were matched.

    PQR's charge and radius are whitespace-delimited and wider than the PDB
    occupancy/B-factor columns, so reading them through PDBParser truncates the
    charges and loses most radii; they are read from the file directly instead."""
    params = {}
    with open(pqr_file) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            fields = line.split()
            try:
                key = (line[21], int(line[22:26]), line[26].strip(), line[12:16].strip())
                params[key] = (float(fields[-2]), float(fields[-1]))
            except (ValueError, IndexError):
                continue
    matched = 0
    for atom in model.get_atoms():
        res = atom.parent
        key = (res.parent.id, res.id[1], res.id[2].strip(), atom.get_id())
        if key in params:
            atom.xtra['charge'], atom.xtra['pqr_radius'] = params[key]
            matched += 1
    return matched

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

def calculate_rin(structure, strict_angle=False, remove_multiples=False, model_num=1, calc_method="sauron", pqr_file=None, calc_energy=False, include_vdw=False, water_bridges=False, vdw_contacts=DEFAULT_VDW_CONTACTS):
    edges = []
    nodes = set()
    clashes = []
    vdw_max_gap = VDW_CONTACT_MAX_GAP[vdw_contacts]

    atoms = list(structure.get_atoms())
    # Whether hydrogen positions are available to check H-bond / X-H...pi geometry.
    has_explicit_h = any(get_element(a) == 'H' for a in atoms)

    energy_cache = {}

    def residue_pair_energy(res_a, res_b):
        """Non-bonded interaction energy (kcal/mol) between two residues, summed
        over ALL their atom pairs (a residue-residue energy, as in protein energy
        networks - e.g. Vijayabaskar & Vishveshwara 2010, Biophys. J. 99:3704 -
        rather than the energy of the single atom pair that defined an edge).
        NaN when either residue lacks PQR charges (e.g. a ligand pdb2pqr could not
        parameterize) or when the two residues are fewer than 3 apart along the
        chain, where bonded terms dominate and are not modelled."""
        key = frozenset((id(res_a), id(res_b)))
        if key in energy_cache:
            return energy_cache[key]
        value = np.nan
        atoms_a, atoms_b = list(res_a), list(res_b)
        if atoms_a and atoms_b and all('charge' in a.xtra for a in atoms_a + atoms_b):
            separation = sequence_separation(res_a, res_b)
            if separation is None or separation >= 3:
                coords_a = np.array([a.coord for a in atoms_a], dtype=float)
                coords_b = np.array([a.coord for a in atoms_b], dtype=float)
                d = np.linalg.norm(coords_a[:, None, :] - coords_b[None, :, :], axis=2)
                d = np.maximum(d, 1e-3)
                q_a = np.array([a.xtra['charge'] for a in atoms_a])
                q_b = np.array([a.xtra['charge'] for a in atoms_b])
                # atom pairs that count: everything except bonded neighbours
                # across a covalent link between the two residues
                keep = np.ones(d.shape, dtype=bool)
                if frozenset((id(res_a), id(res_b))) in linked_residue_pairs:
                    for i, x in enumerate(atoms_a):
                        for j, y in enumerate(atoms_b):
                            if frozenset((id(x), id(y))) in bonded_across_link:
                                keep[i, j] = False
                value = COULOMB_CONSTANT * float(np.sum((np.outer(q_a, q_b) / (4.0 * d * d))[keep]))
                if include_vdw:
                    # Polar hydrogens carry no Lennard-Jones term, as in the AMBER
                    # and CHARMM force fields (zero or negligible LJ parameters):
                    # in an H-bond they sit well inside the acceptor's vdW radius
                    # and their interaction is the electrostatic one.
                    r_a = np.array([get_vdw_radius(get_element(a)) for a in atoms_a])
                    r_b = np.array([get_vdw_radius(get_element(a)) for a in atoms_b])
                    lj_a = np.array([polar_h_parent.get(id(a)) is None for a in atoms_a])
                    lj_b = np.array([polar_h_parent.get(id(a)) is None for a in atoms_b])
                    ratio6 = ((r_a[:, None] + r_b[None, :]) / d) ** 6
                    lj = ratio6 ** 2 - 2 * ratio6
                    value += LJ_WELL_DEPTH * float(np.sum(lj[keep & np.outer(lj_a, lj_b)]))
        energy_cache[key] = value
        return value

    if calc_method == "insty":
        try:
            import prody
            import warnings
            import re
            warnings.filterwarnings('ignore')
            # ProDy's interaction code selects residues with strings like
            # "resid N and chain X", which breaks on crystal numbering: negative
            # residue numbers are parsed as a subtraction (crash) and insertion
            # codes (e.g. 58, 58A, 58B) get merged (crash or mixed-up edges).
            # ProDy's whitespace-based PQR reader also garbles residue numbers
            # that carry an insertion code. So ProDy gets its own PDB file, written
            # from the structure Sauron already parsed (hydrogens included), with
            # residues renumbered 1..N per chain and no insertion codes. Every
            # ProDy residue is then mapped back to its original node id.
            import tempfile
            prody_node_map = {}
            pdb_lines = []
            serial = 0
            for chain in structure:
                counter = 0
                for residue in chain:
                    counter += 1
                    prody_node_map[(chain.id, counter)] = format_node_id(residue)
                    record = "ATOM  " if residue.id[0] == ' ' else "HETATM"
                    for atom in residue:
                        serial += 1
                        x, y, z = atom.coord
                        element = (getattr(atom, 'element', '') or '').strip()
                        pdb_lines.append("%-6s%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n" % (
                            record, serial % 100000, atom.get_fullname()[:4], residue.resname[:3], chain.id[:1],
                            counter, x, y, z, 1.0, 0.0, element[:2]))
            pdb_lines.append("END\n")
            prody_tmp = tempfile.NamedTemporaryFile('w', suffix='.pdb', delete=False)
            prody_tmp.writelines(pdb_lines)
            prody_tmp.close()
            try:
                prody_atoms = prody.parsePDB(prody_tmp.name)
            finally:
                os.remove(prody_tmp.name)
            if not prody_atoms:
                raise ValueError("ProDy failed to parse structure.")

            res = prody.calcProteinInteractions(prody_atoms)
            
            def parse_prody_node(res_str, chid):
                m = re.match(r'^(.*?)(\-?\d+)$', res_str)
                if not m: return None, None, None
                resname = m.group(1)
                node_id = prody_node_map.get((chid, int(m.group(2))))
                if node_id is None: return None, None, None
                return node_id, resname, node_id.split(':')[1]

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
                # Only when requested, like the other methods (--water-bridges is off by default).
                wb_res = prody.calcWaterBridges(prody_atoms) if water_bridges else []
                for bridge in wb_res:
                    if hasattr(bridge, 'proteins') and len(bridge.proteins) == 2:
                        a1, a2 = bridge.proteins
                        n1 = prody_node_map.get((a1.getChid(), int(a1.getResnum())))
                        n2 = prody_node_map.get((a2.getChid(), int(a2.getResnum())))
                        if not n1 or not n2: continue
                        dist = np.linalg.norm(a1.getCoords() - a2.getCoords())
                        ctype = get_prody_contact_type(a1.getName(), a2.getName())
                        add_edge(n1, n2, f'WATERBRIDGE:{ctype}', dist, 'nan', a1.getName(), a2.getName())
            except Exception as e:
                pass
            
            for a in atoms:
                if a.parent.id[0] != 'W':
                    nodes.add(a.parent)
                    

        except Exception as e:
            print(f"ProDy calculation failed: {e}. Falling back to the Sauron method.")
            calc_method = "sauron"

    pairs = []
    if calc_method == "sauron":
        ns = NeighborSearch(atoms)
        pairs = ns.search_all(6.5)
        
    # Sequence separation from chain connectivity, not from residue numbers.
    # Crystal numbering often skips (e.g. papain-family 99 -> 102, or a
    # propeptide 62P -> 1) or repeats with insertion codes (0, 0A, 0B, ...), so
    # |number difference| can call a peptide-bonded neighbor "far" (and report
    # the 1.33 A C-N bond as a VDW:MC_MC contact) or call residues that are far
    # apart in the chain "adjacent" (and drop their real contacts). Peptide
    # bonds (C(i)-N(i+1), or O3'-P for nucleic acids) are found geometrically -
    # file order can't be trusted, e.g. modified residues dropped by pdb2pqr are
    # re-appended at the end of their chain - and each chain is split into
    # covalently connected segments; inside a segment the separation is the
    # number of links between the two residues.
    PEPTIDE_LINK_MAX = 2.0  # A; a real C-N peptide bond is ~1.33 A
    link_atoms = [a for a in atoms if a.parent.id[0] != 'W' and a.get_id() in ('C', 'N', "O3'", 'P')]
    next_res = {}   # id(residue) -> following residue
    prev_res = {}   # id(residue) -> preceding residue
    if link_atoms:
        for la, lb in NeighborSearch(link_atoms).search_all(PEPTIDE_LINK_MAX):
            ra, rb = la.parent, lb.parent
            if ra is rb or ra.parent is not rb.parent:
                continue
            names = (la.get_id(), lb.get_id())
            if names == ('N', 'C') or names == ('P', "O3'"):
                ra, rb = rb, ra
            elif names not in (('C', 'N'), ("O3'", 'P')):
                continue
            # ra's C (or O3') is bonded to rb's N (or P): ra precedes rb.
            if id(ra) not in next_res and id(rb) not in prev_res:
                next_res[id(ra)] = rb
                prev_res[id(rb)] = ra

    chain_position = {}   # id(residue) -> (id(segment head), position in segment)
    for chain in structure:
        for residue in chain:
            if residue.id[0] == 'W' or id(residue) in prev_res:
                continue
            # residue starts a segment (no preceding peptide bond); walk forward
            position = 0
            current = residue
            while current is not None and id(current) not in chain_position:
                chain_position[id(current)] = (id(residue), position)
                position += 1
                current = next_res.get(id(current))
    # Cyclic chains have no head; give any residue left over its own segment.
    for chain in structure:
        for residue in chain:
            if residue.id[0] != 'W' and id(residue) not in chain_position:
                chain_position[id(residue)] = (id(residue), 0)

    def sequence_separation(res_a, res_b):
        """Residues between two residues along the same chain; None if unrelated."""
        if res_a.parent is not res_b.parent:
            return None
        pos_a = chain_position.get(id(res_a))
        pos_b = chain_position.get(id(res_b))
        if pos_a is None or pos_b is None:
            return None
        if pos_a[0] == pos_b[0]:
            # Same covalently connected stretch: count the links between them.
            return abs(pos_a[1] - pos_b[1])
        # A hetero group with no polymer link (ligand, ion) has no sequence
        # position, whatever residue number it was given.
        for res in (res_a, res_b):
            if res.id[0] != ' ' and id(res) not in next_res and id(res) not in prev_res:
                return None
        # Across a chain break (unmodelled residues), numbering is the only
        # estimate of how many residues are missing.
        return abs(res_a.id[1] - res_b.id[1])

    res_pairs = {}
    for a1, a2 in pairs:
        r1, r2 = a1.parent, a2.parent
        if r1 == r2: continue
        
        if r1.id[0] == 'W' or r2.id[0] == 'W': continue

        # Exclude "trivial" local backbone contacts. Any two residues within the
        # same chain closer than 3 apart in sequence are geometrically forced into
        # contact by the fixed peptide-bond geometry alone (e.g. O(i)...N(i+1) sits
        # at ~2.2 A purely from the C=O/C-N bond angle, well inside every distance
        # cutoff below) - this is not a meaningful non-covalent interaction. The
        # minimum sequence separation of 3 is the usual convention for residue
        # interaction networks; without it every peptide bond in the structure
        # was being mis-reported as an HBOND.
        # Separation follows chain connectivity, not residue numbers (see above).
        separation = sequence_separation(r1, r2)
        if separation is not None and separation < 3:
            continue

        key = tuple(sorted([r1, r2], key=lambda x: format_node_id(x)))
        if key not in res_pairs:
            res_pairs[key] = []
        res_pairs[key].append((a1, a2))
        
    # Per-atom lookups used inside the pair loop (computed once).
    hbond_role = {id(a): get_hbond_role(a) for a in atoms}
    polar_h_parent = {id(a): polar_hydrogen_parent(a) for a in atoms if get_element(a) == 'H'}
    use_h_geometry = strict_angle and has_explicit_h

    # Covalent links between non-local residues (disulfides, glycosylation,
    # covalent ligands). As in force fields, atom pairs within three bonds of
    # each other across the link (1-2, 1-3 and 1-4 pairs) are bonded neighbours,
    # not non-bonded contacts: they are excluded from VDW, clashes and energies.
    def atoms_within_bonds(res, start, max_bonds):
        depth = {start: 0}
        frontier = [start]
        for level in range(1, max_bonds + 1):
            nxt = []
            for a in frontier:
                for b in res:
                    if b in depth:
                        continue
                    limit = 1.3 if 'H' in (get_element(a), get_element(b)) else 1.95
                    if np.linalg.norm(a.coord - b.coord) <= limit:
                        depth[b] = level
                        nxt.append(b)
            frontier = nxt
        return depth
    bonded_across_link = set()
    linked_residue_pairs = set()
    heavy_atoms = [a for a in atoms if a.parent.id[0] != 'W' and get_element(a) != 'H'
                   and get_element(a) not in METAL_COORD_DIST]
    if heavy_atoms:
        for la, lb in NeighborSearch(heavy_atoms).search_all(COVALENT_LINK_MAX):
            if la.parent is lb.parent:
                continue
            separation = sequence_separation(la.parent, lb.parent)
            if separation is not None and separation < 3:
                continue
            linked_residue_pairs.add(frozenset((id(la.parent), id(lb.parent))))
            depth_a = atoms_within_bonds(la.parent, la, 2)
            depth_b = atoms_within_bonds(lb.parent, lb, 2)
            for xa, da in depth_a.items():
                for xb, db in depth_b.items():
                    if da + db + 1 <= 3:
                        bonded_across_link.add(frozenset((id(xa), id(xb))))

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
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': 'SG', 'Atom2': 'SG',
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
                    # Atom2 reports the anion oxygen closest to the cation atom.
                    anion_atom = min(
                        (anion_res[a] for a in ANION_ATOMS[anion_res.resname] if a in anion_res),
                        key=lambda a: np.linalg.norm(a.coord - cp)
                    )
                    edges.append({
                        'NodeId1': n1, 'Interaction': 'IONIC:SC_SC', 'NodeId2': n2,
                        'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': cation_atom.name, 'Atom2': anion_atom.name,
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

        # PIHBOND (X-H...pi). With hydrogen positions available (strict angle on
        # and explicit hydrogens) the donor's hydrogen must point at the ring:
        # H...centroid <= PIHBOND_MAX_H_CENTROID and X-H...centroid >=
        # PIHBOND_MIN_ANGLE (Steiner & Koellner 2001). Without them only the
        # heavy-atom distance donor...centroid <= PIHBOND_MAX_DONOR_CENTROID can be
        # applied, which cannot tell whether the hydrogen points at the ring.
        pihbond_found = False
        aromatics = ['PHE', 'TYR', 'TRP', 'HIS']
        if r1.resname in aromatics or r2.resname in aromatics:
            for r_arom, r_donor in [(r1, r2), (r2, r1)]:
                if r_arom.resname not in aromatics:
                    continue
                c_arom = get_centroid(r_arom)
                if c_arom is None:
                    continue
                for a_donor in r_donor:
                    # a genuine H-bond donor (carries a polar hydrogen), not any N/O/S
                    if not hbond_role[id(a_donor)][0]:
                        continue
                    dist = np.linalg.norm(c_arom - a_donor.coord)
                    if dist > PIHBOND_MAX_DONOR_CENTROID:
                        continue
                    angle_str = 'nan'
                    if use_h_geometry:
                        best = None
                        for h in get_attached_hydrogens(a_donor):
                            d_h = np.linalg.norm(c_arom - h.coord)
                            ang = calc_angle(a_donor.coord, h.coord, c_arom)
                            if d_h <= PIHBOND_MAX_H_CENTROID and ang >= PIHBOND_MIN_ANGLE and (best is None or d_h < best[0]):
                                best = (d_h, ang)
                        if best is None:
                            continue
                        dist, ang = best
                        angle_str = f"{ang:.3f}"
                    edges.append({
                        'NodeId1': n1, 'Interaction': f'PIHBOND:SC_{"MC" if a_donor.name in MC_ATOMS else "SC"}', 'NodeId2': n2,
                        'Distance': f"{dist:.3f}", 'Angle': angle_str, 'Energy': 'nan', 'Atom1': '', 'Atom2': a_donor.name,
                        'Donor': format_node_id(r_donor), 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                    })
                    pihbond_found = True
                    break
                if pihbond_found:
                    break

        # METAL/HALOGEN/HBOND are each independent (at most one edge per category
        # and residue pair). VDW is reported only for residue pairs held together
        # by van der Waals contact alone, i.e. with no more specific interaction:
        # a residue pair that already has an H-bond is not re-counted as VDW under
        # a second label. This is a reporting convention (one edge per mechanism),
        # not a statement that H-bonded residues have no vdW contact.
        hbond_found = False
        best_vdw = None   # (surface gap, a1, a2, distance)
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

            # METAL coordination (per-metal distance cutoff, see METAL_COORD_DIST)
            if el1 in metals or el2 in metals:
                if not metal_found:
                    if el1 in metals and el2 in donors_acceptors:
                        metal_el = el1
                    elif el2 in metals and el1 in donors_acceptors:
                        metal_el = el2
                    else:
                        metal_el = None
                    if metal_el is not None and dist <= get_metal_coord_dist(metal_el):
                        edges.append({
                            'NodeId1': format_node_id(a1.parent), 'Interaction': f'METAL:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                            'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                            'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                        })
                        metal_found = True
                # coordination geometry is not a vdW contact nor a clash
                continue
            if frozenset((id(a1), id(a2))) in bonded_across_link:
                continue

            r_sum = get_vdw_radius(el1) + get_vdw_radius(el2)
            if dist > r_sum + vdw_max_gap and dist > 3.5:
                continue   # too far for any atom-pair criterion below

            # HALOGEN bond: X...A <= sum of vdW radii and C-X...A >= HALOGEN_MIN_ANGLE
            special_geometry = False
            if (el1 in halogens and el2 in donors_acceptors) or (el2 in halogens and el1 in donors_acceptors):
                x, acc = (a1, a2) if el1 in halogens else (a2, a1)
                if dist <= r_sum:
                    carbon = min((a for a in x.parent if a is not x and get_element(a) == 'C'
                                  and np.linalg.norm(a.coord - x.coord) <= 2.3),
                                 key=lambda a: np.linalg.norm(a.coord - x.coord), default=None)
                    if carbon is not None:
                        ang = calc_angle(carbon.coord, x.coord, acc.coord)
                        if ang >= HALOGEN_MIN_ANGLE:
                            special_geometry = True
                            if not halogen_found:
                                edges.append({
                                    'NodeId1': format_node_id(a1.parent), 'Interaction': f'HALOGEN:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                                    'Distance': f"{dist:.3f}", 'Angle': f"{ang:.3f}", 'Energy': 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                                    'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                                })
                                halogen_found = True

            # HBOND: a genuine donor and acceptor (by their role in the residue),
            # donor...acceptor <= 3.5 A and, with hydrogen positions, a donor-H...acceptor
            # angle > 120 deg (Baker & Hubbard 1984). Each direction is tried only
            # when the roles allow it, and the Donor column reports the direction
            # that passed.
            if dist <= 3.5:
                d1, acc1 = hbond_role[id(a1)]
                d2, acc2 = hbond_role[id(a2)]
                directions = []
                if d1 and acc2:
                    directions.append((a1, a2))
                if d2 and acc1:
                    directions.append((a2, a1))
                if directions:
                    special_geometry = True
                    chosen = None
                    if strict_angle:
                        for donor, acceptor in directions:
                            ok, ang = check_hbond_angle(donor, acceptor)
                            if ok and (chosen is None or ang > chosen[2]):
                                chosen = (donor, acceptor, ang)
                    else:
                        chosen = (directions[0][0], directions[0][1], 150.0)
                    if chosen is not None:
                        donor, _, angle_val = chosen
                        edges.append({
                            'NodeId1': format_node_id(a1.parent), 'Interaction': f'HBOND:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                            'Distance': f"{dist:.3f}", 'Angle': f"{angle_val:.3f}", 'Energy': 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
                            'Donor': format_node_id(donor.parent), 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                        })
                        hbond_found = True

            # A polar hydrogen and an acceptor overlap by design in an H-bond
            # (H...O ~1.9 A vs 2.72 A radii sum): that geometry belongs to HBOND,
            # not to vdW contact or clash assessment.
            if (polar_h_parent.get(id(a1)) is not None and hbond_role[id(a2)][1]) or \
               (polar_h_parent.get(id(a2)) is not None and hbond_role[id(a1)][1]):
                continue
            # Covalent links between residues (glycosylation, covalent ligands).
            if el1 != 'H' and el2 != 'H' and dist < COVALENT_LINK_MAX:
                continue

            gap = dist - r_sum   # separation of the two vdW surfaces
            if gap <= -CLASH_MIN_OVERLAP:
                continue   # steric clash: never a VDW edge (listed by the clash scan below)
            if gap <= vdw_max_gap and (best_vdw is None or gap < best_vdw[0]):
                best_vdw = (gap, a1, a2, dist)

        # VDW is the fallback contact type (see note above)
        if not hbond_found and not metal_found and not halogen_found and not pihbond_found and best_vdw:
            _, a1, a2, vdw_d = best_vdw
            edges.append({
                'NodeId1': format_node_id(a1.parent), 'Interaction': f'VDW:{get_contact_type(a1, a2)}', 'NodeId2': format_node_id(a2.parent),
                'Distance': f"{vdw_d:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': a1.name, 'Atom2': a2.name,
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
                                'Distance': f"{dist:.3f}", 'Angle': 'nan', 'Energy': 'nan', 'Atom1': pa1.name, 'Atom2': pa2.name,
                                'Donor': '', 'Positive': '', 'Cation': '', 'Orientation': '', 'Model': model_num
                            })
                            nodes.add(r1)
                            nodes.add(r2)
            except Exception as e:
                pass

    # Steric clashes, scanned independently of the network method so the same
    # structure always gives the same list: non-local residue pairs (>= 3 apart
    # along the chain) whose atoms overlap by >= CLASH_MIN_OVERLAP. Not clashes:
    # metal coordination, covalent links, and pairs that are legitimately closer
    # than their radii sum (donor/acceptor atoms within H-bond distance, a polar
    # hydrogen and an acceptor, a halogen and an acceptor).
    scan_atoms = [a for a in atoms if a.parent.id[0] != 'W']
    if scan_atoms:
        for c1, c2 in NeighborSearch(scan_atoms).search_all(3.4):
            if c1.parent is c2.parent:
                continue
            separation = sequence_separation(c1.parent, c2.parent)
            if separation is not None and separation < 3:
                continue
            e1, e2 = get_element(c1), get_element(c2)
            if e1 in METAL_COORD_DIST or e2 in METAL_COORD_DIST:
                continue
            d = np.linalg.norm(c1.coord - c2.coord)
            overlap = get_vdw_radius(e1) + get_vdw_radius(e2) - d
            if overlap < CLASH_MIN_OVERLAP:
                continue
            if e1 != 'H' and e2 != 'H' and d < COVALENT_LINK_MAX:
                continue
            (d1, acc1), (d2, acc2) = hbond_role[id(c1)], hbond_role[id(c2)]
            if d <= 3.5 and ((d1 and acc2) or (d2 and acc1)):
                continue
            if (polar_h_parent.get(id(c1)) is not None and acc2) or (polar_h_parent.get(id(c2)) is not None and acc1):
                continue
            if (e1 in ('F', 'CL', 'BR', 'I') and e2 in ('N', 'O', 'S')) or (e2 in ('F', 'CL', 'BR', 'I') and e1 in ('N', 'O', 'S')):
                continue
            if frozenset((id(c1), id(c2))) in bonded_across_link:
                continue
            if format_node_id(c1.parent) > format_node_id(c2.parent):
                c1, c2 = c2, c1
            clashes.append({
                'NodeId1': format_node_id(c1.parent), 'NodeId2': format_node_id(c2.parent),
                'Atom1': c1.name, 'Atom2': c2.name, 'Distance': f"{d:.3f}",
                'Overlap': f"{overlap:.3f}", 'Model': model_num
            })

    # Residue-residue interaction energy on every edge of a residue pair
    # (--energy). Water bridges are water-mediated, so the direct energy between
    # the two residues does not describe them and stays NaN.
    if calc_energy:
        residue_by_node = {format_node_id(r): r for r in structure.get_residues()}
        for e in edges:
            if e['Interaction'].startswith('WATERBRIDGE'):
                continue
            ra = residue_by_node.get(e['NodeId1'])
            rb = residue_by_node.get(e['NodeId2'])
            if ra is not None and rb is not None and ra is not rb:
                value = residue_pair_energy(ra, rb)
                e['Energy'] = 'nan' if np.isnan(value) else f"{value:.4f}"

    clashes_df = pd.DataFrame(clashes, columns=['NodeId1', 'NodeId2', 'Atom1', 'Atom2', 'Distance', 'Overlap', 'Model'])

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
        
    return edges_df, nodes, clashes_df

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
    parser = argparse.ArgumentParser(description="Sauron: residue interaction network (RIN) calculator. Output follows the RING edge/node file format.")
    parser.add_argument("input", help="Input PDB or CIF file")
    # On by default to stay in sync with SauronGUI.py, which always adds hydrogens:
    # --strict-angle (also on by default) rejects every HBOND without explicit hydrogens.
    parser.add_argument("--add-h", action=argparse.BooleanOptionalAction, default=True, help="Add hydrogens using pdb2pqr. On by default; use --no-add-h to disable (note that --strict-angle then keeps almost no HBOND edges on structures without hydrogens)")
    # On by default to stay in sync with SauronGUI.py, which always enforces both.
    parser.add_argument("--strict-angle", action=argparse.BooleanOptionalAction, default=True, help="Enforce strict angle constraints for Hydrogen Bonds (e.g. >120 deg). On by default; use --no-strict-angle to disable. Needs explicit hydrogens (native or via --add-h), otherwise no HBOND edges are kept")
    parser.add_argument("--remove-multiples", action=argparse.BooleanOptionalAction, default=True, help="Remove multiple interactions of the same type between the same residue pair. On by default; use --no-remove-multiples to disable")
    def calc_method_arg(value):
        value = value.strip().lower()
        if value == 'standard':
            return 'sauron'   # former name of the Sauron method, kept for existing scripts
        if value == 'voronoi':
            raise argparse.ArgumentTypeError(
                "the Voronoi method was removed: with Sauron's distance criteria its networks were a "
                "subset of the Sauron method's (residue-pair Jaccard >= 0.98); use 'sauron'")
        return value
    parser.add_argument("--calc-method", type=calc_method_arg, choices=["sauron", "insty", "rinpy"], default="sauron",
                        help="Method for RIN calculation: 'sauron' (Sauron's own geometric criteria, default), 'insty' (ProDy InSty) or 'rinpy' (elastic network model). 'standard' is accepted as the former name of 'sauron'.")
    parser.add_argument("--chains", type=str, help="Comma-separated list of chains to calculate (e.g. A,B,C)")
    parser.add_argument("--energy", action="store_true", help="Residue-residue interaction energy (kcal/mol) on every edge, summed over all atom pairs of the two residues: Coulomb with PQR partial charges and a distance-dependent dielectric (4r). Requires PQR charges (--add-h, the default, or a PQR input).")
    parser.add_argument("--vdw-energy", action="store_true", help="With --energy, add a generic 12-6 Lennard-Jones term (well depth 0.1 kcal/mol, r_min = sum of Bondi radii). If not set, only electrostatics is computed.")
    parser.add_argument("--vdw-contacts", choices=sorted(VDW_CONTACT_MAX_GAP), default=DEFAULT_VDW_CONTACTS,
                        help="Van der Waals contact classes accepted as VDW edges (Word et al. 1999 / MolProbity): 'close' (default) = vdW surfaces at most 0.25 A apart; 'wide' = up to 0.5 A apart. Overlaps of 0.4 A or more are steric clashes in both modes.")
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
            pqr_file = input_file if input_file.lower().endswith('.pqr') else temp_pdb

        # Partial charges for --energy come only from a real PQR file (pdb2pqr
        # output with --add-h, or a PQR given as input). Occupancy and B-factor of
        # an ordinary PDB/CIF are not charges and are never used as such.
        if args.energy:
            n_param = load_pqr_parameters(struct_to_calc, pqr_file) if pqr_file.lower().endswith('.pqr') else 0
            if n_param == 0:
                print("Warning: --energy needs partial charges (use --add-h or a PQR input); Energy will be NaN.")
            
        if args.chains:
            valid_chains = [c.strip() for c in args.chains.split(',')]
            for chain in list(struct_to_calc):
                if chain.id not in valid_chains:
                    struct_to_calc.detach_child(chain.id)
            
        edges_df, nodes, clashes_df = calculate_rin(struct_to_calc, strict_angle=args.strict_angle, remove_multiples=args.remove_multiples, model_num=model_num, calc_method=args.calc_method, pqr_file=pqr_file, calc_energy=args.energy, include_vdw=args.vdw_energy, water_bridges=args.water_bridges, vdw_contacts=args.vdw_contacts)
        nodes_df = build_nodes_df(nodes, edges_df, os.path.basename(args.input), structure=struct_to_calc, model_num=model_num)
        
        edges_out = f"{prefix}.edges"
        nodes_out = f"{prefix}.nodes"
        
        if not edges_df.empty:
            edges_df.to_csv(edges_out, sep='\t', index=False)
        if not nodes_df.empty:
            nodes_df.to_csv(nodes_out, sep='\t', index=False)
            
        # Steric clashes (vdW overlap >= 0.4 A) are structural defects, not
        # interactions: listed apart and never used as network edges.
        clashes_out = f"{prefix}_clashes.tsv"
        clashes_df.to_csv(clashes_out, sep='\t', index=False)

        print(f"Generated {edges_out} and {nodes_out}")
        print(f"Steric clashes (overlap >= {CLASH_MIN_OVERLAP} A) between non-local residues: {len(clashes_df)} (see {clashes_out})")
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
            if args.calc_method == 'insty':
                logf.write("VDW Contacts: not used (ProDy InSty applies its own criteria)\n")
            else:
                logf.write(f"VDW Contacts: {args.vdw_contacts} (vdW surface gap <= {VDW_CONTACT_MAX_GAP[args.vdw_contacts]} A)\n")
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
        # only the pdb2pqr output generated here; never a PQR given as input
        if args.add_h and pqr_file not in (temp_pdb, input_file) and os.path.exists(pqr_file): os.remove(pqr_file)

if __name__ == "__main__":
    main()
