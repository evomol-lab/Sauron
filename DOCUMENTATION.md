# Sauron — Residue Interaction Network Calculator: Full Documentation

## Table of Contents
1. [Overview](#overview)
2. [Installation](#installation)
3. [Command-Line Interface (`sauron.py`)](#command-line-interface-sauronpy)
4. [Interaction Types and Geometric Definitions](#interaction-types-and-geometric-definitions)
5. [Output Files](#output-files)
6. [Web Interface (`SauronGUI.py`)](#web-interface-sauronguipy)
7. [Validation Against RING](#validation-against-ring)
8. [Known Limitations](#known-limitations)
9. [Container Deployment](#container-deployment)
10. [Third-Party Libraries](#third-party-libraries)
11. [License](#license)
12. [Artificial Intelligence Disclosure](#artificial-intelligence-disclosure)
13. [Citing Sauron](#citing-sauron)
14. [Credits](#credits)

---

## Overview

**Sauron** parses a protein structure file (PDB or PDBx/mmCIF, single-model or
multi-model/NMR ensembles) and computes a **Residue Interaction Network (RIN)**: a
graph where nodes are residues and edges are typed, geometrically-defined
non-covalent (and covalent, for disulfides) contacts between them — hydrogen bonds,
van der Waals contacts, salt bridges, π-stacking, cation-π, metal coordination,
halogen bonds, π-hydrogen bonds, and (optionally) water-mediated bridges.

Its outputs are designed to be compatible with the [RING](https://ring.biocomputingup.it/) tool (Del Conte et al., *NAR*
2024) — see [Comparisons Against other tools](#comparisons-against-other-tools). Sauron is not a
reimplementation of RING's (closed-source) engine, and it is not supposed to produce the same results as RING. It independently computes the
same interaction *types* from documented geometric criteria and published
literature, using available Python packages (e.g. ProDy, RinPy, and NetworkX, etc.).

Sauron has two front ends built on the same [`sauron.py`](sauron.py) calculation
engine:
- **`sauron.py`** — a command-line tool for batch/scripted use.
- **`SauronGUI.py`** — a local Flask web application (`templates/index.html`,
  `static/script.js`) that wraps the same engine with a 3D structure viewer, an
  interactive network graph, and summary tables/charts.

---

## Installation

Requires **Python ≥ 3.12**.

```bash
pip install -r requirements.txt
# or, with uv (the project also ships a pyproject.toml):
uv sync
```

`pdb2pqr` (used by `--add-h` and `--energy`) is installed as a Python dependency and
invoked as a subprocess; no separate system install is needed. See
[Container Deployment](#container-deployment) for a dependency-free alternative.

---

## Command-Line Interface (`sauron.py`)

```bash
python sauron.py <input_file> [options]
```

`<input_file>` is a `.pdb` or `.cif`/`.mmcif` file. Multi-model files (e.g. NMR
ensembles) are processed model by model; see [Output Files](#output-files).

### Options

| Flag | Description |
|---|---|
| `--calc-method {standard,voronoi,insty,rinpy}` | Interaction-calculation methodology. Default: `standard`. See below. |
| `--add-h` | Add explicit hydrogens with `pdb2pqr` (force field `PARSE`) before computing interactions. No-op (with a fallback warning) if `pdb2pqr` is unavailable. Strongly recommended for `--calc-method insty`, and required for `--strict-angle` / `--energy` to have any effect on a structure that doesn't already carry hydrogens (e.g. most X-ray depositions; NMR ensembles usually already have them). |
| `--strict-angle` | Enforce a real donor–H···acceptor angle (> 120°) for `HBOND` edges, using explicit hydrogens already present in the structure (native or added via `--add-h`). Without this flag, `HBOND` edges use a fixed placeholder `Angle` of `150.000` and are accepted on distance alone. If no hydrogens are present at all, every candidate is rejected (no false positives, but also no `HBOND` edges) — use `--add-h` together with this flag on H-less inputs. |
| `--remove-multiples` | Collapse duplicate edges of the *same type* between the same residue pair down to the single shortest-distance one. Multiple edges of *different* types on the same pair (e.g. an `HBOND` and a separate `VDW` contact) are unaffected. |
| `--chains A,B,C` | Restrict the calculation to the given chain(s) (comma-separated, case-sensitive to match the file). Default: all chains. |
| `--energy` | Compute a per-edge electrostatic interaction energy (`Energy` column, kcal/mol, Coulomb's law with a distance-dependent dielectric: `332.0637·q1·q2/(4·r²)`) using atomic partial charges/radii from a PQR file. Requires either `--add-h` (which produces a PQR as a side effect) or an already-charge-annotated input. Without `--energy`, `Energy` is always `nan`. |
| `--vdw-energy` | With `--energy`, also add a generic Lennard-Jones term (using the PQR radii) to `Energy`. Without it, `--energy` reports electrostatics only. No effect without `--energy`. |
| `--water-bridges` | Detect `WATERBRIDGE` edges (see below). Off by default; not computed for `--calc-method insty`. A Sauron-specific extension. |

### `--calc-method` choices

- **`standard`** (default) — Sauron's own distance/angle/role-based geometric
  engine; see [Interaction Types and Geometric Definitions](#interaction-types-and-geometric-definitions).
  This is the method validated against RING.
- **`voronoi`** — Same interaction typing as `standard`, but the candidate
  atom-contact list comes from a 3D Voronoi tessellation of all non-water atoms
  (via `scipy.spatial.Voronoi`) instead of a fixed-radius neighbor search, capped
  at 6.5 Å. This trades a hard distance cutoff for a topological "who's actually
  adjacent in space" contact definition; falls back to `standard` if tessellation
  fails (e.g. too few points).
- **`insty`** — Delegates interaction detection to
  [ProDy](http://www.bahargroup.org/prody/)'s `calcProteinInteractions` /
  `calcWaterBridges`, which independently implements a similar interaction
  taxonomy (HBOND, IONIC/salt bridges, PIPISTACK, PICATION, VDW/hydrophobic,
  SSBOND, WATERBRIDGE). Requires a structure with explicit hydrogens (use
  `--add-h`) and reads coordinates from a PQR file when available. Falls back to
  `standard` if ProDy fails to parse the input.
- **`rinpy`** — Runs the [RinPy](https://pypi.org/project/rinpy/) package instead
  of Sauron's own typed-edge engine. RinPy builds a **residue-affinity network**
  from an Elastic Network Model (not RING-style typed contacts) and computes
  degree/betweenness/closeness centrality plus hinge-residue detection via
  Laplacian normal-mode analysis. Output is written under
  `<prefix>_rinpy_out/`, and a `<prefix>_rinpy_summary.json` (top-10 residues per
  centrality measure, hinge residues per mode) is produced for the GUI. This mode
  does **not** produce `.edges`/`.nodes`/network-metrics files.

### Examples

```bash
# Default RIN, strict H-bond geometry, one edge per type per residue pair
python sauron.py 1AFW.pdb --strict-angle --remove-multiples

# Add hydrogens first, compute electrostatic + VDW energy per edge, chains A and B only
python sauron.py 1AFW.pdb --add-h --energy --vdw-energy --chains A,B

# Include Sauron's water-bridge extension
python sauron.py 4AKE.cif --water-bridges

# Elastic Network Model / hinge analysis via RinPy
python sauron.py 4AKE.cif --calc-method rinpy
```

---

## Interaction Types and Geometric Definitions

All distance thresholds below apply to the **`standard`**/`voronoi` engine
(`calculate_rin()` in [sauron.py](sauron.py)). Candidate atom pairs are first
gathered with a 6.5 Å neighbor search (the widest of all thresholds below), then
grouped by residue pair. **Residue pairs on the same chain closer than 3 apart in
sequence are excluded outright** (`|i − j| < 3`): the fixed peptide-bond geometry
alone puts adjacent residues' backbone atoms well inside every threshold below
(e.g. `O(i)···N(i+1)` sits at ≈2.2 Å), which would otherwise be misreported as a
genuine interaction on *every* residue in every chain.

Each type below is evaluated **independently** — a residue pair can (and often
does) carry edges of several different types simultaneously (e.g. a backbone
`HBOND` and a side-chain `PICATION` on the same pair). Per residue pair, most
types report at most one edge (`SSBOND`, `IONIC`, `PIPISTACK`, `PICATION`,
`METAL`, `HALOGEN`, `PIHBOND`, and the `VDW` fallback described below all keep
only their single best/first-qualifying atom pair); `HBOND` and `WATERBRIDGE` are
the exceptions and can report one edge per qualifying atom pair (e.g. two
side-chain acceptors each hydrogen-bonding to the same donor produce two `HBOND`
rows) — use `--remove-multiples` to collapse those down to one per type per pair.

| Type | Criterion | Notes |
|---|---|---|
| **SSBOND** | `SG···SG` distance ≤ 2.5 Å between two `CYS` residues. | Supersedes all other checks for that pair (a disulfide-bonded pair gets nothing else). |
| **IONIC** | Distance ≤ 4.0 Å between a cation's charge-bearing point and an anion's charge centroid. Cation point: `ARG:CZ` or `LYS:NZ` (the single atom RING itself reports). Anion centroid: mean of `ASP:{OD1,OD2}` or `GLU:{OE1,OE2}`. | `HIS` is **not** treated as a cation (its protonation state is ambiguous without explicit hydrogens; RING's own reference output does not use it as one either). |
| **PIPISTACK** | Ring-centroid-to-ring-centroid distance ≤ 6.5 Å between two aromatic residues (`PHE`, `TYR`, `TRP`, `HIS`). | `Orientation` is classified as `P` (parallel/stacked) or `T` (perpendicular/edge-to-face) from the angle between the two rings' plane normals (computed via SVD of the ring atom coordinates): ≤ 45° → `P`, otherwise `T`. |
| **PICATION** | Distance ≤ 5.0 Å between an aromatic ring centroid and a cation point (same `CZ`/`NZ` convention as `IONIC`), **and** the angle between the cation vector and the ring's plane normal is ≤ 45° (0° = cation directly above the ring face). | The angle check rejects cations sitting off to the side of the ring rather than over its face — without it, this type over-counts by roughly 3× relative to RING. |
| **METAL** | Distance ≤ a per-element coordination cutoff between a recognized metal atom and an `N`/`O`/`S` atom. Cutoffs (Å): Na 3.0, K 3.4, Mg 2.6, Ca 3.2, Mn 2.7, Fe 2.7, Co 2.6, Ni 2.5, Cu 2.6, Zn 2.7, Cd 2.9, Hg 3.0, Li 2.6, Cs 3.6; unlisted metals default to 3.0. | Element identity is read from the file's own element field when present, falling back to a two-letter-aware guess from the atom name otherwise (disambiguating e.g. a calcium ion `CA` from an α-carbon `CA`). |
| **HALOGEN** | Distance ≤ (sum of the two atoms' van der Waals radii) + 0.3 Å between a halogen (`F`,`Cl`,`Br`,`I`) and an `N`/`O`/`S` atom. | Scales per element pair rather than a flat cutoff, so it doesn't under-count for the larger halogens (Br, I). |
| **PIHBOND** | Distance ≤ 4.0 Å between an aromatic ring centroid and a genuine hydrogen-bond donor heavy atom (per the same donor/acceptor table `HBOND` uses — not just any N/O/S). | Heavy-atom-to-centroid distance (no explicit hydrogen assumed); this is intentionally conservative, since without real hydrogen geometry this type is prone to false positives. |
| **HBOND** | Distance ≤ 3.5 Å between a genuine donor and a genuine acceptor atom, matched by their actual role in their residue (see table below) — not simply "any two polar atoms". With `--strict-angle` and explicit hydrogens available, also requires donor–H···acceptor angle > 120°; `Angle` is otherwise a fixed `150.000` placeholder. | Donor/acceptor roles for the 20 standard amino acids: backbone `N` is a donor (except `PRO`, which has none) and backbone `O`/`OXT` is an acceptor; side chains follow standard chemistry (e.g. `SER`/`THR`/`TYR` hydroxyls and `HIS` ring nitrogens can act as either; `LYS:NZ`, `ARG:{NE,NH1,NH2}`, `ASN:ND2`, `GLN:NE2`, `TRP:NE1` are donor-only; `ASP`/`GLU` carboxylate oxygens and `MET:SD` are acceptor-only). Non-standard residues (ligands, nucleic acids) fall back to "any N/O/S can be either", since no chemical-component dictionary is used to type them precisely. |
| **VDW** | Distance ≤ (sum of the two atoms' van der Waals radii) + a tolerance, using each atom's real element (not guessed from its PDB name). Tolerance is **0.5 Å** if the structure has no explicit hydrogens, or **0.0 Å** (the literal touching distance) if it does. | Reported **only as a fallback**: a residue pair only gets a `VDW` edge if it has no `HBOND`/`METAL`/`HALOGEN`/`PIHBOND` edge at all; when it does report one, it's the single closest qualifying atom pair. The tolerance split exists because the +0.5 Å slack was originally meant to compensate for a heavy-atom-only structure's *missing* hydrogens — applying it on top of already-explicit hydrogens (NMR structures, or after `--add-h`) double-counts that slack and was the largest single source of over-counting relative to RING (see [Validation Against RING](#validation-against-ring)). |
| **WATERBRIDGE** (opt-in, `--water-bridges`) | Two residues both within 3.5 Å of the same crystallographic water oxygen, where both contact atoms are genuine H-bond donor/acceptor atoms. Each water considers only its 4 closest such atoms before pairing them up (avoids a combinatorial blow-up in crowded solvent pockets). `Distance` is the **total water-mediated path length** (atom1↔water + water↔atom2), *not* the direct atom1↔atom2 distance (which can otherwise read as a nonsensical 5–7 Å "contact"). | Not computed for `--calc-method insty` (which has its own, ProDy-driven water-bridge detection) or `rinpy` (no typed edges at all). |

`Interaction` values are further suffixed with a contact-type tag describing which
part of each residue is involved: `MC` (main chain), `SC` (side chain), or `LIG`
(a HETATM/ligand residue) — e.g. `HBOND:MC_SC`, `VDW:SC_SC`.

---

## Output Files

For an input `structure.pdb` (single-model) or `structure_model_<N>.pdb`
(multi-model, 1-indexed to match the file's own model numbering), Sauron writes:

### `*.edges` (tab-separated)

| Column | Meaning |
|---|---|
| `NodeId1`, `NodeId2` | `Chain:Position:InsertionCode:ResidueName` (insertion code is `_` when absent), alphabetically ordered. |
| `Interaction` | `TYPE:CONTACT` — see [Interaction Types](#interaction-types-and-geometric-definitions). |
| `Distance` | Å, meaning depends on the type (see the table above — most are direct atom or centroid distances; `WATERBRIDGE`'s is a two-leg path length). |
| `Angle` | Degrees, where applicable (`HBOND` with `--strict-angle`, `PICATION`); `nan` otherwise. |
| `Energy` | kcal/mol, only populated with `--energy`; `nan` otherwise. |
| `Atom1`, `Atom2` | Specific atom names involved (blank for centroid-based types: `PIPISTACK`, `PICATION`). |
| `Donor` | For `HBOND`/`PIHBOND`, the node id of the donor residue. |
| `Positive`, `Cation` | Reserved columns (kept for RING format compatibility; currently unpopulated). |
| `Orientation` | `P`/`T` for `PIPISTACK` only. |
| `Model` | 1-indexed model number. |

### `*.nodes` (tab-separated)

`NodeId`, `Chain`, `Position`, `Residue`, `Type` (`RES` or `LIG`), `Dssp`
(secondary-structure code from [pydssp](https://pypi.org/project/pydssp/), if it
could be computed), `Degree` (in the produced edge graph), `Bfactor_CA`, `x`/`y`/`z`
(Cα, or centroid if no Cα), `pdbFileName`, `Model`.

### Network metrics

- **`*_network_metrics.tsv`** — Degree, Clustering Coefficient, Betweenness
  Centrality, and Eigenvector Centrality for every node (via NetworkX).
- **`*_top25_metrics.tsv`** — the top 25 nodes for each of the four metrics above,
  side by side.
- **`*_degree.pdb`** — the original structure with each atom's B-factor replaced
  by its residue's degree, for quick visualization of hub residues in any
  structure viewer.
- **`*_sauron.log`** — a plain-text run log (timestamp, input file, all flags
  used, node/edge counts).

### Multi-model input

Each model gets its own full set of the files above, named
`<prefix>_model_<N>.*` (`N` starting at 1, matching the file's own model
numbering). A plain `<prefix>.*` is only used for single-model input.

### `--calc-method rinpy`

Produces `<prefix>_rinpy_out/` (RinPy's own CSV/PNG/HTML/PyMOL-script outputs) and
`<prefix>_rinpy_summary.json` (top-10 residues by degree/betweenness/closeness,
and hinge residues per Laplacian mode) instead of `.edges`/`.nodes`/metrics files.

---

## Web Interface (`SauronGUI.py`)

```bash
python SauronGUI.py
```

Opens a Flask server at `http://localhost:5000` (auto-opens a browser window when
run as a packaged executable). Uploaded/fetched files and all generated output are
stored under `uploads/` by default (override with the `SAURON_UPLOAD_DIR`
environment variable).

### Input

- **Upload**: drag-and-drop or browse for a local `.pdb`/`.cif` file.
- **Fetch from database**: enter an RCSB PDB ID (fetched as mmCIF) or a UniProt
  accession (fetched from the AlphaFold DB, model v6).

### Calculation parameters

All of [`sauron.py`](#command-line-interface-sauronpy)'s options are exposed as
toggles: calculation method, chain selection (all vs. a specific comma-separated
list), Protein Energy Network (`--energy`, with a nested van der Waals sub-toggle
for `--vdw-energy`), and Water Bridges (`--water-bridges`, hidden automatically for
methods that don't support it). `--strict-angle` and `--remove-multiples` are
always enforced by the GUI backend (not exposed as toggles).

### Interactive dashboard

Once a calculation finishes:

- **3D structure viewer** ([PDBe Molstar](https://github.com/molstar/pdbe-molstar)):
  cartoon representation of the structure, background synced to the page's
  light/dark theme.
- **Interaction graph** ([3d-force-graph](https://github.com/vasturiano/3d-force-graph)):
  a WebGL force-directed 3D layout of the same network, with a type-colored edge
  legend (`HBOND`, `VDW`, `IONIC`, `PIPISTACK`, `PICATION`, `SSBOND`, `WATERBRIDGE`,
  `METAL`, `HALOGEN`, `PIHBOND`) and four node-coloring modes: default, DSSP
  secondary structure, chemical class (non-polar/polar/basic/acidic), or chain.
- **Bidirectional residue selection**: clicking a residue in *any* of the 3D
  structure viewer, the interaction graph, the primary-sequence bar, or any of the
  data tables below highlights that same residue — in a distinct highlight color,
  independent of whichever node-coloring mode is active — everywhere else at once
  (camera focus in the graph, selection/focus in the 3D viewer, outline + scroll in
  the sequence bar), and refreshes the node interaction table.
- **Primary sequence bar**: one-letter-code sequence per chain, colored by DSSP.
- **Data tables**: node interaction details for the current selection; top-10
  residues by Degree/Betweenness/Clustering/Eigenvector centrality; inter-chain
  interactions; ligand interactions (residues Sauron doesn't recognize as one of
  the 20 standard amino acids); and, for `--calc-method rinpy`, dedicated tabs for
  its centrality rankings and hinge-residue modes.
- **Network complexity analysis**: degree-distribution histogram/scatter (with a
  log-log toggle) plus overall density and average-shortest-path metrics.
- **Download**: a `.zip` bundling every output file described in
  [Output Files](#output-files).

Multi-model input shows a one-time notice that only the first model is
visualized in the dashboard (all models are still fully computed and included in
the downloaded `.zip`); for comparing RINs across an ensemble's conformations, use our other tool
[SlyTheRINs](https://github.com/evomol-lab/slytherins) ([Bradaschia et al. 2026](https://www.biorxiv.org/content/10.64898/2026.02.26.708270v1.abstract)).

---

## Comparisons Against other tools

Sauron's `standard` engine was checked edge-by-edge against real RING output for
two structures with very different characteristics: **4AKE** (X-ray, no
hydrogens) and **2N8N** (a 20-model NMR ensemble, hydrogens present natively).
Current state (aggregate edge-count ratio, Sauron ÷ RING):

| Type | 4AKE (X-ray, no H) | 2N8N (NMR, 20 models, w/ H) |
|---|---|---|
| HBOND | 1.53× | 1.10× |
| IONIC | **1.00×** | 1.02× |
| PICATION | **1.00×** | 1.00× |
| PIPISTACK | **1.00×** | — (none in this structure) |
| VDW | 1.06× | 1.12× |
| **Total** | **1.31×** | **1.11×** |

`IONIC`/`PICATION`/`PIPISTACK` match RING's reference almost exactly. The remaining `HBOND`
gap on hydrogen-less input is an inherent limitation of geometric hydrogen-bond
detection without real donor–H···acceptor geometry (`--add-h --strict-angle`
narrows it substantially, as the 2N8N column shows).

Comparisons with other tools will be added in the future.

---

## Known Limitations

- **No hydrogens ⇒ approximate `HBOND`/`PIHBOND`/`WATERBRIDGE` geometry.** These
  types rely on real donor–H···acceptor angles for full accuracy; without
  `--add-h`/native hydrogens they fall back to distance- and role-only heuristics,
  which typically over-count relative to RING (see above).
- **No chemical-component dictionary.** Ligands, modified residues, and nucleic
  acids are not individually typed for `HBOND` donor/acceptor roles — they use a
  permissive "any N/O/S can be either" fallback, unlike RING 4.0's PDB-CCD-driven
  typing of >35,000 known compounds.
- **`HIS` is never treated as a cation** for `IONIC`/`PICATION` (protonation state
  is ambiguous without explicit hydrogens).
- **Sauron's thresholds are reverse-engineered from published literature** and empirically compared to RING outputs. Sauron's is still under initial development, so **it is not supposed to produce the same results as RING**, which is a mature and widely-tested tool.
- **`WATERBRIDGE` is a Sauron-only heuristic** treat it as indicative rather than authoritative.

---

## Container Deployment

### Docker / Podman

```bash
docker build -t sauron-gui .
docker run -p 5000:5000 -v $(pwd)/uploads:/opt/sauron/uploads sauron-gui
```
*(The `-v` flag mounts a local `uploads` directory so you can access the
downloaded `.zip` output files easily.)*

### Apptainer / Singularity

```bash
apptainer build sauron.sif Apptainer.def

export SAURON_UPLOAD_DIR=/tmp/sauron_uploads
mkdir -p $SAURON_UPLOAD_DIR
apptainer run --bind $SAURON_UPLOAD_DIR:/opt/sauron/uploads sauron.sif
```

---

## Third-Party Libraries

### Python backend

| Library | Used for | License |
|---|---|---|
| [Biopython](https://biopython.org/) | PDB/mmCIF parsing, neighbor search, structure geometry | Biopython License (BSD-3-style) |
| [Flask](https://flask.palletsprojects.com/) | Web application / API for the GUI | BSD-3-Clause |
| [NetworkX](https://networkx.org/) | Graph construction and centrality metrics | BSD-3-Clause |
| [NumPy](https://numpy.org/) | Numerical/vector geometry | BSD-3-Clause |
| [pandas](https://pandas.pydata.org/) | Tabular edge/node output | BSD-3-Clause |
| [SciPy](https://scipy.org/) | Voronoi tessellation (`--calc-method voronoi`) | BSD-3-Clause |
| [Requests](https://requests.readthedocs.io/) | Fetching structures from RCSB/AlphaFold | Apache-2.0 |
| [pdb2pqr](https://pdb2pqr.readthedocs.io/) | `--add-h` / PQR generation for `--energy` | BSD-3-Clause |
| [PROPKA](https://github.com/jensengroup/propka) | pKa/protonation support used by pdb2pqr | LGPL-2.1 (used as an unmodified runtime dependency; not statically linked or redistributed) |
| [pydssp](https://pypi.org/project/pydssp/) | Secondary-structure (DSSP-style) assignment | MIT |
| [ProDy](http://www.bahargroup.org/prody/) | `--calc-method insty` interaction detection | MIT |
| [RinPy](https://pypi.org/project/rinpy/) | `--calc-method rinpy` (Elastic Network Model, hinge detection) | MIT |
| [mmcif-pdbx](https://pypi.org/project/mmcif-pdbx/) | mmCIF I/O support | Apache-2.0 |
| [SymPy](https://www.sympy.org/) | Symbolic math (transitive, ProDy/others) | BSD-3-Clause |
| [PyTorch](https://pytorch.org/) | Transitive dependency (RinPy/ProDy tooling) | BSD-3-Clause |

Werkzeug, Jinja2, MarkupSafe, click, itsdangerous, and blinker (Flask's own stack)
are all BSD-3-Clause under the Pallets Projects.

### JavaScript frontend (loaded from CDN, not vendored)

| Library | Used for | License |
|---|---|---|
| [PDBe Molstar](https://github.com/molstar/pdbe-molstar) | 3D structure viewer | Apache-2.0 |
| [3d-force-graph](https://github.com/vasturiano/3d-force-graph) | 3D interaction graph rendering | MIT |
| [Plotly.js](https://plotly.com/javascript/) | Complexity-analysis charts | MIT |

**License compatibility note**: none of Sauron's direct dependencies are
copyleft in a way that would force Sauron's own source to adopt their license —
PROPKA's LGPL-2.1 only governs its own source and is used here strictly as an
unmodified runtime dependency (imported by pdb2pqr, not redistributed or linked
into Sauron itself). Everything else above is permissive (BSD/MIT/Apache-2.0).
This means Sauron's own code is free to be released under any license the authors
choose, independent of these terms.

---

## License

This project is licensed under the **PolyForm Noncommercial License 1.0.0** — see
the [LICENSE](LICENSE) file for the full, canonical text. In short: you may use,
copy, modify, and distribute Sauron (including modified versions) for any purpose
**other than** commercial use, provided attribution and the license text are
preserved. See <https://polyformproject.org/licenses/noncommercial/1.0.0/> for the
license's own plain-language summary. Commercial use requires a separate
agreement with the authors ([EvoMol-lab](https://evomol-lab.imd.ufrn.br/)).

Sauron was previously released under the MIT License; the switch to PolyForm
Noncommercial 1.0.0 applies going forward and does not retroactively change the
terms under which earlier releases were already distributed.

---

## Artificial Intelligence Disclosure
The authors used generative AI technologies for code auditing and performance
optimization. Throughout this process, the authors
maintained full control over the tool's design and interpretation of its results;
the AI acted solely as a technical and linguistic aid. All AI-generated content was
reviewed and approved by the authors, who take full responsibility for the final
tool.

---

## Citing Sauron

If you use Sauron in your research, please cite the Sauron GitHub repository. A
preprint is forthcoming.

## Credits

Developed by the [EvoMol-lab Team](https://evomol-lab.imd.ufrn.br/) at
[BioMe](https://bioinfo.imd.ufrn.br/), [UFRN](https://www.ufrn.br/).

<div align="center">
  <img src="EvoMol_v2.png" alt="EvoMol Logo" width="400">
</div>
