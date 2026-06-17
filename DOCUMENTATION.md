# Sauron - Residue Interaction Network Calculator Documentation

## Table of Contents
1. [Overview](#overview)
2. [Web Interface (GUI)](#web-interface)
3. [Command-Line Interface (CLI)](#command-line-interface)
4. [Container Deployment](#container-deployment)
5. [Outputs](#outputs)
6. [Required Libraries](#required-libraries)
7. [Credits](#credits)

---

## Overview

**Sauron** is an advanced Residue Interaction Network (RIN) Calculator designed to parse structural files (PDB/CIF) and compute intra-molecular interactions such as Hydrogen Bonds, Salt Bridges, pi-pi stacking, and more. 

Sauron integrates a powerful backend written in Python with a beautiful, synchronized front-end web interface that allows real-time visualization of the 3D protein structure alongside its complex interaction network graph.

---

## Web Interface

The easiest way to interact with Sauron is via the web application, which provides both calculation capabilities and an interactive 3D dashboard.

### Starting the Web App

To launch the web app locally, run:
```bash
python SauronGUI.py
```
Open your browser and navigate to `http://localhost:5000`.

### Data Input

Sauron provides two main ways to input structural data:

1. **Upload File**: Drag and drop any `.pdb` or `.cif` file from your local machine.
2. **Fetch Structure**: Input an ID and Sauron will automatically fetch the structure for you. 
   - **RCSB PDB**: Enter a 4-letter PDB code (e.g., `1AFW`).
   - **AlphaFold DB**: Enter a UniProt accession code (e.g., `P00533`) to fetch the v6 structure from AlphaFold.

### Calculation Parameters

Before running the calculation, you can toggle several parameters:
- **Calculation Method**: 
  - *Standard*: Employs static distance bounds and neighbor-search algorithms. Perfect for strict cutoff network designs.
  - *Voronoi*: Calculates interaction networks via Voronoi Tessellation. This effectively overcomes distance thresholding by mapping 3D space directly between adjoining residues, ensuring dense and topological biological relevance.
  - *ProDy InSty*: Discovers atomic-level chemical interactions (such as Pi-Stacking, Cation-Pi, and Hydrogen Bonds) with dynamic subtyping (e.g. `MC_MC`, `SC_MC`). Automatically enforces `RING` compatibility.
- **Add Hydrogens**: Automatically runs `pdb2pqr` to add missing hydrogens to the structure before computing interactions. Explicitly modeling hydrogens is vital to maximizing the accuracy of ProDy's computations.
- **Strict Angle**: Enforces a strict geometric constraint (>120°) for explicit Hydrogen Bonds.
- **Remove Multiples**: Prunes redundant edges of the same type between two residues, preserving only the one with the shortest distance.
- **Chain Selection**: Choose to compute the network for all chains in the structure, or specify a comma-separated list of chains (e.g., `A,B`).

### Interactive Dashboard Features

Once calculation finishes, the results are visualized seamlessly:

1. **Synchronized 3D Viewers**: 
   - On the left, the structural viewer powered by **PDBe Molstar** shows the 3D protein model.
   - On the right, a 3D force-directed graph (powered by **3d-force-graph**) displays the topological interactions. Selecting a node in one viewer instantly highlights it in the other.
2. **Node Coloring**: Color nodes in the graph based on their default setup, DSSP secondary structure assignment, chemical features, or Chain ID.
3. **Primary Sequence**: A linear, 1-letter sequence viewer dynamically colored by DSSP. Hover or click on residues to select them in the 3D space.
4. **Data Tables**: 
   - **Node Interactions**: Lists every bond and distance connected to the currently selected node.
   - **Top 10 Residues**: A leaderboard of the most central residues (configurable by Degree, Betweenness, Clustering, or Eigenvector Centrality).
   - **Inter-chain Interactions**: A dedicated table showing connections occurring between distinct chains.
   - **Ligand Interactions**: Detects non-standard amino acids (ligands) and maps their interactions.
5. **Network Complexity Analysis**: A dedicated panel mapping the degree distribution of the network (with Log-Log toggle) and providing overarching metrics like Network Density and Average Shortest Path.

---

## Command-Line Interface

For batch processing or server-side automation, use `sauron.py`.

```bash
python sauron.py <input_file> [options]
```

### Options:
- `--calc-method`: Choose the interaction calculation methodology (`standard`, `voronoi`, or `insty`). Default is `standard`.
- `--add-h`: Add explicit hydrogens using `pdb2pqr` prior to computations (Highly recommended when using `--calc-method insty`).
- `--strict-angle`: Enforce strict >120° angle constraints for Hydrogen Bonds.
- `--remove-multiples`: Remove multiple interactions of the same type between the same residue pair.
- `--chains`: Comma-separated list of chains to calculate (e.g., A,B,C).

**Example:**
```bash
python sauron.py 1AFW.pdb --chains A,B --remove-multiples
```

---

## Container Deployment

Sauron can be deployed seamlessly using containers, eliminating the need to configure Python dependencies manually.

### Docker / Podman

A `Dockerfile` is included in the project root.

**Build:**
```bash
docker build -t sauron-gui .
```

**Run:**
```bash
docker run -p 5000:5000 -v $(pwd)/uploads:/opt/sauron/uploads sauron-gui
```
*Note: The `-v` flag mounts a local `uploads` directory so you can access the downloaded `.zip` output files easily.*

### Apptainer / Singularity

An `Apptainer.def` is provided for HPC environments.

**Build:**
```bash
apptainer build sauron.sif Apptainer.def
```

**Run:**
```bash
# Set a writable directory for the uploads (since Apptainer is read-only by default)
export SAURON_UPLOAD_DIR=/tmp/sauron_uploads
mkdir -p $SAURON_UPLOAD_DIR

apptainer run --bind $SAURON_UPLOAD_DIR:/opt/sauron/uploads sauron.sif
```

---

## Outputs

All calculations generate a comprehensive `.zip` package containing:

1. **`*.edges`**: A TSV file defining every network edge (Source, Target, Distance, Angle, Interaction Type).
2. **`*.nodes`**: A TSV file defining every node (Chain, Position, Residue, DSSP, 3D Coordinates).
3. **`*_network_metrics.tsv`**: Full topological metrics for every node in the graph.
4. **`*_top25_metrics.tsv`**: A summarized leaderboard of the top 25 central nodes.

*(Note: Multi-model files like NMR structures will generate separate files for each individual model).*

---

## Required Libraries

Sauron is built on the shoulders of the following exceptional open-source libraries:

**Python Backend:**
- **Flask**: Web application framework.
- **Biopython**: PDB and MMCIF parsing and geometry algorithms.
- **NetworkX**: Complex network generation and centrality metrics.
- **Pandas** & **NumPy**: Data structuring and linear algebra.
- **Requests**: Fetching structures from RCSB and AlphaFold APIs.
- **pdb2pqr**: Explicit hydrogen addition capabilities.
- **pydssp**: Neural-network based secondary structure assignment.

**JavaScript Frontend:**
- **PDBe Molstar**: 3D structural rendering.
- **3d-force-graph**: WebGL-accelerated 3D force-directed network graphs.
- **Plotly.js**: Dynamic and interactive charting for complexity metrics.

---

## Artificial Intelligence Disclosure
The authors used generative AI technologies for code auditing and performance optimization (Gemini® and GitHub-CoPilot®). Throughout this process, the authors maintained full control over the tool's design and interpretation of its results; the AI acted solely as a technical and linguistic aid. All AI-generated content was reviewed and approved by the authors, who take full responsibility for the final tool.

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Credits

Developed by the [EvoMol-lab's Team](https://evomol-lab.imd.ufrn.br/) at [BioMe](https://bioinfo.imd.ufrn.br/), [UFRN](https://www.ufrn.br/).

<div align="center">
  <img src="EvoMol_v2.png" alt="EvoMol Logo" width="400">
</div>
