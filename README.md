<div align="center">
  <img src="static/sauron-rins-logo.png" alt="Sauron Logo" width="300">
</div>

# Sauron

Sauron is a Residue Interaction Network (RIN) Calculator that identifies and maps non-covalent intra-molecular interactions in protein structures. It analyzes PDB and MMCIF files to calculate interactions such as hydrogen bonds, salt bridges, disulfide bonds, and van der Waals contacts, outputting node and edge files perfectly suited for network analysis. 

Sauron is designed to align with and improve upon established standards (like RING 4.0), providing strict geometry controls and graph-based topological metrics.

<div align="center">
  <img src="EvoMol-logo.png" alt="Rictusempra Logo" width="100">
</div>

Developed by the [EvoMol-Lab](github.com/evomol-lab), [BioME](bioinfo.imd.ufrn.br), UFRN, Brazil.

## Features

- **Multiple Interaction Types**: Calculates SSBOND, IONIC (salt bridge), PIPISTACK, PICATION, HBOND, and VDW interactions.
- **Topological Metrics**: Automatically computes network properties including Degree, Clustering Coefficient, Betweenness Centrality, and Eigenvector Centrality for every node.
- **Rigorous Filters**: 
  - Optionally enforce strict >120° angle constraints for Hydrogen Bonds.
  - Remove multiple/redundant interactions of the same type between the same residue pair to prevent inflated network degrees.
- **Hydrogen Addition**: Built-in support to add missing hydrogens using `pdb2pqr`.
- **Dual Interface**: Use Sauron either via a powerful Command-Line Interface (CLI) or a modern, user-friendly Web UI.

## Requirements

The core dependencies are listed in `requirements.txt`. Key packages include:
- `biopython`
- `networkx`
- `pandas`
- `numpy`
- `Flask` (for the web interface)
- `pdb2pqr` (optional, for explicit hydrogen addition)
- `pydssp` (optional, for secondary structure assignment)

To install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

### 1. Command-Line Interface (CLI)

You can run Sauron directly on any PDB or CIF file.

```bash
python sauron.py <input_file> [options]
```

**Options:**
- `--add-h`: Add hydrogens using `pdb2pqr` before calculating interactions.
- `--strict-angle`: Enforce strict angle constraints for Hydrogen Bonds (e.g. >120°).
- `--remove-multiples`: Remove multiple interactions of the same type between the same residue pair.

**Example:**
```bash
python sauron.py 1AFW.pdb --add-h --strict-angle --remove-multiples
```

### 2. Web Interface

Sauron includes a beautiful, interactive web interface for uploading files and downloading ZIP packages of the results.

Start the Flask server:
```bash
python SauronGUI.py
```
Then, open your browser and navigate to `http://localhost:5000`. 
Drag and drop your PDB or CIF file, toggle your desired parameters, and click "Calculate Network" to download your results.

## Output Files

For a given input `structure.pdb`, Sauron generates:

- **`structure.edges`**: A list of all identified interactions (edges) between residues, including interaction type, distance, and angle.
- **`structure.nodes`**: A list of all residues (nodes) in the structure, including chain, position, DSSP secondary structure, and 3D coordinates.
- **`structure_network_metrics.tsv`**: Full topological network metrics for every node in the structure.
- **`structure_top25_metrics.tsv`**: A summary table of the top 25 nodes ranked by different centrality measures.

If the file contains multiple models, Sauron will generate these outputs individually for each model (`structure_model_1.edges`, etc.).

<div align="center">
  <img src="EvoMol_v2.png" alt="Rictusempra Logo" width="400">
</div>
