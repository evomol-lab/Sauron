import sys
import re

with open('sauron.py', 'r') as f:
    content = f.read()

# Modify analyze_network to return degree_dict
content = content.replace(
    "def analyze_network(edges_df, nodes, output_prefix):",
    "def analyze_network(edges_df, nodes, output_prefix):\n    degree_dict_out = {}"
)
content = content.replace(
    "            f.write(f\"{n}\\t{deg}\\t{cc:.4f}\\t{bc:.4f}\\t{ec:.4f}\\t{sp:.4f}\\n\")",
    "            f.write(f\"{n}\\t{deg}\\t{cc:.4f}\\t{bc:.4f}\\t{ec:.4f}\\t{sp:.4f}\\n\")\n            degree_dict_out[n] = deg"
)
content = content.replace(
    "    print(f\"Network metrics saved to {output_prefix}_network_metrics.tsv and {output_prefix}_top25_metrics.tsv\")",
    "    print(f\"Network metrics saved to {output_prefix}_network_metrics.tsv and {output_prefix}_top25_metrics.tsv\")\n    return degree_dict_out"
)

# Call analyze_network and write PDB
new_main_end = """
        degree_dict = analyze_network(edges_df, nodes, output_prefix)
        print(f"Generated {edges_out} and {nodes_out}")
        
        # Write Log
        import datetime
        with open(f"{output_prefix}_sauron.log", 'w') as logf:
            logf.write(f"Sauron RIN Calculation Log\\n")
            logf.write(f"Timestamp: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\\n")
            logf.write(f"Input File: {args.input}\\n")
            logf.write(f"Method: {args.calc_method}\\n")
            logf.write(f"Strict Angle: {args.strict_angle}\\n")
            logf.write(f"Remove Multiples: {args.remove_multiples}\\n")
            logf.write(f"Energy (PEN): {args.energy}\\n")
            logf.write(f"Include VDW: {args.vdw_energy}\\n")
            logf.write(f"Chains: {args.chains if args.chains else 'ALL'}\\n")
            logf.write(f"Total Nodes: {len(nodes)}\\n")
            logf.write(f"Total Edges: {len(edges_df)}\\n")

        # Write Degree PDB
        from Bio.PDB import PDBIO
        for model in orig_struct:
            for chain in model:
                for res in chain:
                    node_id = format_node_id(res)
                    deg = degree_dict.get(node_id, 0)
                    for atom in res:
                        atom.set_bfactor(deg)
        io = PDBIO()
        io.set_structure(orig_struct)
        io.save(f"{output_prefix}_degree.pdb")
"""

content = content.replace(
    "        analyze_network(edges_df, nodes, output_prefix)\n        print(f\"Generated {edges_out} and {nodes_out}\")",
    new_main_end
)

with open('sauron.py', 'w') as f:
    f.write(content)

