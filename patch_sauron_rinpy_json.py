import sys

with open('sauron.py', 'r') as f:
    content = f.read()

json_logic = """
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
                    
"""

# Insert json_logic after log writing and before rm temp_pdb
content = content.replace(
    "                if os.path.exists(temp_pdb): os.remove(temp_pdb)",
    json_logic + "                if os.path.exists(temp_pdb): os.remove(temp_pdb)"
)

with open('sauron.py', 'w') as f:
    f.write(content)

