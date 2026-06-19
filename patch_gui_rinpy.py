import sys

with open('SauronGUI.py', 'r') as f:
    content = f.read()

# Replace the zipping logic
old_zip_logic = """
    out_files = glob.glob(os.path.join(UPLOAD_FOLDER, "*.edges")) + \\
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.nodes")) + \\
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.tsv")) + \\
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.log")) + \\
                glob.glob(os.path.join(UPLOAD_FOLDER, "*_degree.pdb"))
                
    if not out_files:
        return jsonify({'error': 'No output files were generated.'}), 500
        
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for f in out_files:
            zipf.write(f, os.path.basename(f))
            
    # Identify model network files if they exist
    edges_file = None
    nodes_file = None
    metrics_file = None
    for f in out_files:
        basename = os.path.basename(f)
        if basename.endswith(".edges") and (edges_file is None or "model_1." in basename):
            edges_file = basename
        if basename.endswith(".nodes") and (nodes_file is None or "model_1." in basename):
            nodes_file = basename
        if basename.endswith("network_metrics.tsv") and (metrics_file is None or "model_1_" in basename):
            metrics_file = basename
            
    is_multimodel = any("model_2" in f for f in out_files)
"""

new_zip_logic = """
    out_files = glob.glob(os.path.join(UPLOAD_FOLDER, "*.edges")) + \\
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.nodes")) + \\
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.tsv")) + \\
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.log")) + \\
                glob.glob(os.path.join(UPLOAD_FOLDER, "*_degree.pdb"))
                
    prefix = os.path.splitext(filename)[0]
    rinpy_out_dir = os.path.join(UPLOAD_FOLDER, f"{prefix}_rinpy_out")
    
    is_rinpy = calc_method == 'rinpy'
    rinpy_html = None
    
    if is_rinpy:
        if not os.path.exists(rinpy_out_dir):
            return jsonify({'error': 'RinPy output directory not found.'}), 500
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for root, dirs, files in os.walk(rinpy_out_dir):
                for f in files:
                    file_path = os.path.join(root, f)
                    arcname = os.path.relpath(file_path, rinpy_out_dir)
                    zipf.write(file_path, arcname)
                    if f.endswith("interactive_clusters_3d.html") and rinpy_html is None:
                        # Find the first HTML plot to show in GUI
                        rinpy_html = arcname
            # also add the log file
            log_files = glob.glob(os.path.join(UPLOAD_FOLDER, "*.log"))
            for lf in log_files:
                zipf.write(lf, os.path.basename(lf))
    else:
        if not out_files:
            return jsonify({'error': 'No output files were generated.'}), 500
        with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
            for f in out_files:
                zipf.write(f, os.path.basename(f))
                
    # Identify model network files if they exist
    edges_file = None
    nodes_file = None
    metrics_file = None
    if not is_rinpy:
        for f in out_files:
            basename = os.path.basename(f)
            if basename.endswith(".edges") and (edges_file is None or "model_1." in basename):
                edges_file = basename
            if basename.endswith(".nodes") and (nodes_file is None or "model_1." in basename):
                nodes_file = basename
            if basename.endswith("network_metrics.tsv") and (metrics_file is None or "model_1_" in basename):
                metrics_file = basename
            
    is_multimodel = not is_rinpy and any("model_2" in f for f in out_files)
"""

content = content.replace(old_zip_logic, new_zip_logic)

# Replace the return jsonify
old_return = """
    return jsonify({
        'downloadUrl': f'/download/{zip_name}',
        'pdbFile': filename,
        'edgesFile': edges_file,
        'nodesFile': nodes_file,
        'metricsFile': metrics_file,
        'isMultimodel': is_multimodel
    })
"""

new_return = """
    return jsonify({
        'downloadUrl': f'/download/{zip_name}',
        'pdbFile': filename,
        'edgesFile': edges_file,
        'nodesFile': nodes_file,
        'metricsFile': metrics_file,
        'isMultimodel': is_multimodel,
        'isRinpy': is_rinpy,
        'rinpyHtml': rinpy_html,
        'rinpyPrefix': prefix
    })
"""

content = content.replace(old_return, new_return)

# Also need to add a route to serve RinPy HTML files!
new_route = """
@app.route('/rinpy_view/<prefix>/<path:filepath>')
def rinpy_view(prefix, filepath):
    rinpy_dir = os.path.join(UPLOAD_FOLDER, f"{prefix}_rinpy_out")
    return send_from_directory(rinpy_dir, filepath)

"""

if "@app.route('/rinpy_view" not in content:
    content += new_route

with open('SauronGUI.py', 'w') as f:
    f.write(content)

