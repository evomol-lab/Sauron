import os
import sys
import glob
import zipfile
import subprocess
import requests
from flask import Flask, request, render_template, send_file, jsonify, send_from_directory

if getattr(sys, 'frozen', False):
    template_folder = os.path.join(sys._MEIPASS, 'templates')
    static_folder = os.path.join(sys._MEIPASS, 'static')
    app = Flask(__name__, template_folder=template_folder, static_folder=static_folder)
else:
    app = Flask(__name__)

UPLOAD_FOLDER = os.environ.get('SAURON_UPLOAD_DIR', os.path.abspath('uploads'))
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload():
    input_mode = request.form.get('input_mode', 'upload')
    
    if input_mode == 'upload':
        if 'file' not in request.files:
            return jsonify({'error': 'No file uploaded'}), 400
        file = request.files['file']
        if file.filename == '':
            return jsonify({'error': 'No file selected'}), 400
        filename = file.filename
        filepath = os.path.join(UPLOAD_FOLDER, filename)
        file.save(filepath)
    else:
        fetch_db = request.form.get('fetch_db')
        fetch_id = request.form.get('fetch_id', '').strip()
        if not fetch_id:
            return jsonify({'error': 'No ID provided for fetching'}), 400
            
        if fetch_db == 'pdb':
            url = f"https://files.rcsb.org/download/{fetch_id.upper()}.cif"
            filename = f"{fetch_id.upper()}.cif"
        elif fetch_db == 'alphafold':
            url = f"https://alphafold.ebi.ac.uk/files/AF-{fetch_id.upper()}-F1-model_v6.cif"
            filename = f"AF-{fetch_id.upper()}.cif"
        else:
            return jsonify({'error': 'Invalid database selected'}), 400
            
        filepath = os.path.join(UPLOAD_FOLDER, filename)
        try:
            r = requests.get(url, allow_redirects=True)
            r.raise_for_status()
            with open(filepath, 'wb') as f:
                f.write(r.content)
        except Exception as e:
            return jsonify({'error': f"Failed to fetch structure from {fetch_db.upper()}: {str(e)}"}), 400

    calc_method = request.form.get('calc_method', 'standard')
    strict_angle = True
    remove_multiples = True
    chains_opt = request.form.get('chains_opt') # 'all' or 'specific'
    chains_input = request.form.get('chains_input')
    calc_energy = request.form.get('calc_energy') == 'true'
    vdw_energy = request.form.get('vdw_energy') == 'true'
    
    # Clear previous output files in UPLOAD_FOLDER
    for ext in ["*.edges", "*.nodes", "*.tsv", "*.zip", "*.log", "*_degree.pdb"]:
        for f in glob.glob(os.path.join(UPLOAD_FOLDER, ext)):
            try: os.remove(f)
            except: pass
    import shutil
    for d in glob.glob(os.path.join(UPLOAD_FOLDER, "*_rinpy_out")):
        try: shutil.rmtree(d)
        except: pass
    
    # Build command for sauron.py processing
    if getattr(sys, 'frozen', False):
        cmd = [sys.executable, "--run-sauron", filename]
    else:
        cmd = [sys.executable, os.path.abspath(__file__), "--run-sauron", filename]

    if calc_energy or calc_method == 'insty':
        cmd.append('--add-h')
    cmd.extend(['--calc-method', calc_method])
    if strict_angle: cmd.append('--strict-angle')
    if remove_multiples: cmd.append('--remove-multiples')
    if calc_energy: cmd.append('--energy')
    if vdw_energy: cmd.append('--vdw-energy')
    if chains_opt == 'specific' and chains_input and chains_input.strip():
        cmd.extend(['--chains', chains_input.strip().upper()])
    
    # Run the script inside the UPLOAD_FOLDER as CWD to keep output files there
    try:
        subprocess.run(cmd, check=True, cwd=UPLOAD_FOLDER, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        return jsonify({'error': f"Computation failed:\n{e.stderr}"}), 500
        
    # Zip the generated files
    prefix = os.path.splitext(filename)[0]
    zip_name = f"{prefix}_rin_results.zip"
    zip_path = os.path.join(UPLOAD_FOLDER, zip_name)
    
    out_files = glob.glob(os.path.join(UPLOAD_FOLDER, "*.edges")) + \
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.nodes")) + \
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.tsv")) + \
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.log")) + \
                glob.glob(os.path.join(UPLOAD_FOLDER, "*_degree.pdb"))
                
    prefix = os.path.splitext(filename)[0]
    rinpy_out_dir = os.path.join(UPLOAD_FOLDER, f"{prefix}_rinpy_out")
    
    is_rinpy = calc_method == 'rinpy'
    rinpy_html = None
    rinpy_summary_data = None
    
    summary_path = os.path.join(UPLOAD_FOLDER, f"{prefix}_rinpy_summary.json")
    if os.path.exists(summary_path):
        import json
        with open(summary_path, 'r') as f:
            rinpy_summary_data = json.load(f)

    
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
            # add json
            if os.path.exists(summary_path):
                zipf.write(summary_path, os.path.basename(summary_path))

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
    
    return jsonify({
        'downloadUrl': f'/download/{zip_name}',
        'pdbFile': filename,
        'edgesFile': edges_file,
        'nodesFile': nodes_file,
        'metricsFile': metrics_file,
        'isMultimodel': is_multimodel,
        'isRinpy': is_rinpy,
        'rinpyHtml': rinpy_html,
        'rinpyPrefix': prefix,
        'rinpySummary': rinpy_summary_data
    })

@app.route('/download/<filename>')
def download(filename):
    path = os.path.join(UPLOAD_FOLDER, filename)
    if os.path.exists(path):
        return send_file(path, as_attachment=True)
    return "File not found", 404

@app.route('/view/<filename>')
def view_file(filename):
    path = os.path.join(UPLOAD_FOLDER, filename)
    if os.path.exists(path):
        return send_file(path)
    return "File not found", 404


@app.route('/rinpy_view/<prefix>/<path:filepath>')
def rinpy_view(prefix, filepath):
    rinpy_dir = os.path.join(UPLOAD_FOLDER, f"{prefix}_rinpy_out")
    return send_from_directory(rinpy_dir, filepath)

if __name__ == '__main__':
    if len(sys.argv) > 1:
        if sys.argv[1] == "--run-sauron":
            import sauron
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            sauron.main()
            sys.exit(0)
        elif sys.argv[1] == "--run-pdb2pqr":
            from pdb2pqr.main import main as pdb2pqr_main
            sys.argv = [sys.argv[0]] + sys.argv[2:]
            sys.exit(pdb2pqr_main())

    import threading
    import webbrowser
    def open_browser():
        webbrowser.open_new("http://127.0.0.1:5000/")
    
    if getattr(sys, 'frozen', False):
        threading.Timer(1.25, open_browser).start()
        
    app.run(host='0.0.0.0', port=5000, debug=not getattr(sys, 'frozen', False))

