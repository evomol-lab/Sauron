import os
import sys
import glob
import zipfile
import subprocess
from flask import Flask, request, render_template, send_file, jsonify

app = Flask(__name__)
UPLOAD_FOLDER = os.path.abspath('uploads')
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload():
    if 'file' not in request.files:
        return jsonify({'error': 'No file uploaded'}), 400
        
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No file selected'}), 400
        
    add_h = request.form.get('add_h') == 'true'
    strict_angle = request.form.get('strict_angle') == 'true'
    remove_multiples = request.form.get('remove_multiples') == 'true'
    
    # Save the file
    filepath = os.path.join(UPLOAD_FOLDER, file.filename)
    file.save(filepath)
    
    # Clear previous output files in UPLOAD_FOLDER
    for ext in ["*.edges", "*.nodes", "*.tsv", "*.zip"]:
        for f in glob.glob(os.path.join(UPLOAD_FOLDER, ext)):
            try: os.remove(f)
            except: pass
    
    # Build command using relative path to the script from UPLOAD_FOLDER
    script_path = os.path.abspath('sauron.py')
    cmd = [sys.executable, script_path, file.filename]
    if add_h: cmd.append('--add-h')
    if strict_angle: cmd.append('--strict-angle')
    if remove_multiples: cmd.append('--remove-multiples')
    
    # Run the script inside the UPLOAD_FOLDER as CWD to keep output files there
    try:
        subprocess.run(cmd, check=True, cwd=UPLOAD_FOLDER, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        return jsonify({'error': f"Computation failed:\\n{e.stderr}"}), 500
        
    # Zip the generated files
    prefix = os.path.splitext(file.filename)[0]
    zip_name = f"{prefix}_rin_results.zip"
    zip_path = os.path.join(UPLOAD_FOLDER, zip_name)
    
    out_files = glob.glob(os.path.join(UPLOAD_FOLDER, "*.edges")) + \
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.nodes")) + \
                glob.glob(os.path.join(UPLOAD_FOLDER, "*.tsv"))
                
    if not out_files:
        return jsonify({'error': 'No output files were generated.'}), 500
        
    with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for f in out_files:
            zipf.write(f, os.path.basename(f))
            
    return jsonify({'downloadUrl': f'/download/{zip_name}'})

@app.route('/download/<filename>')
def download(filename):
    path = os.path.join(UPLOAD_FOLDER, filename)
    if os.path.exists(path):
        return send_file(path, as_attachment=True)
    return "File not found", 404

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=5000, debug=True)
