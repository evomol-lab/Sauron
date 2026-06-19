import sys

with open('SauronGUI.py', 'r') as f:
    content = f.read()

old_clean = """
    # Clear previous output files in UPLOAD_FOLDER
    for ext in ["*.edges", "*.nodes", "*.tsv", "*.zip", "*.log", "*_degree.pdb"]:
        for f in glob.glob(os.path.join(UPLOAD_FOLDER, ext)):
            try: os.remove(f)
            except: pass
"""

new_clean = """
    # Clear previous output files in UPLOAD_FOLDER
    for ext in ["*.edges", "*.nodes", "*.tsv", "*.zip", "*.log", "*_degree.pdb"]:
        for f in glob.glob(os.path.join(UPLOAD_FOLDER, ext)):
            try: os.remove(f)
            except: pass
    import shutil
    for d in glob.glob(os.path.join(UPLOAD_FOLDER, "*_rinpy_out")):
        try: shutil.rmtree(d)
        except: pass
"""

content = content.replace(old_clean, new_clean)

with open('SauronGUI.py', 'w') as f:
    f.write(content)

