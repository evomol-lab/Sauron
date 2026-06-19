import sys

with open('SauronGUI.py', 'r') as f:
    content = f.read()

# Add rinpy_summary logic
old_is_rinpy = "    is_rinpy = calc_method == 'rinpy'\n    rinpy_html = None"
new_is_rinpy = """    is_rinpy = calc_method == 'rinpy'
    rinpy_html = None
    rinpy_summary_data = None
    
    summary_path = os.path.join(UPLOAD_FOLDER, f"{prefix}_rinpy_summary.json")
    if os.path.exists(summary_path):
        import json
        with open(summary_path, 'r') as f:
            rinpy_summary_data = json.load(f)
"""

content = content.replace(old_is_rinpy, new_is_rinpy)

# Add to the zip!
old_zip_rinpy = """            # also add the log file
            log_files = glob.glob(os.path.join(UPLOAD_FOLDER, "*.log"))
            for lf in log_files:
                zipf.write(lf, os.path.basename(lf))"""

new_zip_rinpy = """            # also add the log file
            log_files = glob.glob(os.path.join(UPLOAD_FOLDER, "*.log"))
            for lf in log_files:
                zipf.write(lf, os.path.basename(lf))
            # add json
            if os.path.exists(summary_path):
                zipf.write(summary_path, os.path.basename(summary_path))
"""

content = content.replace(old_zip_rinpy, new_zip_rinpy)

# Add to jsonify
old_return = """        'isRinpy': is_rinpy,
        'rinpyHtml': rinpy_html,
        'rinpyPrefix': prefix
    })"""

new_return = """        'isRinpy': is_rinpy,
        'rinpyHtml': rinpy_html,
        'rinpyPrefix': prefix,
        'rinpySummary': rinpy_summary_data
    })"""

content = content.replace(old_return, new_return)

with open('SauronGUI.py', 'w') as f:
    f.write(content)

