import sys

with open('SauronGUI.py', 'r') as f:
    content = f.read()

route_code = """
@app.route('/rinpy_view/<prefix>/<path:filepath>')
def rinpy_view(prefix, filepath):
    rinpy_dir = os.path.join(UPLOAD_FOLDER, f"{prefix}_rinpy_out")
    return send_from_directory(rinpy_dir, filepath)
"""

content = content.replace(route_code, "")

content = content.replace("if __name__ == '__main__':", route_code + "\nif __name__ == '__main__':")

with open('SauronGUI.py', 'w') as f:
    f.write(content)

