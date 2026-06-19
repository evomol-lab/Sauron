with open('SauronGUI.py', 'r') as f:
    content = f.read()

content = content.replace(
    "chains_input = request.form.get('chains_input')",
    "chains_input = request.form.get('chains_input')\n    calc_energy = request.form.get('calc_energy') == 'true'\n    vdw_energy = request.form.get('vdw_energy') == 'true'"
)

content = content.replace(
    "if remove_multiples: cmd.append('--remove-multiples')",
    "if remove_multiples: cmd.append('--remove-multiples')\n    if calc_energy: cmd.append('--energy')\n    if vdw_energy: cmd.append('--vdw-energy')"
)

with open('SauronGUI.py', 'w') as f:
    f.write(content)

