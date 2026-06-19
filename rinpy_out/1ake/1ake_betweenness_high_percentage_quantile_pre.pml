load ./1ake_pre.pdb
hide everything
show cartoon
color gray80

select residues, (chain A and resi 7 and resn GLY)+(chain A and resi 88 and resn ARG)+(chain A and resi 158 and resn ASP)+(chain A and resi 167 and resn ARG)+(chain A and resi 168 and resn LEU)+(chain A and resi 169 and resn VAL)+(chain A and resi 171 and resn TYR)+(chain A and resi 172 and resn HIS)+(chain A and resi 173 and resn GLN)+(chain A and resi 174 and resn MET)+(chain A and resi 175 and resn THR)+(chain B and resi 14 and resn GLY)+(chain B and resi 18 and resn GLN)+(chain B and resi 131 and resn ARG)+(chain B and resi 146 and resn ASP)+(chain B and resi 148 and resn VAL)+(chain B and resi 149 and resn THR)+(chain B and resi 150 and resn GLY)+(chain B and resi 153 and resn LEU)+(chain B and resi 154 and resn THR)+(chain B and resi 155 and resn THR)+(chain B and resi 156 and resn ARG)
show cartoon, residues
color warmpink, residues
set stick_radius, 0.25, residues

set label_size, 16
set label_font_id, 7
label residues and name CA, resn + resi
set label_position, [0.0, -0.3, -0.5], residues
set label_color, white


zoom residues, 10
deselect
