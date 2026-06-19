load ./4AKE_temp_laplacian_mode_2.pdb
hide everything
show cartoon
color gray80

select hinge_residues, (chain A and resi 7 and resn GLY) + (chain A and resi 12 and resn GLY) + (chain A and resi 110 and resn ASP) + (chain A and resi 169 and resn VAL) + (chain A and resi 198 and resn GLY) + (chain A and resi 201 and resn PRO) + (chain B and resi 7 and resn GLY) + (chain B and resi 13 and resn LYS) + (chain B and resi 15 and resn THR) + (chain B and resi 44 and resn GLU) + (chain B and resi 109 and resn PHE) + (chain B and resi 171 and resn TYR) + (chain B and resi 173 and resn GLN) + (chain B and resi 195 and resn LYS) + (chain B and resi 196 and resn VAL) + (chain B and resi 202 and resn VAL) + (chain B and resi 204 and resn GLU) + (chain B and resi 205 and resn VAL)
show cartoon, hinge_residues
color yellow, hinge_residues
set stick_radius, 0.25, hinge_residues

select hub_residues, (chain A and resi 29 and resn ILE)+(chain A and resi 30 and resn SER)+(chain A and resi 32 and resn GLY)+(chain A and resi 36 and resn ARG)+(chain A and resi 37 and resn ALA)+(chain A and resi 38 and resn ALA)+(chain A and resi 48 and resn GLN)+(chain A and resi 88 and resn ARG)+(chain A and resi 156 and resn ARG)+(chain B and resi 29 and resn ILE)+(chain B and resi 30 and resn SER)+(chain B and resi 32 and resn GLY)+(chain B and resi 36 and resn ARG)+(chain B and resi 37 and resn ALA)+(chain B and resi 38 and resn ALA)+(chain B and resi 39 and resn VAL)+(chain B and resi 44 and resn GLU)+(chain B and resi 45 and resn LEU)+(chain B and resi 46 and resn GLY)+(chain B and resi 92 and resn GLN)+(chain B and resi 158 and resn ASP)+(chain B and resi 159 and resn ASP)
show cartoon, hub_residues
color warmpink, hub_residues
set stick_radius, 0.25, hub_residues

set label_size, 16
set label_font_id, 7
label hinge_residues and name CA, resn + resi
set label_position, [0.0, -0.3, -0.5], hinge_residues
label hub_residues and name CA, resn + resi
set label_position, [0.0, -0.3, -0.5], hub_residues
set label_color, white


zoom hinge_residues or hub_residues, 10
deselect
