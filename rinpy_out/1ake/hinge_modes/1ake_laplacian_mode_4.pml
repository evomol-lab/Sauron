load ./1ake_laplacian_mode_4.pdb
hide everything
show cartoon
color gray80

select hinge_residues, (chain A and resi 8 and resn ALA) + (chain A and resi 12 and resn GLY) + (chain A and resi 55 and resn ALA) + (chain A and resi 57 and resn LYS) + (chain A and resi 58 and resn LEU) + (chain A and resi 88 and resn ARG) + (chain A and resi 93 and resn ALA) + (chain A and resi 112 and resn PRO) + (chain A and resi 127 and resn ALA) + (chain A and resi 130 and resn GLY) + (chain A and resi 157 and resn LYS) + (chain A and resi 183 and resn SER) + (chain B and resi 9 and resn PRO) + (chain B and resi 12 and resn GLY) + (chain B and resi 14 and resn GLY) + (chain B and resi 31 and resn THR) + (chain B and resi 68 and resn VAL) + (chain B and resi 69 and resn LYS) + (chain B and resi 71 and resn ARG) + (chain B and resi 112 and resn PRO) + (chain B and resi 174 and resn MET) + (chain B and resi 175 and resn THR)
show cartoon, hinge_residues
color yellow, hinge_residues
set stick_radius, 0.25, hinge_residues

select hub_residues, (chain A and resi 7 and resn GLY)+(chain A and resi 88 and resn ARG)+(chain A and resi 158 and resn ASP)+(chain A and resi 167 and resn ARG)+(chain A and resi 168 and resn LEU)+(chain A and resi 169 and resn VAL)+(chain A and resi 171 and resn TYR)+(chain A and resi 172 and resn HIS)+(chain A and resi 173 and resn GLN)+(chain A and resi 174 and resn MET)+(chain A and resi 175 and resn THR)+(chain B and resi 14 and resn GLY)+(chain B and resi 18 and resn GLN)+(chain B and resi 131 and resn ARG)+(chain B and resi 146 and resn ASP)+(chain B and resi 148 and resn VAL)+(chain B and resi 149 and resn THR)+(chain B and resi 150 and resn GLY)+(chain B and resi 153 and resn LEU)+(chain B and resi 154 and resn THR)+(chain B and resi 155 and resn THR)+(chain B and resi 156 and resn ARG)
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
