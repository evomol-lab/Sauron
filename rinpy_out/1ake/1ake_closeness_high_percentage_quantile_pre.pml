load ./1ake_pre.pdb
hide everything
show cartoon
color gray80

select residues, (chain A and resi 7 and resn GLY)+(chain A and resi 8 and resn ALA)+(chain A and resi 9 and resn PRO)+(chain A and resi 88 and resn ARG)+(chain A and resi 167 and resn ARG)+(chain A and resi 168 and resn LEU)+(chain A and resi 169 and resn VAL)+(chain A and resi 170 and resn GLU)+(chain A and resi 171 and resn TYR)+(chain A and resi 172 and resn HIS)+(chain A and resi 173 and resn GLN)+(chain A and resi 174 and resn MET)+(chain A and resi 175 and resn THR)+(chain A and resi 176 and resn ALA)+(chain B and resi 18 and resn GLN)+(chain B and resi 125 and resn VAL)+(chain B and resi 126 and resn HIS)+(chain B and resi 127 and resn ALA)+(chain B and resi 128 and resn PRO)+(chain B and resi 129 and resn SER)+(chain B and resi 130 and resn GLY)+(chain B and resi 131 and resn ARG)+(chain B and resi 132 and resn VAL)+(chain B and resi 133 and resn TYR)+(chain B and resi 134 and resn HIS)+(chain B and resi 138 and resn ASN)+(chain B and resi 139 and resn PRO)+(chain B and resi 140 and resn PRO)+(chain B and resi 141 and resn LYS)+(chain B and resi 142 and resn VAL)+(chain B and resi 143 and resn GLU)+(chain B and resi 144 and resn GLY)+(chain B and resi 145 and resn LYS)+(chain B and resi 146 and resn ASP)+(chain B and resi 147 and resn ASP)+(chain B and resi 148 and resn VAL)+(chain B and resi 149 and resn THR)+(chain B and resi 150 and resn GLY)+(chain B and resi 151 and resn GLU)+(chain B and resi 152 and resn GLU)+(chain B and resi 153 and resn LEU)+(chain B and resi 154 and resn THR)+(chain B and resi 155 and resn THR)
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
