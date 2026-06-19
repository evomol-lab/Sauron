load ./1ake_original_centrality_degree.pdb
hide everything
show cartoon
color gray80

select residues, (chain A and resi 3 and resn ILE)+(chain A and resi 4 and resn ILE)+(chain A and resi 5 and resn LEU)+(chain A and resi 13 and resn LYS)+(chain A and resi 16 and resn GLN)+(chain A and resi 17 and resn ALA)+(chain A and resi 19 and resn PHE)+(chain A and resi 20 and resn ILE)+(chain A and resi 21 and resn MET)+(chain A and resi 26 and resn ILE)+(chain A and resi 29 and resn ILE)+(chain A and resi 31 and resn THR)+(chain A and resi 35 and resn LEU)+(chain A and resi 38 and resn ALA)+(chain A and resi 39 and resn VAL)+(chain A and resi 53 and resn MET)+(chain A and resi 64 and resn VAL)+(chain A and resi 65 and resn ILE)+(chain A and resi 71 and resn ARG)+(chain A and resi 81 and resn PHE)+(chain A and resi 82 and resn LEU)+(chain A and resi 84 and resn ASP)+(chain A and resi 86 and resn PHE)+(chain A and resi 92 and resn GLN)+(chain A and resi 93 and resn ALA)+(chain A and resi 96 and resn MET)+(chain A and resi 105 and resn TYR)+(chain A and resi 106 and resn VAL)+(chain A and resi 109 and resn PHE)+(chain A and resi 116 and resn ILE)+(chain A and resi 119 and resn ARG)+(chain A and resi 120 and resn ILE)+(chain A and resi 123 and resn ARG)+(chain A and resi 124 and resn ARG)+(chain A and resi 126 and resn HIS)+(chain A and resi 133 and resn TYR)+(chain A and resi 163 and resn THR)+(chain A and resi 164 and resn VAL)+(chain A and resi 167 and resn ARG)+(chain A and resi 169 and resn VAL)+(chain A and resi 170 and resn GLU)+(chain A and resi 171 and resn TYR)+(chain A and resi 172 and resn HIS)+(chain A and resi 173 and resn GLN)+(chain A and resi 174 and resn MET)+(chain A and resi 178 and resn LEU)+(chain A and resi 179 and resn ILE)+(chain A and resi 181 and resn TYR)+(chain A and resi 182 and resn TYR)+(chain A and resi 205 and resn VAL)+(chain A and resi 209 and resn LEU)+(chain B and resi 3 and resn ILE)+(chain B and resi 4 and resn ILE)+(chain B and resi 5 and resn LEU)+(chain B and resi 6 and resn LEU)+(chain B and resi 13 and resn LYS)+(chain B and resi 16 and resn GLN)+(chain B and resi 17 and resn ALA)+(chain B and resi 19 and resn PHE)+(chain B and resi 20 and resn ILE)+(chain B and resi 26 and resn ILE)+(chain B and resi 29 and resn ILE)+(chain B and resi 35 and resn LEU)+(chain B and resi 36 and resn ARG)+(chain B and resi 38 and resn ALA)+(chain B and resi 39 and resn VAL)+(chain B and resi 53 and resn MET)+(chain B and resi 64 and resn VAL)+(chain B and resi 71 and resn ARG)+(chain B and resi 81 and resn PHE)+(chain B and resi 84 and resn ASP)+(chain B and resi 86 and resn PHE)+(chain B and resi 92 and resn GLN)+(chain B and resi 93 and resn ALA)+(chain B and resi 96 and resn MET)+(chain B and resi 105 and resn TYR)+(chain B and resi 106 and resn VAL)+(chain B and resi 107 and resn LEU)+(chain B and resi 109 and resn PHE)+(chain B and resi 116 and resn ILE)+(chain B and resi 119 and resn ARG)+(chain B and resi 120 and resn ILE)+(chain B and resi 123 and resn ARG)+(chain B and resi 124 and resn ARG)+(chain B and resi 133 and resn TYR)+(chain B and resi 145 and resn LYS)+(chain B and resi 163 and resn THR)+(chain B and resi 167 and resn ARG)+(chain B and resi 168 and resn LEU)+(chain B and resi 170 and resn GLU)+(chain B and resi 171 and resn TYR)+(chain B and resi 178 and resn LEU)+(chain B and resi 181 and resn TYR)+(chain B and resi 182 and resn TYR)+(chain B and resi 205 and resn VAL)+(chain B and resi 209 and resn LEU)
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
