import re
import os

os.chdir("/home/ptruong/git/computationalPhylogenetics/experiments/20231211_parvoVP_fasta_correction")
os.listdir()

with open('ParvoVP_datamonkey.new', 'r') as f:
    contents = f.read()

#input_string = '(MF177244_1_CANINE_PARVOVIRUS_2C_ISOLATE_ARG33_NON_:0.0004702093650215471,KM457131_1_CANINE_PARVOVIRUS_2C_ISOLATE_UY368_COM:0.0009401978173969662,(MN451675_1_CANINE_PARVOVIRUS_ISOLATE_CPV601_COMPL:0.001439343144448098,(MW732026_1_CANINE_PARVOVIRUS_STRAIN_RVC50_KNA_NON_:0,(MW732025_1_CANINE_PARVOVIRUS_STRAIN_CN10_NON_STRUC:0.0009360756524870998,MW732024_1_CANINE_PARVOVIRUS_STRAIN_CN20_NON_STRUC:0.0004684842564962879):0.0004682017575046694,MW732023_1_CANINE_PARVOVIRUS_STRAIN_CN14_NON_STRUC:0.0009387147202990378):0.003818971743603336):0.0004515118618758366,(EU659121_1_CANINE_PARVOVIRUS_ISOLATE_CPV_411B_US_:0.001410270880106874,AY742936_1_CANINE_PARVOVIRUS_STRAIN_CPV_395_NONSTR:0):0.0009419820228860797,(((KM457128_1_CANINE_PARVOVIRUS_2C_ISOLATE_UY346_COM:0.0009425915386684539,MF177286_1_CANINE_PARVOVIRUS_2C_ISOLATE_UY375_NON_:0.0004704054637333875,KM457123_1_CANINE_PARVOVIRUS_2C_ISOLATE_UY261_COM:0):0.0004698299651243442,KM457109_1_CANINE_PARVOVIRUS_2C_ISOLATE_UY95_COMP:0.001411313894637966):0.0004709225651318982,KM457116_1_CANINE_PARVOVIRUS_2C_ISOLATE_UY185_COM:0.000469684273404215):0.0004700956829562168,(MF177262_1_CANINE_PARVOVIRUS_2C_ISOLATE_PY1_NON_ST:0.0009418192135609029,((MF177283_1_CANINE_PARVOVIRUS_2C_ISOLATE_UY181_NON_:0.0004684618320378695,MF177284_1_CANINE_PARVOVIRUS_2C_ISOLATE_UY196_NON_:0,MF177285_1_CANINE_PARVOVIRUS_2C_ISOLATE_UY269_NON_:0.0009371839382296125):0.0004691121653336132,MF177249_1_CANINE_PARVOVIRUS_2C_ISOLATE_ARG71_NON_:0.0009396850444780578,MF177247_1_CANINE_PARVOVIRUS_2C_ISOLATE_ARG66_NON_:0.0004689415976634255)'

# Define a regular expression pattern to find occurrences of "CANINE" or "CARNIVORE" in the string
canine_pattern = re.compile(r'(CANINE|CARNIVORE)[^:]+')
feline_pattern = re.compile(r'(FELINE)[^:]+')

# Use the sub() function to add {CANINE} before the first colon in each match
output_string = canine_pattern.sub(r'\g<0>{CANINE}', contents)
final_string = feline_pattern.sub(r'\g<0>{FELINE}', output_string)

# Save the modified text to a new file
output_file_path = 'ParvoVP_datamonkey_add_group.nwk'
with open(output_file_path, 'w') as output_file:
    output_file.write(final_string)

print(f"Modified text saved to {output_file_path}")