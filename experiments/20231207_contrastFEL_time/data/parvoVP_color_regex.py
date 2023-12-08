import re

text = '''[&!color=#1900ff]
[&!color=#ff0000]
[&!color=#00ff00]
[&!color=#0000ff]'''

pattern = r'\[&!color=#\w{6}\]'

matches = re.findall(pattern, text)

for match in matches:
    print(match)

def replace_color_code(match):
    return match

modified_text = re.sub(pattern, replace_color_code, text)
modified_text

# hardcode matches
sub_text = "#NEXUS\nbegin trees;\n\ttree tree_1 = [&R] ((KP769859_1_Feline_panleukopenia_virus_isolate_MG13[&!color=#1900ff]:0.003563,((KX434462_1_Feline_panleukopenia_virus_strain_FPV_I[&!color=#1900ff]:0.004624,MN127780_1_Feline_panleukopenia_virus_isolate_19R8[&!color=#1900ff]):0.002663,MN127780_1_Feline_panleukopenia_virus_isolate_19R8[&!color=#1900ff]):0.002025);\nend;\n"
sub_text_2="(MW889106_1_Canine_parvovirus_2c_strain_CPV12_comp[&!color=#ff0000]:5.04E-4,(MW889110_1_Canine_parvovirus_2c_strain_CPV16_comp[&!color=#ff0000]:0.001012,MW889108_1_Canine_parvovirus_2c_strain_CPV14_comp[&!color=#ff0000]:0.0)[&!color=#ff0000]:0.001008)[&!color=#ff0000]:0.0)[&!color=#ff0000]:0.001009,"
text = sub_text + sub_text_2
text
import re

feline_pattern = r'\[&!color=#1900ff\]'
canine_pattern = r'\[&!color=#ff0000\]'

updated_text = re.sub(feline_pattern, "{Feline}", text)
final_text = re.sub(feline_pattern, "{Canine}", updated_text)

# regex on the parvoVP

import os
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20231207_contrastFEL_time")
os.listdir()

with open('data/ParvoVP.nex', 'r') as f:
    contents = f.read()


matches = re.findall(pattern, contents)
unique_matches= list(set(matches))
unique_matches

match_dict = {'[&!color=#1900ff]':"{Feline}", '[&!color=#ff0000]':"{Canine}"}



### final regex


sub_text = "#NEXUS\nbegin trees;\n\ttree tree_1 = [&R] ((KP769859_1_Feline_panleukopenia_virus_isolate_MG13[&!color=#1900ff]:0.003563,((KX434462_1_Feline_panleukopenia_virus_strain_FPV_I[&!color=#1900ff]:0.004624,MN127780_1_Feline_panleukopenia_virus_isolate_19R8[&!color=#1900ff]):0.002663,MN127780_1_Feline_panleukopenia_virus_isolate_19R8[&!color=#1900ff]):0.002025);\nend;\n"
sub_text_2="(MW889106_1_Canine_parvovirus_2c_strain_CPV12_comp[&!color=#ff0000]:5.04E-4,(MW889110_1_Canine_parvovirus_2c_strain_CPV16_comp[&!color=#ff0000]:0.001012,MW889108_1_Canine_parvovirus_2c_strain_CPV14_comp[&!color=#ff0000]:0.0)[&!color=#ff0000]:0.001008)[&!color=#ff0000]:0.0)[&!color=#ff0000]:0.001009,"
text = sub_text + sub_text_2
text
import re

feline_pattern = r'\[&!color=#1900ff\]'
canine_pattern = r'\[&!color=#ff0000\]'

#version 1 without colour
updated_text = re.sub(feline_pattern, "{Feline}", text)
final_text = re.sub(canine_pattern, "{Canine}", updated_text)

#verison 2 with colour
updated_text = re.sub(feline_pattern, "[&!color=#1900ff]{Feline}", text)
final_text = re.sub(canine_pattern, "[&!color=#ff0000]{Canine}", updated_text)


