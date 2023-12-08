import re
import os

os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20231207_contrastFEL_time")
os.listdir()

with open('data/ParvoVP.nex', 'r') as f:
    contents = f.read()

feline_pattern = r'\[&!color=#1900ff\]'
canine_pattern = r'\[&!color=#ff0000\]'

#version 1 without colour
updated_text = re.sub(feline_pattern, "{Feline}", contents)
final_text = re.sub(canine_pattern, "{Canine}", updated_text)

# Save the modified text to a new file
output_file_path = 'data/ParvoVP_modified_no_color.nex'
with open(output_file_path, 'w') as output_file:
    output_file.write(final_text)

print(f"Modified text saved to {output_file_path}")
#verison 2 with colour
updated_text = re.sub(feline_pattern, "[&!color=#1900ff]{Feline}", contents)
final_text = re.sub(canine_pattern, "[&!color=#ff0000]{Canine}", updated_text)
output_file_path = 'data/ParvoVP_modified_color.nex'
with open(output_file_path, 'w') as output_file:
    output_file.write(final_text)

