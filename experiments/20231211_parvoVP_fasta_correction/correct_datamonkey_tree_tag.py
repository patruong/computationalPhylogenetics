import re
import os

os.chdir("/home/ptruong/git/computationalPhylogenetics/experiments/20231211_parvoVP_fasta_correction")
os.listdir()

with open('ParvoVP_datamonkey.new', 'r') as f:
    contents = f.read()

contents = contents.replace("_CANINE_:", "{CANINE}:")
contents = contents.replace("_FELINE_:", "{FELINE}:")

output_file_path = 'ParvoVP_datamonkey_add_group.nwk'
with open(output_file_path, 'w') as output_file:
    output_file.write(contents)

print(f"Modified text saved to {output_file_path}")


