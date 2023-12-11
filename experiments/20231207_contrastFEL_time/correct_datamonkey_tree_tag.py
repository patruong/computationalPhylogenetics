import re
import os

os.chdir("/home/ptruong/git/computationalPhylogenetics/experiments/20231207_contrastFEL_time/data/datamonkey_tree")
os.listdir()

with open('tree_rename.nwk', 'r') as f:
    contents = f.read()

contents = contents.replace("_CANINE_:", "{CANINE}:")
contents = contents.replace("_FELINE_:", "{FELINE}:")

output_file_path = 'tree_tag.nwk'
with open(output_file_path, 'w') as output_file:
    output_file.write(contents)

print(f"Modified text saved to {output_file_path}")


