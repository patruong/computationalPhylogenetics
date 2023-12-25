import os
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/data/20231221_feline_canine_parvoVP_build/selected_sequences_based_on_length_from_full_fasta/tagged")

tree_file = "parvo_panleu_trans_tagged.fasta"

with open(tree_file, "r") as file:
    file_contents = file.read()

updated_contents = file_contents.translate(str.maketrans("ÖÄ", "{}"))

# Specify the path to your output text file
output_file_path = tree_file[:-6] + "_fixed.fasta"

# Open the file in write mode and write the translated contents
with open(output_file_path, "w") as output_file:
    output_file.write(updated_contents)