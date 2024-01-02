import re
import os
import random

os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20231221_removed_duplicated_nexus_difFUBAR")
os.listdir()

with open('ParvoVP.fasta', 'r') as f:
    fasta = f.readlines()

seq_names = []
seqs = []
for line in fasta:
    if line[0] == ">":
        seq_names.append(line)
    else:
        seqs.append(line)
seq_names
seqs

len(set(seqs)) # 320 unique sequences
len(set(seq_names)) # 371 unique sequence names

reduced_seq_names = []
reduced_seqs = []
seq_list = list(zip(seq_names, seqs))
shuffled_seq_list = seq_list.copy()  # or shuffled_seq_list = seq_list[:]
#random.shuffle(shuffled_seq_list)
for seq_name, seq in shuffled_seq_list:
    if seq not in reduced_seqs:
        reduced_seq_names.append(seq_name)
        reduced_seqs.append(seq)
    else:
        pass

reduced_seq_names[0]
reduced_seqs[0]

# fix names ending with -


fixed_seq_names = []
for seq in reduced_seq_names:
    seq = seq.replace("-", "_")
    seq = seq.replace("__", "_")
    if seq[-2] == "_":
        seq = seq[:-2] + "\n"
        fixed_seq_names.append(seq)
    else:
        fixed_seq_names.append(seq)

fixed_seq_names

added_tag_seq_names = []
for seq_name, seq in zip(fixed_seq_names, reduced_seqs):
    print(seq_name)
    if "Canine" in seq_name:
        seq_name = seq_name[:-1] + "{CANINE}\n"
    elif "Carnivore" in seq_name:
        seq_name = seq_name[:-1] + "{CANINE}\n"
    elif "Feline" in seq_name:
        seq_name = seq_name[:-1] + "{FELINE}\n"
    else:
        seq_name = seq_name
    added_tag_seq_names.append(seq_name)
seq_name

added_tag_seq_names
fixed_seq_names


# For datamonkey tree
with open("ParvoVP_removed_duplicates.fasta", "w") as file:
    for seq_name, seq in zip(added_tag_seq_names, reduced_seqs):
        file.write(f"{seq_name.upper()}{seq.upper()}")

# without tag for the contrast-fel local run
with open("ParvoVP_removed_duplicates.fasta", "w") as file:
    for seq_name, seq in zip(fixed_seq_names, reduced_seqs):
        file.write(seq_name.upper())
        file.write(seq.upper())
        #file.write(f"{seq_name.upper()}{seq.upper()}")





