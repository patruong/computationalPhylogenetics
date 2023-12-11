import re
import os

os.chdir("/home/ptruong/git/computationalPhylogenetics/experiments/20231207_contrastFEL_time/data/")
os.listdir()

with open('ParvoVP.fasta', 'r') as f:
    contents = f.read()
