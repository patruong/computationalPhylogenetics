
import os
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20240215_nextclade_pipeline/results/nextclade_output")

with open('nextclade.aligned.fasta', 'r') as file:
    content = file.read()

modified_content = content.replace("__", "{").replace("_-_","}")

with open('nextclade.aligned.fixed.fasta', 'w') as new_file:
    new_file.write(modified_content)

print("Content has been modified and saved to 'nextclade.aligned.fixed.fasta'")





