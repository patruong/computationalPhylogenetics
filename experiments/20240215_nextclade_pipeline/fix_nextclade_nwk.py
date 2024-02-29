
import os
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20240215_nextclade_pipeline/results/nextclade_output")

with open('nextclade.nwk', 'r') as file:
    content = file.read()

modified_content = content.replace("__", "{").replace("_-_","}")

with open('nextclade_fixed.nwk', 'w') as new_file:
    # Write the modified content to the new file
    new_file.write(modified_content)

print("Content has been modified and saved to 'nextclade_fixed.nwk'")


