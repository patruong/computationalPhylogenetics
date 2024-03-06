def modify_sequence_names(fasta_file, fasta_output, suffix = ""):
    modified_lines = []

    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                sequence_name = line.split(' ')[0]  # Extract the desired sequence name
                modified_lines.append(sequence_name + suffix + '\n')  # Append the modified sequence name to the list
            else:
                modified_lines.append(line)  # Append non-header lines as they are

    with open(fasta_output, 'w') as f:
        f.writelines(modified_lines)  # Write the modified lines back to the file


import os
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20240215_nextclade_pipeline/data")
deer_file = "delta_deer.fasta"  
deer_output_file = "delta_deer_renamed.fasta"
modify_sequence_names(deer_file, deer_output_file)

#modify_sequence_names(deer_file, deer_output_file, "__" + "DEER" + "_-_")
omicron_file = "delta_human.fasta"  
omicron_output_file = "delta_human_renamed.fasta"
modify_sequence_names(omicron_file, omicron_output_file)

#modify_sequence_names(omicron_file, omicron_output_file, "__" + "OMICRON" + "_-_")

# Open the input files and read their contents
with open(deer_output_file, "r") as omicron_file:
    omicron_content = omicron_file.read()

with open(omicron_output_file, "r") as deer_covid_file:
    deer_covid_content = deer_covid_file.read()

# Combine the contents
combined_content = omicron_content + deer_covid_content

# Write the combined content to a new file
with open("combined.fasta", "w") as combined_file:
    combined_file.write(combined_content)

####
# combined difFUBAR input


import os
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20240215_nextclade_pipeline/data")
deer_file = "deer_covid.fasta"  
deer_output_file = "deer_covid_difFUBAR.fasta"
modify_sequence_names(deer_file, deer_output_file, "{" + "DEER" + "}")
omicron_file = "omicron.fasta"  
omicron_output_file = "omicron_difFUBAR.fasta"
modify_sequence_names(omicron_file, omicron_output_file, "{" + "OMICRON" + "}")

# Open the input files and read their contents
with open(deer_output_file, "r") as omicron_file:
    omicron_content = omicron_file.read()

with open(omicron_output_file, "r") as deer_covid_file:
    deer_covid_content = deer_covid_file.read()

# Combine the contents
combined_content = omicron_content + deer_covid_content

# Write the combined content to a new file
with open("combined_difFUBAR_input.fasta", "w") as combined_file:
    combined_file.write(combined_content)



