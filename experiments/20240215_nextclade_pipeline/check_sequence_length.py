import os
os.chdir("/home/patrick/git/computationalPhylogenetics/experiments/20240215_nextclade_pipeline/results/nextclade_output")

def read_fasta(file_path):
    sequences = {}
    current_sequence_name = None
    current_sequence = ''

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                # New sequence found
                if current_sequence_name:
                    sequences[current_sequence_name] = current_sequence
                current_sequence_name = line[1:]
                current_sequence = ''
            else:
                current_sequence += line
        # Adding the last sequence after reaching end of file
        if current_sequence_name:
            sequences[current_sequence_name] = current_sequence

    return sequences

# Example usage:
fasta_file_path = "nextclade.cds_translation.S.fasta"  # Provide the path to your FASTA file
sequences = read_fasta(fasta_file_path)


# Print sequences
for name, sequence in sequences.items():
    print(f"Sequence Name: {name}")
    print(f"Sequence: {len(sequence)}")

1274/3
29994/3
29904/3


# CHATGPT

"""
does the number 29903 have any specific significance to covid sequence length? 

ChatGPT
Yes, the number 29903 is significant in the context of the SARS-CoV-2 virus, which is responsible for COVID-19. This number represents the length of the reference genome of the SARS-CoV-2 virus, which is approximately 29,903 base pairs long. The reference genome is the standard sequence used as a basis for comparing and analyzing genetic variations in different strains or isolates of the virus.

The SARS-CoV-2 genome consists of a single-stranded RNA molecule containing several genes encoding various viral proteins, including structural proteins (such as spike, envelope, membrane, and nucleocapsid proteins) and non-structural proteins. Understanding the sequence and organization of the viral genome is crucial for studying the virus's biology, evolution, and pathogenesis, as well as for developing diagnostic tests, vaccines, and treatments.

Researchers have sequenced thousands of SARS-CoV-2 genomes from different geographic locations and time points since the emergence of the virus in late 2019. By comparing these sequences to the reference genome, scientists can identify mutations, deletions, insertions, and other genetic variations that may affect the virus's properties, such as transmissibility, virulence, and susceptibility to antiviral drugs or vaccines.

Therefore, the length of the SARS-CoV-2 genome, represented by the number 29903, serves as a standard reference for genomic studies and plays a crucial role in understanding the molecular biology and evolution of the virus, as well as in the development of strategies for controlling the COVID-19 pandemic.





"""