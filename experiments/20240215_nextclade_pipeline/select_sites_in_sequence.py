import argparse
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

def write_fasta(sequences, output_file):
    with open(output_file, 'w') as file:
        for name, sequence in sequences.items():
            file.write(f'>{name}\n')
            file.write(f'{sequence}-\n')