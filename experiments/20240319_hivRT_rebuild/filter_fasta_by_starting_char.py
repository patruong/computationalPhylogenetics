import sys

def read_fasta(file_path):
    """
    Read FASTA file and return a dictionary with sequence names as keys and sequences as values.
    """
    sequences = {}
    current_seq_name = None
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_seq_name = line[1:]
                sequences[current_seq_name] = ''
            else:
                sequences[current_seq_name] += line
    return sequences

def write_fasta(sequences, output_file):
    """
    Write sequences to a new FASTA file.
    """
    with open(output_file, 'w') as file:
        for seq_name, sequence in sequences.items():
            file.write('>' + seq_name + '\n')
            file.write(sequence + '\n')

def filter_sequences_by_starting_char(sequences, starting_char):
    """
    Filter sequences by sequence names starting with the given character.
    """
    filtered_sequences = {}
    for seq_name, sequence in sequences.items():
        if seq_name.startswith(starting_char):
            filtered_sequences[seq_name] = sequence
    return filtered_sequences

def main(input_file, output_file, starting_char):
    # Read sequences from input file
    sequences = read_fasta(input_file)
    
    # Filter sequences by starting character
    filtered_sequences = filter_sequences_by_starting_char(sequences, starting_char)
    
    # Write filtered sequences to output file
    write_fasta(filtered_sequences, output_file)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py input.fasta output.fasta starting_char")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    starting_char = sys.argv[3]

    main(input_file, output_file, starting_char)