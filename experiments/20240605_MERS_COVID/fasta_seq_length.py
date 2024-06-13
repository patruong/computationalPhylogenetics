import argparse

def read_fasta(file_path):
    sequences = {}
    with open(file_path, 'r') as fasta_file:
        current_sequence_name = None
        current_sequence = ''
        for line in fasta_file:
            line = line.strip()
            if line.startswith('>'):
                if current_sequence_name:
                    sequences[current_sequence_name] = current_sequence
                current_sequence_name = line[1:]
                current_sequence = ''
            else:
                current_sequence += line
        if current_sequence_name:
            sequences[current_sequence_name] = current_sequence
    return sequences

def print_sequence_lengths(file_path):
    sequences = read_fasta(file_path)
    for name, seq in sequences.items():
        print(f"Sequence: {name}, Length: {len(seq)}")

def main():
    parser = argparse.ArgumentParser(description="Check sequence length of all sequences in a FASTA file")
    parser.add_argument("fasta_file", help="Path to the FASTA file")
    args = parser.parse_args()
    
    print_sequence_lengths(args.fasta_file)

if __name__ == "__main__":
    main()