import argparse

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

def main():
    parser = argparse.ArgumentParser(description='Add "-" at the end of each sequence in a FASTA file')
    parser.add_argument('input_file', type=str, help='Input FASTA file path')
    parser.add_argument('output_file', type=str, help='Output FASTA file path with modified sequences')
    args = parser.parse_args()

    sequences = read_fasta(args.input_file)
    write_fasta(sequences, args.output_file)

if __name__ == "__main__":
    main()