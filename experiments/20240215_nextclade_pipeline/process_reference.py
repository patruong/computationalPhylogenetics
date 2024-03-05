# For work code ./experiments/data/20231221_feline_canine_parvoVP_build/selected_sequences_based_on_length_from_full_fasta/process_fasta.py
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Process fasta files with sequence data.")
    parser.add_argument("fasta_file", help="Input fasta file")
    parser.add_argument("--output", default="output.fasta", help="Output fasta file")
    return parser.parse_args()

def process_sequences(fasta_file, output):
    seqs = []       # List to store sequences
    seqnames = []   # List to store sequence names
    seqlens = []    # List to store sequence lengths

    with open(fasta_file, 'r') as file:
        sequence = ""
        for line in file:
            line = line.strip()  # Remove leading and trailing whitespaces
            if line.startswith('>'):
                # This line is a header line, indicating the start of a new sequence
                header = line[1:]  # Exclude the '>' character
                seqnames.append(header)
                if sequence:
                    seqs.append(sequence)  # Append the previous sequence if not empty
                    seqlens.append(len(sequence))

                    sequence = ""  # Reset sequence
            else:
                # This line is part of the sequence
                sequence += line

        # Append the last sequence after the loop
        if sequence:
            seqs.append(sequence)
            seqlens.append(len(sequence))

    # Take sequence of equal length
    #print("Setting uniform sequence length for difFUBAR proportionmap function.")
    df = pd.DataFrame([seqnames, seqs, seqlens], index = ["seqname", "sequence","length"]).T
    
    df_res = df

    print('Transforming special characters "| ,()/.-;:" --> "__________"')
    print('Transforming double underscore "__" --> "_"')
    print("Setting upper case for sequence name.")
    seqnames = list(df_res.seqname)
    seqs = list(df_res.sequence)
    transformed_seqnames = []
    for seqname in seqnames:
        seqname = seqname.upper()
        specie = ""        
        seqname = seqname.translate(str.maketrans("| ,()/.-;:", "__________")).replace("__", "_").upper()[:50] + specie.upper()
        transformed_seqnames.append(seqname)

    print("Generating output fasta.")
    with open(output, 'w') as fasta_file:
        for seqname, seq in zip(transformed_seqnames, seqs):
            fasta_file.write(f'>{seqname}\n{seq}\n')
    print(f"Generated {len(seqs)} sequences in {output}")

if __name__ == "__main__":
    args = parse_args()
    process_sequences(args.fasta_file, args.output)

