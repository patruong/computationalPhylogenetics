#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def prune_sequences(input_file, output_file):
    records = []
    with open(input_file, "r") as infile:
        for record in SeqIO.parse(infile, "fasta"):
            sequence = record.seq
            remainder = len(sequence) % 3
            if remainder != 0:
                # Adjust sequence length to be divisible by 3
                sequence = sequence[:-remainder]
            # Create a new SeqRecord with the pruned sequence
            pruned_record = SeqRecord(Seq(str(sequence)), id=record.id, description=record.description)
            records.append(pruned_record)
    
    # Write the pruned sequences to the output file
    with open(output_file, "w") as outfile:
        SeqIO.write(records, outfile, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prune sequences in a FASTA file to be divisible by 3.")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output FASTA file with pruned sequences")
    
    args = parser.parse_args()
    
    input_file = args.input_file
    output_file = args.output_file
    
    prune_sequences(input_file, output_file)
    
    print(f"FASTA file processed. Pruned sequences saved to '{output_file}'.")
