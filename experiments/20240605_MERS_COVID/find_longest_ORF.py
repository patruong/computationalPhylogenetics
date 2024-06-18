import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def find_longest_orf(sequence):
    """Find the longest ORF in all six reading frames."""
    longest_orf = ""
    longest_length = 0
    
    # Forward frames
    for frame in range(3):
        for i in range(frame, len(sequence), 3):
            codon = sequence[i:i+3]
            if codon in ["TAA", "TAG", "TGA"]:  # Stop codons
                if len(longest_orf) > longest_length:
                    longest_length = len(longest_orf)
                    longest_orf = ""
            else:
                longest_orf += codon
    
    # Reverse complement frames
    rev_sequence = sequence.reverse_complement()
    for frame in range(3):
        for i in range(frame, len(rev_sequence), 3):
            codon = rev_sequence[i:i+3]
            if codon in ["TAA", "TAG", "TGA"]:  # Stop codons
                if len(longest_orf) > longest_length:
                    longest_length = len(longest_orf)
                    longest_orf = ""
            else:
                longest_orf += codon
    
    return longest_orf

def process_fasta(input_fasta, output_fasta):
    """Process a FASTA file to find and write the longest ORF for each sequence."""
    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        seq = record.seq
        longest_orf = find_longest_orf(seq)
        new_record = SeqRecord(Seq(longest_orf), id=record.id, description="Longest ORF")
        records.append(new_record)
    
    SeqIO.write(records, output_fasta, "fasta")
    print(f"Processed {len(records)} sequences. Results saved to {output_fasta}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find the longest ORF for each sequence in a FASTA file.")
    parser.add_argument("input_fasta", help="Input FASTA file containing nucleotide sequences.")
    parser.add_argument("output_fasta", help="Output file to save sequences with longest ORFs.")
    
    args = parser.parse_args()
    
    input_fasta = args.input_fasta
    output_fasta = args.output_fasta
    
    process_fasta(input_fasta, output_fasta)
