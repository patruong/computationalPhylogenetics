import argparse
from pytrimal import Alignment, AutomaticTrimmer

def trim_alignment(input_clw, output_fasta, method="strictplus", remove_duplicates = True):
    ali = Alignment.load(input_clw)
    trimmer = AutomaticTrimmer(method=method)
    trimmed = trimmer.trim(ali)
    if remove_duplicates == False:
        with open(output_fasta, "w") as fasta_file:
            for name, seq in zip(trimmed.names, trimmed.sequences):
                name = name.decode()[:-1] + "]"
                seq = seq.upper()
                fasta_file.write(f">{name}\n{seq}\n")
    elif remove_duplicates == True:
        unique_seqs = set()
        with open(output_fasta, "w") as fasta_file:
            for name, seq in zip(trimmed.names, trimmed.sequences):
                name = name.decode()[:-1] + "]"
                seq = seq.upper()
                if seq not in unique_seqs:
                    unique_seqs.add(seq)
                    fasta_file.write(f">{name}\n{seq}\n")
                else:
                    print(f"Duplicate sequence {name} removed.")

def main():
    parser = argparse.ArgumentParser(description="Trim alignments using pytrimal and output in FASTA format")
    parser.add_argument("input_clw", help="Input CLUSTAL alignment file")
    parser.add_argument("output_fasta", help="Output FASTA file")
    parser.add_argument("--method", default="strictplus", help="Trimming method to use (default: strictplus)")
    parser.add_argument("--remove_duplicates", action="store_true", help="Remove duplicate sequences (default: True)")

    args = parser.parse_args()
    
    trim_alignment(args.input_clw, args.output_fasta, args.method, args.remove_duplicates)

if __name__ == "__main__":
    main()