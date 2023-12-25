import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Process fasta files with sequence data.")
    parser.add_argument("fasta_file", help="Input fasta file")
    parser.add_argument("--tagged", action="store_true", help="Add tags to sequence names")
    parser.add_argument("--group1", default="Canine", help="Name of group 1")
    parser.add_argument("--group2", default="Feline", help="Name of group 2")
    parser.add_argument("--group1_n", type=int, default=324, help="Number of sequences to sample from group 1")
    parser.add_argument("--group2_n", type=int, default=47, help="Number of sequences to sample from group 2")
    parser.add_argument("--all_sequences", action="store_true", help="Use all sequences instead of sampling")
    parser.add_argument("--output", default="output.fasta", help="Output fasta file")
    return parser.parse_args()

def process_sequences(fasta_file, tagged, group1, group2, group1_n, group2_n, all_sequences, output):
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
    df = pd.DataFrame([seqnames, seqs, seqlens], index = ["seqname", "sequence","length"]).T
    df.length.value_counts()
    df_subset = df[df.length == 1755] 
    df_subset = df_subset[~df_subset["seqname"].str.contains("UNVERIFIED")]
    sequence_counts = df_subset.sequence.value_counts()
    unique_sequences = sequence_counts[sequence_counts==1].index
    df_subset = df_subset[df_subset.sequence.isin(unique_sequences)]

    # remove odd characters found by fasttree
    df_subset = df_subset[~df_subset.sequence.str.contains("K|M|N|R|S|W|Y")]

    # subset size 
    if all_sequences:
        df_group_1 = df_subset[df_subset["seqname"].str.contains(group1)]
        df_group_2 = df_subset[df_subset["seqname"].str.contains(group2)]
    else:
        df_group_1 = df_subset[df_subset["seqname"].str.contains(group1)].sample(group1_n)
        df_group_2 = df_subset[df_subset["seqname"].str.contains(group2)].sample(group2_n)
    df_res = pd.concat([df_group_1, df_group_2])

    seqnames = list(df_res.seqname)
    seqs = list(df_res.sequence)
    transformed_seqnames = []
    for seqname in seqnames:
        seqname = seqname.translate(str.maketrans("| ,()/.-;:", "__________")).replace("__", "_").upper()
        transformed_seqnames.append(seqname)


    if tagged:
        # for building tree with mafft and fasttree, #requires fix_tags.py to transform äö to {}
        #with open('parvo_panleu_trans_tagged.fasta', 'w') as fasta_file:
        with open(output, 'w') as fasta_file:
            for seqname, seq in zip(transformed_seqnames, seqs):
                if "CANINE" in seqname:
                    seqname += "ÖG2Ä" # Can change to CANINE 
                elif "FELINE" in seqname:
                    seqname += "ÖG1Ä" # Can change to FELINE
                fasta_file.write(f'>{seqname}\n{seq}\n')
    else:
        # for difFUBAR
        #with open('parvo_panleu_trans.fasta', 'w') as fasta_file:
        with open(output, 'w') as fasta_file:
            for seqname, seq in zip(transformed_seqnames, seqs):
                fasta_file.write(f'>{seqname}\n{seq}\n')
    print(f"Generated {len(seqs)} sequences in {output}")

if __name__ == "__main__":
    args = parse_args()
    process_sequences(args.fasta_file, args.tagged, args.group1, args.group2,
                       args.group1_n, args.group2_n, args.all_sequences, args.output)

