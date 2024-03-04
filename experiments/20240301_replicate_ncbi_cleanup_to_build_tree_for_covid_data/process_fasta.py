# For work code ./experiments/data/20231221_feline_canine_parvoVP_build/selected_sequences_based_on_length_from_full_fasta/process_fasta.py
import pandas as pd
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Process fasta files with sequence data.")
    parser.add_argument("fasta_file", help="Input fasta file")
    parser.add_argument("--group1", default="deer", help="Name of group 1")
    parser.add_argument("--group2", default="human", help="Name of group 2")
    parser.add_argument("--group1_n", type=int, default=50, help="Number of sequences to sample from group 1")
    parser.add_argument("--group2_n", type=int, default=50, help="Number of sequences to sample from group 2")
    parser.add_argument("--all_sequences", action="store_true", help="Use all sequences instead of sampling")
    parser.add_argument("--output", default="output.fasta", help="Output fasta file")
    parser.add_argument("--output_tagged", default="output_tagged.fasta", help="Output tagged fasta file")
    return parser.parse_args()

def process_sequences(fasta_file, group1, group2, group1_n, group2_n, all_sequences, output, output_tagged):
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
    print("Setting uniform sequence length for difFUBAR proportionmap function.")
    df = pd.DataFrame([seqnames, seqs, seqlens], index = ["seqname", "sequence","length"]).T
    
    
    #df.length.value_counts()
    #df_subset = df[df.length == 1755] 


    #print("Removing sequence names containing UNVERIFIED.")
    #df_subset = df_subset[~df_subset["seqname"].str.contains("UNVERIFIED")]
    #sequence_counts = df_subset.sequence.value_counts()
    #unique_sequences = sequence_counts[sequence_counts==1].index
    #df_subset = df_subset[df_subset.sequence.isin(unique_sequences)]

    # remove odd characters found by fasttree
    # DOES THIS WORK IF I REMOVE THIS
    #print("Removing sequences containing ambigous nucleotides K|M|N|R|S|W|Y.")
    #df = df[~df.sequence.str.contains("K|M|N|R|S|W|Y")]

    # sequence name exclusion list
    # this sequence is very far away from the other
    # print(f"{len(df_subset[df_subset.seqname.str.contains('MN2709')])} sequences with id starting MN2709 omitted because they are very far away from other sequences.")
    # print(f"{len(df_subset[df_subset.seqname.str.contains('LY660937')])} sequences with id starting MN2709 omitted because they are very far away from other sequences.")

    #df_subset = df_subset[~df_subset.seqname.str.contains("MN270957")] 
    #df_subset = df_subset[~df_subset.seqname.str.contains("MN2709")] 
    #df_subset = df_subset[~df_subset.seqname.str.contains("MN270944")]
    #df_subset = df_subset[~df_subset.seqname.str.contains("MN270954")] 
    #df_subset = df_subset[~df_subset.seqname.str.contains("MN270954")] 
    #df_subset = df_subset[~df_subset.seqname.str.contains("LY660937")] 
    

    # sequence must not contain stop-codons or contrast-FEL will get Error: The input alignment must have the number of sites that is divisible by 3 and must not contain stop codons in call to assert(fel.codon_filter.sites*3==fel.codon_data.sites, error_msg);
    # Remove stop codons
    # print("Removing stop codons for constrast-FEL run.")
    # df_subset.sequence = df_subset.sequence.map(lambda x:x[:-3])
    # df_subset #1372
    # def is_stop_codon(sequence):
    #     stop_codons = ['TAA', 'TAG', 'TGA']
    #     return sequence[-3:].upper() not in stop_codons

    # df_subset = df_subset[df_subset['sequence'].apply(is_stop_codon)] # this should not remove anything, because that means we have two stop codons in a row

    #print(f"Sampling {group1_n} from {group1} and {group2_n} from {group2}")
    #
    #print(len(df[df["seqname"].str.contains(group2)]))
    #print(len(df[df["seqname"].str.contains(group1)]))
    #print()
    ## subset size 
    #if all_sequences:
    #    df_group_1 = df[df["seqname"].str.contains(group1)]
    #    df_group_2 = df[df["seqname"].str.contains(group2)]
    #else:
    #    df_group_1 = df[df["seqname"].str.contains(group1)].sample(group1_n)
    #    df_group_2 = df[df["seqname"].str.contains(group2)].sample(group2_n)
    #df_res = pd.concat([df_group_1, df_group_2])

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
        if group1.upper() in seqname:
            specie = group1
        elif group2.upper() in seqname:
            specie = group2
        
        seqname = seqname.translate(str.maketrans("| ,()/.-;:", "__________")).replace("__", "_").upper()[:50] + specie.upper()
        transformed_seqnames.append(seqname)

    # for building tree with mafft and fasttree, #requires fix_tags.py to transform äö to {}
    #with open('parvo_panleu_trans_tagged.fasta', 'w') as fasta_file:
    print("Generating output fasta.")
    with open(output_tagged, 'w') as fasta_file:
        for seqname, seq in zip(transformed_seqnames, seqs):
            if group1.upper() in seqname:
                seqname += "-_-" + group1.upper() + "._." # Can change to CANINE 
            elif group2.upper() in seqname:
                seqname += "-_-" + group2.upper() + "._." # Can change to FELINE
            fasta_file.write(f'>{seqname}\n{seq}\n')
    # for difFUBAR
    #with open('parvo_panleu_trans.fasta', 'w') as fasta_file:
    with open(output, 'w') as fasta_file:
        for seqname, seq in zip(transformed_seqnames, seqs):
            fasta_file.write(f'>{seqname}\n{seq}\n')
    print(f"Generated {len(seqs)} sequences in {output}")

if __name__ == "__main__":
    args = parse_args()
    process_sequences(args.fasta_file, args.group1, args.group2,
                       args.group1_n, args.group2_n, args.all_sequences, args.output, args.output_tagged)

