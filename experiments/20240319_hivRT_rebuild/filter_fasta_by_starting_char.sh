#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 input.fasta output.fasta starting_char"
    exit 1
fi

# Assign input arguments to variables
input_file=$1
output_file=$2
starting_char=$3

# Run Python script to filter sequences by starting character
python filter_fasta_by_starting_char.py "$input_file" "$output_file" "$starting_char"
