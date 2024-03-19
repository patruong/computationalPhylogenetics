#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 input.fasta"
    exit 1
fi

# Assign input argument to a variable
input_file=$1

# Run IQ-TREE 2
iqtree2 -s "$input_file" -m MFP -bb 1000

# Optional: Move output files to a desired directory
# iqtree2_output=$(basename "$input_file" .fasta)
# mv "$iqtree2_output".* /path/to/output/directory