#!/bin/bash

# Function to process FASTA file and remove newlines between header and sequence lines
process_fasta() {
    local input_file="$1"
    local output_file="$2"

    # Sed script to process the FASTA file
    sed '/^>/{N;s/\n//}' "$input_file" > "$output_file"
}

# Main script starts here
if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <input.fasta> <output.fasta>"
    exit 1
fi

input_file="$1"
output_file="$2"

# Check if input file exists
if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found!"
    exit 1
fi

# Process the FASTA file
process_fasta "$input_file" "$output_file"

echo "FASTA file processed. Results saved to '$output_file'."