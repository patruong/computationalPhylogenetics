#!/bin/bash

# Check if there are exactly 3 arguments provided
if [ $# -ne 3 ]; then
  echo "Usage: $0 <input_fasta> <tag> <output_fasta>"
  exit 1
fi

# Assign arguments to variables
input_file="$1"
tag="$2"
output_file="$3"

# Perform sed command with input, tag, and output arguments
sed -r "/^>/s/([^ ]+)(.*)/\1${tag}\2/" "$input_file" > "$output_file"

echo "Successfully added tag '${tag}' to FASTA headers in $output_file"