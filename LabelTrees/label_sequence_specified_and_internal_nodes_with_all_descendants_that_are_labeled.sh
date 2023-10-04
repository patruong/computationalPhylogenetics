#!/bin/bash

# Set the input and output file paths
tree_file="data/unlabeled.nwk"
list_file="data/list.txt"
output_file="data/output/all.nwk"

# Run the hyphy command
hyphy label-tree.bf \
--tree "$tree_file" \
--list "$list_file" \
--output "$output_file"

