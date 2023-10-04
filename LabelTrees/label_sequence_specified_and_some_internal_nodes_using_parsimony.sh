#!/bin/bash

# Set the input and output file paths
tree_file="data/unlabeled.nwk"
list_file="data/list.txt"
output_file="data/output/parsimony.nwk"
internal_nodes="Parsimony"

# Run the hyphy command
hyphy label-tree.bf \
--tree "$tree_file" \
--list "$list_file" \
--output "$output_file" \
--internal-nodes "$internal_nodes"

