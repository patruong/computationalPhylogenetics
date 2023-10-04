#!/bin/bash

# Set the input and output file paths
tree_file="data/unlabeled.nwk"
regexp_pattern="sp$"
output_file="data/output/regexp-none.nwk"
internal_nodes="None"

# Run the hyphy command
hyphy label-tree.bf \
--tree "$tree_file" \
--regexp "$regexp_pattern" \
--output "$output_file" \
--internal-nodes "$internal_nodes"

