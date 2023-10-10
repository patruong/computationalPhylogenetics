#!/bin/bash

# Define your variables
POS_THRESH=(0.9 0.95)
REPLICATE=()
SCRIPT="run_difFUBAR.jl"
FASTA_FILE="../../contrastFEL_data/omnibus-multi/datasets/omnibus-multi/sims.24.settings.replicate.1"
TREE_FILE="../../contrastFEL_data/omnibus-multi/datasets/omnibus-multi/sims.24.nwk"
VERBOSITY=1
EXPORTS="true"
ITERS=2500

# Loop over each value in POS_THRESH
for iPos_thresh in "${POS_THRESH[@]}"; do
    OUTPUT_DIR="output/pos_thresh_$iPos_thresh/null"
    
    # Run the Julia command
    julia "$SCRIPT" -f "$FASTA_FILE" -t "$TREE_FILE" -o "$OUTPUT_DIR" \
        --verbosity "$VERBOSITY" --exports "$EXPORTS" --iters "$ITERS" \
        --pos_thresh "$iPos_thresh"
done
