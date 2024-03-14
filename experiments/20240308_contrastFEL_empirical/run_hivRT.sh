#!/bin/bash
output_dir=output/hivRT
#input_fasta=data/rbcl.nex
#input_tree=data/rubisco_C3_vs_C4.nwk

mkdir -p "$output_dir/treesurgery_and_parallel"
#mkdir -p "$output_dir/contrastfel"

(time julia run_hivRT.jl) 2>&1 | tee "$output_dir/treesurgery_and_parallel/difFUBAR_output.txt" "$output_dir/treesurgery_and_parallel/time_output.txt"
#(time hyphy contrast-fel --alignment $input_fasta  --tree $input_tree --branch-set C3 --branch-set C4 --p-value 1.00 --q-value 1.00 --cpu 20 --output $output_dir/contrastfel/contrastfel.FEL.json) 2>&1 | tee "$output_dir/contrastfel/constrastfel_output.txt" "$output_dir/contrastfel/time_output.txt"
