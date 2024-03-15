#!/bin/bash
output_dir=output/hivRT_branchlength_1
input_fasta=data/HIV_RT.fasta
input_tree=data/HIV_RT_branchlength_1.nwk

#mkdir -p "$output_dir/treesurgery_and_parallel"
mkdir -p "$output_dir/contrastfel"

#(time julia run_hivRT.jl) 2>&1 | tee "$output_dir/treesurgery_and_parallel/difFUBAR_output.txt" "$output_dir/treesurgery_and_parallel/time_output.txt"
(time hyphy contrast-fel --alignment $input_fasta  --tree $input_tree --branch-set NAIVE --branch-set TREATED --p-value 1.00 --q-value 1.00 --cpu 20 --output $output_dir/contrastfel/contrastfel.FEL.json) 2>&1 | tee "$output_dir/contrastfel/constrastfel_output.txt" "$output_dir/contrastfel/time_output.txt"
