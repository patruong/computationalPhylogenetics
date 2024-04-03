#!/bin/bash
output_dir=output/hiv_envelope
mkdir -p "$output_dir/treesurgery_and_parallel"
mkdir -p "$output_dir/contrastfel"

input_fasta=data/HIV-env.fasta
input_tree=data/hiv-1_envelope.nwk

(time julia run_hiv_envelope.jl) 2>&1 | tee "$output_dir/treesurgery_and_parallel/difFUBAR_output.txt" "$output_dir/treesurgery_and_parallel/time_output.txt"
(time hyphy contrast-fel --alignment $input_fasta  --tree $input_tree --branch-set HSX --branch-set MSM --p-value 1.00 --q-value 1.00 --cpu 20 --output $output_dir/contrastfel/contrastfel.FEL.json) 2>&1 | tee "$output_dir/contrastfel/constrastfel_output.txt" "$output_dir/contrastfel/time_output.txt"
