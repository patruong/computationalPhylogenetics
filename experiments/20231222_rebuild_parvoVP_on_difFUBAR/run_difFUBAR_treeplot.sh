#!/bin/bash

output_dir="output"

mkdir -p "$output_dir/treeplot"

(time julia run_difFUBAR_treeplot.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned.nwk" --output "$output_dir/treeplot/") 2>&1 | tee "$output_dir/treeplot/difFUBAR_output.txt" "$output_dir/treeplot/time_output.txt"



