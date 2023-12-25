#!/bin/bash

output_dir="output"

mkdir -p "$output_dir/baseline"
mkdir -p "$output_dir/max"
mkdir -p "$output_dir/patrick"
mkdir -p "$output_dir/patrick_max"
mkdir -p "$output_dir/patrick_max_child"
mkdir -p "$output_dir/final"

(time julia run_difFUBAR_baseline_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned.nwk" --output "$output_dir/baseline/") 2>&1 | tee "$output_dir/baseline/difFUBAR_output.txt" "$output_dir/baseline/time_output.txt"
(time julia run_difFUBAR_prune_max_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned.nwk" --output "$output_dir/max/") 2>&1 | tee "$output_dir/max/difFUBAR_output.txt" "$output_dir/max/time_output.txt"
(time julia run_difFUBAR_prune_patrick_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned.nwk" --output "$output_dir/patrick/") 2>&1 | tee "$output_dir/patrick/difFUBAR_output.txt" "$output_dir/patrick/time_output.txt"
(time julia run_difFUBAR_prune_patrick_max_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned.nwk" --output "$output_dir/patrick_max/") 2>&1 | tee "$output_dir/patrick_max/difFUBAR_output.txt" "$output_dir/patrick_max_child/time_output.txt"
(time julia run_difFUBAR_prune_patrick_max_child_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned.nwk" --output "$output_dir/patrick_max_child/") 2>&1 | tee "$output_dir/patrick_max_child/difFUBAR_output.txt" "$output_dir/patrick_max_child/time_output.txt"
(time julia run_difFUBAR_prune_final_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned.nwk" --output "$output_dir/final/") 2>&1 | tee "$output_dir/final/difFUBAR_output.txt" "$output_dir/final/time_output.txt"



