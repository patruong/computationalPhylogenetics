#!/bin/bash

output_dir="output"

mkdir -p "$output_dir/baseline"
mkdir -p "$output_dir/max"
mkdir -p "$output_dir/patrick"
mkdir -p "$output_dir/patrick_max"
mkdir -p "$output_dir/patrick_max_child"
mkdir -p "$output_dir/final"

(time julia run_difFUBAR_baseline.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nex" --output "$output_dir/baseline/") 2>&1 | tee "$output_dir/baseline/difFUBAR_output.txt" "$output_dir/baseline/time_output.txt"
(time julia run_difFUBAR_prune_max.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nex" --output "$output_dir/max/") 2>&1 | tee "$output_dir/max/difFUBAR_output.txt" "$output_dir/max/time_output.txt"
(time julia run_difFUBAR_prune_patrick.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nex" --output "$output_dir/patrick/") 2>&1 | tee "$output_dir/patrick/difFUBAR_output.txt" "$output_dir/patrick/time_output.txt"
(time julia run_difFUBAR_prune_patrick_max.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nex" --output "$output_dir/patrick_max/") 2>&1 | tee "$output_dir/patrick_max/difFUBAR_output.txt" "$output_dir/patrick_max_child/time_output.txt"
(time julia run_difFUBAR_prune_patrick_max_child.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nex" --output "$output_dir/patrick_max_child/") 2>&1 | tee "$output_dir/patrick_max_child/difFUBAR_output.txt" "$output_dir/patrick_max_child/time_output.txt"
(time julia run_difFUBAR_prune_final.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nex" --output "$output_dir/final/") 2>&1 | tee "$output_dir/final/difFUBAR_output.txt" "$output_dir/final/time_output.txt"



