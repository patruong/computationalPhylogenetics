#!/bin/bash

#output_dir="output"
#input_dir="/home/patrick/git/computationalPhylogenetics/contrastFEL_data/omnibus/"
#input_fasta="${input_dir}sims.421.settings.replicate.2"
#input_tree="${input_dir}sims.421.nwk"


input_dir="$1"
input_fasta="$2"
input_tree="$3"
output_dir="$4"

mkdir -p "$output_dir/baseline"
mkdir -p "$output_dir/max"
mkdir -p "$output_dir/patrick"
mkdir -p "$output_dir/patrick_max"
mkdir -p "$output_dir/patrick_max_child"
mkdir -p "$output_dir/final"
mkdir -p "$output_dir/treeplot"
mkdir -p "$output_dir/contrastfel"



# make tree plot
(time julia run_difFUBAR_treeplot.jl --fasta_file $input_fasta --nexus_file $input_tree --output "$output_dir/treeplot/") 2>&1 | tee "$output_dir/treeplot/difFUBAR_output.txt" "$output_dir/treeplot/time_output.txt"

# run difFUBAR
(time julia run_difFUBAR_baseline_parvoVP.jl --fasta_file $input_fasta --nexus_file $input_tree --output "$output_dir/baseline/") 2>&1 | tee "$output_dir/baseline/difFUBAR_output.txt" "$output_dir/baseline/time_output.txt"
(time julia run_difFUBAR_prune_max_parvoVP.jl --fasta_file $input_fasta --nexus_file $input_tree --output "$output_dir/max/") 2>&1 | tee "$output_dir/max/difFUBAR_output.txt" "$output_dir/max/time_output.txt"
(time julia run_difFUBAR_prune_patrick_parvoVP.jl --fasta_file $input_fasta --nexus_file $input_tree --output "$output_dir/patrick/") 2>&1 | tee "$output_dir/patrick/difFUBAR_output.txt" "$output_dir/patrick/time_output.txt"
(time julia run_difFUBAR_prune_patrick_max_parvoVP.jl --fasta_file $input_fasta --nexus_file $input_tree --output "$output_dir/patrick_max/") 2>&1 | tee "$output_dir/patrick_max/difFUBAR_output.txt" "$output_dir/patrick_max/time_output.txt"
(time julia run_difFUBAR_prune_patrick_max_child_parvoVP.jl --fasta_file $input_fasta --nexus_file $input_tree --output "$output_dir/patrick_max_child/") 2>&1 | tee "$output_dir/patrick_max_child/difFUBAR_output.txt" "$output_dir/patrick_max_child/time_output.txt"
(time julia run_difFUBAR_prune_final_parvoVP.jl --fasta_file $input_fasta --nexus_file $input_tree --output "$output_dir/final/") 2>&1 | tee "$output_dir/final/difFUBAR_output.txt" "$output_dir/final/time_output.txt"
#
## run constrast-fel
(time hyphy contrast-fel --alignment $input_fasta  --tree $input_tree --branch-set TEST --branch-set REFERENCE --p-value 1.00 --q-value 1.00 --cpu 20 --output $output_dir/contrastfel/contrastfel.FEL.json) 2>&1 | tee "$output_dir/contrastfel/constrastfel_output.txt" "$output_dir/contrastfel/time_output.txt"

# if above command does not print out contrastfel.FEL.json
#
##hyphy contrast-fel --alignment @input_fasta  --tree $input_tree --branch-set G1 --branch-set G2 --p-value 1.00 --q-value 1.00 --cpu 20 --output $output_dir/contrastfel/contrastfel.FEL.json
#