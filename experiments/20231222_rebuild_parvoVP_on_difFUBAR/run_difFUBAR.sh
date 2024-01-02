#!/bin/bash

output_dir="output"

mkdir -p "$output_dir/baseline"
mkdir -p "$output_dir/max"
mkdir -p "$output_dir/patrick"
mkdir -p "$output_dir/patrick_max"
mkdir -p "$output_dir/patrick_max_child"
mkdir -p "$output_dir/final"
mkdir -p "$output_dir/treeplot"
mkdir -p "$output_dir/contrastfel"

# make tree plot
(time julia run_difFUBAR_treeplot.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped_rerooted_reassigned_groups.nwk" --output "$output_dir/treeplot/") 2>&1 | tee "$output_dir/treeplot/difFUBAR_output.txt" "$output_dir/treeplot/time_output.txt"

# run difFUBAR
(time julia run_difFUBAR_baseline_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped_rerooted_reassigned_groups.nwk" --output "$output_dir/baseline/") 2>&1 | tee "$output_dir/baseline/difFUBAR_output.txt" "$output_dir/baseline/time_output.txt"
(time julia run_difFUBAR_prune_max_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped_rerooted_reassigned_groups.nwk" --output "$output_dir/max/") 2>&1 | tee "$output_dir/max/difFUBAR_output.txt" "$output_dir/max/time_output.txt"
(time julia run_difFUBAR_prune_patrick_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped_rerooted_reassigned_groups.nwk" --output "$output_dir/patrick/") 2>&1 | tee "$output_dir/patrick/difFUBAR_output.txt" "$output_dir/patrick/time_output.txt"
(time julia run_difFUBAR_prune_patrick_max_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped_rerooted_reassigned_groups.nwk" --output "$output_dir/patrick_max/") 2>&1 | tee "$output_dir/patrick_max/difFUBAR_output.txt" "$output_dir/patrick_max_child/time_output.txt"
(time julia run_difFUBAR_prune_patrick_max_child_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped_rerooted_reassigned_groups.nwk" --output "$output_dir/patrick_max_child/") 2>&1 | tee "$output_dir/patrick_max_child/difFUBAR_output.txt" "$output_dir/patrick_max_child/time_output.txt"
(time julia run_difFUBAR_prune_final_parvoVP.jl --fasta_file "parvo_panleu_trans.fasta" --nexus_file "parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped_rerooted_reassigned_groups.nwk" --output "$output_dir/final/") 2>&1 | tee "$output_dir/final/difFUBAR_output.txt" "$output_dir/final/time_output.txt"

# run constrast-fel
(time hyphy contrast-fel --alignment parvo_panleu_trans.fasta  --tree parvo_panleu_trans_tagged_fixed_aligned_internalNode_grouped_rerooted_reassigned_groups.nwk --branch-set G1 --branch-set G2 --p-value 1.00 --q-value 1.00 --cpu 20) 2>&1 | tee "$output_dir/contrastfel/constrastfel_output.txt" "$output_dir/contrastfel/time_output.txt"

