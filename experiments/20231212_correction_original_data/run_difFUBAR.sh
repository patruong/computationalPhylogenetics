#!/bin/bash

mkdir -p output/baseline
mkdir -p output/max
mkdir -p output/patrick
mkdir -p output/patrick_max
mkdir -p output/patrick_max_child
mkdir -p output/final

(time julia run_difFUBAR_baseline_parvoVP.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nwk" --output "output/baseline/") 2>&1 | tee output/baseline/difFUBAR_output.txt output/baseline/time_output.txt
(time julia run_difFUBAR_prune_max_parvoVP.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nwk" --output "output/max/") 2>&1 | tee output/max/difFUBAR_output.txt output/max/time_output.txt
(time julia run_difFUBAR_prune_patrick_parvoVP.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nwk" --output "output/patrick/") 2>&1 | tee output/patrick/difFUBAR_output.txt output/patrick/time_output.txt
(time julia run_difFUBAR_prune_patrick_max_parvoVP.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nwk" --output "output/patrick_max/") 2>&1 | tee output/patrick_max/difFUBAR_output.txt output/patrick_max_child/time_output.txt
(time julia run_difFUBAR_prune_patrick_max_child_parvoVP.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nwk" --output "output/patrick_max_child/") 2>&1 | tee output/patrick_max_child/difFUBAR_output.txt output/patrick_max_child/time_output.txt
(time julia run_difFUBAR_prune_final_parvoVP.jl --fasta_file "ParvoVP_removed_duplicates.fasta" --nexus_file "ParvoVP_removed_duplicates.nwk" --output "output/final/") 2>&1 | tee output/final/difFUBAR_output.txt output/final/time_output.txt



