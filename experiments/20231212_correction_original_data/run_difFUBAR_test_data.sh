#!/bin/bash

mkdir -p output/baseline
mkdir -p output/max
mkdir -p output/patrick
mkdir -p output/patrick_max
mkdir -p output/patrick_max_child
mkdir -p output/final

(julia run_difFUBAR_baseline.jl --fasta_file "data/Ace2_tiny_test.fasta" --nexus_file "data/Ace2_no_background.nex" --output "output/baseline/") 2>&1 | tee output/baseline/time_output.txt output/baseline/difFUBAR_output.txt
(julia run_difFUBAR_prune_max.jl --fasta_file "data/Ace2_tiny_test.fasta" --nexus_file "data/Ace2_no_background.nex" --output "output/max/") 2>&1 | tee output/max/time_output.txt output/max/difFUBAR_output.txt
(julia run_difFUBAR_prune_patrick.jl --fasta_file "data/Ace2_tiny_test.fasta" --nexus_file "data/Ace2_no_background.nex" --output "output/patrick/") 2>&1 | tee output/patrick/time_output.txt output/patrick/difFUBAR_output.txt
(julia run_difFUBAR_prune_patrick_max.jl --fasta_file "data/Ace2_tiny_test.fasta" --nexus_file "data/Ace2_no_background.nex" --output "output/patrick_max/") 2>&1 | tee output/patrick_max/time_output.txt output/patrick_max_child/difFUBAR_output.txt
(julia run_difFUBAR_prune_patrick_max_child.jl --fasta_file "data/Ace2_tiny_test.fasta" --nexus_file "data/Ace2_no_background.nex" --output "output/patrick_max_child/") 2>&1 | tee output/patrick_max_child/time_output.txt output/patrick_max_child/difFUBAR_output.txt
(julia run_difFUBAR_prune_final.jl --fasta_file "data/Ace2_tiny_test.fasta" --nexus_file "data/Ace2_no_background.nex" --output "output/final/") 2>&1 | tee output/final/time_output.txt output/final/difFUBAR_output.txt



