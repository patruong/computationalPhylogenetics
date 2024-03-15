#!/bin/bash
output_dir=output/cytb
#input_fasta=data/Cytb.fasta
#input_tree=data/cytb.nwk

#input_tree_no_branch=data/cytb_no_branchlength.nwk

mkdir -p "$output_dir/treesurgery_and_parallel/$group1_$group2"
mkdir -p "$output_dir/contrastfel"

#(time julia run_cytB_argparse.jl --group1 --group2) 2>&1 | tee "$output_dir/treesurgery_and_parallel/$group1_$group2/difFUBAR_output.txt" "$output_dir/treesurgery_and_parallel/$group1_$group2/time_output.txt"

groups=("Leucocytozoon" "mammals" "birds" "Haemoproteidae")

for group1 in "${groups[@]}"; do
    for group2 in "${groups[@]}"; do
        if [ "$group1" != "$group2" ]; then
            mkdir -p "$output_dir/treesurgery_and_parallel/${group1}_${group2}"
            (time julia run_cytB_argparse.jl --group1 "$group1" --group2 "$group2") 2>&1 | tee "$output_dir/treesurgery_and_parallel/${group1}_${group2}/difFUBAR_output.txt" "$output_dir/treesurgery_and_parallel/${group1}_${group2}/time_output.txt"
        fi
    done
done
#(time hyphy contrast-fel --alignment $input_fasta  --tree $input_tree --branch-set Leucocytozoon --branch-set mammals --branch-set birds --branch-set Haemoproteidae --p-value 1.00 --q-value 1.00 --cpu 20 --output $output_dir/contrastfel/contrastfel.FEL.json) 2>&1 | tee "$output_dir/contrastfel/constrastfel_output.txt" "$output_dir/contrastfel/time_output.txt"
#(time hyphy contrast-fel --alignment $input_fasta  --tree $input_tree_no_branch --branch-set Leucocytozoon --branch-set mammals --branch-set birds --branch-set Haemoproteidae --p-value 1.00 --q-value 1.00 --cpu 20 --output $output_dir/contrastfel/contrastfel.FEL.json) 2>&1 | tee "$output_dir/contrastfel/constrastfel_output.txt" "$output_dir/contrastfel/time_output.txt"
