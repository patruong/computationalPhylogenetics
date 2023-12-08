using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using Cairo
using ArgParse
using CSV
using DataFrames



# Parvo data
cd("/home/patrick/git/computationalPhylogenetics/experiments/20231207_contrastFEL_time/data/")
file_dir = "/home/patrick/git/computationalPhylogenetics/experiments/20231207_contrastFEL_time/data/" # original data
fasta = file_dir * "ParvoVP.fasta"
tree_file = file_dir * "datamonkey_tree" * "/" * "tree.new"
analysis_name = "ParvoVP"
seqnames, seqs = read_fasta(fasta);

tree = MolecularEvolution.read_newick_tree(tree_file)

#treestring, tags, tag_colors = import_colored_figtree_nexus_as_tagged_tree(tree_file)
