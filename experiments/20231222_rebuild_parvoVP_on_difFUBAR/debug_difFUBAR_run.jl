using Revise # for code dev
# dev ~/git/CodonMolecularEvolution.jl/ # in package environment
# dev ~/git/MolecularEvolution.jl/ # in package environment
using MolecularEvolution #read_fasta
using CodonMolecularEvolution
using Compose
using ArgParse
using CSV
using DataFrames
using Tables


cd("/home/patrick/git/computationalPhylogenetics/experiments/20231222_rebuild_parvoVP_on_difFUBAR")

analysis_name = "analysis_name"
fasta_file = "parvo_panleu_trans.fasta"
newick_file = "parvo_panleu_trans_tagged_fixed_aligned.nwk"
exports = true
verbosity = 1
iters = 2500
pos_thresh = 0.95
output_csv = "output"

seqnames, seqs = read_fasta(fasta_file)
treestring_group_labeled, treestring, group_tags, original_tags, tag_colors = import_grouped_label_tree(newick_file)

treestring = read(newick_file, String)
tree = gettreefromnewick(treestring, FelNode)

tree.children

treestring_group_labeled, group_tags, original_tags = replace_newick_tags(treestring)
tag_colors = generate_hex_colors(length(original_tags))